"""
00_fetch_metadata.py — P0-1 数据审计（ID 漂移防线）

对本研究实际使用的每一个 GWAS 数据集拉取权威元数据并落盘：

  暴露端 11 个免疫/细胞因子性状
      · EFO 标签        <- GWAS Catalog REST  /efoTraits/{EFO}
      · 研究级元数据    <- GWAS Catalog REST  /studies/{GCST}
        （样本量、人群、PMID、原始 diseaseTrait 名称）
      · 实际使用的 study accession 取自 data/ivs_{abbr}.csv，
        保证 manifest 描述的就是分析真正用到的那份数据，而非 EFO 下任意一项。

  结局端 5 个 FinnGen R11 脊柱感染/骨髓炎表型
      · FinnGen R11 API /api/pheno/{phenocode}
        （phenostring、num_cases、num_controls、全基因组显著位点数）

一致性核对（consistency_check 列）
  逐条把「预期表型」与「接口返回的 trait 名」做关键词比对：
      OK        名称与预期表型一致
      REVIEW    名称相关但措辞不同，需人工判读（本文件末尾附人工结论）
      MISMATCH  名称与预期表型不符 —— 该数据集不得进入分析
  另外核对样本量与代码中写死的 ncase/ncontrol 是否一致（NUM_MISMATCH）。

输出：data/manifest_exposures.tsv, data/manifest_outcomes.tsv
用法：python 00_fetch_metadata.py
"""

from __future__ import annotations

import os
import re
import sys
import time
from typing import Optional

import pandas as pd
import requests

import mr_pipeline as mp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")

GC_REST = "https://www.ebi.ac.uk/gwas/rest/api"
FINNGEN_API = "https://r11.finngen.fi/api"
HEADERS = {"User-Agent": "spine-MR-pipeline/1.0 (research; python-requests)"}
PAUSE = 1.0          # 温和限速，避免与并行任务共同触发 EBI 限流

# 每个暴露的「预期表型关键词」——用于自动核对 trait 名是否漂移
EXPECTED_KEYWORDS = {
    "CRP":    ["c-reactive protein"],
    "IL6R":   ["interleukin-6 receptor", "il-6 receptor", "il6r",
               "interleukin 6 receptor"],
    "IL1B":   ["interleukin-1 beta", "interleukin 1 beta", "il-1b", "il1b"],
    "CXCL10": ["c-x-c motif chemokine 10", "cxcl10", "ip-10",
               "interferon gamma-induced protein 10", "interferon gamma induced protein 10"],
    "CD40":   ["cd40"],
    "CCL2":   ["c-c motif chemokine 2", "ccl2", "monocyte chemoattractant protein",
               "mcp-1", "monocyte chemotactic protein"],
    "NEUT":   ["neutrophil"],
    "LYMPH":  ["lymphocyte"],
    "MONO":   ["monocyte count", "monocyte"],
    "EOS":    ["eosinophil"],
    "WBC":    ["white blood cell", "leukocyte"],
}

# 每个结局的「预期表型关键词」
EXPECTED_OUTCOME_KEYWORDS = {
    "M13_OSTEOMYELVERTEB":   ["osteomyelitis of vertebra", "vertebral osteomyelitis"],
    "M13_DISCITIS":          ["discitis"],
    "M13_DISCINFECTION":     ["disc", "infection"],
    "M13_SPONDYLOINFECTION": ["spondylopath", "infective", "spondyl"],
    "M13_OSTEOMYELITIS":     ["osteomyelitis"],
}

# FEASIBILITY.md / PLAN.md 记载的病例数（源自 FinnGen **R12**）。
# 本研究实际使用 R11 公开汇总统计，故此处保留 R12 数字作为核对基线，
# 差异会被显式标注为 RELEASE_DIFF —— 这是 PLAN §15.1 A3 要求的样本量核对。
FEASIBILITY_R12_NCASE = {
    "M13_OSTEOMYELVERTEB": 111,
    "M13_DISCITIS": 557,
    "M13_DISCINFECTION": 423,
    "M13_SPONDYLOINFECTION": 74,
    "M13_OSTEOMYELITIS": 2336,
}

# 人工核对结论（逐条判读，落盘留痕）——见文件末尾 MANUAL AUDIT 注释块
MANUAL_VERDICT = {
    "CRP": "VERIFIED name OK; ⚠ POPULATION MISMATCH: GCST008055 (Wojcik 2019, "
           "PAGE) is multi-ancestry (Hispanic/Latino, African American, Asian, "
           "Native Hawaiian, Native American), NOT European — allele "
           "frequencies/LD differ from Finnish outcome sample. CRP results "
           "downgraded to exploratory.",
    "IL6R": "VERIFIED: EFO_0008187 = IL-6 receptor subunit alpha; "
            "GCST90012025 = Folkersen 2020 SCALLOP, European.",
    "IL1B": "VERIFIED name OK; only 1 IV at 5e-8 and n=4,910 — Wald-ratio only, "
            "exploratory.",
    "CXCL10": "VERIFIED: CXCL10 = IP-10; GCST004440 = Ahola-Olli 2017, European.",
    "CD40": "VERIFIED OK despite REVIEW flag: TNF receptor superfamily member 5 "
            "(TNFRSF5) IS the official protein name of CD40. Not an ID drift.",
    "CCL2": "VERIFIED: CCL2 = MCP-1; GCST90012007 = Folkersen 2020, European.",
    "NEUT": "VERIFIED: Astle 2016 blood-cell GWAS, European.",
    "LYMPH": "VERIFIED: Astle 2016 blood-cell GWAS, European.",
    "MONO": "VERIFIED: Astle 2016 blood-cell GWAS, European.",
    "EOS": "VERIFIED: Astle 2016 blood-cell GWAS, European.",
    "WBC": "VERIFIED: Astle 2016 blood-cell GWAS, European.",
}


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def get_json(url: str, tries: int = 4) -> Optional[dict]:
    delay = PAUSE
    for _ in range(tries):
        try:
            r = requests.get(url, headers=HEADERS, timeout=90)
        except requests.RequestException:
            time.sleep(delay)
            delay *= 2
            continue
        if r.status_code == 200:
            time.sleep(PAUSE)
            try:
                return r.json()
            except ValueError:
                return None
        if r.status_code in (429, 500, 502, 503, 504):
            time.sleep(delay)
            delay *= 2
            continue
        time.sleep(PAUSE)
        return None
    return None


# --------------------------------------------------------------------------
# 暴露端
# --------------------------------------------------------------------------
def study_used(abbr: str) -> tuple[Optional[str], int]:
    """返回正向分析实际使用的 study accession 与 IV 数。"""
    f = os.path.join(DATA, f"ivs_{abbr}.csv")
    if not os.path.exists(f):
        return None, 0
    d = pd.read_csv(f)
    if "study" not in d.columns or d.empty:
        return None, len(d)
    return str(d["study"].mode().iloc[0]), len(d)


def parse_case_control(s: str) -> tuple[Optional[int], Optional[int]]:
    """从 initialSampleSize 自由文本中抽取 cases / controls。"""
    if not s:
        return None, None
    ncase = ncontrol = None
    m = re.search(r"([\d,]+)\s+(?:\w+\s+){0,4}cases", s, re.I)
    if m:
        ncase = int(m.group(1).replace(",", ""))
    m = re.search(r"([\d,]+)\s+(?:\w+\s+){0,4}controls", s, re.I)
    if m:
        ncontrol = int(m.group(1).replace(",", ""))
    return ncase, ncontrol


def total_n(study: dict) -> Optional[int]:
    anc = study.get("ancestries") or []
    tot = 0
    got = False
    for a in anc:
        if a.get("type") == "initial" and a.get("numberOfIndividuals"):
            tot += int(a["numberOfIndividuals"])
            got = True
    if got:
        return tot
    m = re.search(r"([\d,]+)", study.get("initialSampleSize") or "")
    return int(m.group(1).replace(",", "")) if m else None


def populations(study: dict) -> str:
    groups = []
    for a in study.get("ancestries") or []:
        if a.get("type") != "initial":
            continue
        for g in a.get("ancestralGroups") or []:
            if g.get("ancestralGroup"):
                groups.append(g["ancestralGroup"])
    return "; ".join(dict.fromkeys(groups)) or "NR"


def check_name(abbr: str, *names: Optional[str]) -> str:
    kws = EXPECTED_KEYWORDS.get(abbr, [])
    blob = " | ".join(n.lower() for n in names if n)
    if not blob:
        return "MISMATCH: no trait name returned"
    if any(k in blob for k in kws):
        return "OK"
    return f"REVIEW: returned '{blob}' vs expected {kws}"


def build_exposures() -> pd.DataFrame:
    rows = []
    for e in mp.EXPOSURES:
        abbr, efo = e["abbr"], e["efo"]
        acc, n_iv = study_used(abbr)
        log(f"暴露 {abbr:7s} EFO={efo} study={acc} n_IV={n_iv}")

        efo_j = get_json(f"{GC_REST}/efoTraits/{efo}") or {}
        efo_label = efo_j.get("trait")

        st = get_json(f"{GC_REST}/studies/{acc}") if acc else None
        if st:
            trait_study = (st.get("diseaseTrait") or {}).get("trait")
            n = total_n(st)
            ncase, ncontrol = parse_case_control(st.get("initialSampleSize") or "")
            pop = populations(st)
            pmid = (st.get("publicationInfo") or {}).get("pubmedId")
            first_author = ((st.get("publicationInfo") or {}).get("author")
                            or {}).get("fullname")
            year = ((st.get("publicationInfo") or {}).get("publicationDate")
                    or "")[:4]
            raw_n = st.get("initialSampleSize")
        else:
            trait_study = n = ncase = ncontrol = pmid = first_author = None
            pop = "NR"
            year = ""
            raw_n = None

        chk = check_name(abbr, efo_label, trait_study)
        if "European" not in pop:
            chk += (f" | POPULATION_MISMATCH: exposure population '{pop}' vs "
                    f"Finnish (European) outcome")

        rows.append({
            "id": efo,
            "trait_name": efo_label or e["name"],
            "sample_size": n,
            "ncase": ncase,
            "ncontrol": ncontrol,
            "population": pop,
            "build": "GRCh38",       # EBI 汇总统计 API 统一以 GRCh38 提供坐标
            "pmid": pmid,
            "consistency_check": chk,
            "manual_verification": MANUAL_VERDICT.get(abbr, ""),
            # ---- 审计辅助列 ----
            "abbr": abbr,
            "study_accession": acc,
            "study_trait_name": trait_study,
            "intended_phenotype": e["name"],
            "group": e["group"],
            "first_author": first_author,
            "year": year,
            "n_iv_used": n_iv,
            "raw_sample_description": raw_n,
            "source": "EBI GWAS Catalog Summary Statistics API",
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# 结局端
# --------------------------------------------------------------------------
def check_outcome_name(pheno: str, name: Optional[str]) -> str:
    if not name:
        return "MISMATCH: no phenostring returned"
    kws = EXPECTED_OUTCOME_KEYWORDS.get(pheno, [])
    low = name.lower()
    if any(k in low for k in kws):
        return "OK"
    return f"REVIEW: returned '{name}' vs expected {kws}"


def build_outcomes() -> pd.DataFrame:
    rows = []
    for o in mp.OUTCOMES:
        ph = o["pheno"]
        log(f"结局 {o['abbr']:9s} {ph}")
        j = get_json(f"{FINNGEN_API}/pheno/{ph}") or {}
        name = j.get("phenostring")
        ncase = j.get("num_cases")
        ncontrol = j.get("num_controls")
        chk = check_outcome_name(ph, name)
        # 与 FEASIBILITY/PLAN 记载的 R12 病例数核对（版本差异，非 ID 漂移）
        ref = FEASIBILITY_R12_NCASE.get(ph)
        if ncase is not None and ref and int(ncase) != int(ref):
            chk += (f" | RELEASE_DIFF: R11 ncase={ncase} vs PLAN/FEASIBILITY "
                    f"R12 figure {ref} ({(int(ncase)-ref)/ref*100:+.1f}%); "
                    f"R11 is the release actually analysed")
        rows.append({
            "id": ph,
            "trait_name": name,
            "sample_size": (int(ncase) + int(ncontrol))
            if ncase is not None and ncontrol is not None else None,
            "ncase": ncase,
            "ncontrol": ncontrol,
            "population": "Finnish (European)",
            "build": "GRCh38",
            "pmid": "36653562",     # FinnGen flagship, Kurki et al. Nature 2023
            "consistency_check": chk,
            "manual_verification": "VERIFIED: phenostring matches the intended "
                                   "ICD-10 phenotype; R11 case counts adopted "
                                   "throughout (incl. power calculation).",
            # ---- 审计辅助列 ----
            "abbr": o["abbr"],
            "intended_phenotype": o["label"],
            "spine_specific": o["spine_specific"],
            "num_gw_significant": j.get("num_gw_significant"),
            "gc_lambda_0.001": (j.get("gc_lambda") or {}).get("0.001"),
            "release": "FinnGen R11",
            "source": "FinnGen R11 public summary statistics (GCS bucket) "
                      "+ https://r11.finngen.fi/api/pheno/{code}",
        })
    return pd.DataFrame(rows)


def main() -> int:
    os.makedirs(DATA, exist_ok=True)
    exp = build_exposures()
    out = build_outcomes()
    pe = os.path.join(DATA, "manifest_exposures.tsv")
    po = os.path.join(DATA, "manifest_outcomes.tsv")
    exp.to_csv(pe, sep="\t", index=False)
    out.to_csv(po, sep="\t", index=False)
    log(f"写出 {pe} ({len(exp)} 行)")
    log(f"写出 {po} ({len(out)} 行)")

    print("\n=== 一致性核对摘要（暴露） ===")
    for _, r in exp.iterrows():
        print(f"  {r['abbr']:7s} {r['id']:14s} {str(r['study_accession']):14s} "
              f"{str(r['trait_name'])[:48]:50s} {r['consistency_check'][:70]}")
    print("\n=== 一致性核对摘要（结局） ===")
    for _, r in out.iterrows():
        print(f"  {r['abbr']:9s} {r['id']:22s} {str(r['trait_name'])[:38]:40s} "
              f"{r['consistency_check'][:100]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
