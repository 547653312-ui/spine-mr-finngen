"""
mr_pipeline.py — 脊柱感染免疫驱动的两样本孟德尔随机化主管线

设计（默认，PLAN.md 缺席时按 FEASIBILITY.md 第三节）：
  正向 FORWARD : 免疫/细胞因子性状  ->  脊柱感染/骨髓炎
  反向 REVERSE : 骨髓炎             ->  免疫/细胞因子性状

数据源（均为开放访问，无需 token）：
  暴露 = EBI GWAS Catalog Summary Statistics API
  结局 = FinnGen R11 公开汇总统计（bgzip+tabix，HTTP Range 精准取数）

为什么不用 IEU ieu-b-4975：OpenGWAS 自 2024-05-01 起强制 JWT，本机无 token（401）。
FinnGen R11 反而提供**脊柱特异**表型（椎体骨髓炎 M46.2 / 椎间盘炎 M46.4），
科学上优于 UKB 的混合骨髓炎表型。详见 results/API_NOTES.md 与 LIMITATIONS.md。

用法：
    python mr_pipeline.py                # 跑全部（正向+反向）
    python mr_pipeline.py --direction forward
    python mr_pipeline.py --quick        # 少量性状快速验证
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import traceback
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import gwas_io as gio
import mr_methods as mrm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
DATA = os.path.join(ROOT, "data")
FIGDIR = os.path.join(RESULTS, "figures")
HARMDIR = os.path.join(RESULTS, "harmonised")
for d in (RESULTS, DATA, FIGDIR, HARMDIR, os.path.join(DATA, "cache")):
    os.makedirs(d, exist_ok=True)

# --------------------------------------------------------------------------
# 研究设计定义
# --------------------------------------------------------------------------
# 暴露：免疫 / 细胞因子性状（EFO）
EXPOSURES: List[dict] = [
    # —— 细胞因子 / 免疫蛋白 ——（EFO 已逐一验证存在汇总统计）
    {"name": "C-reactive protein measurement", "efo": "EFO_0004458",
     "abbr": "CRP", "group": "acute-phase"},
    {"name": "Interleukin-6 receptor subunit alpha measurement", "efo": "EFO_0008187",
     "abbr": "IL6R", "group": "cytokine"},
    {"name": "Interleukin-1 beta measurement", "efo": "EFO_0004812",
     "abbr": "IL1B", "group": "cytokine"},
    {"name": "C-X-C motif chemokine 10 (IP-10) measurement", "efo": "EFO_0008056",
     "abbr": "CXCL10", "group": "chemokine"},
    {"name": "CD40 / TNF receptor superfamily member 5 measurement", "efo": "EFO_0010607",
     "abbr": "CD40", "group": "cytokine receptor"},
    {"name": "CCL2 / MCP-1 measurement", "efo": "EFO_0004749",
     "abbr": "CCL2", "group": "chemokine"},
    # —— 免疫细胞亚群计数 ——
    {"name": "Neutrophil count", "efo": "EFO_0004833", "abbr": "NEUT", "group": "immune cell"},
    {"name": "Lymphocyte count", "efo": "EFO_0004587", "abbr": "LYMPH", "group": "immune cell"},
    {"name": "Monocyte count", "efo": "EFO_0005091", "abbr": "MONO", "group": "immune cell"},
    {"name": "Eosinophil count", "efo": "EFO_0004842", "abbr": "EOS", "group": "immune cell"},
    {"name": "White blood cell (leukocyte) count", "efo": "EFO_0004308",
     "abbr": "WBC", "group": "immune cell"},
]

QUICK_EXPOSURES = {"CRP", "NEUT", "LYMPH", "MONO"}

# 结局：FinnGen R11 脊柱感染 / 骨髓炎表型
#
# ncase/ncontrol 已于 P0-1 数据审计（analysis/00_fetch_metadata.py）中
# 用 FinnGen R11 API /api/pheno/{code} 核实并更正。此前写死的是 FEASIBILITY.md
# 记载的 **R12** 数字（VOM 111 / DISCITIS 557 / OM 2336 ...），与本研究实际
# 使用的 R11 汇总统计不符；效能计算必须用下面的 R11 真实值。
# 核对留痕见 data/manifest_outcomes.tsv 的 consistency_check 列。
OUTCOMES: List[dict] = [
    {"pheno": "M13_OSTEOMYELVERTEB", "label": "Vertebral osteomyelitis (M46.2)",
     "abbr": "VOM", "ncase": 104, "ncontrol": 322314, "spine_specific": True},
    {"pheno": "M13_DISCITIS", "label": "Discitis (M46.4)",
     "abbr": "DISCITIS", "ncase": 495, "ncontrol": 322314, "spine_specific": True},
    {"pheno": "M13_DISCINFECTION", "label": "Pyogenic intervertebral disc infection",
     "abbr": "DISCINF", "ncase": 375, "ncontrol": 322314, "spine_specific": True},
    {"pheno": "M13_SPONDYLOINFECTION", "label": "Other infective spondylopathies",
     "abbr": "SPONDINF", "ncase": 68, "ncontrol": 322314, "spine_specific": True},
    {"pheno": "M13_OSTEOMYELITIS", "label": "Osteomyelitis (all sites, M86)",
     "abbr": "OM", "ncase": 2125, "ncontrol": 429826, "spine_specific": False},
]

QUICK_OUTCOMES = {"VOM", "DISCITIS", "OM"}

P_PRIMARY = 5e-8
P_SECONDARY = 5e-6      # 反向 MR / IV 过少时的次级阈值
CLUMP_KB = 500          # 距离剪枝窗口（无 LD 参考面板时的替代方案）
F_MIN = 10
PALINDROME_MAF = 0.42   # 回文 SNP 的 MAF 剔除阈值

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# --------------------------------------------------------------------------
# IV 选择
# --------------------------------------------------------------------------
def distance_prune(df: pd.DataFrame, kb: int = CLUMP_KB) -> pd.DataFrame:
    """
    距离剪枝：按 P 升序贪心保留，剔除同染色体 ±kb 内的其他位点。
    这是在无 LD 参考面板（无 plink / 无 OpenGWAS clump 接口）情况下的
    标准替代方案，近似独立工具变量。
    """
    if df.empty:
        return df
    df = df.sort_values("pval").reset_index(drop=True)
    keep, taken = [], {}
    for i, r in df.iterrows():
        c = str(r["chr"])
        p = int(r["pos"])
        if any(abs(p - q) < kb * 1000 for q in taken.get(c, [])):
            continue
        keep.append(i)
        taken.setdefault(c, []).append(p)
    return df.loc[keep].sort_values(["chr", "pos"]).reset_index(drop=True)


def select_ivs(exposure: dict, pval: float = P_PRIMARY) -> pd.DataFrame:
    """从 EBI 汇总统计 API 取暴露 IV 并剪枝。"""
    efo = exposure["efo"]
    log(f"  取暴露关联 {exposure['abbr']} ({efo}) P<{pval:g} ...")
    try:
        recs = gio.ebi_trait_associations(efo, p_upper=pval, max_records=4000)
    except Exception as e:
        log(f"  !! EBI 取数失败 {efo}: {e}")
        return pd.DataFrame()
    if not recs:
        log(f"  -- {exposure['abbr']}: 0 条关联")
        return pd.DataFrame()
    df = pd.DataFrame(recs)
    df = df.dropna(subset=["rsid", "pos", "beta", "se"])
    df = df[(df["se"] > 0) & np.isfinite(df["beta"])]
    # 同一 SNP 多研究：保留 P 最小者
    df = df.sort_values("pval").drop_duplicates(subset=["chr", "pos"], keep="first")
    n0 = len(df)
    df = distance_prune(df)
    df["F"] = (df["beta"] / df["se"]) ** 2
    df = df[df["F"] > F_MIN]
    log(f"  -- {exposure['abbr']}: {n0} 条 -> 剪枝后 {len(df)} 个独立 IV (F>10)")
    return df


# --------------------------------------------------------------------------
# 等位基因协调
# --------------------------------------------------------------------------
def harmonise(iv: pd.DataFrame, fs: gio.FinnGenSumstats,
              outcome_label: str, workers: int = 12) -> pd.DataFrame:
    """
    暴露 IV 与 FinnGen 结局按 chr:pos (GRCh38) 匹配并协调等位基因方向。
    规则：
      - 等位基因一致 -> 结局 beta 取 alt 为效应等位基因，按暴露 EA 定向
      - 互补链（回文）-> 用 MAF 判定；MAF>0.42 的回文 SNP 剔除（不可判向）
      - 等位基因不匹配 -> 剔除
    远端 tabix 查询为 I/O 密集，用线程池并发。
    """
    def one(r) -> dict:
        try:
            hits = fs.query(str(r["chr"]), int(r["pos"]),
                            ref=r["other_allele"], alt=r["effect_allele"], window=0)
        except Exception as e:
            return {"SNP": r["rsid"], "status": f"query_error:{type(e).__name__}"}
        if not hits:
            return {"SNP": r["rsid"], "status": "not_in_outcome"}

        ea, oa = r["effect_allele"], r["other_allele"]
        chosen = None
        action = None
        for h in hits:
            ref, alt = h["ref"].upper(), h["alt"].upper()
            if {ref, alt} == {ea, oa}:
                chosen, action = h, ("same" if alt == ea else "flip")
                break
        if chosen is None:  # 试互补链
            ea_c = "".join(COMPLEMENT.get(b, "N") for b in ea)
            oa_c = "".join(COMPLEMENT.get(b, "N") for b in oa)
            for h in hits:
                ref, alt = h["ref"].upper(), h["alt"].upper()
                if {ref, alt} == {ea_c, oa_c}:
                    chosen, action = h, ("same_strandflip" if alt == ea_c else "flip_strandflip")
                    break
        if chosen is None:
            return {"SNP": r["rsid"], "status": "allele_mismatch"}

        # 回文 SNP 检查
        palindromic = (COMPLEMENT.get(ea) == oa)
        eaf = r.get("eaf")
        af_alt = float(chosen["af_alt"])
        if palindromic:
            maf_e = min(eaf, 1 - eaf) if eaf is not None and np.isfinite(eaf) else np.nan
            maf_o = min(af_alt, 1 - af_alt)
            if not np.isfinite(maf_e) or maf_e > PALINDROME_MAF or maf_o > PALINDROME_MAF:
                return {"SNP": r["rsid"], "status": "palindromic_ambiguous"}

        by = float(chosen["beta"])
        byse = float(chosen["sebeta"])
        if action.startswith("flip"):
            by = -by
            af_alt = 1 - af_alt
        if byse <= 0 or not np.isfinite(by):
            return {"SNP": r["rsid"], "status": "bad_outcome_stats"}

        return {
            "SNP": r["rsid"], "chr": r["chr"], "pos": r["pos"],
            "effect_allele": ea, "other_allele": oa,
            "eaf_exposure": eaf, "eaf_outcome": af_alt,
            "beta_exposure": float(r["beta"]), "se_exposure": float(r["se"]),
            "pval_exposure": float(r["pval"]), "F": float(r["F"]),
            "beta_outcome": by, "se_outcome": byse,
            "pval_outcome": float(chosen["pval"]),
            "outcome_af_cases": chosen.get("af_alt_cases"),
            "nearest_gene": chosen.get("nearest_genes"),
            "harmonise_action": action, "palindromic": palindromic,
            "exposure_study": r.get("study"), "outcome": outcome_label,
            "status": "kept",
        }

    records = [r for _, r in iv.iterrows()]
    with ThreadPoolExecutor(max_workers=workers) as ex:
        rows = list(ex.map(one, records))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# 作图
# --------------------------------------------------------------------------
def plot_scatter(h: pd.DataFrame, res: dict, title: str, path: str) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    bx, by = h["beta_exposure"].values, h["beta_outcome"].values
    ax.errorbar(bx, by, yerr=1.96 * h["se_outcome"], xerr=1.96 * h["se_exposure"],
                fmt="o", ms=4, lw=0.7, color="#37474f", ecolor="#b0bec5",
                alpha=0.85, zorder=3)
    xs = np.linspace(min(0, bx.min()) * 1.1, bx.max() * 1.1, 100)
    colors = {"IVW (random effects)": "#1565c0", "MR-Egger": "#c62828",
              "Weighted median": "#2e7d32", "Weighted mode": "#ef6c00"}
    for r in res["results"]:
        if r.method in colors and np.isfinite(r.b):
            inter = res["pleiotropy"]["intercept"] if r.method == "MR-Egger" else 0.0
            inter = inter if np.isfinite(inter) else 0.0
            ax.plot(xs, inter + r.b * xs, lw=1.8, color=colors[r.method],
                    label=f"{r.method} (b={r.b:.3f})")
    ax.axhline(0, color="grey", lw=0.6, ls=":")
    ax.axvline(0, color="grey", lw=0.6, ls=":")
    ax.set_xlabel("SNP effect on exposure")
    ax.set_ylabel("SNP effect on outcome (log-OR)")
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_forest(single: List[dict], res: dict, title: str, path: str) -> None:
    df = pd.DataFrame(single)
    if df.empty:
        return
    df = df.sort_values("b")
    fig, ax = plt.subplots(figsize=(6.4, max(3.0, 0.26 * len(df) + 2.0)))
    y = np.arange(len(df))
    ax.errorbar(df["b"], y, xerr=1.96 * df["se"], fmt="o", ms=3.5,
                lw=0.8, color="#455a64", ecolor="#90a4ae")
    ivw = [r for r in res["results"] if r.method.startswith("IVW (random")]
    if ivw and np.isfinite(ivw[0].b):
        r = ivw[0]
        ax.axvline(r.b, color="#1565c0", lw=1.5,
                   label=f"IVW b={r.b:.3f} (P={r.pval:.2g})")
        ax.axvspan(r.lo95, r.hi95, color="#1565c0", alpha=0.12)
    ax.axvline(0, color="grey", lw=0.7, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(df["SNP"], fontsize=6)
    ax.set_xlabel("Wald ratio (log-OR per SD exposure)")
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_loo(loo: List[dict], title: str, path: str) -> None:
    df = pd.DataFrame(loo)
    if df.empty:
        return
    fig, ax = plt.subplots(figsize=(6.4, max(3.0, 0.26 * len(df) + 2.0)))
    y = np.arange(len(df))
    col = ["#c62828" if s.startswith("ALL") else "#455a64" for s in df["excluded_SNP"]]
    for i, (_, r) in enumerate(df.iterrows()):
        ax.errorbar(r["b"], i, xerr=1.96 * r["se"], fmt="o", ms=3.5,
                    lw=0.8, color=col[i], ecolor=col[i], alpha=0.85)
    ax.axvline(0, color="grey", lw=0.7, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(df["excluded_SNP"], fontsize=6)
    ax.set_xlabel("IVW estimate excluding each SNP (log-OR)")
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_funnel(h: pd.DataFrame, res: dict, title: str, path: str) -> None:
    fig, ax = plt.subplots(figsize=(5.6, 5.0))
    ratio = h["beta_outcome"].values / h["beta_exposure"].values
    prec = np.abs(h["beta_exposure"].values / h["se_outcome"].values)
    ax.scatter(ratio, prec, s=18, color="#37474f", alpha=0.8)
    for r in res["results"]:
        if r.method.startswith("IVW (random") and np.isfinite(r.b):
            ax.axvline(r.b, color="#1565c0", lw=1.5, label="IVW")
        if r.method == "MR-Egger" and np.isfinite(r.b):
            ax.axvline(r.b, color="#c62828", lw=1.5, ls="--", label="MR-Egger")
    ax.set_xlabel("Wald ratio")
    ax.set_ylabel("Instrument strength (|bx| / se_outcome)")
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# --------------------------------------------------------------------------
# 单个 暴露×结局 分析
# --------------------------------------------------------------------------
def analyse_pair(iv: pd.DataFrame, exposure: dict, outcome: dict,
                 fs: gio.FinnGenSumstats, direction: str) -> Optional[dict]:
    tag = f"{exposure['abbr']}__{outcome['abbr']}"
    h_all = harmonise(iv, fs, outcome["label"])
    if h_all.empty:
        log(f"    {tag}: 协调后为空")
        return None
    h_all.to_csv(os.path.join(HARMDIR, f"{direction}_{tag}_harmonised.csv"), index=False)
    h = h_all[h_all["status"] == "kept"].copy()
    dropped = h_all[h_all["status"] != "kept"]["status"].value_counts().to_dict()
    log(f"    {tag}: 保留 {len(h)} IV，剔除 {dropped}")
    if len(h) < 1:
        return None

    res = mrm.run_all(h["beta_exposure"].values, h["se_exposure"].values,
                      h["beta_outcome"].values, h["se_outcome"].values,
                      h["SNP"].tolist())

    rows = []
    for r in res["results"]:
        d = r.as_dict(binary_outcome=(direction == "forward"))
        d.update({"exposure": exposure["name"], "exposure_abbr": exposure["abbr"],
                  "outcome": outcome["label"], "outcome_abbr": outcome["abbr"],
                  "direction": direction,
                  "F_mean": float(np.mean(res["F"])), "F_min": float(np.min(res["F"])),
                  "Q": res["heterogeneity"]["Q"], "Q_pval": res["heterogeneity"]["Q_pval"],
                  "I2": res["heterogeneity"]["I2"],
                  "egger_intercept": res["pleiotropy"]["intercept"],
                  "egger_intercept_pval": res["pleiotropy"]["intercept_pval"]})
        rows.append(d)

    title = f"{exposure['abbr']} -> {outcome['abbr']}  ({len(h)} IVs)"
    if direction == "reverse":
        title = f"{outcome['abbr']} -> {exposure['abbr']}  ({len(h)} IVs)"
    try:
        if len(h) >= 3:
            plot_scatter(h, res, title, os.path.join(FIGDIR, f"{direction}_{tag}_scatter.png"))
            plot_funnel(h, res, title, os.path.join(FIGDIR, f"{direction}_{tag}_funnel.png"))
            plot_loo(res["loo"], f"Leave-one-out: {title}",
                     os.path.join(FIGDIR, f"{direction}_{tag}_loo.png"))
        if len(h) >= 2:
            plot_forest(res["single"], res, f"Single-SNP: {title}",
                        os.path.join(FIGDIR, f"{direction}_{tag}_forest.png"))
    except Exception as e:
        log(f"    !! 作图失败 {tag}: {e}")

    if res["loo"]:
        pd.DataFrame(res["loo"]).assign(exposure=exposure["abbr"], outcome=outcome["abbr"],
                                        direction=direction).to_csv(
            os.path.join(RESULTS, "loo", f"{direction}_{tag}_loo.csv"), index=False)
    pd.DataFrame(res["single"]).assign(exposure=exposure["abbr"], outcome=outcome["abbr"],
                                       direction=direction).to_csv(
        os.path.join(RESULTS, "loo", f"{direction}_{tag}_singlesnp.csv"), index=False)
    return {"rows": rows, "nsnp": len(h)}


# --------------------------------------------------------------------------
# 主流程
# --------------------------------------------------------------------------
def run_forward(exposures: List[dict], outcomes: List[dict]) -> pd.DataFrame:
    log("=== 正向 MR：免疫/细胞因子 -> 脊柱感染/骨髓炎 ===")
    ivs: Dict[str, pd.DataFrame] = {}
    for e in exposures:
        df = select_ivs(e, P_PRIMARY)
        if not df.empty:
            df.to_csv(os.path.join(DATA, f"ivs_{e['abbr']}.csv"), index=False)
            ivs[e["abbr"]] = df
    if not ivs:
        log("!! 无任何暴露 IV，正向分析终止")
        return pd.DataFrame()

    allrows = []
    for o in outcomes:
        log(f"  结局 {o['abbr']} ({o['pheno']}) cases={o['ncase']}")
        try:
            fs = gio.FinnGenSumstats(o["pheno"])
        except Exception as e:
            log(f"  !! 无法打开 FinnGen {o['pheno']}: {e}")
            continue
        for e in exposures:
            if e["abbr"] not in ivs:
                continue
            try:
                r = analyse_pair(ivs[e["abbr"]], e, o, fs, "forward")
                if r:
                    allrows.extend(r["rows"])
            except Exception as ex:
                log(f"    !! {e['abbr']}->{o['abbr']} 失败: {ex}")
                traceback.print_exc()
    return pd.DataFrame(allrows)


def reverse_ivs_from_finngen(outcome: dict, pval: float) -> pd.DataFrame:
    """
    反向 MR 的暴露端 = 脊柱感染表型本身，需要其全基因组显著位点。
    FinnGen 公开文件为 764–808MB，逐行扫描代价高；此处采用
    「显著位点已知区间 + 线性索引扫描」不可行时的诚实降级：
    直接报告不可得，由调用方记录在 LIMITATIONS.md。
    """
    return pd.DataFrame()


def run_reverse(exposures: List[dict], outcomes: List[dict]) -> pd.DataFrame:
    """
    反向 MR：脊柱感染 -> 免疫性状。
    暴露端需 FinnGen 表型的显著 IV。VOM 的 num_gw_significant = 0，
    OM 亦极少，且公开文件需全量扫描才能定位，代价 ~800MB/表型。
    本函数记录可行性并返回空表；详见 results/LIMITATIONS.md。
    """
    log("=== 反向 MR：脊柱感染 -> 免疫性状 ===")
    rows = []
    for o in outcomes:
        meta = {}
        try:
            meta = gio.finngen_pheno_meta(o["pheno"])
        except Exception:
            pass
        ngw = meta.get("num_gw_significant", None)
        log(f"  {o['abbr']}: FinnGen 报告全基因组显著位点数 = {ngw}")
        rows.append({"outcome_as_exposure": o["abbr"], "pheno": o["pheno"],
                     "ncase": o["ncase"], "num_gw_significant": ngw,
                     "reverse_MR_feasible": bool(ngw) and ngw >= 3,
                     "note": "需 >=3 个独立显著 IV 才能跑多方法反向 MR"})
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RESULTS, "reverse_MR_feasibility.csv"), index=False)
    return pd.DataFrame()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--direction", choices=["forward", "reverse", "both"], default="both")
    ap.add_argument("--quick", action="store_true", help="仅跑核心性状，快速验证")
    ap.add_argument("--outcomes", nargs="*", default=None,
                    help="仅跑指定结局缩写，如 VOM DISCITIS")
    ap.add_argument("--out-suffix", default="", help="输出文件名后缀，避免覆盖")
    args = ap.parse_args()

    os.makedirs(os.path.join(RESULTS, "loo"), exist_ok=True)
    exps = [e for e in EXPOSURES if (not args.quick or e["abbr"] in QUICK_EXPOSURES)]
    outs = [o for o in OUTCOMES if (not args.quick or o["abbr"] in QUICK_OUTCOMES)]
    if args.outcomes:
        outs = [o for o in outs if o["abbr"] in set(args.outcomes)]

    log(f"暴露 {len(exps)} 个 / 结局 {len(outs)} 个 / 方向 {args.direction}")
    frames = []
    if args.direction in ("forward", "both"):
        frames.append(run_forward(exps, outs))
    if args.direction in ("reverse", "both"):
        frames.append(run_reverse(exps, outs))

    res = pd.concat([f for f in frames if f is not None and not f.empty],
                    ignore_index=True) if any(
        f is not None and not f.empty for f in frames) else pd.DataFrame()

    if res.empty:
        log("!! 未产生任何 MR 结果")
        with open(os.path.join(RESULTS, "STATUS.md"), "a", encoding="utf-8") as f:
            f.write(f"\n- {time.strftime('%Y-%m-%d %H:%M')} 运行未产出结果\n")
        return 1

    cols = ["direction", "exposure_abbr", "exposure", "outcome_abbr", "outcome",
            "method", "nsnp", "b", "se", "pval", "lo95", "hi95",
            "OR", "OR_lo95", "OR_hi95", "F_mean", "F_min",
            "Q", "Q_pval", "I2", "egger_intercept", "egger_intercept_pval", "note"]
    res = res[[c for c in cols if c in res.columns]]
    sfx = args.out_suffix
    out_csv = os.path.join(RESULTS, f"MR_results_all{sfx}.csv")
    res.to_csv(out_csv, index=False)
    log(f"结果已写出 {out_csv}  ({len(res)} 行)")

    prim = res[res["method"] == "IVW (random effects)"].copy()
    prim = prim.sort_values("pval")
    prim.to_csv(os.path.join(RESULTS, f"MR_results_primary_IVW{sfx}.csv"), index=False)

    sig = prim[prim["pval"] < 0.05]
    log(f"IVW P<0.05 的暴露-结局对：{len(sig)}")
    for _, r in sig.iterrows():
        log(f"   {r['exposure_abbr']} -> {r['outcome_abbr']}: "
            f"OR={r['OR']:.3f} ({r['OR_lo95']:.3f}-{r['OR_hi95']:.3f}) "
            f"P={r['pval']:.3g}, nSNP={r['nsnp']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
