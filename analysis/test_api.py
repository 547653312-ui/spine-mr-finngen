"""
test_api.py — 数据源连通性探测，产出 results/API_NOTES.md

探测对象：
  1. IEU / OpenGWAS（旧域名 gwas.mrcieu.ac.uk + 新域名 api.opengwas.io）
     目标 ID：ieu-b-4975（骨髓炎）、met-b-*（免疫细胞）、prot-a-*（免疫蛋白）
  2. EBI GWAS Catalog REST + Summary Statistics API
  3. FinnGen R11/R12 公开汇总统计（bucket + tabix + browser API）
"""

from __future__ import annotations

import json
import os
import time
from typing import List

import requests

import gwas_io as gio

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
os.makedirs(RESULTS, exist_ok=True)

H = {"User-Agent": "spine-MR-pipeline/1.0 (research)"}
IEU_IDS = ["ieu-b-4975", "met-b-1", "met-b-10", "prot-a-1", "prot-a-670", "ieu-b-49720"]


def probe(url: str, note: str = "", timeout: int = 30, head: bool = False,
          redirects: bool = True) -> dict:
    t = time.time()
    try:
        if head:
            r = requests.head(url, headers=H, timeout=timeout, allow_redirects=redirects)
            body = ""
        else:
            r = requests.get(url, headers=H, timeout=timeout, allow_redirects=redirects)
            body = r.text[:260].replace("\n", " ")
        return {"url": url, "status": r.status_code, "ms": int((time.time() - t) * 1000),
                "note": note, "body": body,
                "size": r.headers.get("Content-Length", "")}
    except Exception as e:
        return {"url": url, "status": None, "ms": int((time.time() - t) * 1000),
                "note": note, "body": f"EXCEPTION {type(e).__name__}: {e}", "size": ""}


def main() -> None:
    rows: List[dict] = []

    # ---- 1. IEU / OpenGWAS ----
    rows.append(probe("https://gwas.mrcieu.ac.uk/api/gwasinfo/ieu-b-4975",
                      "IEU 旧端点（不跟随重定向）", redirects=False))
    rows.append(probe("https://gwas.mrcieu.ac.uk/api/gwasinfo/ieu-b-4975",
                      "IEU 旧端点（跟随重定向）"))
    rows.append(probe("https://api.opengwas.io/api/status", "OpenGWAS 服务状态"))
    for gid in IEU_IDS:
        rows.append(probe(f"https://api.opengwas.io/api/gwasinfo/{gid}",
                          f"OpenGWAS gwasinfo {gid}（无 token）"))
    rows.append(probe("https://api.opengwas.io/api/gwasinfo",
                      "OpenGWAS gwasinfo 列表（无 token）"))
    rows.append(probe("https://api.opengwas.io/api/associations/ieu-b-4975/rs1205",
                      "OpenGWAS 单点关联（无 token）"))

    # ---- 2. EBI GWAS Catalog ----
    rows.append(probe("https://www.ebi.ac.uk/gwas/rest/api/efoTraits/MONDO_0005246",
                      "GWAS Catalog REST：骨髓炎性状"))
    rows.append(probe("https://www.ebi.ac.uk/gwas/summary-statistics/api/",
                      "GWAS Catalog 汇总统计 API 根"))
    rows.append(probe("https://www.ebi.ac.uk/gwas/summary-statistics/api/traits/"
                      "EFO_0004458/associations?p_upper=0.00000005&size=5",
                      "汇总统计：CRP P<5e-8"))

    # ---- 3. FinnGen ----
    rows.append(probe("https://r12.finngen.fi/api/phenos", "FinnGen R12 表型清单"))
    rows.append(probe("https://r11.finngen.fi/api/pheno/M13_OSTEOMYELVERTEB",
                      "FinnGen R11 表型元数据 VOM"))
    for rel, pheno in [("r12", "M13_OSTEOMYELVERTEB"), ("r11", "M13_OSTEOMYELVERTEB"),
                       ("r11", "M13_DISCITIS"), ("r11", "M13_OSTEOMYELITIS")]:
        R = rel.upper()
        rows.append(probe(
            f"https://storage.googleapis.com/finngen-public-data-{rel}/summary_stats/"
            f"finngen_{R}_{pheno}.gz",
            f"FinnGen {R} 汇总统计 {pheno}", head=True))
        rows.append(probe(
            f"https://storage.googleapis.com/finngen-public-data-{rel}/summary_stats/"
            f"finngen_{R}_{pheno}.gz.tbi",
            f"FinnGen {R} tabix 索引 {pheno}", head=True))

    # ---- 4. 实测 tabix 随机访问 ----
    tabix_ok, tabix_msg = False, ""
    try:
        t = time.time()
        fs = gio.FinnGenSumstats("M13_OSTEOMYELVERTEB")
        hits = fs.query("1", 159684665, window=300)
        tabix_ok = True
        tabix_msg = (f"打开 {fs.bgzf.size/1e6:.0f}MB 远端文件，"
                     f"chr1:159684665±300 命中 {len(hits)} 行，耗时 {time.time()-t:.1f}s")
    except Exception as e:
        tabix_msg = f"FAILED {type(e).__name__}: {e}"

    # ---- 输出 ----
    with open(os.path.join(RESULTS, "api_probe_raw.json"), "w", encoding="utf-8") as f:
        json.dump({"probes": rows, "tabix_random_access": {"ok": tabix_ok, "msg": tabix_msg}},
                  f, indent=1, ensure_ascii=False)

    ok = [r for r in rows if r["status"] and 200 <= r["status"] < 300]
    print(f"探测 {len(rows)} 个端点，可达 {len(ok)} 个")
    print(f"tabix 随机访问: {tabix_ok} — {tabix_msg}")
    for r in rows:
        print(f"  [{r['status']}] {r['note']}")


if __name__ == "__main__":
    main()
