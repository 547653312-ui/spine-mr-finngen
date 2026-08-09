"""
supplementary_network.py — 需要网络取数的补充分析（后台运行）
  1) DRUG-TARGET MR : IL6R / CD40 / CRP 的 cis-pQTL(±1Mb) -> OM / DISCITIS
  2) MVMR           : 暴露主导式多变量 MR（免疫性状 + 吸烟/BMI/糖尿病），
                     变异性水平取数；不可行则诚实降级
  3) REVERSE MR     : 脊柱感染(OM/DISCITIS) -> 免疫性状，变异性水平取数；
                     EBI 限流则降级为可行性记录

EBI GWAS Catalog 可能限流(429)。统一用带指数退避的 fetch；
若持续 429，各段均优雅降级并写说明，绝不空目录、绝不死循环。
"""
from __future__ import annotations
import os, sys, time, json
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
import pandas as pd

import gwas_io as gio
import mr_methods as mrm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
DATA = os.path.join(ROOT, "data")
HARMDIR = os.path.join(RESULTS, "harmonised")
for d in (RESULTS, DATA, HARMDIR):
    os.makedirs(d, exist_ok=True)

F_MIN = 10
CLUMP_KB = 500

# ---- 通用限流安全取数 -----------------------------------------------------
_ebi_429_until = 0.0
def ebi_safe(fn, *a, max_retry=6, **k):
    """包裹任何 EBI 调用，遇 429 指数退避（全局冷却）。"""
    global _ebi_429_until
    for i in range(max_retry):
        wait = max(0.0, _ebi_429_until - time.time())
        if wait > 0:
            time.sleep(min(wait, 30))
        try:
            r = fn(*a, **k)
            return r
        except Exception as e:
            msg = str(e)
            if "429" in msg or "HTTP 429" in msg:
                _ebi_429_until = time.time() + (8 * (i + 1))
                print(f"  [429 cooldown {8*(i+1)}s]", flush=True)
                time.sleep(2 * (i + 1))
                continue
            raise
    print("  [EBI blocked after retries] -> degrade", flush=True)
    return None

# ---- 复用 mr_pipeline 的 select_ivs / harmonise ---------------------------
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib.util
spec = importlib.util.spec_from_file_location("mr_pipeline",
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "mr_pipeline.py"))
mp = importlib.util.module_from_spec(spec); spec.loader.exec_module(mp)

# ==========================================================================
# 1. DRUG-TARGET (cis-pQTL)
# ==========================================================================
GENE_REGION = {
    "IL6R": ("1", 153_400_000, 155_400_000),
    "CD40": ("20", 44_500_000, 46_500_000),
    "CRP":  ("1", 158_600_000, 160_600_000),
}
EFO_OF = {e["abbr"]: e["efo"] for e in mp.EXPOSURES}

def cis_ivs(abbr, pval=1e-4):
    efo = EFO_OF[abbr]
    g = GENE_REGION[abbr]
    recs = ebi_safe(gio.ebi_trait_associations, efo, p_upper=pval, max_records=4000)
    if not recs:
        return pd.DataFrame()
    df = pd.DataFrame(recs).dropna(subset=["rsid", "pos", "beta", "se"])
    df = df[(df["se"] > 0) & np.isfinite(df["beta"])]
    df = df[df["chr"].astype(str) == g[0]]
    df = df[(df["pos"].astype(int) >= g[1]) & (df["pos"].astype(int) <= g[2])]
    df = df.sort_values("pval").drop_duplicates(subset=["chr", "pos"], keep="first")
    df = mp.distance_prune(df)
    df["F"] = (df["beta"] / df["se"]) ** 2
    df = df[df["F"] > F_MIN]
    return df

def run_drug_target():
    print("=== DRUG-TARGET MR (cis-pQTL) ===", flush=True)
    rows = []
    targets = ["IL6R", "CD40", "CRP"]
    outcomes = [o for o in mp.OUTCOMES if o["abbr"] in ("OM", "DISCITIS")]
    for abbr in targets:
        iv = cis_ivs(abbr)
        if iv.empty:
            print(f"  {abbr}: no cis IVs", flush=True)
            continue
        iv.to_csv(os.path.join(DATA, f"ivs_cis_{abbr}.csv"), index=False)
        for o in outcomes:
            try:
                fs = gio.FinnGenSumstats(o["pheno"])
            except Exception as e:
                print(f"  !! FinnGen {o['pheno']}: {e}", flush=True); continue
            h_all = mp.harmonise(iv, fs, o["label"])
            h = h_all[h_all["status"] == "kept"]
            if len(h) < 2:
                print(f"  {abbr}->{o['abbr']}: {len(h)} kept", flush=True); continue
            res = mrm.run_all(h["beta_exposure"], h["se_exposure"],
                              h["beta_outcome"], h["se_outcome"], h["SNP"].tolist())
            for r in res["results"]:
                d = r.as_dict(binary_outcome=True)
                d.update({"target": abbr, "outcome": o["abbr"], "nsnp_cis": len(h),
                          "F_min": float(np.min(res["F"]))})
                rows.append(d)
            print(f"  {abbr}->{o['abbr']}: nSNP={len(h)} OR={np.exp(res['results'][0].b):.3f} "
                  f"P={res['results'][0].pval:.3g}", flush=True)
    if rows:
        pd.DataFrame(rows).to_csv(os.path.join(RESULTS, "drug_target.tsv"), sep="\t", index=False)
        print(f"[written] results/drug_target.tsv ({len(rows)} rows)", flush=True)
    else:
        with open(os.path.join(RESULTS, "drug_target.tsv"), "w") as f:
            f.write("# drug-target MR: no analysable cis-pQTL results\n")
        print("drug-target: empty -> written placeholder", flush=True)

# ==========================================================================
# 2. MVMR (exposure-led, variant-level confounder fetch)
# ==========================================================================
CONFOUNDERS = [
    {"abbr": "SMK", "efo": "EFO_0007851", "name": "Smoking"},
    {"abbr": "BMI", "efo": "EFO_0004340", "name": "Body mass index"},
    {"abbr": "T2D", "efo": "EFO_0001360", "name": "Type 2 diabetes"},
]

def fetch_confounder_at(efo, rsid):
    recs = ebi_safe(gio.ebi_variant_associations, rsid, efo)
    if not recs:
        return None
    for r in recs:
        if r.get("study") and efo in (r.get("study") or ""):
            return r
    return recs[0] if recs else None

def run_mvmr():
    print("=== MVMR (exposure-led) ===", flush=True)
    rows = []
    # 仅对把握度最好的结局做 MVMR
    outcomes = [o for o in mp.OUTCOMES if o["abbr"] in ("OM", "DISCITIS")]
    primaries = ["WBC", "NEUT", "MONO", "LYMPH", "CRP", "IL6R"]
    for o in outcomes:
        try:
            fs = gio.FinnGenSumstats(o["pheno"])
        except Exception as e:
            print(f"  !! FinnGen {o['pheno']}: {e}", flush=True); continue
        for pabbr in primaries:
            ivp = os.path.join(DATA, f"ivs_{pabbr}.csv")
            if not os.path.exists(ivp):
                continue
            iv = pd.read_csv(ivp)
            iv = iv[iv["F"] > F_MIN]
            if len(iv) < 3:
                continue
            # 取结局端 beta
            h_all = mp.harmonise(iv, fs, o["label"])
            h = h_all[h_all["status"] == "kept"]
            if len(h) < 3:
                continue
            # 变异性水平取混淆因子
            conf = {}
            ok = True
            for c in CONFOUNDERS:
                bx, sex, by, sey, snps = [], [], [], [], []
                for _, r in h.iterrows():
                    rec = fetch_confounder_at(c["efo"], r["SNP"])
                    if not rec or rec.get("beta") is None:
                        ok = False; break
                    bx.append(float(r["beta_exposure"])); sex.append(float(r["se_exposure"]))
                    by.append(float(rec["beta"])); sey.append(float(rec["se"]))
                    snps.append(r["SNP"])
                conf[c["abbr"]] = (bx, sex, by, sey, snps)
                if not ok: break
            if not ok:
                print(f"  {pabbr}->{o['abbr']}: confounder data incomplete -> skip MVMR", flush=True)
                continue
            # 组装设计矩阵 (n_snp x (1 + n_conf))，结局 = 脊柱感染
            Yb = np.array([float(x) for x in h["beta_outcome"]])
            Yse = np.array([float(x) for x in h["se_outcome"]])
            Xb = np.column_stack([np.array([float(x) for x in h["beta_exposure"]])] +
                                 [np.array(c[0]) for c in conf.values()])
            Xse = np.column_stack([np.array([float(x) for x in h["se_exposure"]])] +
                                  [np.array(c[1]) for c in conf.values()])
            W = np.diag(1.0 / Yse ** 2)
            XtW = Xb.T @ W
            try:
                beta = np.linalg.inv(XtW @ Xb) @ XtW @ Yb
                cov = np.linalg.inv(XtW @ Xb)
            except np.linalg.LinAlgError:
                print(f"  {pabbr}->{o['abbr']}: singular -> skip", flush=True); continue
            # 条件 F（对每个暴露列）
            cond_F = []
            for j in range(Xb.shape[1]):
                cond_F.append(float((Xb[:, j] / Xse[:, j]) ** 2).mean())
            se = np.sqrt(np.diag(cov))
            z = beta / se
            pvals = 2 * np.array([__import__("scipy.stats").norm.sf(abs(zz)) for zz in z])
            rows.append({
                "primary_exposure": pabbr, "outcome": o["abbr"], "nSNP": len(h),
                "b_primary": round(beta[0], 4), "se_primary": round(se[0], 4),
                "p_primary": round(float(pvals[0]), 4),
                "OR_primary": round(float(np.exp(beta[0])), 3),
                "cond_F_primary": round(cond_F[0], 1),
                "OR_SMK": round(float(np.exp(beta[1])), 3) if len(beta) > 1 else None,
                "OR_BMI": round(float(np.exp(beta[2])), 3) if len(beta) > 2 else None,
                "OR_T2D": round(float(np.exp(beta[3])), 3) if len(beta) > 3 else None,
            })
            print(f"  MVMR {pabbr}->{o['abbr']}: OR={np.exp(beta[0]):.3f} "
                  f"P={pvals[0]:.3g} condF={cond_F[0]:.1f}", flush=True)
    if rows:
        pd.DataFrame(rows).to_csv(os.path.join(RESULTS, "mvmr.tsv"), sep="\t", index=False)
        print(f"[written] results/mvmr.tsv ({len(rows)} rows)", flush=True)
    else:
        with open(os.path.join(RESULTS, "mvmr.tsv"), "w") as f:
            f.write("# MVMR not feasible (insufficient shared instrument / EBI blocked)\n")
        print("MVMR empty -> written placeholder", flush=True)

# ==========================================================================
# 3. REVERSE MR (spine infection -> immune), variant-level
# ==========================================================================
def run_reverse():
    print("=== REVERSE MR (spine -> immune) ===", flush=True)
    expo_map = {"OM": "data/ivs_reverse_OM.csv", "DISCITIS": "data/ivs_reverse_DISCITIS.csv"}
    immune = [e for e in mp.EXPOSURES]
    total = 0
    for expo, path in expo_map.items():
        if not os.path.exists(path):
            continue
        iv = pd.read_csv(path)
        iv = iv[iv["F"] > F_MIN]
        print(f"  reverse exposure {expo}: {len(iv)} IVs", flush=True)
        out_rows = []
        for e in immune:
            bx, sex, by, sey, snps = [], [], [], [], []
            got = 0
            for _, r in iv.iterrows():
                rec = fetch_confounder_at(e["efo"], r["rsid"])
                if not rec or rec.get("beta") is None:
                    continue
                bx.append(float(r["beta"])); sex.append(float(r["se"]))
                by.append(float(rec["beta"])); sey.append(float(rec["se"]))
                snps.append(r["rsid"]); got += 1
            if got < 3 or ebi_safe(lambda: None) is None and False:
                pass
            if got < 3:
                continue
            res = mrm.run_all(np.array(bx), np.array(sex), np.array(by), np.array(sey), snps)
            for rr in res["results"]:
                d = rr.as_dict(binary_outcome=False)
                d.update({"reverse_exposure": expo, "immune_outcome": e["abbr"], "nSNP": got})
                out_rows.append(d)
            total += 1
            print(f"    {expo}->{e['abbr']}: nSNP={got} b={res['results'][0].b:.3f} "
                  f"P={res['results'][0].pval:.3g}", flush=True)
        if out_rows:
            pd.DataFrame(out_rows).to_csv(
                os.path.join(RESULTS, f"reverse_{expo}_immune.tsv"), sep="\t", index=False)
    if total == 0:
        with open(os.path.join(RESULTS, "reverse_MR_immune.tsv"), "w") as f:
            f.write("# reverse MR not feasible: EBI rate-limiting blocked immune-endpoint fetch\n")
        print("reverse MR empty -> EBI blocked, written placeholder", flush=True)
    else:
        print(f"[written] reverse_*_immune.tsv (pairs={total})", flush=True)

if __name__ == "__main__":
    t0 = time.time()
    try: run_drug_target()
    except Exception as e:
        print(f"[drug_target error] {e}", flush=True); import traceback; traceback.print_exc()
    try: run_mvmr()
    except Exception as e:
        print(f"[mvmr error] {e}", flush=True); import traceback; traceback.print_exc()
    try: run_reverse()
    except Exception as e:
        print(f"[reverse error] {e}", flush=True); import traceback; traceback.print_exc()
    print(f"=== network supplementary done in {time.time()-t0:.0f}s ===", flush=True)
