"""
supplementary.py — 纯本地补充分析（无需任何网络取数）
  1) POWER   : 各暴露×结局对的 80% 统计效能最小可检 OR (MDE)
  2) STEIGER : 方向性检验（暴露 vs 结局方差解释量，Hemani 2017）
  3) META    : 跨脊柱表型的逆方差随机效应 meta（含/不含最小结局 SPONDINF）

输入：results/MR_results_all.csv （IVW 随机效应 b/se/pval）
      results/harmonised/forward_*_harmonised.csv （per-SNP 用于 Steiger）
输出：results/power.tsv, results/steiger.tsv, results/spine_replication.tsv
"""
from __future__ import annotations
import os, glob, math
import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
HARM = os.path.join(RESULTS, "harmonised")

# --------------------------------------------------------------------------
# 1. POWER  (post-hoc MDE at 80% power, two-sided alpha=0.05)
#    MDE_logOR = (z_{alpha/2} + z_{1-beta}) * se_IVW ;  z=1.95996+0.84162=2.8016
# --------------------------------------------------------------------------
def power_table():
    df = pd.read_csv(os.path.join(RESULTS, "MR_results_all.csv"))
    ivw = df[df["method"] == "IVW (random effects)"].copy()
    z = 1.959964 + 0.841621  # 80% power two-sided
    rows = []
    for _, r in ivw.iterrows():
        se = float(r["se"])
        mde_log = z * se
        rows.append({
            "exposure": r["exposure_abbr"], "outcome": r["outcome_abbr"],
            "nsnp": int(r["nsnp"]), "se_IVW": round(se, 4),
            "MDE_logOR": round(mde_log, 4),
            "MDE_OR": round(math.exp(mde_log), 3),
            "obs_OR": round(float(r["OR"]), 3),
            "obs_p": float(r["pval"]),
        })
    out = pd.DataFrame(rows).sort_values(["outcome", "exposure"])
    out.to_csv(os.path.join(RESULTS, "power.tsv"), sep="\t", index=False)
    return out

# --------------------------------------------------------------------------
# 2. STEIGER directionality (per-SNP variance explained)
#    R2X_i = (bX/seX)^2 ; R2Y_i = (bY/seY)^2
#    stat = (ΣR2X − ΣR2Y) / sqrt(2*(ΣR2X + ΣR2Y)) ; p = 2*pnorm(|stat|)
# --------------------------------------------------------------------------
def steiger_table():
    rows = []
    for fp in sorted(glob.glob(os.path.join(HARM, "forward_*_harmonised.csv"))):
        fn = os.path.basename(fp)
        parts = fn.replace("forward_", "").replace("_harmonised.csv", "")
        abbr, outc = parts.split("__", 1)
        h = pd.read_csv(fp)
        h = h[h["status"] == "kept"].copy()
        if h.empty or len(h) < 2:
            continue
        bx, sex = h["beta_exposure"].values.astype(float), h["se_exposure"].values.astype(float)
        by, sey = h["beta_outcome"].values.astype(float), h["se_outcome"].values.astype(float)
        if np.any(~np.isfinite(sex)) or np.any(sex <= 0):
            continue
        if np.any(~np.isfinite(sey)) or np.any(sey <= 0):
            continue
        r2x = (bx / sex) ** 2
        r2y = (by / sey) ** 2
        Sx, Sy = r2x.sum(), r2y.sum()
        denom = math.sqrt(2 * (Sx + Sy))
        stat = (Sx - Sy) / denom if denom > 0 else float("nan")
        p = 2 * stats.norm.sf(abs(stat)) if np.isfinite(stat) else float("nan")
        supported = "exposure->outcome" if (stat > 0 and p < 0.05) else \
                    ("outcome->exposure" if (stat < 0 and p < 0.05) else "inconclusive")
        rows.append({
            "pair": f"{abbr}->{outc}", "nsnp": len(h),
            "sum_R2_exposure": round(Sx, 2), "sum_R2_outcome": round(Sy, 2),
            "steiger_stat": round(stat, 3), "pval": round(p, 4),
            "directionality": supported,
        })
    out = pd.DataFrame(rows).sort_values("pair")
    out.to_csv(os.path.join(RESULTS, "steiger.tsv"), sep="\t", index=False)
    return out

# --------------------------------------------------------------------------
# 3. META : random-effects (DerSimonian-Laird) across spinal phenotypes
# --------------------------------------------------------------------------
def _dl_meta(b, se):
    b, se = np.asarray(b, float), np.asarray(se, float)
    w = 1.0 / se ** 2
    pooled_f = np.sum(w * b) / np.sum(w)
    Q = float(np.sum(w * (b - pooled_f) ** 2))
    df = len(b) - 1
    # tau^2 DL
    num = max(0.0, Q - df)
    den = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = num / den if den > 0 else 0.0
    wr = 1.0 / (se ** 2 + tau2)
    pooled = float(np.sum(wr * b) / np.sum(wr))
    se_r = math.sqrt(1.0 / np.sum(wr))
    z = pooled / se_r
    p = 2 * stats.norm.sf(abs(z))
    I2 = max(0.0, (Q - df) / Q * 100) if Q > 0 else 0.0
    return pooled, se_r, p, I2, tau2

def meta_table():
    df = pd.read_csv(os.path.join(RESULTS, "MR_results_all.csv"))
    ivw = df[df["method"] == "IVW (random effects)"].copy()
    out_rows, nospond_rows = [], []
    for exp, g in ivw.groupby("exposure_abbr"):
        # 全部 5 个结局
        b = g["b"].values.astype(float); se = g["se"].values.astype(float)
        pooled, ser, p, I2, tau2 = _dl_meta(b, se)
        out_rows.append({
            "exposure": exp, "n_outcomes": len(g),
            "pooled_OR": round(math.exp(pooled), 3),
            "OR_lo95": round(math.exp(pooled - 1.96 * ser), 3),
            "OR_hi95": round(math.exp(pooled + 1.96 * ser), 3),
            "pval": round(p, 4), "I2_pct": round(I2, 1), "tau2": round(tau2, 4),
            "note": "meta across all 5 spinal/OM phenotypes",
        })
        # 剔除最小且不可靠的 SPONDINF
        g2 = g[g["outcome_abbr"] != "SPONDINF"]
        if len(g2) >= 2:
            b2 = g2["b"].values.astype(float); se2 = g2["se"].values.astype(float)
            p2, ser2, p2p, I2b, _ = _dl_meta(b2, se2)
            nospond_rows.append({
                "exposure": exp, "n_outcomes": len(g2),
                "pooled_OR": round(math.exp(p2), 3),
                "OR_lo95": round(math.exp(p2 - 1.96 * ser2), 3),
                "OR_hi95": round(math.exp(p2 + 1.96 * ser2), 3),
                "pval": round(p2p, 4), "I2_pct": round(I2b, 1),
                "note": "meta EXCLUDING SPONDINF",
            })
    out = pd.DataFrame(out_rows); nos = pd.DataFrame(nospond_rows)
    comb = pd.concat([out, nos], ignore_index=True)
    comb.to_csv(os.path.join(RESULTS, "spine_replication.tsv"), sep="\t", index=False)
    return out, nos

if __name__ == "__main__":
    print("=== POWER (80% power MDE) ===")
    pw = power_table()
    print(pw.to_string(index=False))
    print(f"\n[written] results/power.tsv  ({len(pw)} pairs)")
    print("\n=== STEIGER directionality ===")
    st = steiger_table()
    print(st.to_string(index=False))
    print(f"\n[written] results/steiger.tsv  ({len(st)} pairs)")
    print("\n=== META across spinal phenotypes (IVW random effects) ===")
    m_all, m_ns = meta_table()
    print(m_all.to_string(index=False))
    print("\n--- meta EXCLUDING SPONDINF ---")
    print(m_ns.to_string(index=False))
    print(f"\n[written] results/spine_replication.tsv")
