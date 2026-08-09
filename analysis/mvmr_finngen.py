"""
mvmr_finngen.py — 多变量 MR（FinnGen 混杂因子版，绕开 EBI 限流）

设计：暴露主导式 MVMR
  - 工具变量：主暴露（免疫性状）的 IV SNP（来自主管线 harmonised 缓存）
  - 结局：OM / DISCITIS（FinnGen R11，本地缓存的 harmonised beta）
  - 混杂因子：BMI_IRN / T2D / SMOKING（FinnGen R11 全量汇总统计，BGZF 随机访问）
  - 估计：多变量 IVW（固定效应 + 过散布校正 min(1,sqrt(phi))）
  - 报告：校正后 OR、条件 F（近似 mean z^2）、单变量 IVW 对照

等位基因定向：混杂因子与结局同为 FinnGen 文件，同一 chr:pos 的 ref/alt 一致，
故沿用 harmonised 缓存中的 harmonise_action 决定符号（flip* -> 取负）。
"""
from __future__ import annotations
import os, sys, time
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gwas_io as gio

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
HARMDIR = os.path.join(RESULTS, "harmonised")
CACHE = os.path.join(ROOT, "data", "cache", "mvmr_conf")
os.makedirs(CACHE, exist_ok=True)

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

EXPOSURES = ["WBC", "NEUT", "MONO", "LYMPH", "CRP", "IL6R", "EOS"]
OUTCOMES = {"OM": "M13_OSTEOMYELITIS", "DISCITIS": "M13_DISCITIS"}
CONFOUNDERS = [("BMI", "BMI_IRN"), ("T2D", "T2D"), ("SMK", "SMOKING")]

_conf_fs = {}

def conf_query(pheno: str, chrom: str, pos: int):
    """带磁盘缓存的 FinnGen 混杂因子查询（每次新建对象，避免块缓存损坏）。"""
    key = f"{pheno}_{chrom}_{pos}"
    fp = os.path.join(CACHE, key + ".json")
    import json
    if os.path.exists(fp):
        try:
            cached = json.load(open(fp))
            if cached:  # 仅信任非空缓存（空缓存可能是中断运行的毒数据）
                return cached
        except Exception:
            pass
    # 每次新建对象：.tbi 索引磁盘缓存，开销小；避免复用对象的块缓存在大文件多点位随机访问时损坏
    fs = gio.FinnGenSumstats(pheno)
    try:
        # harmonised CSV 中 chr 可能是 float (如 1.0)；必须转为 int 再 str()，否则 tabix 索引查找失败
        hits = fs.query(str(int(float(chrom))), int(float(pos)), window=0)
    except Exception as e:
        hits = []
    with open(fp, "w") as f:
        json.dump(hits, f)
    return hits


def match_row(hits, ea, oa):
    """在 FinnGen 查询结果中找与 {ea,oa}（或其互补）匹配的行。"""
    ea, oa = ea.upper(), oa.upper()
    ea_c = "".join(COMPLEMENT.get(b, "N") for b in ea)
    oa_c = "".join(COMPLEMENT.get(b, "N") for b in oa)
    for h in hits:
        ref, alt = h["ref"].upper(), h["alt"].upper()
        if {ref, alt} == {ea, oa} or {ref, alt} == {ea_c, oa_c}:
            return h
    return None


def mvmr_ivw(Yb, Yse, Xb, Xse):
    """多变量 IVW：W = 1/se_Y^2；返回 beta, se(过散布校正), p, phi。"""
    from scipy import stats
    W = 1.0 / Yse ** 2
    XtW = Xb.T * W
    XtWX = XtW @ Xb
    try:
        cov = np.linalg.inv(XtWX)
    except np.linalg.LinAlgError:
        return None
    beta = cov @ (XtW @ Yb)
    resid = Yb - Xb @ beta
    phi = float(np.sum(W * resid ** 2) / max(1, len(Yb) - Xb.shape[1]))
    rse = max(1.0, np.sqrt(phi))
    se = np.sqrt(np.diag(cov)) * rse
    z = beta / se
    p = 2 * stats.norm.sf(np.abs(z))
    return beta, se, p, phi


def main():
    t0 = time.time()
    rows = []
    for oabbr in OUTCOMES:
        for eabbr in EXPOSURES:
            fp = os.path.join(HARMDIR, f"forward_{eabbr}__{oabbr}_harmonised.csv")
            if not os.path.exists(fp):
                print(f"  {eabbr}->{oabbr}: no cache, skip", flush=True)
                continue
            h = pd.read_csv(fp)
            h = h[h["status"] == "kept"].copy()
            if len(h) < 3:
                print(f"  {eabbr}->{oabbr}: only {len(h)} SNP, skip", flush=True)
                continue
            # 逐 SNP 取混杂因子 beta/se
            conf_b = {c[0]: [] for c in CONFOUNDERS}
            conf_se = {c[0]: [] for c in CONFOUNDERS}
            keep_idx = []
            miss = {c[0]: 0 for c in CONFOUNDERS}
            for i, r in h.iterrows():
                sign = -1.0 if str(r["harmonise_action"]).startswith("flip") else 1.0
                got_all = True
                tmp_b, tmp_se = {}, {}
                for cabbr, pheno in CONFOUNDERS:
                    hits = conf_query(pheno, r["chr"], r["pos"])
                    m = match_row(hits, r["effect_allele"], r["other_allele"]) if hits else None
                    if m is None or not np.isfinite(float(m.get("beta", np.nan))) or float(m.get("sebeta", 0)) <= 0:
                        miss[cabbr] += 1
                        got_all = False
                        break
                    tmp_b[cabbr] = sign * float(m["beta"])
                    tmp_se[cabbr] = float(m["sebeta"])
                if got_all:
                    keep_idx.append(i)
                    for cabbr, _ in CONFOUNDERS:
                        conf_b[cabbr].append(tmp_b[cabbr])
                        conf_se[cabbr].append(tmp_se[cabbr])
            if len(keep_idx) < 4:
                print(f"  {eabbr}->{oabbr}: only {len(keep_idx)} complete SNP "
                      f"(miss={miss}), skip", flush=True)
                continue
            hk = h.loc[keep_idx]
            Yb = hk["beta_outcome"].to_numpy(float)
            Yse = hk["se_outcome"].to_numpy(float)
            Xb = np.column_stack([hk["beta_exposure"].to_numpy(float)] +
                                 [np.array(conf_b[c[0]]) for c in CONFOUNDERS])
            Xse = np.column_stack([hk["se_exposure"].to_numpy(float)] +
                                  [np.array(conf_se[c[0]]) for c in CONFOUNDERS])
            # 单变量 IVW 对照
            w = 1.0 / Yse ** 2
            b_uni = float(np.sum(w * hk["beta_exposure"] * Yb) / np.sum(w * hk["beta_exposure"] ** 2))
            se_uni = float(np.sqrt(1.0 / np.sum(w * hk["beta_exposure"] ** 2)))
            from scipy import stats
            p_uni = float(2 * stats.norm.sf(abs(b_uni / se_uni)))
            # MVMR
            res = mvmr_ivw(Yb, Yse, Xb, Xse)
            if res is None:
                print(f"  {eabbr}->{oabbr}: singular design, skip", flush=True)
                continue
            beta, se, p, phi = res
            condF = [float(np.mean((Xb[:, j] / Xse[:, j]) ** 2)) for j in range(Xb.shape[1])]
            row = {
                "exposure": eabbr, "outcome": oabbr, "nSNP": len(hk),
                "OR_uni": round(float(np.exp(b_uni)), 3), "p_uni": round(p_uni, 4),
                "OR_adj": round(float(np.exp(beta[0])), 3),
                "lo95_adj": round(float(np.exp(beta[0] - 1.96 * se[0])), 3),
                "hi95_adj": round(float(np.exp(beta[0] + 1.96 * se[0])), 3),
                "p_adj": round(float(p[0]), 4),
                "condF_exposure": round(condF[0], 1),
                "condF_BMI": round(condF[1], 1),
                "condF_T2D": round(condF[2], 1),
                "condF_SMK": round(condF[3], 1),
                "phi": round(phi, 3),
            }
            for j, (cabbr, _) in enumerate(CONFOUNDERS, start=1):
                row[f"OR_{cabbr}"] = round(float(np.exp(beta[j])), 3)
                row[f"p_{cabbr}"] = round(float(p[j]), 4)
            rows.append(row)
            print(f"  MVMR {eabbr}->{oabbr}: nSNP={len(hk)} "
                  f"OR_adj={row['OR_adj']} P={row['p_adj']:.3g} "
                  f"condF={row['condF_exposure']}", flush=True)
    if rows:
        out = os.path.join(RESULTS, "mvmr.tsv")
        pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
        print(f"[written] {out} ({len(rows)} rows)", flush=True)
    else:
        print("MVMR: no rows", flush=True)
    print(f"=== done in {time.time()-t0:.0f}s ===", flush=True)


if __name__ == "__main__":
    main()
