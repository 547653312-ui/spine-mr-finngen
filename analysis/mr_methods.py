"""
mr_methods.py — 两样本 MR 估计量与敏感性分析（纯 Python 实现）

说明：PyPI 上不存在 `MendelianRandomization` 包（该名称为 R 包，CRAN 专有），
`genal` 等替代亦不可安装（见 results/API_NOTES.md）。因此本模块按标准公式
自行实现，方法学对齐 R 包 MendelianRandomization / TwoSampleMR：

  - IVW（固定效应 + 乘法随机效应）
  - MR-Egger（加权回归，含截距水平多效性检验）
  - 加权中位数（Bowden 2016，bootstrap SE）
  - 加权众数（Hartwig 2017，bootstrap SE）
  - Cochran Q 异质性、Egger 截距、留一法

约定：bx/bxse 为暴露效应，by/byse 为结局效应，均以同一效应等位基因为参照。
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import List, Optional

import numpy as np
from scipy import stats


@dataclass
class MRResult:
    method: str
    nsnp: int
    b: float
    se: float
    pval: float
    lo95: float
    hi95: float
    or_: Optional[float] = None
    or_lo95: Optional[float] = None
    or_hi95: Optional[float] = None
    note: str = ""

    def as_dict(self, binary_outcome: bool = True) -> dict:
        d = asdict(self)
        if binary_outcome:
            d["or_"] = float(np.exp(self.b))
            d["or_lo95"] = float(np.exp(self.lo95))
            d["or_hi95"] = float(np.exp(self.hi95))
        d["OR"] = d.pop("or_")
        d["OR_lo95"] = d.pop("or_lo95")
        d["OR_hi95"] = d.pop("or_hi95")
        return d


def _finish(method: str, n: int, b: float, se: float, note: str = "") -> MRResult:
    if not np.isfinite(b) or not np.isfinite(se) or se <= 0:
        return MRResult(method, n, float("nan"), float("nan"), float("nan"),
                        float("nan"), float("nan"), note=note or "degenerate")
    z = b / se
    p = 2 * stats.norm.sf(abs(z))
    return MRResult(method, n, float(b), float(se), float(p),
                    float(b - 1.96 * se), float(b + 1.96 * se), note=note)


# --------------------------------------------------------------------------
# IVW
# --------------------------------------------------------------------------
def ivw(bx, bxse, by, byse, random: bool = True) -> MRResult:
    bx, by, byse = np.asarray(bx, float), np.asarray(by, float), np.asarray(byse, float)
    n = len(bx)
    if n < 1:
        return _finish("IVW", 0, np.nan, np.nan, "no SNP")
    if n == 1:
        b = by[0] / bx[0]
        se = abs(byse[0] / bx[0])
        return _finish("Wald ratio", 1, b, se, "single IV")
    w = 1.0 / byse ** 2
    denom = np.sum(w * bx ** 2)
    b = np.sum(w * bx * by) / denom
    se_fixed = np.sqrt(1.0 / denom)
    q = np.sum(w * (by - b * bx) ** 2)
    df = n - 1
    if random:
        phi = max(1.0, np.sqrt(q / df))
        return _finish("IVW (random effects)", n, b, se_fixed * phi,
                       f"Q={q:.3f}, df={df}")
    return _finish("IVW (fixed effects)", n, b, se_fixed, f"Q={q:.3f}, df={df}")


def cochran_q(bx, bxse, by, byse, b_ivw: float) -> dict:
    bx, by, byse = np.asarray(bx, float), np.asarray(by, float), np.asarray(byse, float)
    n = len(bx)
    if n < 2:
        return {"Q": np.nan, "Q_df": 0, "Q_pval": np.nan, "I2": np.nan}
    w = 1.0 / byse ** 2
    q = float(np.sum(w * (by - b_ivw * bx) ** 2))
    df = n - 1
    p = float(stats.chi2.sf(q, df))
    i2 = max(0.0, (q - df) / q) * 100 if q > 0 else 0.0
    return {"Q": q, "Q_df": df, "Q_pval": p, "I2": i2}


# --------------------------------------------------------------------------
# MR-Egger
# --------------------------------------------------------------------------
def mr_egger(bx, bxse, by, byse) -> tuple[MRResult, dict]:
    bx, by, byse = np.asarray(bx, float), np.asarray(by, float), np.asarray(byse, float)
    n = len(bx)
    if n < 3:
        empty = {"intercept": np.nan, "intercept_se": np.nan,
                 "intercept_pval": np.nan, "nsnp": n}
        return _finish("MR-Egger", n, np.nan, np.nan, "need >=3 IV"), empty
    # 定向：使所有暴露效应为正
    sign = np.where(bx >= 0, 1.0, -1.0)
    bx_o, by_o = bx * sign, by * sign
    w = 1.0 / byse ** 2
    X = np.column_stack([np.ones(n), bx_o])
    W = np.diag(w)
    xtwx = X.T @ W @ X
    try:
        inv = np.linalg.inv(xtwx)
    except np.linalg.LinAlgError:
        empty = {"intercept": np.nan, "intercept_se": np.nan,
                 "intercept_pval": np.nan, "nsnp": n}
        return _finish("MR-Egger", n, np.nan, np.nan, "singular"), empty
    coef = inv @ (X.T @ W @ by_o)
    resid = by_o - X @ coef
    dof = n - 2
    # 残差离散度（与 R 包默认一致，下界为 1）
    sigma2 = float(resid.T @ W @ resid) / dof
    disp = max(1.0, sigma2)
    cov = inv * disp
    inter, slope = float(coef[0]), float(coef[1])
    se_inter, se_slope = float(np.sqrt(cov[0, 0])), float(np.sqrt(cov[1, 1]))
    p_inter = 2 * stats.t.sf(abs(inter / se_inter), dof) if se_inter > 0 else np.nan
    pleio = {"intercept": inter, "intercept_se": se_inter,
             "intercept_pval": float(p_inter), "nsnp": n}
    return _finish("MR-Egger", n, slope, se_slope, "dispersion-adjusted"), pleio


# --------------------------------------------------------------------------
# 加权中位数
# --------------------------------------------------------------------------
def _weighted_median(betas: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(betas)
    b, w = betas[order], weights[order]
    s = np.cumsum(w) - 0.5 * w
    s /= np.sum(w)
    below = np.where(s < 0.5)[0]
    if len(below) == 0:
        return float(b[0])
    k = below[-1]
    if k + 1 >= len(b):
        return float(b[-1])
    return float(b[k] + (b[k + 1] - b[k]) * (0.5 - s[k]) / (s[k + 1] - s[k]))


def weighted_median(bx, bxse, by, byse, nboot: int = 1000,
                    seed: int = 20260808) -> MRResult:
    bx, bxse = np.asarray(bx, float), np.asarray(bxse, float)
    by, byse = np.asarray(by, float), np.asarray(byse, float)
    n = len(bx)
    if n < 3:
        return _finish("Weighted median", n, np.nan, np.nan, "need >=3 IV")
    ratio = by / bx
    var_ratio = byse ** 2 / bx ** 2
    w = 1.0 / var_ratio
    b = _weighted_median(ratio, w)
    rng = np.random.default_rng(seed)
    boot = np.empty(nboot)
    for i in range(nboot):
        bxs = rng.normal(bx, bxse)
        bys = rng.normal(by, byse)
        with np.errstate(divide="ignore", invalid="ignore"):
            r = bys / bxs
            v = byse ** 2 / bxs ** 2
        ok = np.isfinite(r) & np.isfinite(v) & (v > 0)
        if ok.sum() < 3:
            boot[i] = np.nan
            continue
        boot[i] = _weighted_median(r[ok], 1.0 / v[ok])
    se = float(np.nanstd(boot, ddof=1))
    return _finish("Weighted median", n, b, se, f"bootstrap n={nboot}")


# --------------------------------------------------------------------------
# 加权众数
# --------------------------------------------------------------------------
def _mode_estimate(ratio: np.ndarray, se_ratio: np.ndarray, phi: float) -> float:
    w = 1.0 / se_ratio ** 2
    w = w / np.sum(w)
    s = np.std(ratio, ddof=1) if len(ratio) > 1 else 1.0
    if not np.isfinite(s) or s <= 0:
        s = 1.0
    # Silverman 经验带宽
    iqr = np.subtract(*np.percentile(ratio, [75, 25]))
    a = min(s, iqr / 1.34) if iqr > 0 else s
    h = phi * 0.9 * a * len(ratio) ** (-0.2)
    if not np.isfinite(h) or h <= 0:
        h = 1e-6
    grid = np.linspace(ratio.min() - 3 * h, ratio.max() + 3 * h, 1024)
    dens = np.sum(w[None, :] * stats.norm.pdf((grid[:, None] - ratio[None, :]) / h), axis=1)
    return float(grid[int(np.argmax(dens))])


def weighted_mode(bx, bxse, by, byse, phi: float = 1.0, nboot: int = 1000,
                  seed: int = 20260808) -> MRResult:
    bx, bxse = np.asarray(bx, float), np.asarray(bxse, float)
    by, byse = np.asarray(by, float), np.asarray(byse, float)
    n = len(bx)
    if n < 3:
        return _finish("Weighted mode", n, np.nan, np.nan, "need >=3 IV")
    ratio = by / bx
    se_ratio = np.sqrt(byse ** 2 / bx ** 2)
    b = _mode_estimate(ratio, se_ratio, phi)
    rng = np.random.default_rng(seed)
    boot = np.empty(nboot)
    for i in range(nboot):
        bys = rng.normal(by, byse)
        bxs = rng.normal(bx, bxse)
        with np.errstate(divide="ignore", invalid="ignore"):
            r = bys / bxs
            sr = np.sqrt(byse ** 2 / bxs ** 2)
        ok = np.isfinite(r) & np.isfinite(sr) & (sr > 0)
        boot[i] = _mode_estimate(r[ok], sr[ok], phi) if ok.sum() >= 3 else np.nan
    se = float(np.nanstd(boot, ddof=1))
    return _finish("Weighted mode", n, b, se, f"phi={phi}, bootstrap n={nboot}")


# --------------------------------------------------------------------------
# 留一法 / 单 SNP
# --------------------------------------------------------------------------
def leave_one_out(bx, bxse, by, byse, snps: List[str]) -> List[dict]:
    bx, bxse = np.asarray(bx, float), np.asarray(bxse, float)
    by, byse = np.asarray(by, float), np.asarray(byse, float)
    out = []
    n = len(bx)
    if n < 3:
        return out
    for i in range(n):
        m = np.ones(n, bool)
        m[i] = False
        r = ivw(bx[m], bxse[m], by[m], byse[m])
        out.append({"excluded_SNP": snps[i], "nsnp": r.nsnp, "b": r.b,
                    "se": r.se, "pval": r.pval, "lo95": r.lo95, "hi95": r.hi95,
                    "OR": float(np.exp(r.b)) if np.isfinite(r.b) else np.nan})
    r_all = ivw(bx, bxse, by, byse)
    out.append({"excluded_SNP": "ALL (none excluded)", "nsnp": r_all.nsnp,
                "b": r_all.b, "se": r_all.se, "pval": r_all.pval,
                "lo95": r_all.lo95, "hi95": r_all.hi95,
                "OR": float(np.exp(r_all.b)) if np.isfinite(r_all.b) else np.nan})
    return out


def single_snp(bx, bxse, by, byse, snps: List[str]) -> List[dict]:
    bx, bxse = np.asarray(bx, float), np.asarray(bxse, float)
    by, byse = np.asarray(by, float), np.asarray(byse, float)
    out = []
    for i in range(len(bx)):
        b = by[i] / bx[i]
        se = abs(byse[i] / bx[i])
        z = b / se if se > 0 else np.nan
        out.append({"SNP": snps[i], "b": float(b), "se": float(se),
                    "pval": float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else np.nan,
                    "OR": float(np.exp(b)),
                    "OR_lo95": float(np.exp(b - 1.96 * se)),
                    "OR_hi95": float(np.exp(b + 1.96 * se))})
    return out


def f_statistic(bx, bxse) -> np.ndarray:
    bx, bxse = np.asarray(bx, float), np.asarray(bxse, float)
    return (bx / bxse) ** 2


def run_all(bx, bxse, by, byse, snps: List[str]) -> dict:
    """跑全套方法 + 敏感性分析。"""
    res: List[MRResult] = []
    r_ivw_re = ivw(bx, bxse, by, byse, random=True)
    r_ivw_fe = ivw(bx, bxse, by, byse, random=False)
    res.append(r_ivw_re)
    if len(bx) > 1:
        res.append(r_ivw_fe)
    r_egg, pleio = mr_egger(bx, bxse, by, byse)
    if len(bx) >= 3:
        res.append(r_egg)
        res.append(weighted_median(bx, bxse, by, byse))
        res.append(weighted_mode(bx, bxse, by, byse))
    het = cochran_q(bx, bxse, by, byse, r_ivw_fe.b if np.isfinite(r_ivw_fe.b) else r_ivw_re.b)
    return {
        "results": res,
        "heterogeneity": het,
        "pleiotropy": pleio,
        "loo": leave_one_out(bx, bxse, by, byse, snps),
        "single": single_snp(bx, bxse, by, byse, snps),
        "F": f_statistic(bx, bxse),
    }
