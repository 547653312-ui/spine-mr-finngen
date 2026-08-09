"""
make_main_figures.py — 生成 4 张投稿主图 (300 DPI PNG)
  Figure 1: 研究设计 + 暴露-结局矩阵 + IV 计数/F 分布
  Figure 2: 主要 IVW 估计森林图 (45 对)
  Figure 3: SPONDINF 信号脆弱性 (方法比较 + 方向不一致 + 留一法)
  Figure 4: 统计功效 + 跨表型 meta
浅色主题，符合 light IDE。
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrow
import matplotlib.gridspec as gridspec

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.edgecolor": "#444", "axes.labelcolor": "#222",
    "xtick.color": "#333", "ytick.color": "#333",
    "axes.titlecolor": "#111", "figure.facecolor": "white",
    "axes.facecolor": "white", "savefig.facecolor": "white",
})

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
FIG = os.path.join(RES, "figures")
os.makedirs(FIG, exist_ok=True)

EXPO_ORDER = ["WBC", "NEUT", "LYMPH", "MONO", "EOS", "CRP", "IL6R", "CCL2", "CD40"]
EXPO_LABEL = {"WBC": "Leukocyte", "NEUT": "Neutrophil", "LYMPH": "Lymphocyte",
              "MONO": "Monocyte", "EOS": "Eosinophil", "CRP": "CRP",
              "IL6R": "IL-6 receptor", "CCL2": "CCL2/MCP-1", "CD40": "CD40"}
OUT_ORDER = ["OM", "DISCITIS", "DISCINF", "VOM", "SPONDINF"]
OUT_LABEL = {"OM": "OM\n(n=2,125)", "DISCITIS": "Discitis\n(n=495)",
             "DISCINF": "Disc infection\n(n=375)", "VOM": "VOM\n(n=104)",
             "SPONDINF": "Spondylopathy\n(n=68)"}
IV_COUNT = {"WBC": 28, "NEUT": 29, "LYMPH": 18, "MONO": 39, "EOS": 13,
            "CRP": 13, "IL6R": 4, "CCL2": 7, "CD40": 3}
# 颜色
C_SPINE = "#1B5E9F"   # 脊柱特异
C_OM = "#5C6B7F"      # OM 结局（中性深灰，区别于显著红色）
C_NULL = "#888"
C_SIG = "#C5283D"
C_KEEP = "#2E7D32"


def load_ivw():
    df = pd.read_csv(os.path.join(RES, "MR_results_all.csv"))
    df = df[(df["direction"] == "forward") &
            (df["method"].str.contains("IVW", na=False)) &
            (df["method"].str.contains("random", na=False))]
    df = df.drop_duplicates(subset=["exposure_abbr", "outcome_abbr"])
    return df


# ============================ Figure 1 ============================
def figure1():
    df = load_ivw()
    fig = plt.figure(figsize=(11, 8.5))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1.05, 1], hspace=0.42, wspace=0.3)

    # (a) 研究设计示意
    axa = fig.add_subplot(gs[0, :])
    axa.set_xlim(0, 10); axa.set_ylim(0, 4); axa.axis("off")
    axa.set_title("(a) Two-sample Mendelian randomization design",
                  loc="left", fontsize=10, fontweight="bold")
    # 暴露框
    for (x, t, c) in [(0.6, "11 immune / cytokine traits\n(EBI GWAS Catalog)", C_SPINE)]:
        axa.add_patch(FancyBboxPatch((x, 2.4), 2.7, 1.2, boxstyle="round,pad=0.08",
                     fc=c, ec="none", alpha=0.92))
        axa.text(x + 1.35, 3.0, t, ha="center", va="center", color="white",
                 fontsize=8.5, fontweight="bold")
    # IV 框
    axa.add_patch(FancyBboxPatch((3.9, 2.4), 2.2, 1.2, boxstyle="round,pad=0.08",
                 fc="#5C6B7F", ec="none", alpha=0.92))
    axa.text(5.0, 3.0, "Genome-wide IVs\n(P<5×10⁻⁸, F>10)\ndistance-pruned",
             ha="center", va="center", color="white", fontsize=8.5, fontweight="bold")
    # 结局框
    axa.add_patch(FancyBboxPatch((6.7, 2.4), 2.7, 1.2, boxstyle="round,pad=0.08",
                 fc=C_OM, ec="none", alpha=0.92))
    axa.text(8.05, 3.0, "5 osteomyelitis phenotypes\n(FinnGen R11)\n4 spine-specific",
             ha="center", va="center", color="white", fontsize=8.5, fontweight="bold")
    # 箭头
    axa.annotate("", xy=(3.9, 3.0), xytext=(3.3, 3.0),
                 arrowprops=dict(arrowstyle="->", lw=1.8, color="#333"))
    axa.annotate("", xy=(6.7, 3.0), xytext=(6.1, 3.0),
                 arrowprops=dict(arrowstyle="->", lw=1.8, color="#333"))
    # 方法行
    axa.text(5.0, 1.5, "IVW · MR-Egger · weighted median · weighted mode",
             ha="center", fontsize=8.5, color="#333", style="italic")
    axa.text(5.0, 0.95, "Sensitivity: Cochran Q · Egger intercept · leave-one-out · "
             "Steiger directionality · power · cross-phenotype meta",
             ha="center", fontsize=7.8, color="#555")
    axa.text(5.0, 0.4, "Extensions: drug-target (cis-pQTL) MR · multivariable MR · reverse MR",
             ha="center", fontsize=7.8, color="#555")

    # (b) 暴露-结局矩阵
    axb = fig.add_subplot(gs[1, 0])
    axb.set_title("(b) Analysable exposure–outcome pairs",
                  loc="left", fontsize=10, fontweight="bold")
    mat = np.full((len(EXPO_ORDER), len(OUT_ORDER)), np.nan)
    pmat = np.full_like(mat, np.nan)
    for _, r in df.iterrows():
        if r["exposure_abbr"] in EXPO_ORDER and r["outcome_abbr"] in OUT_ORDER:
            i = EXPO_ORDER.index(r["exposure_abbr"])
            j = OUT_ORDER.index(r["outcome_abbr"])
            mat[i, j] = 1
            pmat[i, j] = r["pval"]
    for i in range(len(EXPO_ORDER)):
        for j in range(len(OUT_ORDER)):
            if np.isnan(mat[i, j]):
                axb.add_patch(Rectangle((j, len(EXPO_ORDER) - 1 - i), 1, 1,
                             fc="#EEEEEE", ec="white", lw=1.5))
            else:
                p = pmat[i, j]
                col = C_SIG if (p is not np.nan and p < 0.05) else C_KEEP
                axb.add_patch(Rectangle((j, len(EXPO_ORDER) - 1 - i), 1, 1,
                             fc=col, ec="white", lw=1.5, alpha=0.85))
    axb.set_xlim(0, len(OUT_ORDER)); axb.set_ylim(0, len(EXPO_ORDER))
    axb.set_xticks([j + 0.5 for j in range(len(OUT_ORDER))])
    axb.set_xticklabels([OUT_LABEL[o] for o in OUT_ORDER], fontsize=7.5)
    axb.set_yticks([i + 0.5 for i in range(len(EXPO_ORDER))])
    axb.set_yticklabels([EXPO_LABEL[e] for e in EXPO_ORDER][::-1], fontsize=8)
    axb.tick_params(length=0)
    for s in ["top", "right", "left", "bottom"]:
        axb.spines[s].set_visible(False)
    axb.text(0.02, -0.22, "■ significant (P<0.05)", transform=axb.transAxes,
             fontsize=7, color=C_SIG)
    axb.text(0.55, -0.22, "■ null", transform=axb.transAxes,
             fontsize=7, color=C_KEEP)
    axb.text(0.75, -0.22, "■ excluded", transform=axb.transAxes,
             fontsize=7, color="#999")

    # (c) IV 计数与 F
    axc = fig.add_subplot(gs[1, 1])
    axc.set_title("(c) Instrument count per exposure",
                  loc="left", fontsize=10, fontweight="bold")
    names = [EXPO_LABEL[e] for e in EXPO_ORDER]
    counts = [IV_COUNT[e] for e in EXPO_ORDER]
    bars = axc.barh(range(len(names)), counts, color=C_SPINE, alpha=0.85, height=0.65)
    for i, (b, c) in enumerate(zip(bars, counts)):
        axc.text(c + 0.6, i, str(c), va="center", fontsize=8, color="#333")
    axc.set_yticks(range(len(names)))
    axc.set_yticklabels(names, fontsize=8)
    axc.invert_yaxis()
    axc.set_xlabel("Number of independent instruments (F > 10)", fontsize=8.5)
    axc.set_xlim(0, max(counts) * 1.18)
    for s in ["top", "right"]:
        axc.spines[s].set_visible(False)

    fig.suptitle("Figure 1", fontsize=12, fontweight="bold", x=0.02, ha="left", y=0.98)
    fig.savefig(os.path.join(FIG, "Figure1.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("[written] Figure1.png", flush=True)


# ============================ Figure 2 ============================
def figure2():
    df = load_ivw()
    rows = []
    for o in OUT_ORDER:
        for e in EXPO_ORDER:
            sub = df[(df["exposure_abbr"] == e) & (df["outcome_abbr"] == o)]
            if len(sub):
                r = sub.iloc[0]
                rows.append((e, o, r["OR"], r["OR_lo95"], r["OR_hi95"], r["pval"]))
    rows.sort(key=lambda x: (OUT_ORDER.index(x[1]), EXPO_ORDER.index(x[0])))
    fig, ax = plt.subplots(figsize=(8.5, 12))
    y = np.arange(len(rows))
    for i, (e, o, or_, lo, hi, p) in enumerate(rows):
        col = C_SIG if p < 0.05 else C_SPINE  # 红=显著，蓝=空（非显著）
        ax.plot([lo, hi], [i, i], color=col, lw=1.2, alpha=0.85)
        ax.scatter(or_, i, color=col, s=28, zorder=3, edgecolors="white", linewidths=0.4)
    ax.axvline(1, color="#333", ls="--", lw=1)
    ax.set_xscale("log")
    ax.set_xlim(0.02, 50)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{EXPO_LABEL[e]} → {o}" for e, o, *_ in rows], fontsize=7.8)
    # 分组分隔
    prev = None
    for i, (_, o, *_) in enumerate(rows):
        if prev is not None and o != prev:
            ax.axhline(i - 0.5, color="#CCC", lw=0.8)
        prev = o
    # 组标签
    for o in OUT_ORDER:
        idxs = [i for i, r in enumerate(rows) if r[1] == o]
        if idxs:
            ax.text(0.015, np.mean(idxs), OUT_LABEL[o].replace("\n", " "),
                    transform=ax.get_yaxis_transform(), fontsize=8.5,
                    fontweight="bold", va="center", color="#222",
                    bbox=dict(fc="#EEF2F6", ec="none", pad=2))
    ax.set_xlabel("Odds ratio (95% CI), IVW — log scale", fontsize=9)
    ax.set_title("Figure 2. Forest plot of primary IVW estimates (45 pairs)",
                 loc="left", fontsize=10.5, fontweight="bold")
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.text(0.99, 0.01, "● P<0.05 (all in SPONDINF, n=68)", transform=ax.transAxes,
            ha="right", fontsize=7.5, color=C_SIG)
    fig.savefig(os.path.join(FIG, "Figure2.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("[written] Figure2.png", flush=True)


# ============================ Figure 3 ============================
def figure3():
    df = pd.read_csv(os.path.join(RES, "MR_results_all.csv"))
    df = df[df["direction"] == "forward"]
    fig = plt.figure(figsize=(11, 9))
    gs = gridspec.GridSpec(3, 1, hspace=0.5)

    # (a) SPONDINF 三信号方法比较
    axa = fig.add_subplot(gs[0])
    axa.set_title("(a) SPONDINF nominal signals: method comparison",
                  loc="left", fontsize=10, fontweight="bold")
    sigs = [("Leukocyte", "WBC"), ("Neutrophil", "NEUT"), ("CCL2", "CCL2")]
    methods = ["IVW (random effects)", "MR-Egger", "Weighted median", "Weighted mode"]
    mcolors = ["#1B5E9F", "#E65100", "#2E7D32", "#7B1FA2"]
    x = np.arange(len(sigs))
    w = 0.2
    for mi, m in enumerate(methods):
        ors = []
        for _, (lab, e) in enumerate(sigs):
            sub = df[(df["exposure_abbr"] == e) & (df["outcome_abbr"] == "SPONDINF") &
                     (df["method"].str.contains(m.split(" (")[0] if "IVW" not in m else "IVW", na=False, regex=False))]
            if m == "IVW (random effects)":
                sub = sub[sub["method"].str.contains("random", na=False)]
            ors.append(float(sub.iloc[0]["OR"]) if len(sub) else np.nan)
        axa.bar(x + (mi - 1.5) * w, ors, w, label=m, color=mcolors[mi], alpha=0.88)
    axa.axhline(1, color="#333", ls="--", lw=1)
    axa.set_xticks(x)
    axa.set_xticklabels([s[0] for s in sigs], fontsize=8.5)
    axa.set_ylabel("Odds ratio", fontsize=8.5)
    axa.set_yscale("log")
    axa.legend(fontsize=7, ncol=4, loc="upper right", framealpha=0.9)
    for s in ["top", "right"]:
        axa.spines[s].set_visible(False)
    axa.text(0.02, 0.92, "Only IVW reaches P<0.05;\nmedian/mode null",
             transform=axa.transAxes, ha="left", fontsize=7.5, color=C_SIG)

    # (b) WBC 方向不一致
    axb = fig.add_subplot(gs[1])
    axb.set_title("(b) Directional inconsistency: leukocyte count across outcomes",
                  loc="left", fontsize=10, fontweight="bold")
    ors, los, his, labs = [], [], [], []
    for o in OUT_ORDER:
        sub = df[(df["exposure_abbr"] == "WBC") & (df["outcome_abbr"] == o) &
                 (df["method"].str.contains("IVW", na=False)) &
                 (df["method"].str.contains("random", na=False))]
        if len(sub):
            r = sub.iloc[0]
            ors.append(r["OR"]); los.append(r["OR_lo95"]); his.append(r["OR_hi95"])
            labs.append(OUT_LABEL[o].replace("\n", " "))
    yi = np.arange(len(ors))
    cols = [C_SIG if (lo > 1 or hi < 1) else C_NULL for lo, hi in zip(los, his)]
    for i in range(len(ors)):
        axb.plot([los[i], his[i]], [i, i], color=cols[i], lw=1.4, alpha=0.85)
    axb.scatter(ors, yi, c=cols, s=40, zorder=3, edgecolors="white", linewidths=0.5)
    axb.axvline(1, color="#333", ls="--", lw=1)
    axb.set_xscale("log"); axb.set_xlim(0.01, 60)
    axb.set_yticks(yi); axb.set_yticklabels(labs, fontsize=8)
    axb.set_xlabel("Odds ratio (95% CI) — leukocyte count", fontsize=8.5)
    for s in ["top", "right"]:
        axb.spines[s].set_visible(False)
    axb.text(0.98, 0.04,
             "Direction flips: protective in SPONDINF (n=68) → risk-increasing in discitis",
             transform=axb.transAxes, ha="right", fontsize=7.5, color=C_SIG)

    # (c) LOO for WBC -> SPONDINF
    axc = fig.add_subplot(gs[2])
    axc.set_title("(c) Leave-one-out: leukocyte count → SPONDINF",
                  loc="left", fontsize=10, fontweight="bold")
    loo = pd.read_csv(os.path.join(RES, "loo", "forward_WBC__SPONDINF_loo.csv"))
    loo_or_lo = np.exp(loo["lo95"]); loo_or_hi = np.exp(loo["hi95"])
    yi = np.arange(len(loo))
    axc.errorbar(loo["OR"], yi, xerr=[loo["OR"] - loo_or_lo,
                 loo_or_hi - loo["OR"]], fmt="o", color=C_SPINE,
                 ecolor=C_SPINE, alpha=0.7, markersize=3, capsize=0)
    full = df[(df["exposure_abbr"] == "WBC") & (df["outcome_abbr"] == "SPONDINF") &
              (df["method"].str.contains("random", na=False)) &
              (df["method"].str.contains("IVW", na=False))]
    if len(full):
        axc.axvline(float(full.iloc[0]["OR"]), color=C_SIG, ls="--", lw=1.2,
                    label=f"all-SNP OR={float(full.iloc[0]['OR']):.2f}")
    axc.axvline(1, color="#333", ls="-", lw=0.8)
    axc.set_xscale("log"); axc.set_xlim(0.005, 2)
    axc.set_yticks(yi[::2])
    axc.set_yticklabels([str(s) for s in loo["excluded_SNP"][::2]], fontsize=6.5)
    axc.set_xlabel("Odds ratio (95% CI) per SNP excluded", fontsize=8.5)
    axc.legend(fontsize=7.5, loc="lower right")
    for s in ["top", "right"]:
        axc.spines[s].set_visible(False)

    fig.suptitle("Figure 3. Fragility of the SPONDINF nominal signals",
                 fontsize=11, fontweight="bold", x=0.02, ha="left", y=0.985)
    fig.savefig(os.path.join(FIG, "Figure3.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("[written] Figure3.png", flush=True)


# ============================ Figure 4 ============================
def figure4():
    power = pd.read_csv(os.path.join(RES, "power.tsv"), sep="\t")
    meta = pd.read_csv(os.path.join(RES, "spine_replication.tsv"), sep="\t")
    fig = plt.figure(figsize=(13.5, 8.5))
    gs = gridspec.GridSpec(1, 2, wspace=0.38)

    # (a) 功效
    axa = fig.add_subplot(gs[0])
    axa.set_title("(a) Minimum detectable OR (80% power) by outcome",
                  loc="left", fontsize=10, fontweight="bold")
    ncase = {"OM": 2125, "DISCITIS": 495, "DISCINF": 375, "VOM": 104, "SPONDINF": 68}
    out_order_p = ["OM", "DISCITIS", "DISCINF", "VOM", "SPONDINF"]
    sel = ["WBC", "MONO", "CRP", "CD40"]
    x = np.arange(len(out_order_p))
    w = 0.2
    for si, e in enumerate(sel):
        vals = []
        for o in out_order_p:
            sub = power[(power["exposure"] == e) & (power["outcome"] == o)]
            vals.append(float(sub.iloc[0]["MDE_OR"]) if len(sub) else np.nan)
        axa.bar(x + (si - 1.5) * w, vals, w, label=EXPO_LABEL[e],
                color=["#1B5E9F", "#2E7D32", "#E65100", "#7B1FA2"][si], alpha=0.88)
    axa.axhline(1.3, color=C_SIG, ls="--", lw=1)
    axa.text(4.4, 1.35, "OR=1.3", fontsize=7, color=C_SIG)
    axa.set_xticks(x)
    axa.set_xticklabels([OUT_LABEL[o] for o in out_order_p], fontsize=7.5)
    axa.set_ylabel("Minimum detectable OR", fontsize=8.5)
    axa.set_yscale("log")
    axa.legend(fontsize=7, loc="upper left")
    for s in ["top", "right"]:
        axa.spines[s].set_visible(False)

    # (b) meta forest
    axb = fig.add_subplot(gs[1])
    axb.set_title("(b) Cross-phenotype meta-analysis (random effects)",
                  loc="left", fontsize=10, fontweight="bold")
    m_all = meta[meta["note"].str.contains("all 5", na=False)].set_index("exposure")
    m_excl = meta[meta["note"].str.contains("EXCLUDING", na=False)].set_index("exposure")
    yi = np.arange(len(EXPO_ORDER))
    for i, e in enumerate(EXPO_ORDER):
        if e in m_all.index:
            r = m_all.loc[e]
            axb.plot([r["OR_lo95"], r["OR_hi95"]], [i + 0.12, i + 0.12],
                     color="#1B5E9F", lw=1.4, alpha=0.85)
            axb.scatter(r["pooled_OR"], i + 0.12, color="#1B5E9F", s=30, zorder=3,
                        edgecolors="white", linewidths=0.4)
        if e in m_excl.index:
            r = m_excl.loc[e]
            axb.plot([r["OR_lo95"], r["OR_hi95"]], [i - 0.12, i - 0.12],
                     color="#2E7D32", lw=1.4, alpha=0.85)
            axb.scatter(r["pooled_OR"], i - 0.12, color="#2E7D32", s=30, zorder=3,
                        edgecolors="white", linewidths=0.4, marker="D")
    axb.axvline(1, color="#333", ls="--", lw=1)
    axb.set_xscale("log"); axb.set_xlim(0.28, 3.7)
    axb.set_xticks([0.3, 0.5, 1.0, 2.0, 3.5])
    axb.set_xticklabels(["0.3", "0.5", "1.0", "2.0", "3.5"], fontsize=8)
    axb.minorticks_off()
    axb.set_yticks(yi)
    axb.set_yticklabels([EXPO_LABEL[e] for e in EXPO_ORDER], fontsize=8)
    axb.invert_yaxis()
    axb.set_xlabel("Pooled odds ratio (95% CI)", fontsize=8.5)
    axb.text(0.02, 0.03, "● all 5 outcomes", transform=axb.transAxes,
             fontsize=7.5, color="#1B5E9F")
    axb.text(0.55, 0.03, "◆ excl. SPONDINF", transform=axb.transAxes,
             fontsize=7.5, color="#2E7D32")
    for s in ["top", "right"]:
        axb.spines[s].set_visible(False)

    fig.suptitle("Figure 4. Statistical power and cross-phenotype meta-analysis",
                 fontsize=11, fontweight="bold", x=0.02, ha="left", y=0.98)
    fig.savefig(os.path.join(FIG, "Figure4.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("[written] Figure4.png", flush=True)


if __name__ == "__main__":
    figure1()
    figure2()
    figure3()
    figure4()
    print("=== all main figures done ===", flush=True)
