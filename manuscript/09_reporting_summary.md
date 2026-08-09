# Nature Portfolio Reporting Summary — Spine-Infection MR Study

*Generated as a fillable template. The corresponding author will complete this in the Nature submission system at acceptance. The completed PDF will be uploaded as a separate Reporting Summary file per Scientific Reports editorial policy.*

---

## 1. Study design
- **Study type:** Two-sample Mendelian randomization (observational, summary-statistic, in-silico)
- **Objective:** Estimate the causal effect of 11 genetically predicted immune/cytokine traits on 5 osteomyelitis phenotypes (4 spine-specific) from FinnGen R11.
- **Hypothesis:** Pre-specified null framing — no robust causal effect of immune/cytokine traits on osteomyelitis risk.

## 2. Data sources
- **Exposures (11 immune/cytokine traits):** EBI GWAS Catalog Summary Statistics REST API (GRCh38). See Supplementary Table S1.
- **Outcomes (5 osteomyelitis phenotypes):** FinnGen R11 release (GRCh38), accessed by remote BGZF/tabix random access. See Methods table.
- **Confounders for MVMR:** FinnGen R11 type-2 diabetes (T2D, n_cases=71,728); smoking (SMOKING, n_cases=4,271); body-mass index (BMI_IRN, n=321,672).

## 3. Sample sizes
| Resource | n / n_cases |
|----------|-------------|
| Immune/cytokine exposures | 4,910 – 172,435 (per trait, see Supplementary Table S1) |
| FinnGen VOM (M46.2) | 104 cases / 322,314 controls |
| FinnGen DISCITIS (M46.4) | 495 / 322,314 |
| FinnGen DISCINF (pyogenic disc infection) | 375 / 322,314 |
| FinnGen SPONDINF (other infective spondylopathies) | 68 / 322,314 |
| FinnGen OM (all-site osteomyelitis, M86) | 2,125 / 429,826 |

## 4. Statistical methods
- **Primary:** Inverse-variance-weighted (IVW) under multiplicative random effects.
- **Secondary:** MR-Egger, weighted median, weighted mode.
- **Sensitivity:** Cochran Q + I², MR-Egger intercept, leave-one-out, Steiger directionality, Benjamini-Hochberg FDR + Bonferroni (45 tests).
- **Multivariable:** MVMR-IVW adjusting for T2D and smoking (and BMI for the IL6R sensitivity).
- **Meta-analysis:** DerSimonian-Laird random effects across the five phenotypes, with pre-specified exclusion of the n=68 SPONDINF outcome.
- **Power:** Post-hoc minimum detectable OR at 80% power (two-sided α=0.05).
- **Software:** Python 3.13 with numpy, pandas, scipy, matplotlib. No R, no OpenGWAS token. FinnGen R11 + EBI GWAS Catalog as fully open equivalents.

## 5. Reporting guidelines followed
- **STROBE-MR** (Skrivankova et al., JAMA 2021) — completed checklist provided as a Supplementary file.

## 6. Reproducibility
- **Code:** Public GitHub repository: https://github.com/547653312-ui/spine-mr-finngen (MIT licence). Archived version with persistent identifier: https://doi.org/10.5281/zenodo.21863390.
- **Data:** All summary statistics are publicly available (no controlled-access data used).
- **Per-SNP harmonised tables:** Released alongside the code.

## 7. Conflicts of interest
None declared.

## 8. Ethics
This study used only publicly available, de-identified GWAS summary statistics. No individual-level data were accessed. The work is therefore exempt from institutional ethical approval and informed-consent requirements under the policies of the source consortia (EBI GWAS Catalog, FinnGen).

## 9. Funding
[To be added at acceptance.]

## 10. Key results summary
- **45** exposure–outcome pairs analysed (9 exposures × 5 outcomes).
- **0** associations survived FDR/Bonferroni correction or pleiotropy-robust estimators.
- **3** nominally significant IVW signals, all in the smallest outcome (SPONDINF, n=68), all non-robust, all directionally inconsistent across phenotypes.
- **Drug-target MR (IL6R, CD40):** null for OM and DISCITIS.
- **MVMR (T2D + smoking adjusted):** null, consistent with univariable.
- **Reverse MR:** not feasible (insufficient genome-wide instruments in spine-infection GWAS).

## 11. Conclusions
No robust causal evidence that genetically predicted circulating immune effector levels cause osteomyelitis. The findings argue against systemic immune modulation as a prevention strategy and support triangulation with clinical cohorts before any therapeutic inference.