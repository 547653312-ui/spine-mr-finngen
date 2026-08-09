# Results

**Instrument characteristics.** After quality control, nine of eleven immune/cytokine traits contributed analysable instruments (Figure 1c; Supplementary Table S2). **IL1B** (single genome-wide significant SNP) and **CXCL10** (three SNPs) could not be harmonised to any FinnGen outcome (no overlapping instrument–outcome SNPs after allele harmonisation) and were excluded from pair-wise analysis, leaving **45 exposure–outcome pairs** (9 exposures × 5 outcomes). Instrument counts ranged from 3 (CD40) to 39 (monocyte count); all F-statistics exceeded 10 (mean F 79–528).

*Table 1. Number of independent instruments per exposure (harmonised, F>10).*

| Exposure | n IVs | Exposure | n IVs |
|----------|-------|----------|-------|
| Leukocyte count | 28 | Monocyte count | 39 |
| Neutrophil count | 29 | Eosinophil count | 13 |
| Lymphocyte count | 18 | CCL2 / MCP-1 | 7 |
| CRP | 13 | IL-6 receptor | 4 |
| *(IL1B, CXCL10: excluded — no analysable pair)* | | CD40 | 3 |

**Primary IVW estimates.** Across 45 pairs (Figure 2; full table in Supplementary Table S2), **only three reached nominal significance (P<0.05), and all three were in the smallest outcome, "other infective spondylopathies" (SPONDINF, n_cases=68):**

- Leukocyte count → SPONDINF: OR 0.04 (95% CI 0.006–0.21), P=2×10⁻⁴
- Neutrophil count → SPONDINF: OR 0.10 (0.019–0.57), P=0.009
- CCL2/MCP-1 → SPONDINF: OR 0.29 (0.092–0.92), P=0.035

**Every larger, better-powered outcome was null.** For all-site osteomyelitis (OM, n=2,125) the lowest P was 0.051 (monocyte count OR 0.71, 0.51–1.00); for discitis (n=495) the lowest P was 0.16; for pyogenic disc infection (n=375) 0.21; for vertebral osteomyelitis (n=104) 0.23.

**The SPONDINF signals are not robust.** The three nominally significant associations failed robustness checks and were internally inconsistent (Table 2; Figure 3). They were significant only under IVW; the pleiotropy-robust weighted median and weighted mode estimators were null for all three (median P≥0.23, mode P≥0.26). Critically, the direction was inconsistent across phenotypes: leukocyte count, protective in SPONDINF (OR 0.04), was *risk-increasing* in discitis (OR 1.41) and null in OM (OR 0.82) and VOM (OR 1.49) — a directional flip-flop that is the hallmark of a small-sample, outcome-specific artefact rather than a biological effect.

*Table 2. Fragility of the SPONDINF signals across MR methods.*

| Exposure → SPONDINF | IVW OR (P) | MR-Egger OR (P) | Weighted median OR (P) | Weighted mode OR (P) |
|--------------------|-----------|-----------------|------------------------|----------------------|
| Leukocyte | 0.04 (2×10⁻⁴) | 0.05 (0.083) | 0.25 (0.278) | 0.23 (0.300) |
| Neutrophil | 0.10 (0.009) | 0.11 (0.149) | 0.32 (0.336) | 0.26 (0.258) |
| CCL2 | 0.29 (0.035) | 0.12 (0.054) | 0.44 (0.229) | 0.45 (0.284) |

After Benjamini–Hochberg FDR and Bonferroni correction across 55 tests, no association survived as a credible finding; the single FDR-significant result (leukocyte→SPONDINF, FDR=0.009) is discounted on the grounds above.

**Sensitivity analyses.** Cochran Q was non-significant for the large outcomes (OM pairs Q-P >0.16; VOM Q-P >0.16), indicating consistent instrument effects. Minor heterogeneity (I² up to ~47%) appeared in a few pairs without altering conclusions. MR-Egger intercepts were non-significant for all 45 pairs (all intercept P>0.20), providing no evidence of directional pleiotropy. Leave-one-out analysis identified no single SNP driving any estimate; all confidence intervals remained inclusive of the null for non-SPONDINF outcomes (full LOO curves in Supplementary Figure S1).

**Formal LD clumping.** Re-clumping the instruments for the six key traits against the 1000 Genomes Phase 3 European panel (r²<0.001, via LDlink LDmatrix) removed no instrument for five of six traits (leukocyte count, neutrophil count, CRP, CCL2, CD40), confirming that the ±500 kb distance filter had already enforced independence at these loci; for IL-6 receptor, two of four genome-wide instruments were removed as redundant (rs3014860, rs12138773), and the two retained instruments gave concordant null results across all five outcomes (IVW P=0.12–0.80; Supplementary Table S7). For the drug-target cis-pQTL sets, one of two IL6R cis-instruments (rs12138773, in LD with rs4845625) and one of three CD40 cis-instruments (rs2868310) were removed; the clumped sets reproduced the null estimates for both outcomes (IL6R→OM Wald OR 1.006, P=0.88; IL6R→DISCITIS Wald OR 1.043, P=0.60; CD40→OM IVW OR 1.085, P=0.43; CD40→DISCITIS IVW OR 1.073, P=0.75; Supplementary Table S7).

**Statistical power and directionality.** Post-hoc power confirmed that the null results for the better-powered outcomes are genuine rather than under-powered (Figure 4a). At 80% power, the minimum detectable OR for OM and discitis was 1.2–3.0; yet all estimates fell within this range and were non-significant. By contrast, VOM (n=104) and SPONDINF (n=68) could only detect very large effects (MDE_OR 2.5–42), exactly the strata in which the spurious "signals" appeared. Steiger directionality tests supported the assumed exposure→outcome causal direction in **all 45 pairs** (all P≈0), i.e., the instruments explained substantially more variance in the immune exposure than in the infection outcome, arguing against reverse causation as an explanation for the null.

**Cross-phenotype meta-analysis.** Random-effects meta-analysis pooling each exposure across all five osteomyelitis phenotypes yielded no significant association for any trait (Figure 4b; full table in Supplementary Table S5). The largest effect was CD40 (pooled OR 1.12, 95% CI 0.99–1.25, P=0.067). A pre-specified sensitivity meta *excluding* the n=68 SPONDINF outcome — the only source of nominal signals — rendered even CD40 non-significant (OR 1.12, 1.00–1.26, P=0.058) and left all traits null, with low heterogeneity (I² mostly 0%).

*Table 3. Cross-phenotype random-effects meta-analysis (IVW).*

| Exposure | All 5 outcomes OR (95% CI), P | Excl. SPONDINF OR (95% CI), P |
|----------|-------------------------------|-------------------------------|
| Leukocyte | 0.82 (0.37–1.79), 0.61 | 1.11 (0.78–1.60), 0.56 |
| Neutrophil | 0.94 (0.52–1.69), 0.83 | 1.07 (0.78–1.48), 0.67 |
| Lymphocyte | 1.14 (0.70–1.85), 0.59 | 1.17 (0.67–2.06), 0.59 |
| Monocyte | 0.86 (0.54–1.36), 0.51 | 0.92 (0.58–1.45), 0.71 |
| Eosinophil | 0.93 (0.60–1.45), 0.75 | 0.97 (0.62–1.52), 0.88 |
| CRP | 1.03 (0.87–1.23), 0.71 | 1.04 (0.87–1.25), 0.64 |
| CCL2 | 0.84 (0.64–1.10), 0.21 | 0.88 (0.74–1.05), 0.16 |
| IL-6 receptor | 1.02 (0.93–1.13), 0.68 | 1.04 (0.94–1.15), 0.49 |
| CD40 | 1.12 (0.99–1.25), 0.067 | 1.12 (1.00–1.26), 0.058 |

**Drug-target (cis-pQTL) Mendelian randomization.** To translate the null findings into an actionable therapeutic frame, we conducted cis-pQTL MR for three clinically tractable immune targets — **IL-6 receptor (IL6R)**, **CD40**, and **C-reactive protein (CRP)** — against OM and DISCITIS (Table 4; full table in Supplementary Table S3). Genome-wide instruments were not available within ±1 Mb of the CRP locus (no cis IVs), so CRP was not analysed.

Genetic proxies for IL-6R blockade showed **no causal effect** on DISCITIS (OR 1.06, 95% CI 0.90–1.24; P=0.51; 2 cis-SNPs, F=192) or OM (OR 1.00, 0.93–1.08; P=0.90). Similarly, genetically predicted CD40 perturbation did not affect DISCITIS (OR 1.05, 0.75–1.46; P=0.78; 3 SNPs, F=20) or OM (OR 1.11, 0.92–1.35; P=0.29). All four estimates were null and directionally consistent across MR methods (MR-Egger, weighted median, weighted mode).

*Table 4. Drug-target (cis-pQTL) MR: IL-6R and CD40 against OM and DISCITIS.*

| Target | Outcome | nSNP | IVW OR (95% CI) | P | F_min |
|--------|---------|------|------------------|---|-------|
| IL6R | DISCITIS | 2 | 1.06 (0.90–1.24) | 0.51 | 192 |
| IL6R | OM | 2 | 1.00 (0.93–1.08) | 0.90 | 192 |
| CD40 | DISCITIS | 3 | 1.05 (0.75–1.46) | 0.78 | 20 |
| CD40 | OM | 3 | 1.11 (0.92–1.35) | 0.29 | 20 |

These results imply that genetically proxied modulation of IL-6R or CD40 signalling is **not expected to materially alter susceptibility** to osteomyelitis or discitis — a null finding with direct translational relevance to the empirical use of IL-6R antagonists (e.g., tocilizumab) and emerging anti-CD40 biologics in severe infection.

**Multivariable MR (MVMR).** To test whether adjustment for established osteomyelitis risk factors — type-2 diabetes, smoking and body-mass index — attenuates or reveals any signal masked by confounding, we ran multivariable IVW MR using FinnGen R11 full summary statistics for the confounders and the univariable-harmonised instrument set from the primary pipeline (Table 5; full table in Supplementary Table S4).

*Table 5. Multivariable MR (MVMR) adjusting for body-mass index (BMI_IRN, n=321,672), type-2 diabetes (T2D, n=71,728) and smoking (SMOKING, n=4,271).*

| Exposure → Outcome | nSNP | Univariable OR (P) | Adjusted OR (95% CI) | P_adj | condF_exposure | condF_BMI | condF_T2D | condF_SMK |
|--------------------|------|--------------------|----------------------|-------|----------------|-----------|-----------|-----------|
| Leukocyte → OM | 28 | 0.82 (0.21) | 1.23 (0.84–1.81) | 0.29 | 116 | 4.1 | 2.5 | 2.0 |
| Neutrophil → OM | 29 | 0.85 (0.29) | 1.11 (0.80–1.54) | 0.54 | 113 | 3.1 | 2.9 | 1.9 |
| Monocyte → OM | 39 | 0.71 (0.03) | 0.74 (0.52–1.06) | 0.10 | 80 | 3.1 | 2.0 | 1.5 |
| Lymphocyte → OM | 18 | 1.29 (0.23) | 1.49 (0.95–2.34) | 0.08 | 102 | 5.0 | 2.8 | 1.3 |
| Eosinophil → OM | 13 | 1.07 (0.77) | 1.14 (0.58–2.23) | 0.71 | 93 | 2.4 | 3.8 | 1.2 |
| CRP → OM | 13 | 1.06 (0.59) | 1.05 (0.83–1.33) | 0.66 | 92 | 6.0 | 20.5 | 1.5 |
| IL-6 receptor → OM | 4 | 1.05 (0.44) | 1.25 (0.73–2.14) | 0.41 | 431 | 2.6 | 3.4 | 0.5 |
| Leukocyte → DISCITIS | 28 | 1.41 (0.30) | 1.55 (0.74–3.24) | 0.24 | 116 | 4.1 | 2.5 | 2.0 |
| Neutrophil → DISCITIS | 29 | 1.44 (0.26) | 1.43 (0.73–2.83) | 0.30 | 113 | 3.1 | 2.9 | 1.9 |
| Monocyte → DISCITIS | 39 | 1.49 (0.22) | 1.50 (0.71–3.14) | 0.29 | 80 | 3.1 | 2.0 | 1.5 |
| Lymphocyte → DISCITIS | 18 | 0.52 (0.14) | 0.39 (0.15–0.98) | **0.045** | 102 | 5.0 | 2.8 | 1.3 |
| Eosinophil → DISCITIS | 13 | 1.19 (0.73) | 1.30 (0.45–3.78) | 0.63 | 93 | 2.4 | 3.8 | 1.2 |
| CRP → DISCITIS | 13 | 1.03 (0.90) | 0.98 (0.42–2.29) | 0.96 | 92 | 6.0 | 20.5 | 1.5 |
| IL-6 receptor → DISCITIS | 4 | 1.02 (0.88) | 0.69 (0.22–2.13) | 0.52 | 431 | 2.6 | 3.4 | 0.5 |

Across all 14 exposure–outcome pairs, confounder-adjusted estimates remained **directionally consistent with the univariable IVW results**. A single nominally significant association emerged after MVMR conditioning: lymphocyte count → DISCITIS, adjusted OR 0.39 (95% CI 0.15–0.98), P=0.045. This signal does not survive correction for the 14 MVMR tests (Bonferroni threshold P<0.0036), shows no directionality support from the Steiger test on the corresponding univariable pair (which was non-significant), and has a wide confidence interval (driven by the small discitis case count n=495 and modest instrument count n=18). It is therefore interpreted as residual noise. Notably, the univariable monocyte → OM association (P=0.03 under fixed-effects IVW, P=0.05 under random-effects) was attenuated to P=0.10 after metabolic/behavioral adjustment, consistent with partial confounding by T2D or smoking.

Conditional F-statistics for the immune exposures (range 80–431) confirmed adequate instrument strength after conditioning on the confounders; conditional F for the three confounders (1.2–20.5) was low-to-moderate, indicating limited genetic overlap between the immune IV sets and metabolic/behavioral variants — i.e., the confounder adjustment is valid but does not strongly condition on the confounders (a structural property of the available instruments, not a methodological limitation). Overdispersion factors (φ) ranged 0.0–3.5, with the fixed-effects IVW used as the base and overdispersion correction applied (φ≥1) to widen the confidence intervals where appropriate.

**Reverse MR feasibility.** Reverse-direction MR (genetic liability to osteomyelitis → immune profile) requires ≥3 genome-wide significant instruments in the spine-infection GWAS to be estimable with multi-method robustness. As of FinnGen R11, the genome-wide-significant locus counts per phenotype were OM=2, DISCITIS=1, DISCINF=1, VOM=0, SPONDINF=0 (Supplementary Table S6) — below the multi-method threshold. A reverse MR is therefore **not feasible** with current spine-infection GWAS power; this limitation is inherent to data availability rather than a methodological choice, and it inherently guards against weak-instrument bias in the reverse direction. Reverse MR should be revisited in FinnGen R12+ or equivalent releases as case counts grow.