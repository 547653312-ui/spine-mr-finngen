# Methods

**Study design.** We conducted a two-sample Mendelian randomization (MR) study estimating the causal effect of genetically predicted immune and cytokine traits (exposures) on osteomyelitis risk (outcomes). The study is reported in accordance with the STROBE-MR guideline¹³. All summary statistics were obtained from publicly available consortia; no individual-level data were accessed, and the study required no ethical approval or informed consent (consistent with the fully in-silico, summary-statistic scope). Exposure and outcome GWAS were drawn from independent consortia (European and multi-ancestry exposure panels versus the Finnish FinnGen cohort), so the potential for sample overlap is minimal; any residual overlap would bias estimates towards the null and therefore does not threaten the negative conclusion³⁰.

**Data sources: exposures.** Summary statistics for eleven exposure traits were obtained from the **EBI GWAS Catalog Summary Statistics REST API** (GRCh38). Each trait's EFO identifier was manually verified for trait identity (ID-drift check) before use; the curated exposure manifest is provided in Supplementary Table S1.

| Abbrev | Trait | EFO | n (range) | Population | Notes |
|--------|-------|-----|-----------|------------|-------|
| CRP | C-reactive protein | EFO_0004458 | 8,349–15,912 | Multi-ancestry | *Exploratory only — population mismatch vs Finnish outcome* |
| IL6R | IL-6 receptor subunit α | EFO_0008187 | 21,758 | European (SCALLOP) | |
| IL1B | Interleukin-1β | EFO_0004812 | 4,910 | European | 1 IV only (Wald ratio) |
| CXCL10 | CXCL10 / IP-10 | EFO_0008056 | 3,685 | European | |
| CD40 | CD40 / TNFRSF5 | EFO_0010607 | 21,758 | European (SCALLOP) | |
| CCL2 | CCL2 / MCP-1 | EFO_0004749 | 21,758 | European (SCALLOP) | |
| NEUT | Neutrophil count | EFO_0004833 | 170,384 | European (Astle 2016) | |
| LYMPH | Lymphocyte count | EFO_0004587 | 171,643 | European (Astle 2016) | |
| MONO | Monocyte count | EFO_0005091 | 170,721 | European (Astle 2016) | |
| EOS | Eosinophil count | EFO_0004842 | 172,275 | European (Astle 2016) | |
| WBC | Leukocyte count | EFO_0004308 | 172,435 | European (Astle 2016) | |

**Data sources: outcomes.** Outcome summary statistics were obtained from the **FinnGen R11** public release (GRCh38), accessed by remote BGZF/tabix range queries without full-file download. Five phenotypes were analysed (case counts after reconciliation with the R11 API):

| Abbrev | FinnGen code | Phenotype (ICD-10) | n_cases | n_controls | Spine-specific |
|--------|--------------|--------------------|---------|------------|----------------|
| VOM | M13_OSTEOMYELVERTEB | Vertebral osteomyelitis (M46.2) | 104 | 322,314 | Yes |
| DISCITIS | M13_DISCITIS | Discitis (M46.4) | 495 | 322,314 | Yes |
| DISCINF | M13_DISCINFECTION | Pyogenic intervertebral-disc infection | 375 | 322,314 | Yes |
| SPONDINF | M13_SPONDYLOINFECTION | Other infective spondylopathies | 68 | 322,314 | Yes |
| OM | M13_OSTEOMYELITIS | Osteomyelitis, all sites (M86) | 2,125 | 429,826 | No |

**Instrument selection.** For each exposure, SNPs associated at genome-wide significance (P<5×10⁻⁸) were extracted and pruned to approximate linkage independence using a ±500 kb distance filter. Because no local LD reference panel was available in the execution environment, distance pruning substituted for standard LD clumping in the primary analysis. Palindromic SNPs with minor-allele frequency >0.42 were excluded as non-orientable. Instruments with F-statistic <10 (F = (β/se)²) were removed to limit weak-instrument bias¹⁹,²⁰.

**Formal LD-clumping sensitivity analysis.** To verify that the distance-pruned instrument sets were not affected by residual linkage, the instruments for the six traits central to the primary and drug-target analyses (IL-6 receptor, CD40, C-reactive protein, leukocyte count, neutrophil count, CCL2/MCP-1) were re-clumped against the 1000 Genomes Phase 3 European reference panel (r² threshold 0.001, GRCh38) using pairwise LD retrieved from the LDlink LDmatrix service. Stepwise clumping retained the most significant variant per LD block; drug-target cis-pQTL sets were clumped identically. Clumped sets were carried forward to re-estimate the corresponding MR analyses as a sensitivity check.

**Harmonisation.** Exposure and outcome alleles were aligned by chromosome:position (GRCh38). Strand flips and complementary-allele resolution were applied; SNPs with unresolvable allele mismatch were dropped. Effect alleles were oriented to the exposure's effect allele for all ratio estimates; the harmonisation action (same / flip / strand-flip variants) was recorded per SNP for downstream sign-consistent reuse in supplementary analyses (Supplementary Methods).

**MR estimation.** The primary estimate used the inverse-variance-weighted (IVW) method under a multiplicative random-effects model. Secondary estimators were MR-Egger (with intercept test for directional pleiotropy), the weighted median and the weighted mode — both robust to <50% invalid instruments²¹⁻²⁴. Exposure effects on binary outcomes are reported as odds ratios (OR) per genetically predicted standard-deviation increase in the trait. The primary IVW estimates for all 45 exposure–outcome pairs are tabulated in Supplementary Table S2. Analyses followed current MR reporting guidance⁴⁶.

**Sensitivity analyses.** Heterogeneity was assessed by Cochran Q (significant at P<0.05) and I² across instruments. Pleiotropy was assessed by the MR-Egger intercept and its P-value. Leave-one-out analysis re-estimated the effect excluding each SNP sequentially. Directionality was assessed by the Steiger test comparing variance in exposure vs outcome explained by the instruments²⁵; a significant exposure→outcome directionality supports the assumed causal direction. Multiple-testing control applied Benjamini–Hochberg FDR and Bonferroni correction across the 45 exposure–outcome pairs.

**Power and meta-analysis.** Post-hoc power was quantified as the minimum detectable OR (MDE_OR) at 80% power (α=0.05, two-sided), equal to exp[(z₀.₉₇₅ + z₀.₈₀)·SE_IVW]²⁶. Cross-phenotype random-effects meta-analysis (DerSimonian–Laird) pooled each exposure across the five outcomes, with a pre-specified sensitivity meta excluding the smallest outcome (SPONDINF, n=68).

**Drug-target (cis-pQTL) MR.** For three clinically tractable immune targets — IL-6 receptor (IL6R), CD40 and C-reactive protein (CRP) — variants within ±1 Mb of the gene (cis-pQTL window) were used as instruments against the two best-powered outcomes (OM and DISCITIS), linking population-scale findings to actionable therapeutic hypotheses (e.g., IL-6R antagonists, anti-CD40 biologics). Cis-variants were selected at P<1×10⁻⁴ within the cis window (a standard relaxed threshold for cis-pQTL studies, where cis-acting instruments are typically fewer and weaker than genome-wide hits), pruned for independence and retained at F>10. The full cis-IV table and per-method estimates are in Supplementary Table S3.

**Multivariable MR.** To test whether the univariable null estimates were robust to adjustment for established osteomyelitis risk factors, we performed multivariable MR (MVMR) adjusting for body-mass index (BMI_IRN, n=321,672; FinnGen R11), type-2 diabetes (T2D, n_cases=71,728; FinnGen T2D) and smoking (SMOKING, n_cases=4,271; FinnGen SMOKING). Confounder summary statistics were retrieved by remote BGZF random access to the full FinnGen summary files; the custom reader created a fresh FinnGenSumstats object per query to ensure stable variant retrieval across scattered genome-wide instrument positions. The univariable-harmonised instrument set from the primary pipeline was used, with complete-case exclusion per SNP when a confounder record was missing. A fixed-effect IVW estimator was used as the base, with multiplicative overdispersion (φ) widening the confidence intervals where φ≥1. The full confounder-adjusted MVMR estimates for all 14 exposure–outcome pairs are in Supplementary Table S4.

**Reverse MR.** Spine-infection GWAS summary statistics (from FinnGen R11) were used as the exposure and immune/cytokine traits as the outcome. Genome-wide-significant SNP counts per spine-infection phenotype were tabulated from the FinnGen release to assess instrument feasibility before conducting reverse-direction estimation.

**Software and reproducibility.** Analyses used a custom Python pipeline (`analysis/`: `gwas_io.py`, `mr_methods.py`, `mr_pipeline.py`, `supplementary.py`, `supplementary_network.py`, `mvmr_finngen.py`, `make_main_figures.py`); Python 3.13 with numpy, pandas, scipy and matplotlib. No R or OpenGWAS token was required; the IEU OpenGWAS API has required a JWT authentication token since May 2024, and FinnGen R11 + EBI GWAS Catalog served as fully open equivalents. The full pipeline and per-SNP harmonised tables are released with this manuscript (see Code and Data Availability).