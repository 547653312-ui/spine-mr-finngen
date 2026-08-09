# Manuscript Outline — Spine Infection MR Project

**Working title:** Genetically predicted immune-cell and cytokine traits and risk of osteomyelitis, including spine-specific subtypes: a Mendelian randomization study

**Journal target (2–4 IF):** Frontiers in Immunology / Scientific Reports / Journal of Translational Medicine / BMC Genomics / Infectious Diseases and Therapy

**Core message (honest negative framing):**
> Using a spine-specific, FinnGen-based Mendelian randomization design, we find no robust causal effect of genetically predicted immune-cell counts or circulating cytokine/chemokine levels on osteomyelitis risk, including the vertebral (VOM) and discitis (DISCITIS) subtypes. The only nominally significant signals (WBC/neutrophil/eosinophil → "other infective spondylopathies", n=68) failed robustness checks and were directionally inconsistent across phenotypes, and therefore are reported as exploratory, unreliable findings rather than evidence of protection.

---

## Section plan

| File | Content |
|------|---------|
| 01_title_abstract.md | Title, running title, structured abstract, keywords |
| 02_introduction.md | Burden of spine infection; immune/cell mediators hypothesized; gap (no causal MR on spine-specific osteomyelitis); objectives + hypotheses |
| 03_methods.md | Data sources (EBI GWAS Catalog exposure; FinnGen R11 outcomes), IV selection, harmonization, MR methods (IVW/Egger/median/mode), sensitivity (Cochran Q, Egger intercept, LOO), Steiger directionality, power, meta, supplementary (drug-target/MVMR/reverse if available) |
| 04_results.md | Characteristics; forward MR main table; SPONDINF signals + fragility; sensitivity; power; Steiger; meta (±SPONDINF) |
| 05_discussion.md | Interpretation of null; comparison with literature; biological plausibility; why SPONDINF signals are unreliable; strengths |
| 06_limitations.md | Distance-prune vs LD clump; VOM/SPONDINF power; reverse-arm feasibility; CRP multi-ancestry; population (European); absence of experimental validation; triangulation recommendation |

## Flagged novelty points (embedding in Intro/Discussion)
- **N1** First MR to use *spine-specific* osteomyelitis phenotypes (VOM M46.2, DISCITIS M46.4) rather than mixed all-site OM.
- **N2** Comprehensive immune landscape: 11 traits spanning adaptive/microbicidal axes (cytokines IL6R/IL1B, chemokines CXCL10/CCL2, cell counts WBC/NEUT/LYMPH/MONO/EOS, acute-phase CRP, receptor CD40).
- **N3** Full sensitivity + directionality + power + cross-phenotype meta framework, preventing over-interpretation of small-outcome noise.
