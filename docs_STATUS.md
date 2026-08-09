# Spine Infection MR — Project STATUS (top-level)

**Goal:** Pure in-silico Mendelian randomization of immune/cytokine traits → osteomyelitis (incl. spine-specific VOM/DISCITIS), 2–4 IF SCI.

**Roles (per user request):** kimi-K3 = planner/reviewer (plan + audits); Hy3 = executor (pipeline + analyses + manuscript). Final audit = kimi-K3 (`review/REVIEW_FINAL.md` → ACCEPT, submission-ready core).

## Completed (core, sufficient for submission)
- [x] Data feasibility: GWAS availability confirmed (FinnGen R11 spine phenotypes + EBI immune traits). `FEASIBILITY.md`
- [x] Forward MR: 9 immune/cytokine traits × 5 outcomes = 45 pairs. `results/MR_results_*.csv`, `results/figures/` (196 plots)
- [x] Data audit / manifest: trait-ID verification, population-mismatch & power flags. `data/manifest_*.tsv`
- [x] Sensitivity: IVW + Egger + median + mode + Cochran Q + Egger intercept + LOO (all per pair)
- [x] Steiger directionality: 45/45 exposure→outcome supported. `results/steiger.tsv`
- [x] Power (MDE @80%): OM/DISCITIS well-powered → genuine null. `results/power.tsv`
- [x] Cross-phenotype meta (±SPONDINF): all null. `results/spine_replication.tsv`
- [x] Manuscript (7 files): `manuscript/00_outline … 06_limitations`. Honest negative framing; SPONDINF signals flagged unreliable.
- [x] kimi-K3 final review: ACCEPT, submission-ready. `review/REVIEW_FINAL.md`

## Supplementary analyses (network, task zyaDt2 + MVMR re-run juE04k — DONE)
- [x] Drug-target (cis-pQTL) MR: **COMPLETED, all null.** IL6R→DISCITIS OR 1.06 (P=0.51), IL6R→OM OR 1.00 (P=0.90); CD40→DISCITIS OR 1.05 (P=0.78), CD40→OM OR 1.11 (P=0.29); CRP no cis IVs. Integrated into Results §7 + Methods §9. `results/drug_target.tsv`
- [x] MVMR: **NOT FEASIBLE** — confounder summary stats (smoking/BMI/T2D) unretrievable from public server (rate-limited). Honestly documented, not reported. `results/mvmr.tsv`
- [x] Reverse MR: **NOT FEASIBLE** — <3 genome-wide-significant spine loci per outcome. Documented as a limitation. `results/reverse_MR_immune.tsv`
- Both non-feasible items are reported transparently as limitations, not missing analyses; core paper stands.

## Headline finding
No robust causal effect of genetically predicted immune-cell/cytokine traits on osteomyelitis, including spine-specific VOM/DISCITIS. Only nominally significant signals (leukocyte/neutrophil/CCL2 → "other infective spondylopathies", n=68) failed robustness and were directionally inconsistent → reported as unreliable, not as a protective effect.

## How to re-run
```
cd spine_MR/analysis
./venv/Scripts/python.exe supplementary.py           # power/steiger/meta (local)
./venv/Scripts/python.exe mr_pipeline.py --direction forward   # forward MR
./venv/Scripts/python.exe supplementary_network.py   # drug-target/MVMR/reverse (network)
```
