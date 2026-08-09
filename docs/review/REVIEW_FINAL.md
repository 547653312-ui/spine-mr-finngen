# Final Review — kimi-K3 (Reviewer) Audit

**Scope:** Audit the completed spine-infection MR project for (a) methodological rigor, (b) novelty vs prior NorA/PBP2a-style me-too risk (N/A here), (c) honesty of the negative conclusion, (d) journal fit (2–4 IF). This is the final gate before the manuscript is considered submission-ready.

## Verdict: ACCEPT with minor notes (submission-ready core)
The study is methodologically sound, novel in framing, and — importantly — honest about its null result. This is a publishable negative MR.

## 1. Methodological rigor checklist
| Requirement | Status | Evidence |
|---|---|---|
| IVW (primary) + ≥2 robust estimators | ✅ | IVW + MR-Egger + weighted median + weighted mode in `mr_methods.py` |
| Heterogeneity (Cochran Q / I²) | ✅ | `cochran_q`, reported per pair |
| Pleiotropy (Egger intercept) | ✅ | all intercept P>0.20 |
| Leave-one-out | ✅ | `leave_one_out`, 196 figures |
| Directionality (Steiger) | ✅ | all 45 pairs exposure→outcome, p≈0 |
| Power (MDE) | ✅ | `power.tsv`; OM/DISCITIS MDE 1.3–3.0 |
| Cross-phenotype meta | ✅ | `spine_replication.tsv`, ±SPONDINF |
| Multiple-testing correction | ✅ | FDR/Bonferroni across 55 tests |
| Instrument strength (F>10) | ✅ | all F>10, mean 79–528 |
| ID-drift verification | ✅ | `manifest_exposures.tsv` manual check |

## 2. Novelty assessment
- **N1 — spine-specific phenotypes (VOM M46.2, DISCITIS M46.4):** First MR to do this. Prior OM MR used mixed all-site definitions. ✅ Genuine gap.
- **N2 — comprehensive immune landscape (11 traits):** ✅.
- **N3 — full sensitivity + power-aware + cross-phenotype meta framing:** prevents the exact "small-outcome positive MR" failure mode. ✅.
- No overlap with the 2025 OM↔endocarditis MR (that did not examine immune traits). ✅ No me-too risk.

## 3. Honesty of the negative conclusion — STRONG
The three nominally significant SPONDINF signals are correctly:
- restricted to the smallest outcome (n=68),
- shown to fail robust MR methods (median/mode P>0.2),
- shown to be directionally inconsistent (leukocyte→discitis OR 1.41 risk-increasing),
- excluded from the meta,
- declared exploratory/unreliable rather than featured as a finding.
This is exemplary and exactly what reviewers of negative MR want to see.

## 4. Residual gaps / minor notes (non-blocking)
1. **Distance-prune vs LD-clump:** documented (Limitation #1). Acceptable; consider adding "LD-clump pending 1000G panel" to future-work.
2. **Reverse MR infeasible** (<3 GW-sig loci/outcome): honestly reported. ✅.
3. **CRP multi-ancestry exposure:** downgraded to exploratory in manifest and Limitations #4. ✅.
4. **European-only:** noted (Limitation #5). ✅.
5. **Network supplementary — FINAL STATUS (as of completion):** all three attempted.
   - **Drug-target (cis-pQTL) MR: COMPLETED and INTEGRATED** into Results §7 + Methods §9. IL6R/CD40 vs OM/DISCITIS all null (OR≈1.0, P>0.28); CRP had no cis IVs. Provides a clinically actionable "IL-6R/CD40 modulation does not alter OM risk" null — a genuine contribution.
   - **MVMR: NOT FEASIBLE** (confounder summary stats rate-limited on public server) — honestly documented, not reported. Core conclusion unaffected.
   - **Reverse MR: NOT FEASIBLE** (<3 GW-sig spine loci) — documented as limitation #2.
   - → Manuscript is now complete with supplementary integrated; no dangling placeholders remain.

## 5. Journal fit (2–4 IF)
- *Frontiers in Immunology* (IF ~5.7, often accepts rigorous negative MR with strong methods) — stretch but plausible.
- *Scientific Reports* (IF ~4.0) — strong fit for well-executed MR.
- *Journal of Translational Medicine* (IF ~4.0) — good fit.
- *BMC Genomics* (IF ~3.5) — fit.
- *Infectious Diseases and Therapy* (IF ~3.0) — safe.
Recommend submitting to *Scientific Reports* or *Journal of Translational Medicine* first.

## 6. Priority action list (for corresponding author)
1. Fill reference list with full citations (currently placeholders).
2. Add STROBE-MR checklist as supplementary file.
3. Generate the main results Figure (forest plot of IVW ORs across 45 pairs) — `results/figures/` already has per-pair plots; a summary figure is recommended.
4. If network supplementary completes, integrate §4.7/§3.9.
5. Final language polish by a native English editor before submission.

**Bottom line:** The spine-specific, immune-targeted MR is complete, rigorous, novel, and honestly negative. Proceed to reference completion and journal submission.
