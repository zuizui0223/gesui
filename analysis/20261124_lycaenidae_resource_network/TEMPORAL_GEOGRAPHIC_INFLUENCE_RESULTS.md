# Geographic influence on temporal plant–consumer asymmetry (2026-11-24)

## Scope
Read-only post-hoc diagnostic of the existing strict training-early H6-heldout loss table, pinned to PR #153 commit 8d1073414f222146abbf04291d7893c2415d6e25, input Git blob 9cac456111e4660b2d3dd4ecc7b1cf0e8d7421be. No model or target change, independent observations, extra catalogue selection, or new hypothesis-test p-values. H6 groups are equally weighted in the saved losses. Rank within each sample and decomposition by absolute reduction in plant-side MSE, remove 0/1/5/10 highest-contributing H6, and remove the SAME groups from the consumer endpoint. Recompute 100*(sum MSE0 - sum MSE1)/sum MSE0; these are fixed-prediction sensitivity diagnostics, NOT retraining or a prespecified test. Ranking on heldout losses is explicitly outcome-dependent and must not be used for sample exclusion in a headline model.

| Sample | decomposition | removed H6 | plant gain % | consumer gain % | contrast pp |
|---|---|---:|---:|---:|---:|
| quality | symmetric | 0 | 7.328 | -0.864 | 8.192 |
| quality | symmetric | 1 | 6.011 | -0.887 | 6.898 |
| quality | symmetric | 5 | 1.876 | -0.893 | 2.769 |
| quality | symmetric | 10 | -1.905 | -0.616 | -1.289 |
| quality | early fixed | 0 | 5.802 | -1.473 | 7.274 |
| quality | early fixed | 1 | 3.892 | -1.476 | 5.368 |
| quality | early fixed | 5 | -1.567 | -1.326 | -0.241 |
| quality | early fixed | 10 | -6.817 | -1.309 | -5.508 |
| all | symmetric | 0 | 8.228 | -0.243 | 8.472 |
| all | symmetric | 1 | 7.627 | -0.243 | 7.871 |
| all | symmetric | 5 | 4.041 | -0.144 | 4.185 |
| all | symmetric | 10 | 0.964 | -0.409 | 1.373 |
| all | early fixed | 0 | 9.867 | -0.356 | 10.223 |
| all | early fixed | 1 | 9.122 | -0.350 | 9.472 |
| all | early fixed | 5 | 4.637 | -0.381 | 5.018 |
| all | early fixed | 10 | -0.972 | -0.219 | -0.753 |

Symmetric plant endpoint positive in all five original folds in both cohorts, but only 38/77 and 55/89 individual H6 improve; top five contribute 50.0% and 39.9% of gross positive loss reductions. Early-fixed plant endpoint positive in 4/5 and 5/5 folds. Removing ten high-contribution H6 largely erases the signal; do not infer ecological heterogeneity from predictive errors alone. No proof of spatial universality, resource increase, abundance, sewage or time lag. No uncertainty intervals recomputed. Numerical results calculated from pinned CSV in a separate JS implementation; included Python script is a reproducible implementation, not executed here.

## Next independent test (proposal only)
Do not tune the existing Ferrer-Paris/Clarke endpoints. Existing alternate catalogue failed to reproduce the temporal association (expanded all200 gain +2.62%, interval -4.72 to +9.79). The present source panel has only two periods (2017–19/2020–22); it cannot identify plant-before-butterfly lag. First establish a genuinely new third period (candidate 2023–25 GBIF and compatible river monitoring) without opening biological outcomes, verify same H7 and method, effort and catalogue coverage, and retain all five butterfly families on common support. No third-period data have been acquired in this diagnostic. Freeze periods, source lineage, geographic blocks and common catalogue rules before fitting. On eligible three-period panel, compare early EQR and t1 plant-composition change predicting t2 butterfly response beyond t0 and t1 butterfly baselines; reverse direction and richness-only predictions as controls. Use separately fixed, evidence-graded host catalogues; report all, with no threshold tuning or replacement of negative results. If insufficient common support, STOP. For source-specific evidence, report association rather than wastewater effect. Preserve original main and manuscript.
