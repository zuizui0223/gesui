# 2026-11-24: strict training-only early geographic reference

Exploratory reuse of the original source-locked 222 repeated H7 (89 H6); primary source coverage subset 173 H7 (77 H6). No additional field observations. Prior full/early reference results preserved unchanged.

The null graph is defined using the early plant occurrence profiles of the training paired H7 only; no target H6 or late occurrence profile is used. Within the original plant family, frequency-bin and six geographic-cluster strata, use the exact mean adjacency of uniform permutations (not a Monte Carlo approximation). Both training and held-out endpoints are evaluated against the fold-specific reference, so they are not a common fixed target across folds. Static original host catalogue remains known. The fold-assigned H6 heldout Ridge10, covariates, early network baseline, group-equal losses and independent sklearn check are retained. Country test removes all shared H6 from training. Fixed-loss bootstrap does not propagate reference construction/model uncertainty.

Plant/consumer terms are symmetric decomposition, or exact early-fixed components plus interaction. All 454 baseline source reconstruction and additive identities verified.

| Endpoint | 173 H7 gain (95% fixed-loss interval) | 222 H7 gain (95% fixed-loss interval) |
|---|---:|---:|
| Plant symmetric | 7.328 [1.160,13.444] | 8.228 [2.650,14.267] |
| Consumer symmetric | -0.864 [-2.865,0.284] | -0.243 [-1.316,0.402] |
| Plant, early consumer fixed | 5.802 [-2.549,12.313] | 9.867 [1.629,17.673] |
| Consumer, early plant fixed | -1.473 [-4.956,-0.533] | -0.356 [-1.549,0.081] |
| Interaction | -1.510 [-5.077,-0.450] | -1.008 [-3.867,0.032] |

Plant symmetric minus consumer symmetric paired contrast Bonferroni-3 intervals +0.266 to +16.041 points (quality), +1.645 to +16.323 (all). Early-fixed plant minus consumer is uncertain in the quality subset, positive in all paired H7.

Strict reference country-equal plant symmetric gains: quality +10.40% (FR +7.64, IT -5.57, LT +14.30), all +8.30% (FR +7.54, IT +7.97, LT +10.09, NO +7.41); early-consumer-fixed plant gain +18.72% and +12.82%, respectively. No claim of universal country transfer; the country sample is small and differs across coverage regimes.

The EQR association is not a sewage or nutrient effect. GBIF records and potential host catalogue are not observed feeding or plant biomass. No new plant trait outcomes were fitted.
