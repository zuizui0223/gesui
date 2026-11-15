# Wastewater geography, aquatic EQR, plant S50 and butterfly HNS: same-H7 multiresponse diagnostics

**15 November 2026 — source-original, exploratory POSTHOC analysis.** Same 454 original European H7 basin × historical period source rows / 232 distinct H7 / 90 authentic H6 parent basins / 15 countries. These are not 454 independent basins. **NO causal structural equation model or ecological mediation is identified or fitted.**

## Why it was done

The proposed diagram links wastewater exposure, aquatic biota, riparian vegetation and a butterfly resource trait. Those arrows must be separated: (a) the same basin-period may exhibit correlated, different responses; (b) observed mediating nutrient transfer, insect emergence, adult feeding and larval host use are NOT available in this data. This diagnostic tests *co-response and additional prediction at the identical hydrological grain*, while retaining all nulls and country/source-quality limitations.

## Frozen material and data meaning

- Re-used original `results/20261115_sync/hns_macro/H7_454_HNS_panel.csv`, SHA256 `1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca` (full 454 source rows). Zero duplicate `(h,period)` keys; 232 H7; 90 H6; 222 H7 measured in both historical periods.
- **W**: prior HydroWASTE **modelled, static basin-outlet source geography**, NOT observed N/P/Na/dissolved-organic-water dose or actual effluent fraction reaching a butterfly patch. No source W changes in the two periods.
- **EQR**: source-reported method-specific aquatic benthic ecological-quality ratio, not aquatic biomass, measured emergence, adult prey production or species diversity.
- **Plant S50**: expected distinct GBIF-indexed plant speciesKey at 50 source occurrence records, NOT census abundance or host abundance at a given riparian patch.
- **Butterfly HNS**: record-weighted average habitat N affinity of historically known larval host plants among numerically scoreable butterfly species (186 accepted keys), NOT plant tissue N, local larval feeding, adult puddling, reproduction, true individual butterfly abundance or original WallisDeVries community N index.
- Controls: source recording effort (EEA, Animalia, Plantae, five butterflies), river flow, natural outlet flow, area, upstream crop and urban cover, butterfly source coverage and catalogue completeness, four independent family shares for HNS outcomes. Country and assessment method are represented by **one joint country×method stratum**, to avoid the exact four-column collinearity detected in an initial rejected exploratory fit. Both periods of an H7 are always in the same H6-heldout fold.
- Method: fixed five H6-disjoint folds, training-only numeric standardization and category vocabulary, group-H6-equal Ridge alpha 10 prediction, H6-equal MSE; 2,000 fixed heldout H6 error resampling yields *descriptive* 95% uncertainty. Observational H6-equal WLS cluster-robust descriptive95 partial slopes are separate from OOS predictiveness. No causal identification from either.

## Original same-grain results (all 454 rows)

| Response / predictor added last | H6 equal out-of-parent MSE gain % [descriptive95] | Within-source adjusted coefficient per +1 SD predictor [H6 cluster descriptive95] |
|---|---:|---:|
| Aquatic EQR, add W after controls | **+6.725 [−2.453, +13.968]** | **−0.0750 [−0.1133, −0.0367]** |
| Plant GBIF S50, add W | **−1.358 [−3.647, +1.630]** | −0.1702 [−0.7135, +0.3732] |
| Butterfly HNS, add W | **+10.014 [−3.446, +20.979]** | **+0.1448 [+0.0650, +0.2247]** |
| Butterfly HNS, add plant S50 *after W* | +4.306 [−2.242, +10.739] | −0.0647 [−0.1035, −0.0260] |
| Butterfly HNS, add EQR *after W* | +5.079 [−0.228, +10.524] | −0.0538 [−0.0920, −0.0157] |
| Butterfly HNS, add S50+EQR jointly after W | +7.941 [−0.811, +16.318] | S50 **−0.0612 [−0.0991,−0.0232]**; EQR **−0.0473 [−0.0836,−0.0110]** |
| Butterfly HNS, add W after S50+EQR | +5.697 [−9.032, +17.002] | +0.1174 [+0.0379,+0.1969] |
| Plant S50, add EQR after W | +0.482 [−1.952, +2.608] | +0.2367 [−0.0897,+0.5630] |

**Important:** MSE gains are predictive information, not sign or percentage change in organisms or ecological services; adjusted coefficients are descriptive conditional associations, not causal path coefficients. With all 454 rows, the descriptive predictive intervals for ALL EIGHT added-value tests include zero; do not choose one favored +10% model as confirmation. The sign of S50/HNS and EQR/HNS conditional slopes is **negative**, contradicting a simple always-positive bottom-up pathway. S50 has no clearly stable W added-value at this scale.

## Existing France and quality gates retained as post-result checks

| Cohort | n H7×period / H7 / H6 | W→EQR H6 MSE gain | W→HNS H6 MSE gain | W→plant S50 H6 MSE gain |
|---|---|---:|---:|---:|
| France only | 299 / 152 / 36 | +0.55% [−7.55,+8.20] | −5.41% [−10.88,+0.72] | −0.75% [−8.50,+3.91] |
| Outside France | 155 / 82 / 56 | −1.33% [−14.41,+8.83] | +6.07% [−14.67,+20.35] | −1.27% [−4.11,+1.91] |
| Six-group source species coverage ≥80% | 390 / 217 / 83 | +5.63% [−3.40,+13.13] | **+16.91% [+9.00,+26.58]** | −1.57% [−3.28,+1.12] |

The source quality subset and France split were **examined after observing the original spatial HNS signal**; the positive source80 HNS gain is selection-sensitive and does not erase full-cohort/country uncertainty. Country×method stratum and W/urban correlation cannot remove all confounding.

## Within-H7 temporal cross-layer cochange

- Two recorded periods for **222 distinct H7** in 89 H6; only **172 H7** have the same reported *majority* EQR assessment method. Even these 172 are NOT all confirmed identical original physical river monitoring stations. The other 50 have a changed majority/unknown method and must not be labelled matched biomonitoring sites.
- Across all 222, Spearman `ΔEQR↔ΔHNS` **ρ≈−0.00003 [−0.124,+0.147]**, `ΔPlant S50↔ΔHNS` **ρ≈+0.061 [−0.082,+0.201]**. The brackets are H6-resampled descriptive95, not experimental effects.
- Same-majority-method 172: `ΔEQR↔ΔHNS` **ρ≈−0.032 [−0.166,+0.118]**; `ΔPlant S50↔ΔHNS` **ρ≈+0.130 [−0.051,+0.289]**. None establishes a recurrent temporal cross-layer association.
- W is identical in each period. **No treated effluent dose change was observed**, thus no wastewater-mediated change or difference-in-differences was estimated.

## Independent numerical audit and falsifiers

- Initial naive separate `country + method` fixed effects produced rank 30/34 (4 perfectly collinear dummies) because reporting method and country were coupled. It was rejected and replaced by the full-rank `country×method` stratum; all eight final models rank equals column count.
- Independent sklearn StandardScaler/OneHotEncoder + Ridge recomputed four full source H6 cross-validation results without invoking the original model matrix or Ridge functions. Maximum absolute MSE gain discrepancy **<1.4×10⁻¹² percentage points**; this verifies numerical implementation on **the same ecological dataset**, not independent field replication.
- Five original source/model grain unit tests pass. Reproduction gate requires one original source SHA256, unique 454 H7×period, complete fixed source values, 90 H6 parent groups, 222 repeated H7, 172 same-majority-method pairs, H6-purged folds, eight full-rank fitted models, 12 posthoc cohort outputs, and `causal_SEM_or_mediation_fitted: false`.
- Multiple post-result models and selected source-quality cuts mean these intervals are exploratory and have NOT been adjusted for multiple comparisons.

## Scientific conclusion / required next evidence

A *source-geographic mosaic* — lower aquatic EQR together with higher record-weighted butterfly host N affinity at some locations — is consistent with the current data, but **not identifiable as a sewage-induced water→plant→butterfly trophic cascade**. No stable change–change coupling was found. The next field-informed SEM requires verified discharge source and dilution, measured nutrient/Na/DON/DO at the same reach/time, local host plant biomass/leaf chemistry, larval feeding and abundance, aquatic adult emergence dry mass, riparian predator diet/consumption, adult butterfly species×sex×puddling substrate/effort, matched naturally wet/stormwater controls, independent discharge systems, season and detection models.

This is a new draft analysis stacked on existing unmerged PR #87. It must not be auto-merged into `main` or `chocho`. No eBird, no TinyFish.
