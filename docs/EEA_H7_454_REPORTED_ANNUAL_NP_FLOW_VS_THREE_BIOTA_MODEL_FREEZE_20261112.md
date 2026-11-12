# Year-specific reported WWTP N/P as an additional environmental predictor: identical 454 H7 periods

This analysis was FIXED on 2026-11-12 after the original reported-NP support-only gate [Actions 43222864099](https://github.com/zuizui0223/gesui/actions/runs/43222864099) passed, but *before* re-fitting the new annual N/P source models. It is **posthoc**, since the previously registered static HydroWASTE model and associated EQR sign were already seen.

## Physical and observation units

Exact same 454 H7 basin-periods (232 distinct subbasins) in PR #65; 2017–19 EEA2018 n227 and 2020–22 EEA2022 n227. Same 454 original ecological response rows and same original GBIF plant S50, five-family Lycaenidae family log odds, EEA EQR values. Fixed original source-only annual measured WWTP mass-to-station and river-model flow normalized manifest SHA256 `012ad79feabc8c11e0d74f803730b5b7913b96a64a6af26fb0c97263d9aa7c26`. Original source positive **56/454 basin-periods**, **38 H7 basins**, including **28 early and 28 late** and **13 countries**, **132/3525 stations** with positive matched reported upstream N/P. Other 398 are zero in matched report index, NOT known WWTP-free sites.

## Frozen source exposure unit, transform and model

Reported source P intensity: original station-equal average of `P_source_reported_t/year / DIS_AV_CMS` times `1e6 /(365.25 × 86400)` mg/L-equivalent, transformed `log1p(1000 × P_mgL_equivalent)` (µg/L-equivalent numerical transform, not actual observed TP). Reported N same but `log1p(N_mgL_equivalent)` (mg/L-equivalent). Both same time bin, zero only denotes no strictly matched reported source, not real absent discharge. Separate source presence flag `1[P>0 OR N>0]`. The conversion assumes full delivery/no retention and historical river DIS_AV_CMS is modeled flow, not true dilution at that sampling date.

Common original 3 outcomes:
1. EEA EQR actual country/method assessment;
2. GBIF plant speciesKey 50-record rarefied S50;
3. Lycaenidae log odds among GBIF five butterfly-family indexed counts.

**Baseline B** is the exact original PR #65 common country, EEA BQE method, year-bin, basin area, original river flow, native original EQR station count, and original Animalia/Plantae/five-family GBIF recording-opportunity covariates PLUS existing static HydroWASTE estimated outfall flow per original upstream area. **New source models:** B + `reported_source_detected`; B + `log1p(P report-flux equivalent µg/L)`; B + `log1p(N report-flux equivalent mg/L)`; B + both log1p P and N. N and P correlated; neither explains independent N/P causality, and no initial target sign is chosen.

All models: identical original 454 rows, 232 basin IDs, hashed 5fold (same H7 two periods share fold) ridge alpha10 with training-only feature standardization; outcome-specific signed MSE gains relative to identical B; 2000 basin bootstrap descriptive 95% gains. Source-only ordinal P/N report indicator is a metadata/underlying emissions mixture, cannot be interpreted as treatment/no-treatment.

Report adjusted OLS H7 clustered descriptive 95 intervals for NEW features under country/method/year controls, but no warranted causal N/P separation. Do original country-equal total weighted sensitivity and early/late diagnostic if support exists, without optimizing model scale or posthoc family subsets. **France dominance may invert evidence**; never hide negative/null results. Source-specific annual upstream reported P/N remains only in a minority of H7 periods.

**Hard anti-claim**: There is no real measured on-site river TP/TN, no exact nutrient delivery/travel/retention modeling or operation-status QA per facility, and no measured adult emergence/flower/host/seed set in this European macro panel. Even if EQR predicts but terrestrial GBIF outcomes do not, this suggests a difference in INDEX sensitivity only, not demonstrated ecosystem-service decoupling.

No merging, no eBird, preserve all previous PRs and source snapshots, report outcome signs honestly.
