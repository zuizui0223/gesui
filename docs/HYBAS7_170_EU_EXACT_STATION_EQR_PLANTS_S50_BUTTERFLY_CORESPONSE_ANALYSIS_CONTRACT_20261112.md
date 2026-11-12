# Source-gated co-response model: EQR × plant S50 × five-family Lycaenidae composition across repeated EU L7 basins

**2026-11-12, first ecological outcome analysis on SOURCE-FROZEN new plant S50 and 5-family butterfly endpoints.** This is explicitly **POSTHOC / exploratory**, since the 2018/22 EQR source and prior GBIF source-count aggregates were known before this model was proposed.

## Fixed source-reconciled cohort and permission

Previous strictly source-only action [43212724530](https://github.com/zuizui0223/gesui/actions/runs/43212724530) passed originally frozen source gate without reading numeric EQR: 234 H7 basins with measured EQR in 2018 and 2022 (not selected by ecological signs), 468 native basin-periods, 455 period rows with plant 50-rarefied S50 and ≥20 original five-family butterfly records, 223 basins with both valid periods, **170 basins with BOTH periods and original exact EEA station/assessment-method pairs**. Exactly **1,047 matched native river stations** span only **174 H7 basins and 7 European countries**, with **1,004 /1,047 stations in France**. Thus no general pan-European ecosystem-service or nutrient causal estimate can be claimed regardless of p value.

All new model data MUST be those original 170 common source-supported H7 basin pairs (no choosing based on effect sign). Original candidate cross-source SHA256 `07fa9744cd4e9574587a7b25c4ed4fa7f3183e377cac61c28723a4f24a86c875`, source exact-station/method SHA256 `373a3a96ea790a7592cc4d10b0e4f9244fbdddab2c1af4901118792674874930`. Original EU2018 EQR source Git blob `0f30bfbc3dac13f6a182c8a910e9899c04778322`, 2022 `325f199e985d9862ffee5b0edb205c61a7be3833`. Check exact original station IDs and country and assessment-method codes, never cross-join different kinds of station or different country/year.

## Outcomes and predictors defined before opening newly joined actual EQR values

**Biological response `ΔEQR_h7`**: mean over original same physical station and SAME original assessment method `normalised_EQR_2022 − numeric_EQR_2018` within each H7 basin. Count each original matched station once, then each H7 basin once in regression. Do not average different 2018/2022 station inventories. Flag impossible/nonfinite EQR, never impute. Do not rename EQR as species richness, butterfly counts or measured ecosystem service.

**Primary shared-basin co-responses**: (1) `ΔS50 = plant speciesKey rarefied expected richness at 50 indexed plant records, 2020–22 − 2017–19`; (2) `ΔLycaenidae_logodds= log((L2022+0.5)/(five_butterfly2022−L2022+0.5)) − log((L2018+0.5)/(five_butterfly2018−L2018+0.5))`. Five butterfly families are Lycaenidae, Nymphalidae, Pieridae, Papilionidae and Hesperiidae. Their index observations are **NOT butterfly field abundance or actual pollination**.

**Recording opportunity / confound diagnostics fixed for every model**: `Δlog1p(Animalia GBIF counts)`, `Δlog1p(Plantae GBIF counts)`, `Δlog1p(five-butterfly-family count)`, original matched EQR station count log1p, log L7 area, original static `HydroWASTE_upstream_to_BASIN_OUTLET_modelled_flow /(UP_AREA+1)` log1p. These are descriptive covariates, not verified discharge at station. Original country (majority of matched native station source IDs, deterministic tied alphabetically) fixed effects.

**MODEL B** predicts `ΔEQR` from recording opportunity + country + area + flow-model pressure + matched site count. **B+P** adds `ΔS50`; **B+L** adds `ΔLycaenidae_logodds`; **B+P+L** adds both. Training-only standardized predictor columns, original basin hash 5fold, same 170 source rows, Ridge alpha=10 fixed for all models; equal basin weights within each heldout fold. Report ALL four heldout MSEs and signed percentage gain relative to B, no rescue by changing alpha, years, family/count denominator, H7 scale or sample thresholds. Source light/dark timing not yet measured.

Partial OLS with H7-sample HC3 descriptive intervals and 2,000 basin-level bootstrap resamples for both co-response covariates; identical source conditioning. Country-balanced sensitivity weighting and France-only vs non-France subgroup reporting without presenting n<20 groups as generalizable results. Spearman ΔEQR vs plant ΔS50 and vs butterfly Δshare and Δplant vs Δbutterfly to show whether an actual multi-layer co-response exists; no p for 170 repeated station rows treated as independent.

**Source first hard checks**: 170 source basins with two valid rows, 7 exact-station countries; EQR numeric finite and both years in same assessment method, and all required recorded-effort terms. If EQR/GBIF missing shrinks to <120 basins or ≤5 country groups, STOP and do not hunt favorable subgroup signs. Source static wastewater flow is repeated between periods and therefore cannot explain a **change** as a time-varying intervention.

**Interpretation**: a correlation between two basin changes is not WWTP→benthos→plant→butterfly mediation. The times are coarse observational windows and GBIF may have temporal recording biases. Even a positive predictive gain indicates shared geography/observations, not measured pollination/biomass, and an apparent negative/positive effect may not generalize beyond France.

No eBird, no other original ecosystem data replaced, do not merge to `gesui/main` or `chocho` or overwrite old PRs. Save new original-source report and caveats.
