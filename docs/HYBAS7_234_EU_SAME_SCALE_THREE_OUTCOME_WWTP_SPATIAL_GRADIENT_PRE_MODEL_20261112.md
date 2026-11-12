# Same-watershed three-outcome responses to the ONE static modeled WWTP spatial gradient

2026-11-12. This is an exploratory POSTHOC macro analysis designed after the original EEA/GBIF/HydroWASTE source audit and after the n169 temporal EQR–S50–Lycaenidae co-change result was seen. Run regardless of signs.

## Original, fixed candidate source

234 original EEA river-monitoring HydroBASINS L7 basins observed in both 2018 and 2022; 468 original two-period source rows, 455 original rows support genuine rarefied plants S50 and >=20 butterfly-family records. Candidate source SHA256 07fa9744cd4e9574587a7b25c4ed4fa7f3183e377cac61c28723a4f24a86c875. Original EU2018 EQR Git blob 0f30bfbc3dac13f6a182c8a910e9899c04778322, EU2022 325f199e985d9862ffee5b0edb205c61a7be3833, H7 source monitor assignment SHA256 a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89.

A priori numeric gate: at least 320 complete basin-period rows, 170 unique H7 basins, 12 source countries. If it fails, stop; no alternative favorable scale/year/basin threshold. EQR per basin-period is the original equal-physical-station average within the appropriate original biological source year. This cross-sectional design explicitly DOES NOT assert station and method repeated over time.

## Identical basin-period unit, one source predictor, three distinct responses

Predictor W: log1p(original static HydroWASTE upstream-to-H7-basin-outlet modelled wastewater discharge divided by upstream catchment area plus 1). Static, mostly estimated, many WWTP operational statuses unknown, NOT actual annual measured TN/P, NOT delivered pressure to an interior monitoring station.

Responses:
1. Y_EQR = site-equal mean of genuine EU2018 numeric_EQR_2018 or genuine EU2022 normalised_EQR (method effects must be adjusted);
2. Y_Plant = S50 expected indexed original accepted speciesKey richness at 50 original Plant species-indexed GBIF records (NOT field flora);
3. Y_Lycae = log((original Lycaenidae indexed count + 0.5)/(sum of exact five butterfly families count minus Lycaenidae plus 0.5)) (NOT insect abundance, visitation or pollination).

Same fixed baseline across responses: country one-hot, source EQR majority assessment method (ties alphabetical), source year-bin, log source HYBAS7 local area, log1p original station count, log1p original Animalia, Plantae and all-five-butterfly GBIF indexed count opportunities. Original monitoring site river DIS_AV_CMS can be used only if complete for all 3 outcomes; otherwise report missingness and omit it in ALL models, no per-outcome rescue. Some baseline recording counts overlap a response's own taxon numerator, making observational adjustment potentially endogenous: no causal language.

B = identical baseline source predictors; BW = baseline plus W. Train-only feature standardization, Ridge alpha 10 fixed, deterministic 5fold H7-ID-hash heldout (both original periods from a basin stay together). Report signed percent heldout MSE gains, and standardized outcome beta per 1SD W with H7-clustered descriptive 95% from OLS. Bootstrap out-of-fold MSE contributions over whole basins 2000 times. Do not describe basin hash folds as independent geographic-region holdout. Report raw Spearman and optional FR versus non-FR/leave-France-out diagnostics; no small subgroup inference.

The broad scientific test is whether the *same* original modeled wastewater source geography adds different predictive information for benthic EQR, GBIF source-rarefied plants and Lycaenidae record composition. These are separate response families. An apparent sign discrepancy is not a verified water-to-land trophic mechanism, and a null predictive increment is a valid result. The entire analysis is POSTHOC and missing original agricultural/urban landcover, station-specific delivered dose, actual plant hosts, pollination and emergence. Do not merge original source/cohort branches, main or chocho; no eBird.
