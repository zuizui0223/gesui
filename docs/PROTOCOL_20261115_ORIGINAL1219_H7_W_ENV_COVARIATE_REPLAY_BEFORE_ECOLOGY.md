# Source-only freeze: original 1219 H7 HydroATLAS and HydroWASTE predictors for new source-qualified H7

**2026-11-15. This is a post-result geographic expansion feasibility design, not a new wastewater experiment or preregistration. The original 454 H7×period/232 H7/90 H6 source-only and ecological outcomes are not changed.**

## Factual opportunity
Previous source-only 30 additional H7 ×2 periods ×6 GBIF taxa true speciesKey facets are now **360/360 acquired**, after selectively retrying an originally missing Papilionidae facet (true 0 indexed records, not an API-missing zero). Exactly **29/30** new H7 meet **both periods** ≥10 numeric EIV plant taxa, ≥40% original all-Plantae speciesKey numeric-EIV coverage, ≥2 numerical HNS butterfly species, ≥10 scored butterfly indexed records and ≥70% of all five-family raw records scored; all **10 additional genuine singleton H6** are represented. This is source/trait coverage only. Source authenticity preserved by prior runs #43725053634 and #43725154865.

## What must be reproduced before opening new ecological regression
Old 454 source H7-period CSV SHA256 `1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca`. Extract only original geospatial predictors `W, area, crop, urban, naturalout` to verify source value identities, **not** original butterfly/plant outcome signs.

Original 1219 H7 native HydroATLAS `HYDROATLAS_V1_L07_1219_H7_BASIN_RAW_LANDUSE_ATTRIBUTES_SOURCE_ONLY.json` with raw SHA256 `ce6780e01d2b98ed8e852574adb827cb3125be792a164e8be9e2ebb917b3a686`, source run #43225986249. Source-backed fields `crp_pc_use`, `urb_pc_use`, `dis_m3_pyr`. Exact original 454 definitions:
- `crop = crp_pc_use`
- `urban = urb_pc_use`
- `naturalout = log1p(dis_m3_pyr)`

Original 1219 H7×2-period HydroWASTE/GBIF master `HYBAS7_1219x2_GBIF_PLANTAE_LYCAENIDAE_ANIMALIA_HYDROWASTE_SOURCE_ONLY_MASTER.csv` from source run #43202839926. Source field `HYDROWASTE_STATIC_UPSTREAM_TO_BASIN_OUTLET_MODELED_WASTE_DIS_M3_PER_DAY`, `UPSTREAM_TOTAL_AREA_KM2` and `LOCAL_SUBBASIN_AREA_KM2`. Exact previously published source formula in `scripts/gesui_HYBAS7_234EU_EEA_EQR_Plantae_S50_Lycaenidae_static_WWTP_three_response.py`:
- `W=log1p(max(0, static_outlet_modelled_waste_dis_m3_day)/(up_area_km2+1))`
- `area=log(local_H7_area_km2)`.

### Strict source replay
1. Read original 1219 native source rows (2,438 unique H7×period, 1,219 H7), original published HydroATLAS covariate map (1,219 H7). No new local observer data or GBIF queries.
2. For **ALL original 232 H7**, compare reconstructed W/area/crop/urban/naturalout to original 454 rows at original source `atol <=1e-10`, and require those covariates remain exactly identical across the two original periods for each H7. Fail all downstream work if any mismatch.
3. Extend **only** to the 29 newly sourced H7 in 10 original singleton H6, using original fixed official H7→H6 source. Check all 29 H7 are in both 1219 sources, valid finite predictors, bounds crop/urban [0,100], modeled pressure≥0, areas>0, W static over time; no NUTRI or HNS response fits.
4. Inspect available country/EEA station/source-visit and original GBIF recording opportunity fields for all 29 but do **not** assert 2018/2022 measured same-station EQR in single-era EEA sites. Some original original (nEEA), `riverflow` and original `country×method` were conditional on EQR availability and **cannot be blindly carried over**. Separate terrestrial-only comparison shall be defined explicitly as a different model if compatible covariates exist.
5. On reproducing 232/232, freeze source-only one H7 per new singleton H6 from the original immutable ranks in #43723429521 and QA selecting only those passing both-period source EIV/HNS quality, not response size. Minimum 10 H6, 40 total multi-H7 H6 remains. Never lower gate.
6. If original predictors disagree with old source at any H7, STOP and report mismatches, do not tune source formula to get favorable outcome. Inference is posthoc observational, not actual delivered effluent N, measured leaf chemistry or ecological food-web link.

**Not yet measured:** local outfall chemical dose; species populations, larval feeding, biomass, ecosystem service. Leave all earlier Draft PRs open/unmerged. Main and chocho unchanged.
