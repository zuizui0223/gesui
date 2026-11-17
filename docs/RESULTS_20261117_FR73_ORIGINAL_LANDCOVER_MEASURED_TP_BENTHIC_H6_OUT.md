# 2026-11-17 | FR73 original native H7: land-use → measured river TP and benthic EQR

**New numerical descriptive reanalysis** on previously known original FR73 native H7 × same physical 263 TP/I2M2 field stations in 2018 and 2022. This is **posthoc exploratory**; it is neither new field samples nor proof of sewage river→riparian nutrient transfer.

## Original exact-source data

- Original extracted FR73 CSV from original GitHub Actions #43373639391: FR73_HYBAS7_NATIVE_2018_2022_REAL_TP_I2M2_EQR_AND_GBIF_INDEX_DELTA.csv
- Source CSV bytes SHA256: **f2ff4e56f245e1aa8f84f5ac099d68ec1c7e27a786c1067439db388b186dd90b**.
- n=73 unique native H7, 24 genuine native H6 parent groups. Same exact cohort in 2018 and 2022, not two independent studies.
- Pre-existing HydroATLAS upstream cropland% and urban% are catchment land-use **covariates**. Do not interpret as particular local riverbank vegetation or effluent exposure.

## All four previously fixed descriptive rank correlations per year

| Original source columns paired at H7 | 2018 Spearman rho | 2022 Spearman rho |
|---|---:|---:|
| upstream cropland vs annual measured river logTP | +0.409 | +0.411 |
| upstream urban vs annual measured river logTP | +0.397 | +0.318 |
| upstream cropland vs genuine I2M2 benthic EQR | -0.458 | -0.383 |
| upstream urban vs genuine I2M2 benthic EQR | -0.416 | -0.424 |

## Same H6-parent-disjoint 5-fold Ridge predictive test

- Feature standardization fitted on training H7 only, Ridge alpha=10. All H7 within a native parent H6 group in the same heldout fold.
- Geomorphological baseline, unchanged: log H7 area, modeled outlet natural discharge, modeled mean measured site river Q, log(1+paired TP/I2M2 station count).
- Add land predictors **together**: original upstream cropland percentage and urban percentage.
- For EQR, fit geomorphological baseline (B0), plus cropland+urban (B1), plus annual measured river logTP (B2).
- MSE improvement = 100*(1 - added model H6-out mean squared error / baseline model H6-out mean squared error). True responses are raw logTP and raw EQR, not standardized.
- Descriptive 95 intervals use 3,000 native-H6-cluster resamples of fixed fold squared errors, **not** a confirmatory confidence interval or independent resampling of field studies.

| Contrast | 2018 improvement % [descriptive H6 95] | 2022 improvement % [descriptive H6 95] |
|---|---|---|
| Add cropland+urban to baseline → predict river logTP | +8.34 [−2.66,+18.63] | +7.05 [−6.64,+16.92] |
| Add cropland+urban to baseline → predict benthic EQR | **+23.97 [+9.52,+39.50]** | **+15.95 [+4.54,+25.92]** |
| Add measured logTP after cropland+urban to predict EQR | **+28.32 [+5.38,+40.01]** | **+22.94 [+7.79,+37.24]** |

Values are not causal land→river estimates, the two years are correlated same H7, and heldout MSE gains of different outcomes cannot be compared as common effect sizes. These source data were chosen in earlier exploratory work, so posthoc 95 ranges require restrained interpretation.

**Scientific meaning:** There is a reproducible *co-observed macro spatial land-use–aquatic condition gradient* on genuine shared H7×year supports. The positive land-use→EQR prediction gain and additional true TP signal imply that regional land use/flow can explain both water quality and biological condition without ever invoking water→bank transport. They are a strong **confounding-control motivation** for the upcoming real 20–100m vs 500–1000m riparian vegetation extraction. This is not a wastewater-specific effect or demonstrated land uptake.

## Linked new source-only extension (pending real imagery)
- Original 263 genuine field points, same-reach 262 HYRIV IDs: results/20261117_FR263_FIXED_ORIGINAL_EEA2022_STATION_COORDINATES_SOURCE_ONLY.csv
- Geography pre-raster QA: results/20261117_FR263_GEO_DISTANCE_PRE_RASTER_QA.json
- HydroRIVERS original geometry source audit and H6-balanced 2018/22 public S2-STAC metadata: analysis/20261117_FR263_riverland; corresponding workflows on branch analysis/gesui-20261117-onepaper-pathway-identification-decision.
- Status of source-only remote GIS/S2 new runs not yet verified as successful; original HR-VPP 2018/22 QFLAG TPROD raw pixel acquisition not completed, so no river→terrestrial production coefficient exists.

**Source/tool trace:** Reanalysis script generated from exact authentic original GitHub Actions #43373639391 source ZIP, CSV SHA verified and saved as fr73_land_aquatic_predictive_audit.py in current local artifact working directory. Repeated live model fits were evaluated on original existing 73 H7 observations; no new third-party raw datasets were used.
