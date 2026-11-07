# Source-predeclared European river TP/TN–wastewater observational macro pilot (2026-11-07)

## Strict hierarchy of claims

1. **Source geometry:** HydroWASTE WASTE_DIS (m³/day) and HydroRIVERS directional NEXT_DOWN provide *modeled wastewater-water-volume exposure*, never N/P load. Test source coverage without reading measured concentrations.
2. **Observed river chemistry:** Peters et al. Figshare v5 `TN.2_-_Harmonised.zip` and `TP.2_-_Harmonised.zip` are the *original measured-concentration station series* (unit mg/L). Do **not** use `TN.4_estimated` / `TP.4_estimated`, WRTDS modeled series, or source `C_mean/C_median` predictive target columns.
3. **Ecological link (separate):** Existing 48 GBIF Plantae sites are a different, fixed sampling design. Do **not** merge source-selected new nutrient sites with old plant outcomes or call measured nutrient concentration an effluent-specific causal effect.

## Original source, covariate and sample lock (before fitting outcomes)

- HydroWASTE original 58,502 facilities, blob `0506f873f6541c5f324478894a9731290d722bce`. Exclude explicit nonoperational/planned statuses and accept model flow/status unknown with quality ≤2 facility location and ≤3 wastewater flow quality; only graph-linked positive wastewater volumes. Report operational-status-known sensitivity.
- HydroRIVERS Europe v1-derived GeoParquet SHA256 `5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f`, not yet independently byte/topology verified against official release.
- Station geography inner lat45–63°, lon−11–25°; river source outer lat44–64°, lon−12–26°. Station needs ≤500 m projected offset (EPSG:3035) to river line; original Peters v5 metadata `WQ.available=1`, harmonized station file present, 2011–2024 metadata/year-index overlap, original total count ≥12 and HydroRIVERS mean Q ≥0.5 m³/s.
- Source response time window frozen to **2011-01-01 through 2024-12-31**. Later outcome gate requires at least **12 valid positive observed concentrations in ≥2 separate years** per station; use **median of measured observations** (mg/L). Log median is the analysis target; report exact n and years. No source metadata value substitution.
- For each station river reach `r`, primary exposure is `log1p(total modeled WWTP effluent m³/day from facilities on hydrologically upstream reaches with ≤50 km directed graph path to r)`. **Same-reach facilities excluded** because discharge point order within reach is uncertain. Report exposure as continuous, plus a derived *potential water fraction* relative to HydroRIVERS mean river Q. Zero source exposure means *no recorded HydroWASTE upstream source*, not genuine sewage absence.
- Negative control is `log1p(total modeled WWTP effluent volume on reaches downstream of monitoring station within 50 km)` (not a physically delivered upstream dose); compare only on the same complete-case station sample.
- Outcome-independent baseline covariates: geodetic lat/lon, river mean log Q, upstream drainage area; all scales fitted within training. Secondary richer sample adds `crop2` (cropland cover) and `pop_dens` (population density), plus `N_deposition` when present, from public `TN.4_land.csv` / `TP.4_land.csv`; this is only a **source-selected model-input subset**, not all worldwide stations. Column names `C_mean, C_median, L, Y` never enter models or source feature tables; they describe outcome/model-derived properties and must not leak target information.
- One monitoring station = one row; never use daily rows as independent sites. Predefined spatial block is floor(lat/3) × floor(lon/3), all stations in a block always in same evaluation fold. Need at least 5 spatial blocks and ≥40 stations/parameter for an out-of-region prediction contrast; otherwise report insufficient source support.

## Fixed analyses when source gate and measured-value gate pass

All predictors are fixed from source metadata before inspecting measured TN/TP values.

1. **Primary descriptive** Spearman rho between log observed median nutrient concentration and log1p upstream WWTP source volume. Report n, number of geographic blocks, discharge-source-status composition.
2. **Fixed five-fold grouped heldout prediction**, using `GroupKFold(n_splits=5)` (not random site CV). Standardize predictors only inside training folds and fit `Ridge(alpha=10.0)` (no outcome-driven hyperparameter search). Primary baseline M0 = geographic lat/lon + log(1+river mean Q) + log(1+drainage area); M1=M0+upstream WWTP modeled volume; M2=M0+downstream WWTP placebo volume; M3=M0+both. Report identical withheld station set, MAE/MSE for all models, per-block residuals and directional increment signs. No manufactured significance from model coefficient p-values.
3. **Confounder-gated sensitivity**, only if independent cropland+population covariates are actually matched for ≥40 source sites across ≥5 spatial blocks: M0c=M0+cropland cover+log1p(pop_density) (+deposition when present), M1c=M0c+upstream WWTP volume; same GroupKFold and identical denominator. This contrasts urban/agriculture with wastewater-volume infrastructure, **not** chemical treatment effects.
4. **Source reporting quality sensitivity**, if ≥40 stations: replace all-status modeled upstream volume with volume restricted to officially `Operational` facilities only. Not Reported statuses may dominate; do not call them operational.
5. Never claim a physical upstream signal if only downstream placebo improves, or if source support/land-use sensitivity fails. A positive prediction increment without independent nutrient exposure attribution is **an observational association**; the outcome is river chemistry, not verified nutrient loading attributable to WWTP.
6. Keep TP and TN separate. Any new plant linkage must have independent river-site overlap and adequate true upstream/downstream controls. No dataset tuning based on the original 16-outfall null.

## Risk and audit

Because the network graph is truncated to a European bbox, report station distance from domain boundary and incomplete upstream paths. River crossing a national border and wastewater catalog coverage may create spatial bias. HydroRIVERS is third-party transformed rather than independently reproduced from official GIS. HydroWASTE historical operations are retrospective proxies. Original Peters land-use covariate set was assembled for published models: restrict to nonoutcome columns and disclose its selection. Report results as observed pattern only, not sewage-caused river pollution.

**Primary publication gate:** independent measured WWTP N/P **mass flux**, matching river hydrology, urban/cropland controls and robust spatial validation needed before any strong causal manuscript claim.
