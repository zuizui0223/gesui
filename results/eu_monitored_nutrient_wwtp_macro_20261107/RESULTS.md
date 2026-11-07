# Actual monitored river TP/TN vs directed wastewater source — European pilot

**Purpose:** assess whether upstream (not downstream) source-pinned WWTP treated-water volume helps explain independently measured river TN/TP concentration, without confusing WWTP volume with actual nutrient load. Original real-source runs: https://github.com/zuizui0223/gesui/actions/runs/42371629424 and https://github.com/zuizui0223/gesui/actions/runs/42372235775 . New independent zero-network source/metric audit in INDEPENDENT_AUDIT.json.

## Actual source support
- Original HydroWASTE **58,502 facilities** re-verified by Git blob; 6,354 location/flow-quality and graph-compatible source facilities in European network study domain.
- Original Peters 2026 v5 observed TN/TP records; each station measured concentration response is 2011–2024 median with at least 12 positive observations in at least 2 distinct years. No source estimated TN/TP concentrations or response-derived covariates used as predictors.
- **TP: 1,418** eligible stations (28 geographic 3-degree blocks), **85** with at least one upstream WWTP source within 50 km, 510 stations with independent crop-cover and population-density predictors, **73** of those with upstream WWTP support.
- **TN: 236** eligible stations, only **13** with upstream registered WWTPs; land-cover/pop/deposition-source subset is 20 sites with **one exposed**, therefore low exposure support, not robust nutrient-generalization evidence.

## Frozen source-only spatial-block CV
Fixed Ridge(alpha=10) with 5 folds held out by 3-degree geographic blocks, same site denominator throughout; no spatial leave-one-site-out leakage. Response = log observed station-median nutrient concentration in mg/L; baseline includes latitude, longitude, river mean flow and upstream drainage area. Source predictors: log(1+50km hydrologically upstream WWTP model discharge), and as negative control, modeled WWTP discharge on downstream river reaches.

- TP full 1,418: Spearman raw rank TP versus upstream source **rho −0.00192 (p 0.942, descriptive)**; baseline MSE 0.68881, with upstream source **0.71635 (4.00% worse)**, downstream placebo **0.67181 (2.47% better)**.
- TP 510 crop/pop-complete: baseline MSE **0.50821**, with upstream WWTP **0.51374 (1.09% worse)**; downstream placebo **0.50765 (0.11% better)**. This schema repair uses exactly the two source-published land-use covariates for TP; **TP.4_land contains no N_deposition**, so requiring it would have incorrectly discarded all 510 available sites. Previous zero-subset source result was a **code/data-schema mismatch**, not evidence.
- TN full 236: baseline MSE 0.69413, plus upstream 0.69873 (0.66% worse). Only 13 exposed and the source rich subset contains just one exposed site; TN cannot support a meaningful land-adjusted hypothesis test.

**Interpretation:** This macro test finds no reproducible positive predictive effect of the current simplistic upstream WWTP-water-volume indicator. The source graph and source-status proxies do not quantify delivered TN/TP mass, treatment efficacy, non-WWTP sanitation, agriculture, urbanization or seasonal loads. Downstream-infrastructure placebo outperforming upstream in the large TP panel is a serious sign that spatial confounding or source-model mismatch remains. Neither this negative predictive result nor a local positive case disproves all wastewater effects. Original chocho positive butterfly-host resource and port trends remain independently valid.

## Provenance and limitations
All underlying raw site series/coordinates and site-specific source features remain available in original 30-day GitHub Actions research artifacts and pinned public data DOIs; this public report retains only aggregated scientific outputs and file-level hashes to avoid unnecessary redistribution of underlying monitoring observations.
