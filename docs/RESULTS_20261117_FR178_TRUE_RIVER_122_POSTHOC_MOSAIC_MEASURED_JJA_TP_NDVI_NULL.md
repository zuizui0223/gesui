# 2026-11-17 | FR178 physical riverline × original measured TP × 2018→22 cloud-mosaicked riparian NDVI: independently verified results

## Scientific status and scope

**Real original river-water + terrestrial optical observations at the SAME source-authenticated field sites and matched calendar months are now joined and numerically fitted, but they do NOT support a macro-scale river TP → near-versus-far bank NDVI response.**

This is a *posthoc technical sensitivity*, not a preregistered independent confirmatory analysis. Preserve the original **single-scene 116/178 below >=120-site SOURCE STOP** permanently. Technical mosaicking was introduced AFTER that original source failure; the new source support does not retroactively turn the initial test into a success.

### Truly shared observational unit and method

- Fixed source universe: original **178 genuine river monitoring stations /58 H7 native catchments /23 H6 parents**, chosen using authentic original 263-EEA field station ID, fixed HydroRIVERS v1 native `HYRIV_ID` source reach, source station-to-river <=250m and original 900m available riverline, and truly observed river TP summer monitoring dates, all **prior** to satellite original pixel NDVI outcomes. HydroRIVERS geometry is from SHA-pinned **third-party version-1 GeoParquet**, NOT the original official ZIP (original ZIP fetch failed HTTP403).
- Two physical freshwater years **2018 and 2022**, same original individual station/WQ ID and matching June/July/August calendar months. TP from original truly measured daily-date harmonised Peters `TP.2` raw real monitoring CSV (not 2022 reported WWTP N/P plant mass).
- Fixed 10m Sentinel-2 L2A B04+B08 radiometry with original scene-specific `scale, offset`, true SCL4/5 clear land pixel, 2021 ESA WorldCover fixed natural classes 10,20,30,90. Riverline segments ±500m along actual original reach; original lateral bands **20–100m, 100–250m, 250–500m, 500–1000m**, >=100 genuine clear/valid natural pixels in each of all four bands per year and month; >=2 original JJA same-calendar months with both water and land. NDVI is greenness, not standing biomass, foliar nutrient concentration, insects, subsidy, or service.
- New **posthoc** technique: at most 2 **already source-frozen** S2 scenes for each original month and year combined on a common 10m EPSG3035 lattice, arithmetic *mean of valid pixel NDVI* at each grid cell; never choose dates to maximize NDVI, never change distance bands, classes, thresholds, or selected original river stations based on outcome correlation.
- Original native H7 land/flow and 2018 optical baseline; aggregate 122 source-supported stations to equal H7 samples; H6 true parent-block bootstrap and original five H6-disjoint spatial folds for Ridge10 OOS analysis.
- Confirmatory main source STOP is unaffected; this is hypothesis-generating observational analysis AFTER the source technical rescue.

### Original source gates, true satellite results

| Source-quality step | n valid | H7 | H6 | decision |
|---|---:|---:|---:|---|
| Original pre-pixel source true reach + TP JJA monitoring opportunities | 178 | 58 | 23 | source-only |
| Pre-specified original per-month **single-scene** true 10m pixel condition | **116** | 42 | 22 | **STOP: 116 below 120** |
| New posthoc exactly-two-preselected-scene cloud-gap mosaic | **122** | **44** | **22** | posthoc source-only passes same 120/40/15 denominator thresholds |

All 116 old pass sites remain, **6 previously failed sites regained**, none lost. Valid same year/month source pixels for **66 sites with 2 matched months**, **56 sites with 3**. Remaining **56 sites still lack 2 valid months** (41 have 1, 15 have 0). Four original source shards all true no-input-duplication n=12/93/38/35 and posthoc source QA n=6/70/27/19, original 178 unique site keys confirmed. Source NDVI-only `Δ[(NDVI20–100 − NDVI500–1000)]` unweighted 122-site mean **+0.019315**, median +0.020559, 76 positive/46 negative. H6-equal mean +0.01584, 13 H6 positive/9 negative. **Do not treat sign or mean as an effect of river TP or wastewater.**

### Matched summer true TP versus riparian NDVI change: independent numeric verification

The model was computed by [GitHub Actions #44020568960 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/44020568960), original artifact ID **14301175909** (artifact ZIP SHA256 `555f7b63d0f4ce0a8f15804185d2cefeabe78ab7e021dd17874fbbde3af4e18d`; internal JSON SHA256 `d2a66c80527187ea6fcd7dbb6ec1cad8464740644ba17f4fbdd44bb6c1d44b2e`). Real original TP is from source locked Figshare item `70896288` original raw ZIP SHA256 `26f2af11f4972987ee6f7c3f497ad902820f145851f71919d4f45b589e1d809e`.

For each matched site/month, original `dTP=log(median true measured TP2022_m)−log(median TP2018_m)`; `dNDVI=[near−far]2022_m−[near−far]2018_m`. Average over only equal original same-month comparisons for each station (2 or 3); then aggregate equal stations within each H7, leaving 44 original H7 in 22 H6. No year/month imputation.

| Analysis | Real outcome |
|---|---|
| 2018→22 `dTP` versus `dNDVI`, H7 native ecological units | Spearman **ρ=−0.027625** |
| 95% descriptive bootstrap resampling native H6 parent blocks (n2000) | **[−0.435997,+0.315800]** |
| Leave-one-H6-out descriptive rho range | **[−0.11894,+0.06896]** |
| Original H6-disjoint 5fold Ridge10 baseline includes cropland, urban, flow, area, 2018 near−far NDVI; add true summer `dTP` to predict `dNDVI` | **H6-equal MSE gain −3.73059%**, fixed H6 error resampling descriptive95 **[−9.19257,+2.25054]%** |
| Site/H7-unit weighted alternate predictive gain | **−2.63391%** |

**Interpretation:** no detectable covariation at this particular macro spatial/temporal support, and no regional out-of-H6 predictive gain from added true TP. Interval covers moderate positive and negative associations; this is NOT evidence of physical absence of source→terrestrial transfer. Neither a universal null, nor proof that sewage has no benefits to riverbank ecology.

### A single specified aquatic biological sensitivity, not a new fishing exercise

Original same-physical-station I2M2 benthic EQR **2018→2022 numeric change** is present in independently signed original FR263 cohort; same 122 stations aggregated into same 44 H7:
- `dI2M2 EQR ↔ dNDVI` Spearman **ρ=+0.006202**, H6 descriptive95 **[−0.28197,+0.26104]**.
- `dTP ↔ dI2M2 EQR` Spearman **ρ=−0.017900** on this same 44 H7.
- EQR observation date is NOT same exact calendar month as the true TP+S2 remote-sensing comparisons; therefore EQR secondary is not a time-matched mechanistic trophic mediator and does not overturn existing original seasonal difference null.

Independent local reimplementation crosschecked **all 122 exact NDVI site values/month identifiers**, 44 H7 aggregation, original Spearman and both original fivefold MSE-gain estimates against the GitHub Action artifact. No p-hacking / no later site exclusion.

### Severe informative missingness remains (posthoc-pass122 vs fail56)

| Observational original full site annual summary | Pass122 median | Fail56 median |
|---|---:|---:|
| Measured TP 2018 mg/L | 0.1100 | 0.07575 |
| Measured TP 2022 mg/L | 0.1200 | 0.08175 |
| Original annual `ΔlogTP` | +0.07880 | −0.06974 |

Pass site concentration differs; 122 points not random 178-site sample, which further limits any inferred continental or wastewater policy effect. Furthermore the two largest native H6 groups contribute ~36.1% of the 122. Account for H6 grouping; do not treat sites or 2–3 months as independent trial observations. The fixed 2021 WorldCover mask also carries time-varying class selection uncertainty.

### Main manuscript consequence

This genuine **same river, same month, 2018→22** macro panel no longer has the original basic criticism that the freshwater and riparian layers are entirely different spatial/time units. However, on this physically matched panel the expected **water TP↔riverbank greenness change is not empirically supported**. Thus do NOT illustrate a positive wastewater→stream TP→riparian biomass→butterfly/foodweb chain with current data.

A credible one-paper core can compare: (1) original wastewater-reported source geography and within-water benthos status, (2) landscape context predictive of stream EQR, (3) stringent same-site source-compatible water–land near-far optical test with supported **non-detection**, and (4) stronger limitation of source-selection and cross-country generalization. If a strict positive physical transfer is the scientific requirement, the presently available macro observation cannot establish it.

## Code, artifacts, and scope

- Source only 4 shards: [Actions #44018628546](https://github.com/zuizui0223/gesui/actions/runs/44018628546), artifacts 14299413733,14300991114,14300196692,14299953285; SHA256 verified.
- Posthoc S2 source method `analysis/20261117_FR263_riverland/fr178_2018_2022_frozen_2scene_cloud_mosaic_QA_POSTHOC.py`.
- TP↔NDVI model `analysis/20261117_FR263_riverland/fr178_POSTHOC_2scene_mosaic_vs_real_JJA_TP_H6.py`.
- Original independent local numeric QA `FR178_posthoc_2scene_realTP_NDVI_independent_validation_20261117.json`, source QA `FR178_POSTHOC_TWO_FROZEN_SCENE_SOURCE_QA_20261117.json`; original user-visible local downloadable source CSV contains 122 site/year near-minus-far NDVI but not third-party raw pixels.
- All files are Draft PR #120 branch; no `main` merge or automatic redefinition of pre-specified main effect.
