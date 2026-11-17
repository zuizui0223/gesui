# 2026-11-17｜フランス263河川：水陸同地点マクロ解析と衛星二時点検証の判定

Status: **FR109 actual fixed-model numeric results done; FR263 original riverline source QA done with pinned third-party conversion; FR24 H6 2018/22 real Sentinel2 scene metadata PASS; 8-H6 true pixel QA submitted but not yet numerically certified**. Review-stage exploratory study. No new evidence of sewage-derived N/P flux into land, no claims of ecosystem services.

## 1. Genuine common physical site: freshwater chemistry, benthos, riparian mapped classes

The ORIGINAL 263 physically identical French WQ (real measured TP) + I2M2 EQR stations repeated 2018/2022 have been joined with independent official Copernicus Riparian Zones2018 classification polygons by exact `site_key` in original EEA2018 source. All 263 ID joins valid, of which genuine >=70% mapped polygons in both station-centred 0–250m and 250–500m radial zones yield **109/263**. Source input before opening TP/EQR: 49 H7 native catchments, 20 H6 parents, 8 geographical 3° blocks, **only 1/109 original positive upstream reported 2018 WWTP P source**. Both rings have >10% non-water mapped area in all 109. This passes exploratory preoutcome minimum 80 sites /30 H7 /12 H6 /5 blocks.

Across these 109 true river monitoring stations, original observed water logTP 2018↔EQR2018 raw Spearman rho **−0.56511**, but real water logTP2018 vs near−far mapped natural-cover composition rho **+0.06009**. Same near−far 2018 natural-cover response vs 2022 river logTP rho +0.06201 (the static land map is NOT a repeated 2022 measurement). Near natural-land fraction 0.57989, far 0.58233, near−far **−0.00244**.

Fixed original 5 source fold indices are **0..4**, H6-disjoint; Ridge10 train-only standardized, H7-equal station fit weights, original source geography/flow/urban/agricultural landscape/far classes as baseline. Heldout error grouped by 20 actual H6 parents, 2000 resamples of **fixed parent losses** are descriptive uncertainty:
- Added original observed 2018 logTP → predict **2018 near−far natural habitat fraction**: original H6-equal MSE improvement **−0.776%**, descriptive95 **[−5.686,+3.110]**.
- Added near−far terrestrial natural-fraction → predict **2018 actual benthic EQR after same-site measured TP**, original H6-equal MSE improvement **−4.181%**, descriptive95 **[−11.284,+0.707]**.
- Added near−far terrestrial natural-fraction → predict **2018 actual measured logTP**, original H6-equal improvement **−1.020%**, descriptive95 **[−6.001,+4.225]**.
- Source/coverage bias: among all 263 real aquatic stations, original observed 2018 TP medians: RZ map-qualified 109=0.083 mg/L, nonqualified154=0.0975 mg/L; original benthic EQR medians 0.72658 vs 0.70235. No outcome-based sampling, but lack of mapping is not random.

Interpretation: an identifiable matched physical water-monitoring-site aquatic TP/benthos signal, with **no additional out-of-region predictive evidence from mapped 2018 bank land COVER fraction**. Neither material transfer nor lack of transfer can be concluded. This remains posthoc and nonconfirmatory; do not change widths/classes or folds to make an affirmative result.

Primary real-source file authentication from original GitHub Actions FR263 #43373639391 / pre-outcome #43373061409; original RZ source #44004357868. Local SHA-verifying code and JSON in active chat artifact `fr263_RZ2018_realTP_EQR_same_station.py` and `gesui_FR263_RZ2018_TP_EQR_same_station_20261117.json`; separate full numeric account [FR109 report](https://github.com/zuizui0223/gesui/blob/analysis/gesui-20261117-onepaper-pathway-identification-decision/docs/RESULTS_20261117_FR263_RZ2018_REAL_TP_EQR_SAME_STATION_109.md).

## 2. Authentic riverline geometry, independent of biology outcomes
The **official** HydroRIVERS original zipped archive `https://data.hydrosheds.org/file/HydroRIVERS/HydroRIVERS_v10_eu_shp.zip` failed GitHub runner download with **HTTP 403** [source run #44008867522](https://github.com/zuizui0223/gesui/actions/runs/44008867522). This failure is retained. A **different original-source-verified third-party conversion**, fixed SHA256 `5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f` from `cassini-team-todo/hydro-rivers-europe`, preserves HydroRIVERS v1 `HYRIV_ID` and true line geometry; newly executed source-only [GitHub Actions #44009618859](https://github.com/zuizui0223/gesui/actions/runs/44009618859), artifact 14296271408.

- All **262 original unique HYRIV_ID** found; 263 true station coordinates matched by original ID (263 physical sites, one reach shares two stations).
- True original riverline station offset <=250m: **226/263**.
- Real original same-reach line segment station-centered available length >=900m: **224/263**.
- Both original fixed gates: **192/263 stations, 61 H7, 23 H6**, exceeds preoutcome geometry thresholds 120 stations/40 H7/15 H6. Original station-to-true-line offset median118.53m, max490.8m. This is geography-only, not habitat pixel support.
- **Intersection with the earlier 109 station-centred RZ2018 class source**: **66 stations, 37 H7, 20 H6, 7 geographical regions, only 1 upstream WWTP-report-positive**. Because 66 **falls below the original 80-site minimum for the primary RZ2018 ecological comparison**, do not refit a new favored 66-site model with a relaxed gate. Treat 2018 station-centred RZ result as station-proximity analysis and this strict intersection as a separate source feasibility STOP.
- Third-party source transformation provenance is explicit; do NOT say direct original official ZIP was downloaded or byte-identical verified.

## 3. Genuine 2018+2022 public Sentinel-2 imagery is now source-supported, not yet full pixel-eligible
Previously locked 24 real original source H6 parent groups, one station chosen by preoutcome SHA256(field station id) per H6, source year2018 and 2022 / Jun-Jul-Aug, each month up to 2 original scene options with red/NIR/SCL metadata, original reflectance scale+offset individually. Original [source GitHub Actions #44008867555](https://github.com/zuizui0223/gesui/actions/runs/44008867555) original artifact14296265400, JSON SHA256 `22414dec230bc9bf9c308d823f166303cc843143c6ade922af6cb44b1f2c3f3c`: **24/24 original H6 each have >=2 summer months with true original 2018 and 2022 coverage**, exceeding metadata-only threshold18/24. This is not the cloud-free pixel gate.

Next independent *new* signed 8-of-24 H6 smoke selected by SHA256(H6) **before looking at satellite reflectance, water chemistry/EQR or outcome signs**, [workflow run #44009758637](https://github.com/zuizui0223/gesui/actions/runs/44009758637). For each 2018 &2022 month, sample true source B04/B08/SCL original COG windows + ESA WorldCover2021 class-matched natural cover and *source metadata radiometric offset*; require >=2 paired same-calendar-month acceptable scenes per fixed H6 and >=6/8 passes before opening full24 true pixel opportunity. The pixel smoke uses **station-centred rings only** and is not itself source original true riverline-distance habitat proof. Preserve all failures and prior independent AT67 S2 natural-cover 5/10 STOP unchanged.

## 4. Publication/inference decision
**Supported:** original independent same-site aquatic nutrient–benthic co-gradient and landcover description; a 2-year same-station satellite source cohort can now be attempted across distinct European river regions.

**Not supported:** source-specific municipal sewage water→plant nutrient flux, same-date productivity or biomass, butterfly larval feeding/puddling, emergent aquatic insects' terrestrial trophic transfer, or measured ecological services.

**Next step once genuine pixel gate has actually passed:** build near–far 10m vegetation **NDVI** at same original physical river stations for both 2018 and 2022 (or source-authenticated independent HR-VPP TPROD), predeclare landcover stability and riverline coordinate checks, then compare real 2018→22 optical index vs real measured river TP/EQR in H6-group-heldout tests. Never recast new NDVI as the preexisting HR-VPP TPROD product, and never interpret static 2018 RZ map as 2022 land productivity.
