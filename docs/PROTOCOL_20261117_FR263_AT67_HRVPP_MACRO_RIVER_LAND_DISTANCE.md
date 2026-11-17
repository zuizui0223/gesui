# 2026-11-17 | 公開データだけで水域・陸域を結ぶマクロ解析 — 河岸距離勾配

**Status: NEW independent exploratory analysis protocol / source-first plan. No new raster has been retrieved, no new ecology coefficient fitted, not a successful dataset gate, and not a confirmatory preregistration.** This replaces a proposed field survey as the immediate research priority. Preserve all earlier outcome/technical STOP decisions. No changes to main and no assertion of sewage-specific river-to-land transfer.

## 1. The one ecological question

> Do observed river nutrient / aquatic community gradients covary with terrestrial vegetation productivity **preferentially immediately beside the same river reach**, rather than simply with regional terrestrial vegetation?

Only after this spatially explicit land response can a sewage-origin link be tested in a cohort with actual WWTP source/Q; a river TP index itself is **not** sewage N/P. Direct movement of nutrients or of emergent aquatic insects is *not* observed by vegetation indices.

The importance of wastewater is its combination of source-identifiable directional input, water/baseflow subsidy, N/P/organic resources and contaminants/stress. These pathways might have opposite ecological effects; a positive association cannot be recast as biodiversity gain or ecosystem service.

## 2. Actual source-verified starting populations: no fusion of mismatched frames

**Primary France panel: original fixed 263 physically repeated aquatic monitoring stations, 2018 & 2022**, native 73 H7 / 24 H6. Each has authentic datewise measured TP and matched-method I2M2 benthic EQR; original TP has >=4 real positive measurements in >=3 months for both years; 263/263 have >=3 overlapping sampled calendar months (250/263 >=6). Use the original *month-balanced* TP where possible; do not newly select stations on outcome. Existing spatial TP vs EQR strong in both years, but within-site 2018→22 \u0394TP vs \u0394EQR **not robust** after month matching. Source documentation: `docs/RESULTS_20261113_FR263_SAME_SITE_TP_2018_2022_COMMON_MONTHS_I2M2_THREE_LAYERS.md`; Actions #43378917089; source and official TP.2 checksum retained.

**Separate Austria external panel: original 67 genuine 2022 aquatic monitoring stations**, five geographical region blocks; 2022 real receiving-water TP, official benthic EQR, 2022 upstream *reported* WWTP load / Q and WorldCover local landscape. Original 68 minus one source-coordinate disagreement >950m, excluded before satellite analysis. The existing Sentinel-2 masked-reflectance 10-site source-smoke was 5/10 against >=6/10 threshold: **technical STOP remains unchanged**. The new Copernicus HR-VPP product is an independent source opportunity, not a rescored old test. Sources: draft PR #118 and #120. No claim that annual source N/P was measured within the bank pixel, or that N and P are separately identified (N/Q vs P/Q original rho ~0.945).

**Do not** collapse the France, Austria and EU-wide 447 H7×period species-trait panels into one mediation sample. The previously verified plant MOIST→butterfly host MOIST OOS result and plant NUTRI→butterfly HNS result are within-terrestrial, broad-community associations, not a river-bank consumer link.

## 3. New public, spatially continuous terrestrial measurement

Official Copernicus Land Monitoring Service **High Resolution Vegetation Phenology and Productivity (HR-VPP)**:
- 10 m European raster, 2017–present; `Total Productivity` = integral of Plant Phenology Index in growing season; annual season1 `TPROD` is the first source product to verify.
- Independently obtain the **same year's** HR-VPP quality flags for 2018 and 2022; optionally add `Seasonal Productivity` and 10-daily PPI only as prespecified sensitivity. Values measure remote-sensed greenness/productivity *indices*, not measured dry biomass, nutrient concentration, or production in carbon units.
- Documentation: https://land.copernicus.eu/en/products/vegetation/high-resolution-total-productivity ; https://land.copernicus.eu/en/products/vegetation/high-resolution-vegetation-phenology-parameters-quality-flag ; https://eea.github.io/clms-api-docs/download.html
- Land-use context: Copernicus Riparian Zones 2018 55-class layer (only selected rivers, so coverage must be audited), WorldCover 2021 as contextual classification (not a pre-2018 control), original HydroRIVERS line and nearest waterbody/tributary geometry. Do not use Google Maps screenshot pixels as numerical remote sensing.

## 4. Identical physical reach, same year, explicit distance

At each original FR263 aquatic station, verify snapped **true HydroRIVERS line geometry** and native same-reach identity. Construct a fixed +/-500 m *along-channel* analysis segment, both banks. For each land pixel compute the **perpendicular / shortest distance to the correct riverline**, masking water and invalid pixels, and flagging nearby other streams. The fixed lateral bands:

- B1: 20–100 m (primary near-bank)
- B2: 100–250 m (attenuation)
- B3: 250–500 m (attenuation)
- B4: 500–1000 m (primary far reference)

Boundaries are **design values**, not observed effects or biological transfer distances. Match masks/coordinates across 2018 and 2022 and record cover class and quality; report all source denominator counts by band/year. Avoid mixing crops vs natural riparian vegetation under a single physiological baseline; primary contrast within fixed predeclared vegetated cover strata, plus separate all-vegetation sensitivity, reporting the classes excluded. In years where HR-VPP season1/quality was not measurable, mark NA, not zero.

New primary outcome per physical station/year:

`R_it = mean_valid_TPROD(B1_it) − mean_valid_TPROD(B4_it)` (standardize product only within legitimately comparable source grid/season; record units/scale factors).

This is a **lateral vegetation gradient** beside the actual river, not total H7 greenness. Keep B2/B3 as fixed attenuation test, not replacements selected by signs. Examine whether greater river TP or better measured benthic EQR accompanies a distinctive near-bank response.

## 5. Falsifiable, prespecified conditional analysis

**Source eligibility (run before opening HR-VPP values):** all original 263 stations and both years must be attempted. Validate geolocation/reach snap from source independent of water quality and vegetation values. Prespecify technical gates for first extraction: >=100 valid 10-m vegetation/quality pixels per station-year *each primary band*, continuous same cover stratum across both years, >=120 original paired stations, >=40 H7 / >=15 H6 parents; report n by gate/reason (if fail, **STOP**, not biological null). These thresholds are new feasibility proposals, not previously validated or inherited confirmation gates.

**Primary cross-sectional/replicated-year test:** partial river TP–`R_it` association controlling predeclared land-cover mix, river reach size, upstream cropland and urbanization, flow proxy, season/climate, elevation, and geography. Study river EQR→`R_it` **separately**, rather than using EQR as a proven cause/mediator. The key interaction equivalent is: `TP × bank-distance`; only the difference of slopes B1 vs B4 supports a *river-localized statistical spatial signature*. Report the contrast and its uncertainty, not "B1 p<0.05 / B4 p>0.05".

**Independent spatial validation:** hold out entire H6 parent regions (24 groups) with identical folds and honest out-of-H6 MSE improvement against an environment/land-cover-only baseline; additionally report source H7 or parent-cluster uncertainty. Stations/pixels within the same native H6 must never enter different train/test folds. Compare `R_it` model vs near vegetation alone; the far-band reference should control regional greenness.

**Secondary temporal test:** at exact same 263 source sites, `Delta R_i = R_i,2022 − R_i,2018` vs original **calendar-month matched** `Delta logTP_i` and, in a separate analysis, `Delta I2M2 EQR_i`, with changed local climate/observation support. Existing month-matched \u0394TP→\u0394EQR null must be retained and cited. Two years are not a BACI treatment and do not resolve lag.

**Negative controls:** a far belt with genuinely similar cover, station/reach label permutation **within H6** preserving sampling effort, and a no-hydrologic-connection neighboring drainage where verifiable. Baseline rainfall, agricultural/fertilizer landscape, urbanity, floodplain width, canopy, topography and satellite compositing can explain near-river vegetation and water TP together. A near-band coefficient is not proof of N transfer.

**Austria separate replication:** apply identical HR-VPP bands to all original 67 physical 2022 stations, and compare `R_i` against original TP plus reported upstream P/Q, landcover, Q and the five genuine region splits. Primary source term choose P/Q; N/Q only a sensitivity because of extreme collinearity. Report each of the five region-heldout errors. This is an external *source pressure association*, not a causal WWTP-specific pulse or a formal meta-analysis with FR effects.

## 6. When/how butterfly becomes part of the macro chain

Only **after river-bank terrestrial response source support**:
- Search actual standardized geo-referenced eBMS transect/section-year 2018–2022 data within the same river corridors. eBMS full scientific records generally require access request/approval; these are not automatically public unrestricted microdata. Official info: https://butterfly-monitoring.net/ebms-data%20access
- If access not available, test *source-only* original GBIF coordinate + date + uncertainty + observation-effort eligibility along these exact same reaches. Prior 16 WWTP 250m bank source returned 0/16 passing 20 plant records on both sides in either 3-year window; this failure cannot be rescued by arbitrary widening. It says nothing about the station-based FR263 sample, which must be audited anew without assuming success.
- Test terrestrial **community composition/habitat traits** conditional on river distance and vegetation, not adult feeding, larval host consumption, abundance or pollination. Spiders/actual emergence constitute another *separate* trophic-channel hypothesis, and EQR does not measure emergence biomass.

## 7. Interpret possible results before seeing them

1. **TP/EQR linked to B1 more strongly than B4, decay via B2/B3, robust H6 heldout:** observational evidence of **river-localized cross-ecosystem coupling**; wastewater attribution remains Austria/source-only and ambiguous.
2. **Comparable B1 and B4 slopes:** likely broad landscape/shared watershed environmental gradient, **not** a localized water-to-land transfer.
3. **Neither:** no detected river–land signal under this product/source/time scale, **not** absence of ecological subsidies; give power and source-opportunity diagnostics.
4. **Spatial but not 2018→22 temporal:** stable co-gradient with no tested temporal transmission; no causal perturbation inference.

**Next executable operation:** recover original source-only *physical 263 station coordinates and HydroRIVERS IDs*, then query CLMS HR-VPP metadata+quality coverage at exactly those coordinates for 2018/2022 before looking at TPROD values. Record full source hashes, downloads, excluded pixel counts, sample gate and STOP/GO. No new ecology analysis has been run by this protocol commit.

## 8. Source trace
- FR263 source/report: `docs/RESULTS_20261113_FR263_SAME_SITE_TP_2018_2022_COMMON_MONTHS_I2M2_THREE_LAYERS.md`; PR #71; source Actions #43378917089.
- Water-land H7 mismatch FR73: draft PR #114.
- Same WWTP spatial source failed at 250m: draft PR #115.
- Austria67 2022 aquatic/satellite source and independent original S2 5/10 STOP: draft PR #118.
- Ongoing one-paper negative/positive evidence freeze: draft PR #120.
