# 2026-11-16 EEA2022 true annual measured N/P outfalls → true river-path riparian Sentinel2 source study

**Source-only and technical sensor QA, NOT ecological wastewater effect.** Main untouched; new exploratory Draft PR #116. All source values originate from exact signed original EEA2022 annual measured N/P reporting and original public remote sensing. No biological/river function coefficients, independent causal experiments, biomass, adult butterfly resource use or ecosystem services have been estimated here.

## 1. Why change the sampling system?

The original n16 HydroWASTE source-blind outfall ±1500m true river-site design (PR #115) cannot receive annual 2022 reported EEA UWWTD 2022 N/P by assumption: among 8 original DE/AT HydroWASTE outlets, none is within 1km of any public original DE/AT annual both-positive reported source from PR #23. Assigning each source nearest EEA facility's mass would be false. Original local 250m GBIF plant field/year support 0/16 both sides, 1km 2020–22 only 10/16 and two-era only4/16. Thus a new original outfall-annual report source is necessary.

## 2. Exact annual source, real flow direction, true route coordinates

- Source: original published 2022 EEA UWWTD genuine **1,156 facilities with measured reported annual N and P mass and one direct FW outfall**, DE630 /AT526. Original Git blob SHA1 1b4fbdfd8516e4040615a6cc1737febb1b6a33f0, originally source-verified PR #23. Source annual N and P tonnes/year are not actual seasonal discharge concentrations at roots/biota.
- Network: original true HydroRIVERS v1 Europe third-party converted GeoParquet SHA256 5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f, native NEXT_DOWN graph, distinct U2/U1/outlet/D1 linked reaches, 100m nearest river-source snap, runner-up gap≥25m, 100m connecting reach endpoints, exact true route ±1500m up/down and nonoverlapping 250m circles.
- [Original SHA verified Actions #43850677284 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43850677284). **282/1,156 true source/river outfalls valid** — DE137/AT145. Original invalid/uncertain counts: 667 source outlets >100m from line, 51 no mapped local river, 128 no two upstream graph links, 18 short path, 8 ambiguous two reaches, 2 missing mapped downstream. ALL source sites retained as denominators; no imputing those into an ecological treatment. The third-party conversion is not independently byte/topology-verified against official original release.
- Before satellite values, selected fixed independent/5km-space-separated **48** geographic pilot original reported source IDs, **24 DE/24 AT**, stratified across coarse 2° source cells and deterministic SHA tie order; no satellite green pixels or species response selected these.

## 3. Annual 2018 and 2022 satellite source availability

- Original Element84 Earth Search Sentinel-2 C1 L2A Collection1 STAC 2018 scene source valid ≥3 dates for all48, **2022 source scene n0/48**. Source [Actions #43851296445 SUCCESS/STOP](https://github.com/zuizui0223/gesui/actions/runs/43851296445). This is a documented 2022 C1 archive backprocessing gap (Element84 Earth Search issue #29 and #70), NOT biological or satellite observation absence.
- A separately pre-frozen alternative used **SAME legacy Sentinel2-L2A collection in BOTH 2018 and 2022**, no mixing old and C1 reflectance archives. Only original May–Sep, source scene cloud≤40%, true source patch coverage across BOTH 250m circles, source red/nir/scl public COG assets. **48/48 original source outfalls** have ≥3 distinct genuine acquisition dates in EACH year, both countries 24/24. Signed [Actions #43851791936 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43851791936). This inventory selected 917 low-cloud monthly candidate scene assets. The original 2018 red/nir reflectance scale0.0001 offset0, whereas sample 2022 assets have scale0.0001 and offset −0.1. Never compare uncorrected 2022 DN with corrected 2018 reflectance: that would create an artificial spectral trend.

## 4. Actual sensor pixel local cloud/land quality gate

Before reading the first Sentinel COG pixel, froze 12 geographic response-blind original source systems (6 each DE/AT sorted SHA source ID) and June/July/August each year, up to2 lowest cloud scenes per month, 250m circle in actual source 10m CRS, 20m original SCL resampled nearest, accept only SCL4+5 (nominal vegetation+bare land), reject cloud/shadow/water/snow and invalid original DN→reflectance using genuine source metadata per scene. A paired local scene/month passes only if BOTH upstream/downstream circles have ≥100 physical clear land pixels (out of roughly1960 site pixels). Two-era season comparison needs ≥2 common months in 2018 and2022 for ≥8/12 independent source systems and ≥3 from each country.

[Correctly executed real COG/radiometry/SCL smoke #43853111036 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43853111036). The original first attempted run #43852993116 failed on a rasterio integer-vs-string GDAL cache option before any actual pixels; fixed WITHOUT changing biotic/sensor science thresholds. Then real COG:
- Original scene attempts with technical clear-land result: 55 pass, 37 fail local quality, no 404/COG read failures.
- Original **2018** >=2 months per source: **12/12**.
- Original **2022** >=2 months per source: **6/12** (AT3/DE3); at least one month 10/12.
- Original paired 2018+2022 ≥2 SAME Jun/Jul/Aug months in BOTH years: **6/12**, below fixed required8 -> **STOP_SENTINEL_COG_12SOURCE_LOCAL_PIXEL_QUALITY**, NOT a negative sewage/vegetation effect.
- Source original 12-pilot STOP remains valid. Do not retrospectively lower the threshold to6, select greener months, switch back to wider GBIF 14-year windows, or claim a satellite experiment.

## 5. Post-12-source exploratory check of original full48 (not the same confirmatory sample)

A separate AFTER-smoke source-only protocol was frozen to apply the **identical COG/sensor/cloud quality** to ALL original 48 geographically selected measured N/P one-outfall source systems. Primary source opportunity is SAME original source year2022 only, ≥2 locally clear Jun/Jul/Aug months upstream+downstream, minimum24/48 systems and >=8 from each original DE/AT. Secondary paired 2018+2022 has the same count criteria and true same months in both years. The previous 12 pilot failure does not disappear, and this full48 is not an independent confirmation, only broader technical source feasibility. [Full48 original raster QA Actions #43854317398](https://github.com/zuizui0223/gesui/actions/runs/43854317398) initiated; no success/result is asserted in this memo until the run reaches completion.

## 6. Why this still does not identify the wastewater food-web causal chain

- We now HAVE source-specific **annual 2022** outfall reported N/P mass and real ±1500m upstream/downstream river sites, plus real same-season satellite images for all48 metadata source and demonstrably readable original land pixels in the fixed 12 technical pilot.
- We DO NOT YET have 2022 locally-clear seasonal-paired satellite source for a sufficiently large new independent outfall frame, stable 2021 WorldCover tree/grass/wetland-vs-crop mask, independent same-river reference without actual WWTP interposition, receiving-water measured N/P at these identical local sites, or measured aquatic benthic/water-to-land flux and real larval feeding/plant abundance.
- Current SCL land classes 4+5 do not mean ground-verified riverbank plants. A 250m circle can include gardens/fields/streets; satellite NDVI measures optical greenness not biomass or biodiversity. Source annual 2022 mass does not imply local concentration or actual effluent delivery in June-Aug, and there is no comparable facility-specific original 2018 measured mass to identify a treatment change.
- Original 2022 EEA BQE source has 0 genuine German benthic field EQR observation station and only68 Austrian field observations, so a DE+AT 2022 UWWTD-source×Sentinel field experiment cannot be asserted to have balanced matching same-location measured aquatic benthic responses. The French 73 H7/24 H6 real 2018/22 receiving-water TP–benthic EQR and original plant NUTRI–butterfly HNS are OTHER location/biota cohorts; not a water→land causal mediator for this DE/AT source study.

## 7. Revised single-manuscript priority

A defensible one-paper ecological contribution must use truly common WWTP source populations and match a *specific realized* transfer channel, such as environmental receiving-water TN/P plus directly connected site-level riparian primary productivity (with source/year/river controls) or actual emergence mass transferring to spiders. Only then can plant/animal trait exposure/access patterns test response variation. One can write a scientifically honest **methods/data-source paper** on observation-grain accessibility and scale mismatch with current QA; it is **not** a verified ecological response to sewage yet.

Main/chocho unchanged, all audits on Draft PR #116, with original third-party raster data not rehosted.
