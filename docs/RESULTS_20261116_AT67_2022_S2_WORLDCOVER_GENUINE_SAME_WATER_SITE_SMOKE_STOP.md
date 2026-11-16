# 2026-11-16: Same 2022 Austrian water-site TP/EQR/WWTP source vs true riparian Sentinel2 + WorldCover land-class matching: original source result and STOP

**Status: Original official 2022 Sentinel2 L2A metadata availability PASS on all 67 true field aquatic sites; predeclared class-matched riparian and reference-land pixel SOURCE SMOKE gate FAIL (5/10 versus required >=6/10). Therefore full67 pixels, cross-compartment NDVI source–water BQE correlations and any pollution effect are NOT estimated.** A failed source feasibility gate is neither a biological absence nor proof of no wastewater effect.

## Original fixed same-site source design

- Start with 68 actual Austrian 2022 EEA benthic field monitoring sites, exact original BMLUK FW1 unfiltered total river P annual means at all68 and original EEA2022 source modeled upstream eligible 2022 annual WWTP P/N source P/Q exposure at 53/68. National 2022 true GZÜV official PKINT↔EEA 68/68 station-ID validation gives one real coordinate discrepancy approx **954.614 m** at ATFW54110117, excluded BEFORE any Sentinel2 source return. The physically unambiguous primary = **67 original actual water sites**, fixed 3° geographic region counts 15_3:9,15_4:14,15_5:21,16_4:13,16_5:10.
- Real land use is not necessarily confirmed riverbank habitat. Pre-spec inner **250 m radius around the genuine EEA2022 physical river station** and outer **250–500 m annulus** in the same 2022 Sentinel scene. This is materially more defensible geographically than EU H7 whole-basin GBIF plants/larval-host scores joined to a separated river chemistry monitor, but a buffered scene may still include other land and isn't actual water→soil transfer.

## STAC metadata stage: successful genuine 2022 source, no optical pixel effects

[2022 Earth Search public legacy Sentinel-2 L2A metadata original Actions #43884504952 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43884504952).

- All original **67/67** same-source aquatic sites have 2022 June–August genuine Sentinel2 B04 10m/red, B08 10m/NIR and SCL 20m original COG assets with full 500m source box coverage in at least two months. Max two lower scene-cloud original scenes per source month selected by metadata only, **394 source image candidates** stored with original IDs/asset URLs/scale and offset/scene timestamps (not image pixels).
- The original Collection1 2022 archive was previously incomplete for the separate AT/DE WWTP 48-outfall pilot; the same legacy L2A collection is used at all Austrian sites to avoid mixing sensor processing source archives. 2022 actual reflectance offset must be applied from original STAC, e.g. scale0.0001 offset−0.1 where native data so specify.
- This PASS establishes only *public satellite acquisition opportunities*, not cloud-free local pixels, protected natural vegetation or species abundance.

## Primary WorldCover natural-land × S2 pixel source smoke: STOP by prior fixed design

Original source-fixed [pixel QA protocol](https://github.com/zuizui0223/gesui/blob/analysis/gesui-20261116-AT67-same-riverbank-S2-landgreenness/docs/PROTOCOL_20261116_AT67_S2_2022_WORLD_COVER_LAND_PIXELS_SMOKE.md) required exactly **2 original true river stations per each of 5 regions** chosen by station-code SHA256 before opening any pixels, n10, same radius/season/WorldCover mask/classes and SCL criteria as intended full67. Pilot must have **>=6/10 true sites spanning >=4 regions** with >=2 of June/July/August 2022 valid SCL clear natural-land pixels in BOTH a 250m near river circle and a 250–500m outer annulus. Technical per-scene minimum: SCL4/5 clear-land >=100 near />=200 outer, ESA WorldCover v200/2021 stable natural landcover classes tree10, shrub20, grass30, wetlands90 >=40 near/ >=100 outer and >=10 same-class clear pixels in BOTH. Exclude source water6, cloud/shadow, urban/cropland/other classes. Original land source ESA WorldCover2021 v200 public 3°×3° COG tiles, CC-BY4.0, DOI 10.5281/zenodo.7254221.

[2022 original pixel + WorldCover smoke Actions #43885473411 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43885473411):

| Source fact | Result |
|---|---:|
| True same-source water stations attempted | **10**, two per original 5 geography blocks |
| Sites with >=2 monthly matches of proper natural-cover same class near/outer | **5** |
| Source pixel success regional count | 15_3:1, 15_5:1, 16_4:2, 16_5:1, 15_4:0 |
| Source precommitted SMOKE threshold | **6/10**, >=4 geographic groups |
| Site-month scene attempts failing source natural-clear pixel support | **31** |
| Remote original Sentinel2/WorldCover file download error requiring a fake alternative | **0 reported** |
| Decision | **STOP_SOURCE_NATURAL_LAND_RIPARIAN_LOCAL_PIXEL_COVERAGE** |

All ten original source stations and their per-month accepted/failed status are preserved in source artifact, no EQR or TP used for site/image/pixel acceptance. Source sites with >=2 valid 2022 natural class months: ATFW80218017 (block 15_3); ATFW61400237 (15_5); ATFW40505037 and ATFW40506036 (16_4); ATFW31000177 (16_5). The other five original sites fail the fixed source season/class pairing without claiming vegetation absent.

### Interpretation

The source catalog has plenty of real 2022 images, but the **actual 2022 river-centered 250m versus background 250–500m bands often lack sufficient source-clear NATURAL WorldCover classes on both sides** to calculate a class-matched optical land-greenness contrast over >=2 months. Expanding the study after seeing where vegetation is available would select different land use and subject the source study to post-result selection bias. An NDVI from all pixels, including cropland and urban gardens, would not be a comparable natural riparian vegetation response and was explicitly excluded before viewing pixels.

Thus, despite strong same river field source evidence upstream reported P/Q→measured TP→benthic EQR, no same-source terrestrial vegetation landscape phenotype has passed the fixed source gate and **NO TP/EQR/WWTP-P→NDVI numeric correlation was fitted**. That is a source/observation design limit, NOT an ecological absence of terrestrial subsidy.

## One paper decision and safe next empirical path

Supported within aquatic compartment:
1. Austria68 2022 actual upstream annual registered WWTP P/Q source alone added 16.3% unseen region TP prediction, and 19.3% unseen region EQR prediction; actual 2022 stream TP added 33.9% EQR prediction. Same-Q normalized downstream placebo source failed geographically, but may be unstable; no causal mediation identified.
2. Separate original France73 matched 2018/22 real river TP↔I2M2 benthic condition also robust, non-identical sample/flow exposure populations.

Terrestrial component remains an ecological hypothesis:
- EU451 whole-basin plant EIV NUTRI ↔ butterfly host nutrient habitat HNS are correlated geographical *trait affinities* but no local actual larval feeding and no measured receiver local sewage P. H7 measured TP→these terrestrial indices out-of-region predictions did not improve. Neither study captures true source river→plant P uptake or aquatic-emergent insect→terrestrial consumer biomass.
- Same-source Austria67 2022 raw satellite scene availability is **not** enough to claim local riparian plant productivity; genuine stable WorldCover natural-class matched pixel evidence at predetermined sites remains insufficient at the smoke gate. An alternative observational land-response study could choose a **new prospective** field sample of riverbank vegetation plots or nationally standardized riparian LULC habitats across a comparable source cohort; it cannot retroactively rewrite original smoke study or treat all cropland NDVI as unpolluted habitat.

Therefore one paper claiming **wastewater biologically enriches or impoverishes aquatic AND terrestrial communities through a measured nutrient food-web chain remains NO-GO**. A descriptive within-water 2022 source→chemistry→condition paper can be framed more defensibly, and source availability limitations can motivate a prospective coupled river–bank design. Do not average French/Austrian out-of-region prediction percentages as if they are a pooled causal treatment. Original main, chocho, and PR117 untouched by the stacked review-only new source branch.

## Exact scripts and source receipts

- Scene source: analysis/20261116_AT67_riparian_satellite/source_2022_same67_waterstation_stac.py, public source 2022 original run #43884504952.
- Pixel source: analysis/20261116_AT67_riparian_satellite/S2_2022_WorldCover_near250_outer500_pixel_source_QA.py, public source original run #43885473411.
- All site IDs, original scene IDs, source cloud metadata, SCL class counts and the explicit source STOP are archived in the Actions artifact; raw source COGs and publisher WorldCover landcover tiles are not redistributed.
