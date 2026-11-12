# HydroBASINS L7 × 2017–2019/2020–2022: genuine same-scale Plantae, Lycaenidae, reporting opportunity and HydroWASTE

**Source contract fixed 2026-11-12 before any new native basin polygon GBIF count is examined.** Work is posthoc relative to previous EQR/EPT and port-plant positives. This phase measures source availability and may not produce a biotic association if sampling support fails.

## Population, geography and time: fixed prior source

- Exactly **1,219 real Europe/USA HYBAS_ID_L07 river basins** with ≥1 true EEA2018, EEA2022, EPA2018–19 monitoring event, as defined in pre-outcome source manifest SHA256 `a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89`. Do not discard those lacking GBIF butterfly occurrence or rerank by known EQR sign.
- Exactly **2 bins**, `2017,2019` and `2020,2022`; each basin queried in both windows, including structurally unobserved US P2 benthic biology (retains biological NA, not zero). 2,438 basin-periods × 3 taxa = **7,314 intended count requests**, adjusted only for documented MULTIPOLYGON components without changing geography.
- Original exact official HydroBASINS v1c Europe and NA Level07 ZIP SHA256 EU `b9b810f4db271f0fa6b39fecf01784041b948944b3c3f296a393766388e86b44`, NA `4e5150f0fe3a2c3d17b21938ec5c6ca1959d136497f572a0203af6b1842e501d`. Polygons from same original ZIP/attributes source as basin assignment. Source polygons in WGS84. Preserve outer and inner rings and all MultiPolygon components; **never substitute a centroid/radius circle** as a native basin.
- Stage0 generate true polygon WKT with topology-preserving simplification to meet query length ≤3,500 characters and ≤5% relative original planar polygon area error, with maximum tolerance sequence {0, 0.002, 0.004, 0.008, 0.015, 0.03, 0.06} degree as predeclared approximate simplification. Record original geom SHA, WKT SHA, area retention, components, WKT max. If a component remains unqueryable, flag basin *geometry source missing*, not butterfly absence. L7 geometry eligibility requires ≥85% of the initial 1,219 original source basins. Stage0 must be performed **before** new GBIF results. Store source-only pinned per-basin WKT in short-lived Actions artifact (no outcome responses).
- API geometry unsupported classes or original invalid geometries: no silent `make_valid` or largest-part-only rescue. If multiple components, query every disjoint component and sum indexed occurrence counts (beware geometry boundary duplicates; uncommon but note).
- Source-limited sampled pilot subset: **all 1,219 basins** is the intended full macro population. Four immutable shards via SHA256(HYBAS_ID) mod 4 can execute source queries at low GBIF concurrency; their union must equal original source manifest. Run GBIF index requests on whole source frame, not a selected 234-basin ecological subset.

## Original GBIF count-only source gate

GBIF occurrence **search** `/v1/occurrence/search?limit=0` with exact polygon `geometry=WKT`, original `year=2017,2019` or `year=2020,2022`; `hasCoordinate=true`, `hasGeospatialIssue=false`, `occurrenceStatus=PRESENT`; counts **not occupancy or species richness**.

Three source-fixed taxon categories per polygon×period:
1. Plantae: GBIF `taxonKey=6` (verify current accepted taxonomy and kingdom).
2. Lycaenidae: **GBIF legacy Backbone family ID 5473**, explicitly `checklistKey=d7dddbf4-2cf0-4f39-9b2a-bb099caae36c` (do NOT use iNaturalist 47923).
3. Animalia: GBIF `taxonKey=1` (recording-opportunity comparator, not unaffected background group).
Preserve official API exact response SHA256, request URL, HTTP status/error, timestamp, count, original polygon component, and original GBIF taxon lookup response identity. **HTTP errors are missing, never 0**. Stop/retry with bounded exponential backoff for 429; 4 or fewer concurrent source requests per shard (max 16 across matrix jobs); no unspecified adaptive radius or period alteration.

Sample support source gate, evaluated before any cross-taxon ecological model:
- ≥90% of all 2,438 original basin-periods with all 3 valid original source requests (and full polygon components, not just one part),
- ≥200 supported basin-periods for **Plantae≥50 AND Lycaenidae≥10 AND Animalia≥50**, covering ≥8 source countries / ≥50 distinct independent HYBAS7 IDs,
- ≥100 distinct Level7 basins with the two original periods both *complete* for the three taxa (support only, not amount).
If fail, report source geometry/record coverage and keep original model outcome columns **NULL**. Do not loosen thresholds or choose a coarser level to rescue desired positive coefficients. If pass, freeze exact source eligible basin-period rows before opening any ecological association, then perform species faceting and comparison to aquatic EQR/EPT on the same rows.

## WWTP exposure

Use original publisher [HydroWASTE v1](https://www.hydrosheds.org/products/hydrowaste) 58,502 global treatment plants, official CC BY 4.0. Official linked [Figshare source dataset](https://figshare.com/ndownloader/files/31910714) may be inaccessible from runner; first verify response byte source checksum, file structure, row grain and `HYRIV_ID` / estimated outfall coordinates, pressure type. **Do not substitute old WPI port-model WWp**. HydroWASTE is static/mixed estimated flow, not 2018/22 observed annual N/P mass, and a WWTP inside a basin does not necessarily affect every monitoring site therein.
- First reported variable should be registered **estimated outfall pressure within a local L7 subbasin** or explicitly upstream-routed outfalls by `NEXT_DOWN` and `HYRIV_ID`. Keep two variables distinct: local outfall versus cumulative upstream network pressure; cap loops, conservation check and attribution.
- Gate first source: reported 58,502 unique global plants, valid geocoordinates, ≥85% row outfall assignment, and documented consistent units / valid source flow; if not then source-only STOP without coefficients.

## Genuine common-unit biological test (conditional; no results yet)

On **same** `HYBAS_ID_L07 × PERIOD` rows, keep underlying method-specific responses separate:
- EU 2018/2022 standardized EQR by country/assessment methodology, actual WWTP upstream reported N/P as source-specific signal.
- US 2018–19 stream TP/TN, benthic EPT fraction, riparian structure; no US 2020–22 falsely imputed.
- Plantae original accepted species and candidate host richness requires **separate facet/source-quality stage**, not the above raw index counts. No inferential model on raw GBIF counts pretending plant diversity.
- Lycaenidae record proportion relative to other butterflies requires actual same-basin year Papilionoidea denominator and/or survey effort, not the broad Animalia count alone. The present count-only audit is just availability.
- Control flow, cropland, urban cover, river length, basin area, country, source and observation effort. Use basin holdout / 234 EU repeated basins for time stability. No serial SEM until actual matching upstream effluent, water chemistry, aquatic emergence, riparian plant hosts and butterfly visit/seed set.

**Hard exclusions:** no eBird; no merger of coastal WPI port flora with inland river outcomes, no accidental annualizing of old 2011–24 GBIF facet, no release of unrelated unpublished files, no artificial butterfly zeros, no merge to `gesui/main` or `chocho`. All original confirmed positives and nulls remain untouched.
