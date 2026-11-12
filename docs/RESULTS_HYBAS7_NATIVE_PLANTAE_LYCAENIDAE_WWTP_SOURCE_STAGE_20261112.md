# Same-scale river macroecology: NEW original H7-WWTP source completion and GBIF collection (2026-11-12)

**Original monitored EU+US basin source sample fixed before new GBIF results**: 1,219 distinct HYBAS_ID_L07 across two three-year bins 2017–19 and 2020–22 → **2,438 basin-period** rows, from earlier validated PR #62. This is an EU+US macro frame, NOT full worldwide river biota sampling.

## Stage 0: real original basin polygons, no small case-study radius substitutions

Exact official original HydroBASINS v1c 20.99MB EU and 18.30MB NA ZIP source identities reverified to original SHA256. **1,219 /1,219 original monitored subbasins** transformed into true **1,340 polygon components**, retaining holes/islands and every piece of multipolygons. WKT topology-preserving simplification preserved area within ±5%; maximum 3,500 source WKT characters (some polygons simplified). **Original pre-response WKT SHA `47f597e0a798cd6155ee9fb69ce8916ed5629951d046d114feb9f667963310e4`**, four deterministic fixed geographic shards (308/311/337/263 basins). Original native source geometry gate SUCCESS [Actions #43197894183](https://github.com/zuizui0223/gesui/actions/runs/43197894183).

**Technical HTTP400 anomaly, transparent:** first fixed 12-basin /72-query source preflight passed 62 requests, **10 failed HTTP400** due to long WKT GET URLs / intermittent transport. Source rejection log [#43198697091](https://github.com/zuizui0223/gesui/actions/runs/43198697091). Without reselecting a basin, radius, period or species, URL punctuation percent-encoding was corrected uniformly (same WKT coordinates/geometry SHA) and technical HTTP400 retried; **72/72** fixed identical original basin×period×taxon queries then SUCCESS [#43199075513](https://github.com/zuizui0223/gesui/actions/runs/43199075513). These are **index counts only**; success is not butterfly/plant population or a positive biological association.

**Full-source GBIF**: exactly 1,219 H7 basins ×2 fixed periods ×3 original taxa (Plantae legacy key6, Lycaenidae GBIF legacy family5473 with explicit backbone checklistKey, Animalia key1 opportunity comparator) = **7,314 source tasks**. Original native polygon component requests number ≥7,314. Four original source shards request all basins, source-preserved SHA/URL/status; **do not infer full-source results until all four are completed and combined**. Started [GitHub Actions #43199340731](https://github.com/zuizui0223/gesui/actions/runs/43199340731). Same sample unchanged; HTTP failures are missing, not zeros. GBIF source primary gate fixed **before** results: ≥90% complete basin-periods, ≥200 support rows with plant≥50/lycae≥10/Animalia≥50, ≥50 unique basins across ≥8 countries, ≥100 two-period complete native basins.

## Stage 0B: complete official worldwide HydroWASTE genuine original

Authentic publisher HydroSHEDS-linked Figshare `HydroWASTE_v10.zip` at `https://figshare.com/ndownloader/files/31910714`. Original publisher ZIP SHA256 **`6e19d8f31d95d533b2b29056b823eeb915d933fa22ea090f234dc38f6427d044`**; CSV SHA256 **`af040d9b4fdc93bab2eb8d1b52985f5bacc6b0db9f16795c29e31c10fb615dd4`**, 8.81 MB original CSV (some non-UTF8 names, read with reversible Latin-1 raw byte preserving original file). 58,502 unique original registered treatment facilities across 188 countries/territories; **58,502 original outlet coordinates valid**, 58,123 have original HydroRIVERS ID, 58,502 original modeled wastewater discharge numeric. Three treatment levels: Advanced22,855; Secondary34,766; Primary881. Registered original facility status **Not Reported n53,033 (90.65%)**, only 3,767 explicitly Operational; thus status cannot support a census of 58,502 confirmed active facilities.

[Original all-global HydroWASTE 58,502 source and field-quality QA SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43198783733).

## Stage 1: NEW source-LINKED HydroWASTE outlet to all Europe/North America L7

Original publisher 58,502 estimated WWTP outlets mapped through **the authentic 14,758 EU+NA Level7 subbasins** (14,055 valid original polygons; 703 invalid original polygons skipped without invisible repair). Exact genuine `NEXT_DOWN` topological hydrography traversed all 14,758 original basin IDs without cycles.

- **42,328** of the 58,502 global outfalls map uniquely to original valid EU+NA subbasins (the rest lie outside this two-region pilot / original geometry gaps).
- **1,079/1,219 monitored H7 basins (88.5%)** contain ≥1 original registered outfall; **14,695** registered facilities lie directly inside those 1,079 monitored subbasins (including statuses Not Reported and non-operational).
- **1,113/1,219 monitored basins (91.3%)** have ≥1 registered facility either locally or upstream in the HydroBASINS `NEXT_DOWN` graph, counting the full 14,758 node context.
- Source-time invariant measures: count and `WASTE_DIS` original modeled wastewater discharge sum **local to basin**, and cumulative upstream **to that basin outlet**. Note **a monitoring site inside the basin is not necessarily downstream of all local outlets**. This is not true biological site-specific WWTP dose or monitored annual N/P tonnes.
- Since 88–91% already have at least one original registered facility, **binary exposed/unexposed is a poor discriminating macro predictor**; choose source-modelled **load intensity** with area/upstream flow/dilution and registration/effluent-quality flags. Do not interpret Not Reported status as operational.
- No ecological association or response was examined during this source stage.

[Real HydroWASTE 58K to monitored 1,219 H7 true-network original source SUCCESS and CSV](https://github.com/zuizui0223/gesui/actions/runs/43199920991). Original 1,219 pressure CSV SHA256 **`a311127b04e0a081a32a764b183240c8642aa8e4e71f81016bde9a56e537fcac`**.

## Source-correct one-table comparison, not a made-up ecological cascade

The common `HYBAS_ID_L07, PERIOD` master must hold:
- Static original `HYDROWASTE_UPSTREAM_TO_BASIN_OUTLET_MODELLED_WASTE_DIS` **distinct** from temporally observed 2018/2022 EU reported WWTP N/P and genuine water TN/TP measured US2018–19.
- Original EEA real benthic EQR 2018 /2022 on identical native basin-period, US real EPT fraction in P1 only; do not equate EQR and EPT measurement scales.
- TRUE basin-polygon, time-filtered GBIF Plantae and Lycaenidae source counts **not** old 1,055-port 5km and 2011–24 aggregated GBIF. The new GBIF count is recording opportunity only; true plant richness requires speciesKey facet with rarefaction and Lycaenidae biological occurrence requires true butterfly family denominator, not Animalia as "butterfly" substitute.
- Country/assessment method, river size and local basin area, cropland, urbanization and recording opportunity controls before interpreting environmental slopes. HydroATLAS global full data ~2.5–4GB; do not falsely claim downloaded landcover/streamflow.
- Only compare **separate layer associations** with common exposure where real same-basin/period support exists. No full chain from WWTP → larvae → adult emergence → predator → pollination until measured shared-source and temporal alignment.

**Status at report authoring:** Stage0 original 1,219 geometry, worldwide HydroWASTE source and Stage1 HydroWASTE L7 outlet pressures DONE. Entire real GBIF 7,314 original-source full collection NOT YET VERIFIED; keep status visible without assuming future success. Subsequent same-period macro ecological fits and true species richness unestimated.

All changes remain Draft stacked on PR #62. No `gesui/main` or `chocho` merge, no eBird, no downstream pollutant/food-web benefit causality claim.
