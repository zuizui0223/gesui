# EEA/UWWTD 2022 public N/P discharge source gate — frozen before recipient chemistry (2026-11-07)

**Scientific motive:** The earlier fixed European river macro pilot uses `HydroWASTE.WASTE_DIS` (modeled treated water volume) as an upstream source proxy and did **not** predict real TP concentrations robustly. This new independent source uses EEA's **reported nutrient mass discharge** directly; it is NOT a posthoc search for a positive model coefficient.

## Public and authentic source

- EEA [UWWTD Waterbase 2022–2023 final public data, published 2026-03-03](https://www.eea.europa.eu/en/datahub/datahubitem-view/6244937d-1c2c-47f5-bdf1-33ca01ff1715).
- Official ArcGIS REST [TreatmentPlants_2022 public](https://water.discomap.eea.europa.eu/arcgis/rest/services/WISE_UWWTD/TreatmentPlants_2022/MapServer/0) and [DischargePoints_2022 public](https://water.discomap.eea.europa.eu/arcgis/rest/services/WISE_UWWTD/DischargePoints_2022/MapServer/0).
- EEA 2026 [UWWTD Article 15 data dictionary](https://dd.eionet.europa.eu/datasets/latest/UWWTDArt15) defines `uwwNDischargeMeasured`, `uwwPDischargeMeasured`, calculated/estimated alternatives, treatment and plant identifiers. **Only reported measured annual nutrient loads** enter the primary source table; never sum calculated and estimated as measured.
- The EEA public feature layer provides **22,554** treatment plants, **5,515** with positive reported measured N discharge, **6,529** positive P, **5,492** positive **both**, and **24,416** listed discharge points. These are *API source coverage*, not linked monitored river concentrations. [Source receipt](https://github.com/zuizui0223/gesui/actions/runs/42373101790).
- Confirm with dictionary the actual nutrient mass units before derived concentration approximations; a field name alone does not justify units. Distinguish *reported at plant* from *measured at recipient monitoring station*.

## Pre-outcome data linkage

1. Obtain exact REST paging with `OBJECTID` order, total count gate, request receipts and no outcome-dependent filtering. Standardize `repCode` such as `AT_UWWTD_2022` and source `uwwCode`; restrict to unique country+plant IDs to avoid joining unrelated countries.
2. Link 5,492 both-positive-measured-load facilities to discharge points by `repCode/uwwCode`. Main gate only one unambiguous discharge point per facility. Multiple points must be reported as ambiguous; do not allocate the whole WWTP N/P mass to each.
3. Snap coordinates to source-hash-locked HydroRIVERS v1-derived EU GeoParquet within projected **500m** in EPSG:3035. Exclude sea/lake/surface-water incompatibilities where officially coded. The source river mirror is not independently verified against the official original.
4. Directed network reach-to-reach source accumulation: **≤50 km** upstream of fixed original TP/TN monitoring station reach; exclude same reach because station and discharge ordering within a reach is unknown. Distinguish N and P *measured source mass*; explicitly compare downstream WWTP-load negative control.
5. Freeze which existing Peters 2026 v5 monitoring stations have source support before using *2022 measured in-river concentrations*. Keep 2022 source–response temporal alignment separate from 2011–2024 historical observational averages. The 48 earlier plant sites remain untouched.

## Next scientific test if source support passes

For each nutrient, evaluate 2022 station median of original observed river concentrations only when a station has ≥4 positive source observations in at least three distinct months in **2022**. Correctly convert measured EEA source mass to a conservative potential fully mixed concentration based on HydroRIVERS mean Q, with the understanding that removal, retention, tributaries and mixing make it only a proxy. Predefine a spatial block heldout model with river flow, drainage, geography, available cropland/population, **actual upstream reported load**, and a downstream-source negative control; report same-denominator errors and exact source-quality failures. This is a correlation of 2022 reported source load with 2022 river monitoring, never proof the effluent produced the observed river concentration.

**Future publication gate:** validate EEA facility/discharge historic availability and coordinate precision, join to actual monitored receiving water reaches, and quantify crop/urban and other sources. Do not claim N/P discharged by plant if public API only reports calculated or estimated versions.
