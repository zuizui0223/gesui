# Source-gate v0.1 — pre-plant-response freezing

## Scientific unit
One WWTP estimated outfall's **directed HydroRIVERS chain**, selected without any GBIF or vegetation response. Three planned independent plant-observation sites: **U2 → U1 → [outfall-associated reach] → D1**, where U2 and U1 are upstream river reaches, and D1 is the directly downstream reach.

HydroWASTE HYRIV_ID supplies the known outfall-associated river reach; HydroRIVERS HYRIV_ID / NEXT_DOWN define the actual directed network. Upstream is chosen deterministically by greater upstream drainage area. The wastewater facility coordinates and **estimated** outfall coordinates remain different and are not treated as field-confirmed points.

## Inclusion, exposure, spatial controls
- Pinned original 58,502 HydroWASTE CSV, exact Git blob 0506f873f6541c5f324478894a9731290d722bce.
- Official HydroRIVERS version 1 Europe/Middle East GIS shapefile, source ZIP SHA256 recorded at run time.
- Source geography BBOX −12° to 26° E, 44° to 64° N; 11 listed European countries.
- Active or not-reported status (explicitly retired/planned facilities excluded); QUAL_LOC ≤2; QUAL_WASTE ≤3; effluent water flow >0, model receiving river flow ≥0.5m³/s, coastline distance >10km.
- Four distinct consecutive reaches required (U2, U1, outfall reach, D1). Centers for three biological sites >2.2km apart (1km radius buffers cannot overlap) and within 35km maximum separation. A **reaches-to-GBIF** workflow must separately assess exact bank intersection, location uncertainty and sampling effort.
- Potential effluent water fraction = WWTP WASTE_DIS/(WASTE_DIS+86400×RIVER_DIS). **No N/P flux, actual effluent concentration, hydrologic seasonality, water quality or plant outcome is implied**.
- Select at most two sites per high/low **within-country** source-fraction quartile from each country with ≥8 candidate facilities, maintaining ≥45km distance from other selected outfalls. Selection is **before querying plant records** and not optimized on biological response. Pilot sample is an opportunity inventory, not a prevalence estimator.

## Biologically meaningful comparisons (later stage, not run in this gate)
1. U1 vs D1 comparison across the outfall-associated reach, versus the upstream U2 vs U1 background turnover, using the same plant species index.
2. Contrast above incremental turnover between high/low source-only water-fraction strata **with spatial, catchment, urban, river-order and botanical effort controls**. Modelled flow fractions are year-agnostic; do not infer seasonality or stable dry-season flow.
3. GBIF records are detection/cataloguing not abundance. Plant groups need taxonomic validation to distinguish vascular plants, macrophytes and terrestrial plants. If source coverage is weak, report failure and do not fit response models.
4. Upstream may have other WWTPs; a proper cumulative-upstream-outfall loading audit is required before identifying any sole-facility effect.

## Artifacts
The source run outputs a JSON source quality/coverage table and selected pre-outcome outfalls with coordinates and reach IDs. Large GIS input ZIP and raw HydroWASTE CSV are not vendored. This gate performs **real source geometry**, but **no plant ecological correlation or causal inference**.
