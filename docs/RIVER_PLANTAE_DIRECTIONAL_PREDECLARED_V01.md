# Directional river Plantae pilot — frozen before inspecting GBIF responses (2026-11-07)

## Geometry and response units
The real HydroWASTE/HydroRIVERS source run selected 16 outfalls (four each AUT/CZE/DEU/DNK; high/low water fraction) without querying biological outcomes. Each has three different river reaches ordered **U2 → U1 → outfall-associated reach → D1** by the NEXT_DOWN river network. The new 48 biological sites are fixed at the U2/U1/D1 center coordinates, radius 1 km.

## Fixed biological data source
GBIF kingdom Plantae (key 6), occurrence status PRESENT, coordinate available, no geospatial issue, 2011–2024 inclusive; paginate speciesKey facets for every site. Save original speciesKey×site counts as well as all-Plantae record counts. These are indexed plant observations, not abundance or true absence, and may include non-vascular taxa.

## Predeclared responses and gates
- Primary: all three sites require >=30 GBIF records resolved to speciesKey and >=6 distinct plant species; max/min species-key-record counts across the three sites must be <=10 to avoid extreme effort imbalance.
- At each site calculate exact expected recorded species after rarefaction to 30 species-key occurrences. Let Y = [log(1+S_D1)−log(1+S_U1)] − [log(1+S_U1)−log(1+S_U2)]. The first change spans the outfall-associated reach; the second is a pre-outfall upstream comparison. Y is a spatial observation index, not a causal difference-in-differences.
- If >=8 outfalls from >=3 countries qualify, describe the Spearman rank correlation of Y with original modeled effluent-water fraction, otherwise mark support insufficient and do not change the predeclared sample.
- Sensitivity A: 30-record rarefaction without the <=10 effort-imbalance limit. Sensitivity B: 10-record rarefaction with >=10 species-key records and >=4 distinct species at each site. The low-depth estimate must be labelled a weakly sampled species index.
- Companion: unrarefied GBIF speciesKey Jaccard dissimilarity(U1,D1) minus dissimilarity(U2,U1), with explicit effort and unequal segment-distance caveats.
- Source high/low contrasts are descriptive, not causal or statistically exchangeable. No post-hoc replacement of sites or thresholds is allowed.

## Limitations
River segment centers are not equidistant; for several outfalls D1 lies much farther from U1 than U1 from U2. Urbanization, geology, plant-recording effort and other untreated upstream WWTPs have not been controlled. There is no independently measured effluent N/P, downstream uptake, discharge seasonality or verified taxon-specific plant dispersal. Geographic connectivity does not prove pollutant or propagule transport.

The underlying Europe HydroRIVERS v1-derived GeoParquet is a source-hash-verified third-party conversion, NOT yet independently matched to an official release. Its 4,713/4,713 original HydroWASTE outfall HYRIV_ID matches, and <=47 m observed source-line offsets at the 16 selected sites, are positive source-compatibility checks but are not a full independent topology/field validation.

If data are sparse, publish the complete source quality table and null/insufficient-support result without inventing positive patterns.
