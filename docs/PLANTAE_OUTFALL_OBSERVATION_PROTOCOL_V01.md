# Plant-first source-gated pilot v0.1 — noncausal GBIF recording index

## Hypothesis and observation grain
Do areas at HydroWASTE **estimated wastewater outfalls** have different GBIF-indexed Plantae recorded-species diversity compared with nearby **WWTP coordinates**, and does the outfall-minus-WWTP contrast vary across source-only modeled wastewater/river-flow fractions?

Important: the WWTP is **not an upstream hydrological control**. This is exploratory record-index quality/feasibility work, not evidence for downstream effect, wastewater nutrient enrichment or plant production.

## Source-locked selection made before any plant data
- HydroWASTE 58,502 original observations, git blob 0506f873f6541c5f324478894a9731290d722bce.
- Germany, France, UK, Benelux, Alpine, CEE and Sweden/Denmark, same 11 fixed country ISO3 codes as the directed pilot. Known status active or missing; remove retired/proposed/under construction; location quality ≤2, wastewater flow quality ≤3, modeled WWTP effluent >0, receiving river natural flow ≥0.5 m³/s, HYRIV_ID present, outfalls >10km from coast.
- Facility-to-outfall straight-line distance **4.25–35km**, so a 2-km radius around each coordinate has no overlap. This selects on source geography, not on GBIF site count.
- Nominal fraction = WASTE_DIS (m³/day) / (WASTE_DIS + RIVER_DIS (m³/s) × 86,400). It is not dry-season flow fraction or nitrogen load.
- Within each country require ≥6 source-eligible candidate facilities. Choose top/bottom quartile source fraction, then **one matching high–low pair** prioritized by similarity of receiving natural river discharge (ratio ≤6) and within 30–350km, deterministic ID tie-break. Biology not queried until after source pair selection is locked in this code and receipt.

## Frozen GBIF survey
- 2011–2024 **all Plantae** (kingdomKey 6), accepted speciesKey facets, 2km from estimated outfall and from WWTP. No butterfly-host selection. GBIF indexed observations, not plant abundance or truly complete vascular flora.
- coordinate present, no geospatial issue, occurrence status present. Source 2011–2024 is deliberately modern retrospective because HydroWASTE has no historical facility-operating year.
- Facet pagination for all species keys; verify species record sum ≤ all indexed plant records. Flag records without species ID, rather than pretending they're absent.
- For comparability, expected species richness in a fixed **30 species-keyed GBIF occurrence records** via exact hypergeometric rarefaction, requiring ≥30 species-keyed records and ≥6 different species at each site.

## Exploratory outcomes and limitations
- Each facility: log(1+rarefied outfall richness) minus log(1+rarefied WWTP richness).
- Source high/low matched within country: difference of those differences. Report denominators, missingness and full negative/positive results.
- If ≥8 qualified facilities, report descriptive Spearman rank correlation between source-fraction and the outfall-minus-WWTP recording contrast; otherwise explicitly mark insufficient.
- Do **not** regard high fraction correlations as proof of N loading, as outfalls/WWTPs differ in riverine vs urban habitats and source record effort.
- Unlike a proper H1 test, this does not yet include upstream/downstream hydrologic topology, independent agriculture/urban/census controls or biological abundance.

Data acquisition is one small, bounded exploratory GBIF facet pass (≤2 coordinates per selected WWTP), **not a repeat of the 1,596 old chocho port host requests**. Artifacts retain original site source-facet counts and analytic records for a 30-day window.
