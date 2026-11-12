# First actual comparable biodiversity endpoints on repeat Euro-river catchments: prespecified native source cohort

**Frozen scope created 2026-11-12, after prior source-only 1219×2 GBIF COUNT collection and its support-tier audit, but BEFORE looking at newly collected plant speciesKey facets or four other butterfly-family count values.** This is source-first *posthoc* exploratory research, not an independently preregistered discovery.

## The repeat source cohort: 234 actual HydroBASINS level-7 basins, NOT 218 from new GBIF success tier

Use original [EEA2018/EEA2022/EPA source-only Level7 assignment](https://github.com/zuizui0223/gesui/actions/runs/43165326244), immutable SHA256 `a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89` and the original 1,219 basin WKT source SHA256 `47f597e0a798cd6155ee9fb69ce8916ed5629951d046d114feb9f667963310e4`.

Select **all** H7 IDs that have ≥1 EEA original benthic-EQR monitoring station in **both** 2018 and 2022, exactly **234 unique basins**. No requirements about butterfly/plant record count at cohort selection, no choice based on EQR trend. Each selected basin has two fixed source periods (2017–19 and 2020–22) and same genuine full polygon components; **468 basin×periods**.

## Stage 0, before new macro biodiversity response

### A. Original GBIF taxonomy identities

Use original verified legacy GBIF Lycaenidae **5473**; verify **Nymphalidae, Pieridae, Papilionidae, Hesperiidae** accepted FAMILY keys live from GBIF species name matching, checking `rank=FAMILY`, exact accepted name and original GBIF Backbone `checklistKey=d7dddbf4-2cf0-4f39-9b2a-bb099caae36c`. No guessed or iNaturalist family IDs. Preserve exact taxonomy response URL and SHA before any new family occurrence count.

### B. Original same-basin-period plants (true species facets)

Original Plantae `taxonKey=6`, WKT original polygons, fixed source dates, `limit=0`, `facet=speciesKey`, `facetLimit=1000`/facetOffset pagination, true original per-species counts; maximum **20,000 facet species** per basin×period. If facet reaches 20,000 before source exhaustion, source result is **coverage-incomplete NA**, not an artificially capped richness value. Count index species-key facet coverage vs original total occurrence index count; require ≥70% facet record coverage, ≥50 accepted species-indexed occurrences, ≥10 original species keys before 50-draw rarefaction. No speciesKey included == coverage unknown, not 0 plant species. Never call the GBIF indexed speciesKeys field field-observed abundance.

Rarefied plant richness for **exactly 50 original species-assigned indexed records** computed without pseudo-duplicates using log Gamma:

[
 S_{50} = sum_i\left(1-\frac{\binom{N-n_i}{50}}{\binom{N}{50}}\right)
]

where `N=sum speciesKey facet count`, `n_i` per unique indexed speciesKey, and `N≥50`. Numerical zero when N≥50 is not expected. Record observed species-count and coverage as diagnostics.

### C. Original 5-butterfly-family source denominator

Query genuine fixed H7×period occurrence index counts for **Lycaenidae + Nymphalidae + Pieridae + Papilionidae + Hesperiidae** under a single documented GBIF Backbone checklist and same source flags (coordinate, presence, geospatial issue false). Original Lycaenidae source count from prior 7,314 requests retained as independent check, no source filters altered. Sum five family raw count totals as **recorded butterfly opportunity denominator**; use Lycaenidae share of five families, and original animal-source count to control for broad observer intensity. Require ≥20 butterfly records across five families and original source no HTTP failures before attempting a same-period butterfly-family composition model.

Five families are a **proxy** for butterfly guild, not a perfect taxon universe, abundance, realized pollination or butterfly-host use. Do not claim butterfly species richness from family counts.

### D. Gate before ecological EQR cross-layer models

- ≥300/468 original H7×periods with original successful GBIF true species facets (Plantae 50-rare) AND all five family source counts;
- ≥150 distinct H7 basins with both original periods valid;
- ≥12 countries in valid records and ≥80 original EU repeated EQR physical stations/method-comparable original site/time pairs inside candidate catchments (must be verified via original EU 2018/22 measured station files);
- original per-area land use, original country/method biological EQR, native survey report and other covariates separated. Do not use same response data to select positive study years, region boundaries or taxon family denominator.

If any source gate fails, record explicit STOP without EQR associations or claims of water→plant→butterfly service. If passes, source freeze complete biological candidate panel, and only then exploratory response models.

## First macro biological model (conditional)

One row per **H7 basin×period**, with real country/method/site biological EQR measured in 2018 and 2022. Restrict primary within-basin temporal analysis to repeated **original same monitoring stations and compatible biological assessment methods**, equal-mean across physical station within H7, not a new station mix. Distinct methods within same country require method effect. Plant `S50` and five-family Lycaenidae share are observational *biodiversity proxies*, not ecosystem services. Direct H7 static HydroWASTE (local and upstream outlet) cannot explain changes over time without environmental interaction — fit as spatial effect/context, never time-varying real effluent.

Exploratory model: EQR measured-site change vs change in effort-standardized plant richness and Lycaenidae family share with country FE and changes in original Animalia recording effort; robust H7 block bootstrap, country holdout. **Do not treat three correlated original source observations as three independent biological experiments**, nor infer serial causal mediation. True pollution-to-outlet→invertebrate→flower→pollinator service remains unobserved.

Do not use eBird, don't overwrite earlier PRs, no `main` merge, no false global biodiversity claim — this initial repeat panel is EU only within global shared geography.
