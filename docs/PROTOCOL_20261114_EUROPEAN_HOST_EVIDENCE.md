# European evidence-qualified host opportunities — exploratory continuation of PR86

The original 454 H7×period /232 H7 /90 H6 cohort is unchanged. Previous PR86 results are known; this is POSTHOC, not preregistered or independent ecological replication. New evidence-specific numeric outcomes have not yet been fitted when this note is committed.

## Authentic source and material limitation
Clarke 2024, DOI 10.1002/ece3.10834, original Appendix S1 (5589 rows). The Dryad API returned401 and public route failed; no source table was manufactured from those responses. Official NLM 2026 anonymous S3 objects were successfully retrieved, article identity checked, every object MD5 verified, SHA256 archived. Source run43600676744, artifact14115848892, ZIP SHA256 6c58ca11a4f1c5c9f138b33990f9b4ce3d4eee4adfe22ae671bc9a31129ca8a5.

The appendix has accepted butterfly/plant, evidence, countries, questionable_countries and references. The article explicitly defines evidence as the strongest evidence across ALL references for that taxon pair. Therefore a country mention plus a global Larvae/Proven rating DOES NOT establish observed larvae in that country. Raw reference×country×stage rows from Dryad remain unavailable. No country-confirmed-larval tier is fitted.

## Four definitions, same original biodiversity source
1. EU_reported: Proven/Larvae/Eggs/Claim only and at least one nonquestionable country/region report.
2. EU_larval: Proven/Larvae with a nonquestionable country/region report somewhere.
3. country_reported: first definition plus EXACT original panel country in the source country list. Broad regions, Europe and islands are not expanded to countries. Questionable country overrides a conflicting positive mention conservatively.
4. country_larval_global: second definition plus same-country report; this is NOT evidence localized to that country.

Species-rank gate: exact accepted/canonical butterfly name takes priority; a unique alias is accepted only without ambiguity. Plant names projected to species binomials; genus-only, aggregates and malformed names are not expanded. Exact GBIF SPECIES/Plantae/ACCEPTED-or-SYNONYM matching only, no fuzzy/higher-rank acceptance. Infraspecific-to-species projection flagged and removed in sensitivity. Unresolved names and unavailable catalogue coverage never become biological absences.

## Model and network checks
Use identical original GBIF indexed records, not field individuals, for 20-butterfly-record and50-plant-record hypergeometric expected links/support. All four tiers share the same resolved European plant universe. Catalogue coverage is an explicit predictor; when no recorded butterflies have known hosts, support/link outcomes are missing, not zero. Source-record identity coverage80 retains the originally defined390 rows where estimable, with all454 sensitivity. Report any further missingness.

For each tier use499 plant-column permutations within plant family×source occurrence-frequency bin, plus the fixed six-cluster H6 plant-distribution grouping as a stronger spatial sensitivity. Retain butterfly degrees within each country matrix and report movable columns. Both are potential co-recording nulls, not observed feeding or causality tests.

Ridge10 uses original observation effort, land cover, geography, species richness and catalogue coverage; train-only scaling/categories, H6-equal fitting, fixed H6-hash 5fold. Country holdout excludes country features and purges ALL overlapping H6. Report signed W increments, fixed-error3000 H6 bootstrap, country-equal/row-weighted diagnostics; W is modelled static outlet source geography, not delivered effluent chemistry. No optimizing tiers for positive signs.

Independent SciPy row probability arithmetic and sklearn weighted-Ridge verification; source-only, genus, ambiguity, questionable geography and unknown-host tests written before functions. Accepted-only butterflies and species-rank-only source sensitivities retained. No merge to main/chocho, no eBird/TinyFish. Numerical analysis is local unless a separate successful numerical workflow is actually verified.
