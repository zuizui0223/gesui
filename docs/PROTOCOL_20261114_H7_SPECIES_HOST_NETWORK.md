# H7 species-level candidate larval-host network: source and analysis contract

2026-11-14. Continuation of Draft #85. Exploratory posthoc relative to earlier known wastewater/species-richness results; not a preregistration or observed food-web experiment.

## Question and scope
Join butterfly and plant identities in the original 234 H7 catchments and two periods (2017–19 and 2020–22), then test whether modelled wastewater geography adds predictive information about catalogue-supported local host opportunities, beyond species richness, recording effort and land use. Do not shrink to one butterfly family; retain all five original families. Do not infer feeding, reproduction, pollination, biomass or sewage chemistry.

## Source identities
Use exactly recovered artifact 14113695873, run 43593227162; ZIP SHA256 4b37c741d1bd28c178b0cef61c9457af0dd30797c173f9384ebe1d1d0036d07c. It contains 2808 H7×period×group source rows, 234 H7, 468 period units, 303 butterfly speciesKeys and 11613 plant speciesKeys. The original 454-row ecology/exposure cohort stays fixed; source coverage ≥80% in all six groups is the pre-existing post-source primary diagnostic, full cohort and ≥90% are sensitivities.

Ferrer-Paris 2014 host catalogue, DOI 10.6084/m9.figshare.1168861.v1, exact Git blob ea93a5c0fdb17e361a9b6a9637a497a12e90493c and 43645 original rows, is the reproducible initial host source. It shares historical HOSTS origins with chocho and is not an independent replicate of chocho. European host checklist DOI 10.5061/dryad.1vhhmgr12 / article 10.1002/ece3.10834 is a regional cross-check candidate; availability and actual use must be reported separately.

GBIF legacy numeric speciesKey semantics are retained, not silently mixed with Catalogue of Life string identifiers. Retrieve each recorded butterfly's canonical name and species-level synonyms from official GBIF. Match host catalogue butterfly names exactly to these unambiguous species aliases, with source family agreement. Strict binomial host names are resolved through official GBIF /v1/species/match, requiring EXACT match, SPECIES rank, Plantae and unambiguous accepted numeric key. Fuzzy, higher-rank, unknown and ambiguous matches stay excluded/unknown, not no-host evidence. Preserve raw JSON responses, URLs, acquisition timestamps and SHA256. Source syntax/taxonomy filtering must not use wastewater slopes.

## Outcomes before viewing new effects
1. Catalogue coverage by family: fraction of detected butterfly species and source records with at least one verified catalogue host; report unresolved hosts and synonyms. No all-traits-complete taxon requirement.
2. Candidate links at identical H7×period: unique verified catalogue butterfly–plant key pairs with both keys indexed locally. Co-recorded links are opportunities, not realized interactions.
3. Recording-effort-standardized expected candidate links: sum A_bp × Pr(b in 20 butterfly records) × Pr(p in 50 plant records), computed with exact hypergeometric inclusion probabilities. All counts are GBIF records, not censused individuals.
4. Candidate host support: among catalogue-mapped butterfly species standardized to 20 records, expected fraction with ≥1 of their catalogue hosts in 50 local plant records. For each butterfly use the summed local counts of its known hosts in the hypergeometric inclusion formula. Denominator excludes catalogue-unknown butterflies and reports that coverage explicitly.

## Richness and generic host availability controls
Primary baseline retains original country, period, EEA method, geography/area, recording opportunities, stream flow, upstream crop/urban land use; add butterfly S20, plant S50 and catalogue coverage before adding static modelled wastewater W. Compare the same eligible rows in every model. No model covariate selection based on effect sign.

Null sensitivity for expected links: permute plant identities in the host catalogue within plant family and broad occurrence-frequency strata, holding each butterfly's host count, the locally recorded communities, and rarefaction probabilities fixed. Report singleton strata and effective randomization. This tests taxon-specific matching beyond those margins, not a causal sewage pathway.

Use original H6-parent-disjoint 5fold partition, H6-equal Ridge alpha10 training and error evaluation, training-only standardization/categories, and fixed-error H6 bootstrap descriptive intervals. Country holdout uses no country effects and keeps a catchment's periods together. Report effect directions separately from signed prediction improvement. Endpoints and sensitivities are dependent exploratory checks, not independent confirmations or meta-analysis.

## Implementation and verification
Create scripts/gesui_h7_host_source_20261114.py for source-only collection and taxonomy gates; preserve full input sources as an Actions artifact. Create offline scripts and tests for hypergeometric inclusion, alias ambiguity, source identity joins, link construction, null margins and H6 leakage. Recompute primary model predictions with a separate numerical implementation. Save a reproducible archive and source-provenance manifest. Keep main/chocho untouched and use a new unmerged draft branch.
