# Macro research: nutrient enrichment, periphyton biomass, benthic biodiversity, and service proxies

Status: **source-only** audit. New empirical coefficients are NOT present until all QC gates pass.
New sources: EPA National Rivers and Streams Assessment 2018–2019 and 2023–2024 official CSVs.
Pre-existing EPT-vs-TN/TP evidence: [gesui #58](https://github.com/zuizui0223/gesui/pull/58).
Do not merge with earlier European EEA I2M2 river EQR or conflate nutrient with wastewater effluent.

## Frozen analysis decisions (before inspecting biomass outcomes)
- **Primary analysis grain:** EPA native UID (visit), not GBIF photo or city/WWTP buffer. Preserve SOURCE-only uniqueness; separate visit and physical SITE_ID. Avoid post-outcome site selection and repeated-visit pseudoreplication.
- **Primary biological responses, separately measured:** periphyton ash-free dry mass (AFDM) if the exact publisher label, QC flag and units justify this identification; periphyton chlorophyll-a if identified; benthic EPT counted-individual fraction from authentic count/taxonomy; benthic richness/MMI as a distinct composition endpoint.
- **Primary gradient:** same-UID measured water total N and total P, units verified against original metadata. Compare log positive concentrations. No substitutions (DIN for TN, PO4 for TP, chlorophyll for biomass), no zeros imputed for absent measurements.
- **Common covariance structure:** state/eco region, coordinates, river order/size, year, conductivity, upstream cropland and local physical-habitat attributes *only if truly observed for both study periods*. HUC8 grouped holdouts; compare nutrient-added MSE separately for AFDM, chlorophyll-a, EPT share and richness.
- **Survey history:** do not interpret 2018/19 vs 2023/24 as paired temporal changes without explicit repeated native field SITE_ID and method equivalence; first analyze eras separately, then test external time transfer.
- **Candidate conflict test:** does measured nutrient increase correlate with greater biomass but lower EPT fraction/richness? Joint signed effects are the question, not a guaranteed positive result. Pre-fix separate coefficients and uncertainty; avoid choosing nutrients/response transformations based on effect signs.
- **Actual wastewater exposure:** no WWTP claim without year-valid measured load and source-to-visited-reach directed routing, verified upstream placement and dilution. Unverified presence or upstream-to-basin outlet is not delivered nutrient mass.
- **Ecosystem services:** AFDM ≠ primary productivity rate; benthic EPT share ≠ adult emergence flux; GBIF butterflies ≠ pollination. Keep nutrient recycling, insect export biomass, trophic feeding and realised pollination as additional *unmeasured* endpoints until relevant source-compatible datasets are joined.
- Source QA must report all original table identifiers, row grain, counts, SHA256, field dictionaries, core-UID overlap, flag/unit availability before any numerical association. Source failure is a STOP not a biological null result.

## Data source
Official EPA [NARS data listing](https://www.epa.gov/national-aquatic-resource-surveys/data-national-aquatic-resource-surveys), published source files under 2023-01, 2021-04 and 2026-06.
No eBird. No non-public sample source. Original CSVs remain with EPA; workflow artifacts contain only schema + source identity/availability metadata.
