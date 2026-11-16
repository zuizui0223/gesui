# NRSA 2023–2024 measured physical habitat × biofilm nutrient macro study

Status: **source-only, before new physical environment response coefficients**. This extends [PR #75](https://github.com/zuizui0223/gesui/pull/75) and keeps PR #74 wastewater N/P source gate failures unchanged.

Use official EPA NRSA2324 Physical Habitat Metrics and metadata linked by [EPA NARS Data](https://www.epa.gov/national-aquatic-resource-surveys/data-national-aquatic-resource-surveys): `nrsa2324_physicalhabitat_metrics.csv`, `nrsa2324_physicalhabitat_metrics_metadata.csv`. Previous NRSA 2018–19 PHab publisher metadata specifies existing field candidates: XCDENMID canopy density at mid-channel; XPCAN canopy presence; XGVEG riparian vegetated ground layer; XFC_OHV overhanging vegetation; XEMBED streambed embeddedness; XCEMBED center-channel embeddedness; LitSB1_SF/FN fine sediment; XWIDTH mean wetted width. Source-only audit checks current genuine headers and definitions. No field will be swapped based on nutrient/biofilm sign.

**Outcome and comparison originally fixed:** the same 2023–24 pre-PHab n2,129 original native UID site/water/AFDM/Chla opportunity sample; strict n1,828 original physiological quality panel [PR #75](https://github.com/zuizui0223/gesui/pull/75), 946 HUC8, 48 states. Respondents are original log10 areal AFDM, log10 areal Chla, and their algebraic log difference Chla/AFDM. No MMI or EPT from unavailable publisher 2023/24 benthic files, no emergent biomass or ecosystem service inferred.

**Primary environmental falsification:** compare the original baseline geographic+state+eco9+stream order+year+sampling season+conductivity, plus river log TN/TP, with adding predefined field canopy+sediment+width. Crucially evaluate incremental nutrients **conditional on same field environment**, and signed TP slopes, full HUC8 heldout folds, and original n changes due to physical metric missingness. Do not assume canopy is direct RSDS measurement or fine substrate is causal toxic stress. Distinguish confounding from mechanism.

**Original-source gate:** PHab file genuine, source SHA archived, unique UID or uniquely mapped exact site-date, >=400 original source IDs, >=100 HUC8, >=20 states; source-first measurement availability and complete case transparency. If missingness or wrong units cause failure, STOP without changing thresholds. Numeric environmental effects only after source interpretation/measurement QA.

Keep all old positive and null findings, no changed main branch, no eBird.
