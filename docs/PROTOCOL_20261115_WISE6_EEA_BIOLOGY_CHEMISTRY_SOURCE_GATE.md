# Source-only EEA 2026 WISE6 water-quality × EEA biological site matching

**Frozen 15 November 2026 before inspecting any new concentration values; source-feasibility experiment, NOT a sewage effect estimate.**

## Strict original grains
- Same 454 existing H7×2017–19/2020–22 biological cohort (232 H7); source 3,525 actual EEA station-year records (2,485 distinct native station IDs) from original GitHub Actions #43165326244.
- Original physical station matched on BOTH monitoringSiteIdentifier and monitoringSiteIdentifierScheme, and same calendar year (2018 or 2022). Do not replace a missing station measurement with an averaged H7 centroid, nearest neighbor, modelled waterbody chemistry, or another year.
- EEA Waterbase Water Quality ICM, 2026 published July 7, 2026, coverage through 2025; annual site AggregatedData and raw DisaggregatedData are distinct populations. Prefer usable disaggregated observations if available; do not double-count the aggregate and disaggregated sources.
- Only confirmed actual **TOTAL N (TN)** or **TOTAL P (TP)** in water; report determinand code, native units, site ID, sample count, LOQ flags and year. NO substitution with NO3, NO2 or PO4 fractions.
- Minimum for considering a preliminary source-load to measured-receptor-concentration analysis: at least 30 source-positive H7×period units with >=20 independent H7, >=10 distinct H6 parents, and >=3 countries containing a same-year same-station measured nutrient. This is source-support, not statistical power for mediation.

## Existing prior-source eligibility problem
- Reported upstream positive WWTP N/P original 56/454 basin-periods, 38 H7 and 31 true H6. The other 398 original periods are **no matched positive reported source**, not verified sewage-free controls.
- The 2018 and 2022 annual outfall-mass reports do not exactly cover all 2017–19 and 2020–22 GBIF source dates.
- Mass/mean modelled flow is a hypothetical full-delivery mg/L equivalent, not measured TN/TP or an effluent fraction; it cannot independently attribute agricultural sources or river retention.

## Source-only initial preflight
Call the official EEA Discodata SQL endpoint for actual WISE6 2026 schema/version. Download no 4.3GB remote file. Save response source SHA, exact query and errors. No biological response value is read. If no trusted table is accessible, stop and report access failure without interpreting it as zero measured concentration.

**No causal SEM and no invented source measurements:** even an exact same-river WQ record does not directly show delivered WWTP-specific N/P fraction in individual butterfly feeding sites. That would require outfall operation, dilution, actual discharge locations, matched controls and riverbank surveys.
