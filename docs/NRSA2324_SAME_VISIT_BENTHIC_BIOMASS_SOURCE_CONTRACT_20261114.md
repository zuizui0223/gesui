# Same-visit NRSA2023–24 biodiversity × biofilm mass — source-first extension

**2026-11-14. Draft, no results examined before the original source census.** This branch extends [gesui #73](https://github.com/zuizui0223/gesui/pull/73) and [#74](https://github.com/zuizui0223/gesui/pull/74) without changing their original coefficients, signs, support limits or unmerged status.

## Scientific question

Do *measured river* nutrient concentrations predict both areal periphyton biomass (AFDM; mg/cm²), areal chlorophyll-a (µg/cm²), and standardized **benthic community quality** at the **same original 2023–24 EPA NRSA native survey visit**? If available, can benthic raw counts supply genuine EPT fractional composition and taxon richness instead of merely an MMI? Where do quality indices, mass and pigment give contrasting responses?

These are physical river nutrient comparisons, NOT effluent-specific causal impacts. EPA historical fiscal DMR same-source bilateral N/P failed predeclared gates in #74, so no retrospective reclassification of NRSA river TN/TP as wastewater-originated.

## Original source choices before inspecting ecological responses

- **Cohort**: EXACT original 2023–24 source-only 2,129 common native SiteInfo, WaterChem, AFDM, periphyton-Chla UIDs frozen in original [PR #73 source QA #43417864266](https://github.com/zuizui0223/gesui/actions/runs/43417864266). All attempted official benthic files are inspected on this complete set, not outcome-selected.
- **Official benthic sources**: EPA NRSA 2324 Benthic Count, Benthic Metrics, Benthic MMI, Benthic Taxa List and companion metadata, *from the EPA official table's exact URLs* rather than guessed alternatives. Earlier candidate Count and Metrics links returned HTTP404; do not mislabel missing as zero. The source gate examines all endpoints' HTTP status, publisher bytes SHA256, native header and one-to-one site ID overlap. No positive/negative ecological effect can be inferred from successful downloads alone.
- **ID match**: native `UID` first. If missing or inconsistent, require unique `SITE_ID + VISIT_NO + DATE_COL` site-date identity, independently checked against official SiteInfo and frozen source UID manifest. Never arbitrarily pair two dates or duplicate taxonomic counts.
- **Benthic source grain**: distinguish one `UID` MMI row, repeated taxon rows, substrate-specific sample type and separate lab duplicates. Never sum records across `SAMPLE_TYPE` without original method definition. `MMI` is multi-metric biological integrity *index*, not Shannon, richness, individual count, biomass or comparable to EU EQR.
- **Minimum for an explanatory model**: at least 400 verified 2023–24 same visit joined benthic+AFDM+Chla+TN/TP, 100 distinct HUC8, 20 US states, and unambiguous one-to-one biological score or standard benthic-count definition. If below, stop without changing eligibility. Counts may be analyzed additionally only when authentic counts and taxonomic grain are independently verified.
- **Outcome choices after source gate (distinct models)**: areal AFDM, areal periphyton-Chla, original published benthic MMI/quality index, EPT individual fraction only if original exact count field and aquatic orders are verified, and/or richness from source counts. Use common covariates/state/eco/latitude/longitude/river order/season/conductivity, HUC8-group heldout Ridge + HUC8 cluster intervals and source-count-only null controls. No unobserved outcomes are fabricated.
- **Out-of-sample**: same HUC8 5fold assignments deterministic and stable for all outcomes; compare baseline→measured ambient TN+TP per outcome. Detect label leakage, aggregate by unique physical SITE_ID for pseudoreplication sensitivity, and do not treat 2023/24 as repeated 2018/19 identical sites without native site ID proof.
- **Ecosystem-services boundary**: MMI or EPT fraction ≠ emergent aquatic insect biomass, AFDM ≠ primary productivity rate, Chla ≠ pollination; independent consumer data are required for any chain.
- **Status**: preserve positive EPT/TN and negative TP/Chla original spatial results from #73, plus WWTP reporting-support failure from #74. The next benthic file QA is independent new *source access* work, not evidence those findings were false or confirmed.

Official public files: [EPA NARS original download table](https://www.epa.gov/national-aquatic-resource-surveys/data-national-aquatic-resource-surveys), NRSA2023–24 section.
