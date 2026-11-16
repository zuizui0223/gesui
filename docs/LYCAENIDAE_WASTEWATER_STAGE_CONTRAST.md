# Wastewater × Lycaenidae: competing adult, larval, and ant-mediated pathways (2026-11-08)

**Status: mechanistic research plan, NOT new measured sewage–butterfly association.** This note is independent of the main aquatic BQE study and does not change the original `chocho` GEB manuscript. The original `gesui` all-Plantae response remains the primary flora estimand.

## Motivating tension

In existing **chocho** frozen potential-resource geography, **22/25 Lycaenidae (88%)** show added regions when anthropogenically introduced larval host distributions are included; the family median is **+73 geographic resource units**. This is the *geography of potential larval host availability*, not wastewater exposure, adult habitat use, oviposition, reproduction, biomass, or abundance. Major butterfly families are similarly positive; Lycaenidae is not statistically exceptional on this metric. The frozen source is retained at [original reference commit](https://github.com/zuizui0223/chocho/tree/e8b8f86c7a889b6a3224a4ef69ffe5e0c7738660).

Specific examples from the original records: *Lycaena phlaeas* 140→224 (+60% geographic opportunity), *Polyommatus icarus* 128→207 (+61.7%), *Brephidium exilis* 242→267 (+10.3%). *Pseudozizeeria maha* 16→275 (+1618.75%) is **not acceptable as a primary showcase** because the original `host_taxonomy_lower_bound_adequate` flag is **False** and its native-only reference is only 16.

## Three separable mechanisms (do not collapse into sewage nitrogen)

**A — Adult nutrient acquisition: dissolved organic nitrogen/protein-rich moisture versus Na+ salts.** Tropical field bait experiments (Beck, Mühlenberg & Fiedler 1999, *Oecologia* 119:140–148, DOI 10.1007/s004420050770, PMID 28308154) found many male puddling Lycaenidae disproportionately attracted to protein-rich albumin, while Papilionidae/Pieridae were more often associated with NaCl. This contrast is taxon-, setting-, and bait-specific; **neither an increase in bulk treated effluent N load nor an N:P ratio shows that bioavailable albumin/protein or odor cues are present**. The 2025 puddling review (DOI 10.1093/aesa/saaf007) emphasizes micronutrients and unresolved sensory/fitness mechanisms. Wastewater may supply perennial moisture even without protein; sewage treatment can remove organic substrates. Direct adult puddling is a behavior and cannot be inferred from GBIF occurrence alone.

**B — Larval resource geography and nutritional quality can decouple.** Butterflies may acquire a larger geographical portfolio of introduced larval host plants while their local realized fitness declines under intense enrichment or contaminants. A direct fertilization experiment including *Lycaena phlaeas* and *L. tityrus* found increased nitrogen in host vegetation and reduced survival across six Lepidoptera taxa (PMID 30288608); those agricultural fertilization doses are **not equivalent to river water, soil N, or WWTP reports**. Host presence and plant Ellenberg N affinity do not measure larval survival, oviposition or local nitrogen dose. A hump-shaped or sign-reversing local response is a testable **hypothesis**, not an observed wastewater effect.

**C — Conditional ant mediation.** Some Lycaenidae larvae have ant-associated mutualisms, whereas others are ant parasites or ant-independent. The biological N economy of ant-tended caterpillars and host-plant N affinity have published comparative support (Pellissier et al. 2012, DOI 10.1111/j.1420-9101.2012.02555.x), but not a monotonic within-clade nitrogen response; some socially parasitic species rely on ant brood. Ant host identity and ant site occupancy must be independently observed; do not assume all Lycaenidae are obligatorily myrmecophilous.

## Core macro-ecological question

**Does wastewater-associated water and nutrient redistribution alter an adult butterfly's resource geography differently from its larval host-resource geography, and does this discordance depend on aridity, nutrient form, and the requirement for ants?**

Distinct observation units:
- **Source/exposure**: verified original EEA WWTP discharge point × year × downstream river reach. `uwwNDischargeMeasured` and `uwwPDischargeMeasured` in t/yr are WWTP *reported effluent masses*, not dissolved organic N, sodium or concentrations actually delivered to mud-puddling substrate.
- **Adult**: verified photo/observation annotation of puddling or liquid feeding by butterfly species × site × period; independent observer-time controls. An occurrence near sewage is merely *occurrence*, not puddling. iNaturalist user-created observation field 17428 mixes **aggregation, roosting and mud-puddling**; `Yes` is an *imperfect candidate pool*, never a directly validated puddling label (https://www.inaturalist.org/observation_fields/17428).
- **Larval plant resource**: current verified WCVP/GBIF known host identity and present occurrence × same region/season; no inferred feeding from co-occurrence. In matched regions, estimate host resource opportunity, not actual biological production.
- **Ant channel**: original species-level ant ecology (independent review/species natural-history register) and ant presence if source-validated; do not substitute generic ant species count for ant partner.

## Prespecified sequence before looking for positive wastewater associations

1. **Phase 0 — independent source and biology data QC:** pull frozen `chocho` Lycaenidae host-resource summary with its Git blob SHA; preserve 25/22 legacy positivity. Separately source 2022 EEA measured WWTP N/P mass × *actual* freshwater discharge positions and HydroRIVERS directed reach topology from `gesui`; 2020 year-compatibility pilot had only 17 BQE source-exposed stations, but this does *not* preclude a distinct butterfly photo-based cohort if it passes a new independent biological support gate. Keep 2022 source year and 2018/2020 exposure separate.
2. **Phase 1 — detection/behavior feasibility, BLIND to effluent gradients:** predeclare 1–5 km buffers on true water-connected downstream river reaches, physically upstream reaches and no-registered-source/control river corridors within comparable country/urban/road/land-cover/cropland and observation coverage strata. Freeze selected coordinates and inclusion before querying Lycaenidae outcomes. Exclude river/road crossings that break upstream/distance definition. Record no-observation as NA and never as zero abundance. Observer IDs, coordinates/uncertainty, year and target-group butterfly effort are required; screen geolocation jitter/obscured observations.
3. **Behavior label audit:** inspect photo and text for a visible contact with damp substrate and/or extended proboscis. Distinguish **confirmed puddling / possibly drinking / simply resting / grouping**, with blinded duplicate labeling for a random subsample. Use iNat field 17428 only to find *candidates*, never as a gold standard. Compare Lycaenidae and Papilionidae/Pieridae (taxon/area/season effects); separate sex when actually annotated, never infer male from family-level habits.
4. **Exposures with chemical specificity:** W=`dry-season wastewater-supported moisture` (requires independent surface flow evidence); N=`measured organic-N/amino acid/protein fraction in water or damp substrate` (requires independent monitoring); S=`Na+/electrical conductivity`; T=`toxicity/DO` (requires its own monitoring). WWTP mass + river reach may only test *source association* and *direction*, NOT direct A/N/S/T mechanisms.
5. **Independent response tests:** (i) adult confirmed-puddling odds relative to butterfly photo opportunity, (ii) larval host **resource** opportunity after equal-effort species inventory, (iii) actual ant-tending/larval survival **only if observed**. Analyze all plants as background, with Lycaenidae known hosts a targeted secondary analysis. Keep Lycaenidae vs Pieridae/Papilionidae comparisons as effect *interactions*, not subjective significance rankings.
6. **Falsifiers/negative controls:** same-distance upstream/outside hydrological graph; unrelated river source in same municipality; urbanization/cropland covariates; non-puddling photo negative control; equal target-group observation records; source-without-flow control; leakage across spatial blocked CV and observer holdouts.
7. **Stop gate:** if too few verified puddling photos, compatible exposure years, independent water chemistry, or ant-site data, report *which layer fails*; allow descriptive GBIF occurrence/resource geography but no inferred sewage puddling, fertilizer-type effect or ecological trap.

## Hypothesis matrix (all untested for WWTP)

| Hypothesis | Predicted interaction under measured source evidence | Critical alternative |
|---|---|---|
| **H-A adult nutrition** | Organic-N source availability × Lycaenidae yields higher confirmed puddling opportunity than sodium-only conditions, relative to Papilionidae/Pieridae | Generic wet-site/flower/observer concentration |
| **H-B seasonal water subsidy** | Dry-season connected effluent is associated with puddling records especially in intermittently dry rivers | Other permanent water/irrigation and season effects |
| **H-C larval filtering** | After controlling host availability, *measured* high local N shifts larval performance nonlinearly rather than uniformly improving it | Habitat management, pesticides and host species turnover |
| **H-D ant bottleneck** | Conditional response of ant-associated lycaenids depends on confirmed suitable host ant occupancy | Global ant richness is not specific partner occupancy |

**Evidence hierarchy:** chocho source-confirmed 22/25 potential resource expansion = supported but **not sewage**; adult/protein preference and excessive host N effects = published experiments, not direct wastewater; EEA measured discharge = physical source data; **actual wastewater × Lycaenidae effect remains untested**.

## Sources

- [chocho Lycaenidae 25-species summary](https://github.com/zuizui0223/chocho/blob/e8b8f86c7a889b6a3224a4ef69ffe5e0c7738660/provenance/reviewer_defenses/results/butterfly_taxonomic_family_expansion_v0.1.json)
- [chocho frozen 239-species resource opportunities](https://github.com/zuizui0223/chocho/blob/e8b8f86c7a889b6a3224a4ef69ffe5e0c7738660/data/frozen/figure_sources/anthropogenic_species_metrics.csv)
- [gesui source-reported WWTP N/P and aquatic biology (2022)](https://github.com/zuizui0223/gesui/pull/14)
- [gesui 2020 independent source-year feasibility stop](https://github.com/zuizui0223/gesui/pull/19)
- [Beck et al. 1999](https://pubmed.ncbi.nlm.nih.gov/28308154/) • [Pellissier et al. 2012](https://pubmed.ncbi.nlm.nih.gov/22725645/) • [Fertilization/larval survival experiment](https://pubmed.ncbi.nlm.nih.gov/30288608/) • [2025 puddling review](https://academic.oup.com/aesa/article/118/2/110/8010399)
