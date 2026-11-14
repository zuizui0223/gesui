# Wastewater macro/meta integration: source quality and independent-study gate (14 Nov 2026)

**Scope:** source-verified and literature-screened methodological ledger, not a completed meta-analysis. Stacked Draft PR #85 from Draft #84. No merger to main or chocho.

## Actual original-species data
- Fixed upstream geometry: 234 European H7 catchments × 2 historical windows (2017–2019, 2020–2022) × six biological source groups (Plantae and five butterfly families) = 2,808 full source-request opportunities, NOT independent ecological samples.
- Original pilot [GitHub Actions run 43591984062](https://github.com/zuizui0223/gesui/actions/runs/43591984062): 8 H7 × 2 periods × 6 groups = **96/96 source identity facets fetched successfully; zero API failures**. Artifact 14112154992 ZIP SHA256 d8e5edda581f109e4013bfe79d278d831bdd57c24896a060f38b89da117cd2f7.
- Inspecting actual 96-row ZIP speciesKey record coverage: **14/16 H7×periods** have at least 80% identified species-key record coverage in ALL six groups. Remaining two periods are one same H7: Lycaenidae 0.640/0.543 and Pieridae 0.713/0.637. This ≥0.80 threshold is a **post-source QA diagnostic**, not a prespecified outcome exclusion.
- Four-shard whole-population original GBIF identity retrieval launched as [run 43592211212](https://github.com/zuizui0223/gesui/actions/runs/43592211212). Do not assert completion until all four job source audit JSON files and SHA-checked artifact ZIP files are inspected.
- Real speciesKey histograms allow plant/butterfly potential co-occurrence and taxonomic crosswalk at the SAME H7×period scale. They do NOT show true site presence or actual plant–larva, pollination, puddling or ant mutualist interactions. Missing facet identities are not biological absences.

## Research intervention registry (published empirical studies, unextracted effect-size status)
| Publication / identifier | Intervention and sampling design | Nonindependence unit | Eligible response and caution |
|---|---|---|---|
| [de Guzman et al. 2023](https://doi.org/10.1016/j.jenvman.2023.118510) | BACI adding highly diluted treated WWTP effluent to an unpolluted stream | EHU experimental stream | Reported diversity reductions, higher production/herbivory, and lower terrestrial food-web energy flux |
| [de Guzman et al. 2026](https://doi.org/10.1007/s10750-025-06103-w) | Related analysis of same effluent diversion experimental setup, paired before/after reaches | **Same EHU stream as 2023; not an independent source** | Reported biofilm standing stock ratio ~2.93, change to algal pathway, unchanged food-chain length; extract covariance before fitting |
| [Di Prinzio et al. 2024](https://doi.org/10.1016/j.scitotenv.2024.170237) | Up/downstream WWTP contrasts in six Patagonia urban streams | Stream/source ID, full text crosscheck | Fish biomass, density, abundance higher in controls; observational, not a manipulated effluent addition |
| [Bucol et al. 2025](https://doi.org/10.1016/j.aspen.2025.102432) | Before/after WWTP construction that REMOVED much raw wastewater, Dongda Creek (2019–20 vs 2021–23) | Same Dongda stream | Sensitive insect community/richness increase; **opposite intervention direction**, keep separate from added treated effluent |
| [Environmental Pollution 2025](https://doi.org/10.1016/j.envpol.2025.125948) | Upstream/downstream three WWTP systems, insect + riparian spider isotopes and 16S | Individual WWTP/system, within-study dependence | Direct cross-boundary isotope/microbiome source evidence, not food subsidy biomass |
| [USGS 2025](https://doi.org/10.1021/acsenvironau.5c00021) | Effluent-dominated Iowa stream, insects, fish, spiders tissue chemicals | Same stream | Neonicotinoid chemical trophic transfer, NOT emergent insect biomass or predator abundance |
| [Beck et al. 1999](https://doi.org/10.1007/s004420050770) | Butterfly bait comparisons for protein and sodium | Butterfly bait experiment | Mechanism only: EXCLUDE from wastewater effect-size pool |
| [Lamie et al. 2025 review](https://doi.org/10.1093/aesa/saaf007) | Butterfly puddling literature synthesis | Secondary source | Mechanism/background, deduplicate cited original experiments |

## Same-scale inference: source action and response must be aligned
1. Define study × independent discharge system × within-system contrast × period × response metric, retain shared controls and study/system IDs. Do not treat upstream/downstream or repeated papers as independent experiments.
2. Separate exposure regimes: (A) treated WWTP effluent ADDED to natural river, (B) plant installation/remediation REDUCING raw sewage, (C) unmanipulated point source upstream–downstream contrasts, (D) generic nutrient/mineral bait only. They answer different potential-outcome questions.
3. Meta pool within actual compatible biological outcomes only: biomass/density (log response ratio), measured BACI stock contrast (log ratio-of-ratios when observed), diversity (comparable richness/rarefaction or EQR instrument strata), butterfly behavior-specific feeding rate (observed count offset), isotope/contaminant transfer (its own evidence channel). **Predictive MSE improvement is NOT a biological effect size.**
4. Nested macro grain: local feeding patch → measured reach/outfall → H7 catchment×period → independent H6 geographical parent. A basin outlet MODELLED HydroWASTE source index must not be assigned as a chemical measurement at butterflies. Record flow fraction, dilution, NH4, Na, DOC/DON, DO and actual treatment source, plus wet non-sewage controls.
5. Exact plant and butterfly speciesKey matching can identify **candidate host links only** after accepted-name and source-vetted host crosswalk. Grounded larval observation, ant attendance, flower visits, aquatic emergence biomass and predator diet are needed to establish observed interaction edges.

## Status and falsifiers
- Current EPA three-period PHab→EPT predictive gain is authentic but does not identify WWTP chemistry; EPT taxon richness fails one era, clinger composition/tolerance correlate with actual riverbed sediment.
- Current five-family GBIF geography result is variable by country/period, does NOT establish Lycaenidae wastewater increase, or that five-family assemblage Hill-Shannon diversity is functional diversity.
- The central ecosystem hypothesis is an interaction mosaic: a source may increase adult puddling or biofilm biomass while sensitive benthos and larval recruitment decline. All positive/negative/null outcomes survive preregistered source and independent-region tests.

Literature source abstracts, publisher records and study methods were individually checked; numeric meta-effects and sampling covariances have NOT been extracted.
