# Sewage-associated assemblage re-analysis: five butterfly families, temporal turnover, aquatic–plant coupling
**14 November 2026 · new source-grounded exploratory POSTHOC results · Draft PR #84, not independently replicated or peer-reviewed**

## Why the estimand changed
Earlier five-family macro comparison found no reliable Lycaenidae-specific wastewater association. A more suitable macro proxy for *distinctive assemblage structure* is **within-five-family observed compositional diversity**. Separate analyses must distinguish geographic assemblage structure, temporal change, and measured ecological interactions. HydroWASTE source model at the basin outlet is a static geographic index, not field-delivered effluent, nitrogen/protein/sodium or discharge at observed butterflies.

## Original data used (all authentic same source)
- Original 454 European HydroBASINS L7 basin×period records, 232 H7, 90 genuine H6 parents: exact previous PR #84 source reassembly of the four independent five-family GBIF source shards + EEA/plant/source geographic panel.
- Original full 5-family GBIF recording counts: 4,400,070 *source index records*, not butterfly census individuals.
- Added **source-locked original 169 exact same physical EEA field station/method basin pairs**, from original Actions run #43214039321. All 169 original historical Lycaenidae Δlog-odds are independently reconstituted and matched to current original family source counts without row reassignment.
- Original 222 H7 basins have both 2017–19 and 2020–22 family panels.
- Computation input SHA256 for full original-source reassembled CSV: `45eb30d9c04122e5a1cefdfd5cd1566754625e2fd8f92df7704f815772fde9af`. 169 basin-pair original Actions ZIP: `804476aca8ec811b79a5dd649952233d8c66459e486fd8f0eac859808a67b42a`.

## Prespecified-to-this-new-test outcomes (before numerical fit)
T1: observed original 5-family Shannon effective number `exp(-sum(p_f*log p_f))`, same 454 source rows. Add static modelled W after original country/EEA method/period/recording efforts/HydroATLAS crop/urban/flow. 90 authentic H6 group-disjoint 5fold, group-equal training weights; H6 group-equally averaged MSE, 3000 H6 bootstrap error resampling. Repeat selection thresholds at original all n≥20 butterfly records, n≥100, n≥500.

T2: source count half-increment Jensen–Shannon distance between paired family compositions in **222 source-repeat H7**; add W as *context* after observation effort changes, source country and HydroATLAS controls, 89 H6 parents. This is NOT temporal change in discharge.

T3: In exactly same **169 physically matched river-station/method H7 pairs** and 52 H6 parents, add five-family Jensen–Shannon turnover to source-original background, GBIF effort deltas, plant S50 change and HydroATLAS environment to predict actual EEA **ΔEQR** (2018 to 2022).

T4: In the same 169, add five-family turnover to predict **Δplant indexed S50**, adjusting same-source environment, GBIF effort and measured ΔEQR; observational only.

## Signed numerical results
| Endpoint; frozen source sample | Added source/turnover true H6-group equal heldout MSE improvement | H6 group bootstrap descriptive 95 |
|---|---:|---|
| T1 observed five-family Shannon Hill¹ effective number, **454 rows/90 H6** | **+10.486%** | **[+1.559,+18.516]%** |
| T1 source >=100 five-family indexed records, **439/89** | +10.945% | [+1.667,+19.265]% |
| T1 source >=500 five-family indexed records, **384/80** | +10.443% | [+0.934,+17.662]% |
| T2 5-family Jensen–Shannon historical turnover, **222/89** | **−2.137%** | [−7.772,+3.671]% |
| T3 same river monitoring site/method ΔEQR, **169/52**; add butterfly turnover | **−1.366%** | [−16.571,+10.265]% |
| T4 same H7×period Δplant record-index S50, **169/52**; add butterfly turnover | **−11.725%** | [−27.841,+3.132]% |

All numeric bootstrap values are uncertainty diagnostics for fixed predictive errors, **not formal causal effects**, not independent external study replications, and these posthoc exploratory model choices were made after earlier WWTP/EQR, S50, Lycaenidae models were known.

### Post-result sensitivity checks (never label as confirmation)
- Original T1 only 2017–19: +16.117% H6 MSE gain, descriptive95 [+6.174,+24.695]. 2020–22: +1.883%, [−7.289,+9.584].
- France-only (299 original rows): **−0.229%**, [−3.190,+2.730]; non-France (155 rows) +17.818%, [+5.249,+31.501] but mixed country composition. **Whole country-out row-weighted MSE WORSE by −7.374%**; six well-supported countries equal-weighted gain +8.703%, hence weight-dependent cross-country transfer FAIL.
- Original source W vs upstream `urban_up`: Spearman **ρ=+0.781**; W vs Hill¹ unadjusted ρ=+0.120, per +1SD W fitted conditional Hill¹ beta **+0.0473** (H6 descriptive95 **[−0.0237,+0.1184]**, spans zero).
- Alternative source family index, *post-result* sensitivity: Simpson effective family number H6 +8.524%, descriptive95 [−0.066,+16.243]; 20-record exact hypergeometric expected 5-family richness H6 +8.047%, [−1.152,+16.863]. Thus H6 signed positive geographical *predictive signal* depends on the definition of recorded family diversity; not uniform robustness.
- 222 original paired-H7 median count-based five-family JS distance **0.04664** (not actual specimen species turnover).

## Scientific interpretation
1. **New useful result:** Modelled source geography adds H6-region-local predictive information for a *particular* five-family GBIF recording-diversity index at a fixed space-time grain. It does NOT prove sewage increases butterfly family diversity. The adjusted sign is uncertain and country transfer fails.
2. **Directly unsupported:** increased Lycaenidae family share, universal five-family source effect, five-family temporal turnover tracking static source, or shared aquatic→plant→butterfly temporal co-change. Negative predictive results are preserved.
3. **Not measured:** chemistry at the butterflies, active puddling (vs rest/nectar), sex or species, confirmed larval host, ant interactions, emergence mass, predators/diet, fish-mortality mechanisms and ecosystem service. Evenness of five FAMILY index counts is **not species-, trait- or interaction-network functional diversity**. No AFDM was measured at the same 454 ecology source rows.
4. Next source priority: verified discharge-outlet/receptor position and independent matched local natural humid, stormwater and mineral environments; count **behavior-specific butterflies per encounter opportunity**, then separately measure flower/host/larva/ant and aquatic insect emerging biomass and consumer diet. Separate local outfall ecological comparisons from independent-study hierarchical meta-analysis.

## Relevant independent mechanism literature (not empirical tests of five-family WWTP ratio)
- Beck et al. (1999), comparative puddling resources, [DOI 10.1007/s004420050770](https://doi.org/10.1007/s004420050770); **not a sewage experiment**.
- Boggs & Dau (2004), source substrate and sodium specialization, [DOI 10.1603/0046-225X-33.4.1020](https://doi.org/10.1603/0046-225X-33.4.1020).
- A 2026 treated-effluent whole-stream BACI/paired experiment reports biofilm standing biomass effect size about 2.93 and a diet shift towards algal food, with no food-chain-length increase: [DOI 10.1007/s10750-025-06103-w](https://doi.org/10.1007/s10750-025-06103-w). **Different biological layer/site source**, not additive support for the 454 butterfly count profile.
- Three wastewater-treatment outfalls with aquatic insect/riverbank spider isotopic and microbiome evidence: [DOI 10.1016/j.envpol.2025.125948](https://doi.org/10.1016/j.envpol.2025.125948).
- Effluent-associated neonicotinoid transfer into riparian spiders: [DOI 10.1021/acsenvironau.5c00021](https://doi.org/10.1021/acsenvironau.5c00021).

**Integrity:** All local output files and executable scripts retained as conversation source files. New official GH workflow successful source run for *this* augmentation is **not yet asserted**; remain Draft and no main merge.