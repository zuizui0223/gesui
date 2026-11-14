# Original-source species-level butterfly–plant macro test: does a family signal reach species diversity?

**2026-11-14 · posthoc exploratory source-grounded macro test, NOT a causal sewage result or a meta-analysis.** Stacked Draft PR #85. Uses actual zipped original GH source data, not synthetic ecological rows.

## Primary source population and fidelity
- Original 234 EU H7 × two GBIF periods (2017–19, 2020–22) × six groups (Plantae and Lycaenidae, Nymphalidae, Pieridae, Papilionidae, Hesperiidae): all 2,808 unique source requests found across four [original full acquisition run #43592211212](https://github.com/zuizui0223/gesui/actions/runs/43592211212) artifacts. 234 unique H7, 468 H7×period, **2,800 successful first-pass facets plus 8/8 HTTP429 failures recovered by selective genuine-source retry: 2,808/2,808 successful final source facets**; original successful facets untouched.
- Original artifact ZIP ids and digest SHA256: shard0 ID 14113028121 b550d2195755dbf169e81232db1917c48dcce6969828d2ece651352dee661220; shard1 ID 14113381217 a40feaff73ec4263b11eaed9f61caefab1d80ee9dd1e57ed3fa0187b45859d14; shard2 ID 14112564855 71c9aed9afad851409bf7bed62ddb4a84663df32fab37ba4f686c6079e468775; shard3 ID 14112454382 04670247202167aabb41e24228208d80f9bfe32269d96cec848a85b3921da902.
- Exactly **397/468** H7×period rows have actual successful source identity facets and **≥80% of occurrence records assigned speciesKey across ALL SIX groups**; at ≥90% **300/468**. These are provisional post-source QA thresholds, not an experimentally fixed outcome-free exclusion.
- The historical authentic original 454-row EEA/BQE × plant × family × HydroWASTE panel ([run #43215754868](https://github.com/zuizui0223/gesui/actions/runs/43215754868)) aligns with **454/454** six-group source-complete rows after the eight source-only 429 repairs. Exact original family occurrence total reproduced **403/454**, and maximum count drift among the other 51 source requery snapshots was six records. Reconstructed historical speciesKey plant rarefaction S50 has Pearson r **0.9999995932** with original S50, median difference ~1.4×10^-12, so species-count computation can be independently reconstructed. **Not independent ecological sampling**.
- Strict six-group source-record coverage ≥80% within these 454 yields **n390 actual H7×period records, 217 distinct H7 catchments, 83 distinct authentic H6 hydrological parents, 14 countries**. Unrestricted source-success sensitivity n454, 232 H7, 90 H6, 15 countries. Do not call these 383 independent catchments.

## New response, common space-time resolution
- For five butterfly families combined, merge real speciesKey integer histograms in each H7×period (do not infer occurrence absence) and compute exact hypergeometric expected butterfly **speciesKey taxon-ID richness** at **20 sampled SOURCE records**: sum over detected speciesKey j of 1 − choose(N−n_j,20)/choose(N,20). This is not individual butterfly abundance, breeding, realized local diversity or functional trait diversity.
- Reference plant response: analogously reconstructed **50-record rarefied plant speciesKey taxon-ID richness** from actual Plantae source facets. Avoid treating the real 20 and 50 record depths as 20 or 50 field-collected individuals; these are GBIF indexed observation records.
- The predictor W is old signed HydroWASTE **static basin-outlet modelled flow per upstream area**, NOT field effluent dose at the butterfly location. Same source country, monitoring method, period, area, EEA sampling effort, GBIF source effort for Animalia/Plantae/five butterfly families, stream flow, upstream crop/urban and HydroATLAS natural flow in baseline.
- Posthoc method: Ridge alpha 10 with training-only standardization and categorical encoding; fixed H6-hash 5fold heldout, all H7 periods with same H6 kept together, inverse H6 visit frequency weights in training, H6-equal heldout MSE, 2,000 fixed heldout H6 error resampling descriptive 95 intervals. Leave-country-out transfer models exclude country label as feature; all source zeros/missingness and row reuse preserved.

## H6-parent disjoint predictions (strict source coverage ≥80%, n390/83 H6)
| Outcome | New predictor added to baseline | H6-equal heldout MSE gain % | Fixed holdout H6 descriptive95 interval % | Leave-country-out equal-country gain % |
|---|---|---:|---|---:|
| Butterfly S20 | Static modelled wastewater W | **−2.099** | **[−4.279, +0.159]** | −1.381 |
| Butterfly S20 | Plant source S50 | +7.137 | [−2.807, +18.244] | +10.529 |
| Butterfly S20 | W after adding plant S50 | −1.181 | [−3.099, +0.796] | −0.402 |
| Butterfly S20 | Plant S50 after adding W | +7.972 | [−2.476, +19.301] | +11.394 |
| Plant source S50 | Static W | −1.423 | [−3.485, +1.977] | −4.394 |
| Plant source S50 | Butterfly S20 | +6.669 | [−4.427, +13.508] | +8.662 |

**Sensitivity to source identity coverage selection:** all 454 source-complete panel rows without the ≥80% threshold yield W→butterfly S20 **−0.242% [−2.600,+2.323]**, plant→butterfly S20 **+2.227% [−4.294,+11.228]**, W→plant S50 **−1.537% [−3.977,+2.176]**. Do NOT claim that marginally significant negative C80 effect confirms lower butterfly diversity: this is posthoc, model-selection/multiple-test sensitive, with country transfer failures and taxa/source sampling effort.

## Independent numerical-code check
- Offline original-source extraction and source-grain QA using raw zipped artifact JSON and official original 454/H6/Atlas. All 12 test combinations for 2 source-eligibility cohorts × 6 models were independently refitted with **manual NumPy weighted Ridge normal equations**, train-only numeric scale, train-only categorical vocabulary. Maximum absolute pointwise squared-error discrepancy from sklearn **<1.7×10^-11**, consistent within numerical precision. This is independent numerical implementation only, not a second set of ecological field observations.

## Biological reading and explicit falsification
- Earlier source-reassembled **five-family Hill-Shannon effective FAMILY number** (PR #84) had H6-local modelled-source geography gain +10.486% but failed whole-country transfer; **new actual source SPECIES richness S20 does not show the same general source added-value. Family-level relative record evenness cannot be relabelled as species or interaction diversity.**
- There is weak context-dependent signal that plant source species richness may aid predicted butterfly species richness at same H7-period; descriptive H6 bootstrap still crosses zero, source effort and geography may cause common trends. There is no demonstrated local host feeding interaction or ecological mediation.
- The outstanding direct mechanism test still requires actual observed adult puddling by species/sex/substrate, known source chemistry (Na, NH4, DON, conductivity, DO and flow), larval host/ant evidence, emergent insect biomass and riparian predator diets, plus non-sewage naturally wet and carrion comparators.

**Source identity 429 repair:** selective requery and original-row-preserving merge launched by [workflow](https://github.com/zuizui0223/gesui/blob/analysis/gesui-original-h7-specieskey-pilot-trigger-20261114/.github/workflows/gesui-H7-specieskey-429-selective-recovery.yml). Completed successfully: [run 43593227162](https://github.com/zuizui0223/gesui/actions/runs/43593227162), all **8/8 recovered**, artifact 14113695873, ZIP SHA256 **4b37c741d1bd28c178b0cef61c9457af0dd30797c173f9384ebe1d1d0036d07c**; initial 2,800 valid row sources preserved unchanged. No eBird, no TinyFish, all PRs Draft/unmerged.
