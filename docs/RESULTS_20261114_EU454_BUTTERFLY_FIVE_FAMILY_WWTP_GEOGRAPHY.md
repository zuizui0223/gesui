# Actual original-source five butterfly family comparison — 14 November 2026

Status: **new exploratory local original-artifact reanalysis; unreviewed Draft PR #84**. This is an original-data cross-family investigation, not a literature meta-analysis, source-specific sewage chemical attribution or butterfly feeding observation. All original source artifacts were downloaded through connected GitHub and extracted locally. The source-frozen data/denominator were validated before running the new family contrasts. This local calculation is **not** claimed to be a completed GitHub Actions run; the branch contains an Actions recipe for further reproducibility.

## Original-source audited study population

- Same original 454 European H7 basin×2017–19 or 2020–22 rows, 232 unique H7, 90 actual level-6 hydrological parent groups; entire same original source sample, NO sign-based species/family filtering.
- Original source H7-454 EQR/plant/butterfly CSV from Actions #43215754868, unmodified original CSV SHA256 \`c73bedbc4742c182565673867dd7dbd6b6c2842b5813b6906e7c5e986d929eca\`.
- Five butterfly family count source (Lycaenidae/Nymphalidae/Pieridae/Papilionidae/Hesperiidae) from four authentic original query shards Actions #43204891739, not newly queried data.
- Original level6 geography source SHA256 \`372b800b85abe2f395241987cea2000ff901cf5aff31bf1d1fb825304f2495fc\`.
- HydroATLAS source SHA256 \`ce6780e01d2b98ed8e852574adb827cb3125be792a164e8be9e2ebb917b3a686\`.
- Every joined H7×period has all five source-family GBIF query status success, original summed total matches \`round(expm1(log1p_GBIF_5family_butterfly_source_records))\`, and original Lycaenidae share/log odds are reproduced exactly from original numeric counts. **No missing source count is assigned zero.**
- Total original observed family-index records **4,400,070** across 454 H7-periods (not independent butterflies): Lycaenidae 761,024, Nymphalidae 2,353,609, Pieridae 921,651, Papilionidae 133,412, Hesperiidae 230,374.
- Source-country groups by stable per-H7 country: FR 299, IT 48, LT 30, NO 22, HR 16, NL 12, SE 8, IE 6, BE 6, AT 6, HU 1. Only 11 stable country groups (not 15 independent national replications).

## Same original 454 response rows; new actual observed comparisons

Baseline fixed original model: original catchment area, site count, original indexed Animalia/Plantae/butterfly recording opportunity, real EEA river flow, original EQR biological method, period, country, authentic upstream HydroATLAS crop/urban/flow. Added predictor is original modelled HydroWASTE flow routed to *basin outlet*, NOT measured source-delivered wastewater dose. H6-parent-disjoint fivefold Ridge(alpha10), training-only standardization and H6-equal training weights; primary signed MSE improvement is H6-equal heldout error. Fixed H6 parent 2000-resample error descriptive intervals. Country-left-out excludes country FE. No sign cherry-picking, no WWTP N/P/Na lab results, no species abundance or puddling behavior.

| Original GBIF five-family composition endpoint | TRUE H6 equal MSE gain from modelled source % [H6 bootstrap descriptive95] | Country-left-out 6-supported-country equal gain % | Adjusted source slope per +1SD [H6 descriptive95] |
|---|---:|---:|---:|
| **PRIMARY** Lycaenidae / Pieridae original count log-ratio | **+0.891** [-5.189,+5.170] | **-5.679** | -0.057 [-0.141,+0.027] |
| **PRIMARY** Lycaenidae / Nymphalidae original count log-ratio | **-1.610** [-5.333,+1.700] | **-1.523** | +0.006 [-0.072,+0.084] |
| Lycaenidae / other four log-odds | **-2.454** [-4.870,-0.291] | -3.439 | -0.012 [-0.084,+0.060] |
| Nymphalidae / other four log-odds | +3.867 [-3.408,+10.876] | +3.407 | -0.032 [-0.084,+0.020] |
| Pieridae / other four log-odds | **+6.361** [+1.361,+11.208] | **-3.210** | +0.060 [+0.007,+0.112] |
| Papilionidae / other four log-odds | -3.271 [-9.405,+3.149] | +5.826 | +0.013 [-0.186,+0.212] |
| Hesperiidae / other four log-odds | +0.151 [-7.253,+5.461] | -1.322 | +0.078 [-0.048,+0.204] |

### Scientific decisions
1. **No robust evidence of Lycaenidae-specific enhancement at the existing basin-period spatial grain**. Both predeclared-to-this-new-analysis primary Lycaenidae pairwise comparisons have H6 prediction bootstrap intervals crossing zero and fail country transfer.
2. Pieridae family composition has a **local-geography predictive positive** after upstream crop/urban adjustment in 90-H6 analysis but **sign reverses** when leaving countries out. This is NOT a robust ecological or causal result and is a secondary multiple comparison.
3. Declines in original Lycaenidae-versus-other-four prediction when adding WWTP are retained (do not describe a universal Lycaenidae decline; it is *prediction degradation*). The coefficient still spans zero.
4. The observations do **not** test mud puddling, N-rich carrion use, Na preference, insect reproductive fitness, ant–caterpillar interactions, plant–pollinator interactions, aquatic emergence or predator consumption. Natural wet sites and wastewater influence cannot be distinguished from recorded source geography. A share change may come from numerator or denominator.
5. The stronger research target becomes *feeding- and host-plant-specific wastewater resource use*, not global abundance or family-level occurrence shares. Direct paired field data across **independently verified actual outfalls**, natural groundwater seepage / riverbank controls and stormwater-only drains must record species, sex, behaviour, observed substrate, freshwater moisture, dissolved Na/N/DON, floral/host availability, ant activity, plus hydrology, and ultimately aquatic insect and riparian predator interactions.

**No positive effect was protected by redefining endpoints after inspecting the signs.** Prior gesui source-linked EQR positive prediction and all null butterfly/plant outcomes are preserved.

## Numerical audit status / reproducibility
A separate local Python/NumPy original-artifact reimplementation generated 7 outcomes under the unchanged H6 fold assignment and full country exclusions and archived its local CSV/JSON output. This is **a single ecological dataset reanalyzed with source checks**, not independent ecological replication. Source extraction of known full original five-family records and original 454 cohort is completely separate from new ecological outcome tests. The GitHub branch contains the audited input contract and an independent implementation/workflow, but no new workflow-success claim should be made before checking Actions run status.

No main merge or overwriting older PR evidence; no eBird/TinyFish.
