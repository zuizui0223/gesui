# 2025 Denmark EBMS NEW399: positive Lycaenidae three-family share reverses with five-family denominator and geography

**2026-11-10, exploratory POSTHOC. Not a sewage cause or wastewater nutrient dose.** This is a newly source-locked 2025 eBMS survey-site cohort **disjoint from the already analyzed 107 Denmark sites**, with original accepted individualCount data read only after a no-butterfly source/observer event-opportunity gate passed. We already knew the 2022 DE/AT iNaturalist photo-positive and the older Denmark 2023–2025 results when this new sample was defined.

## Data provenance and pre-2025-family-count design

- Original licensed Denmark eBMS v1.4 (GBIF dataset caa38f94-ddeb-43e7-9a27-1ac46d7ddc5c; CC-BY-NC 4.0), ZIP SHA256 `2cc64ec352754223db06a92eef5b01803968eac94bae34057c17f1abdb6da79b`; official complete source **15,606 Event** and **21,617 Occurrence** rows, accepted GBIF taxonomic species-to-family mapping frozen before original effect analyses. Personal observer data, photo pixels, licensed original occurrence rows and raw archives NOT redistributed.
- Original source-only Event site census has 506 quality-qualified public 2025 site coordinates (>=2 Event days, uncertainty <=500m and original coordinate scatter <=2m). Remove the entire set of **107 previously studied** 2023/24 repeated sites BEFORE looking at new 2025 family counts, leaving **399 entirely new site-coordinate groups**.
- Original EEA DK2022 unique facility + direct freshwater both N/P-positive reporting 2022 source contains **203** qualifying outfalls nationally. Distances are geodesic to the nearest registered qualifying FW outfall, not confirmed receiving river water/butterfly larval habitat or 2025 chemistry. Near **<=5km**, remote **>=10km**, intermediate 5–10km. A 'remote' site is NOT pollution-free.
- Of original new399 source/visits: **100 near, 127 remote, 172 intermediate**. Exactly **20** separate WWTP outfall IDs fall nearest to the near 100 coordinates. 5km outfall-to-outfall spatial clustering leaves **19** groups before documentation screening.
- The new original Event-to-ANY-published-Occurrence **coreid** publication audit was completed and Git-committed before new butterfly counts. Require >=2 original 2025 dates/site with any Occurrence row. **286 of 399** pass (near **70**, remote **92**, intermediate **124**); near **18 unique registered outfall IDs** at those sites, **17** outfall connected spatial components at 5km, **14** outfalls with >=2 documented source-near sites, six original country×1° blocks with both near and remote. Source-only gate for outcome access PASSED; no tuning by direction of family counts.
- NOTE: any-Occurrence-linked days are an imperfect observation/reporting-opportunity proxy, not verified complete zero-butterfly events or fully unbiased survey effort. If an original Event has no linked Occurrence, count it as unknown documentation rather than confirmed absence. Counts are not standardized to exact route length, weather, flowering/flight dates or plant habitat.

## New independent-site 2025 measured family individualCount

| eBMS 2025 accepted source-locked counts and reporting opportunities | Registered WWTP near <=5km | >=10km remote |
|---|---:|---:|
| Original new eligible sites (previous 107 excluded) | **70** | **92** |
| Original dates with published ANY Occurrence | **539** | **757** |
| Lycaenidae | **382** | **644** |
| Pieridae | 680 | 2,602 |
| Papilionidae | 0 | 0 |
| Nymphalidae | 3,822 | 3,544 |
| Hesperiidae | 372 | 192 |
| Lycaenidae per documented survey date | **0.709** | **0.851** |
| Lycaenidae SHARE using 3-family denominator only | **36.0%** | **19.8%** |
| Lycaenidae SHARE using 5-family denominator | **7.27%** | **9.22%** |

**Critical contradiction of an assumed beneficial sewage effect:** 3-family pooled relative Lycaenidae fraction is higher near reported WWTPs (pooled near/remote odds ratio **+2.270**), while **the actual accepted Lycaenidae recording rate per documented survey date is LOWER near** (near/remote **0.833**). This is not a butterfly abundance increase. The specific *reference-family composition* drives the sign: Pieridae per documented survey day is 1.262 near vs 3.437 remote (near/remote ratio **0.367**); Nymphalidae instead 7.091 near vs 4.682 remote (ratio **1.515**), and Hesperiidae 0.690 vs 0.254 (ratio **2.721**).

Adding Nymphalidae and Hesperiidae to the 3-family denominator changes the Lycaenidae fraction **from near-high to near-low** (pooled five-family near/remote OR **0.772**). Papilionidae zero counts in both groups, so the five-family alternative effectively contains four recorded families. This is an explicit AFTER-OUTCOME denominator sensitivity; neither denominator choice can establish cause.

## Same-region and genuine source-unit stress tests

| Original country×1° regions with near/remote 2025 published-data opportunity | 3-family near SHARE higher? | 5-family near SHARE higher? |
|---|---|---|
| 55°N 10°E | No | No |
| 55°N 11°E | Yes | No |
| 55°N 12°E | Yes | Yes |
| 55°N 9°E | No | No |
| 56°N 10°E | No | No |
| 56°N 9°E | No | No |

Out of six within-region source comparisons, the original **3-family near-high share occurs in only 2/6**, while the **5-family near-high share occurs in 1/6**. When each original 1° region is weighted equally, the geometric mean regional near/remote **3-family OR=0.908**, despite pooled OR=2.270; the region-equal site-equal logit contrast gives OR=**0.547**. This is a **Simpson-like aggregation sensitivity**, NOT the strict version where every region has the opposite sign.

The five-family region-equal geometric mean near/remote OR is **0.404**. Posthoc 10,000 resamples of the **six** 1° regions yield a descriptive 95% interval for the five-family logOR **[−1.727, −0.061]** (OR approximately 0.18–0.94); three-family logOR interval **[−1.308, +1.313]** (OR approximately 0.27–3.72). These are NOT confirmatory exposure-randomized effect intervals: only six matched regions, nonrandom terrestrial sites, overlapping geographic source catchments and known outcome-driven alternate denominators.

Leave one original 1° region out: 3-family region-equal OR range **0.520–1.494** and even its sign is unstable; five-family region-equal OR remains lower than one (**0.297–0.553**). Source-near documented 70 sites correspond to **17** registered outfall point spatial components at 5km; remote 92 documented sites cluster into **19** connected geographic components at 5km, not 92 independent river ecosystems.

**Observer/region leverage:** the FIVE highest-count remote Pieridae site records contain **1,321/2,602 = 50.8%** of all remote recorded Pieridae individuals. Removing those five sites according to their observed results (posthoc, NOT valid source selection) reduces the pooled three-family near/remote OR from **2.270 to 1.260**; removing ten highest Pieridae remote sites yields OR **0.815**. These outcome-ranked exclusions diagnose leverage, not a confirmatory alternative analysis.

## Benchmark: previously examined 107 repeated sites in the same 2025 source year

The older 107-site panel, excluded from the new399 study, showed 2025 3-family Lycaenidae proportion near **95/(95+148)=39.1%**, remote **585/(585+324)=64.4%** (near-low). With all five families its proportions are near **95/1,331=7.14%**, remote **585/4,011=14.58%** (also near-low). The pooled THREE-family sign therefore **changes** when the standardized 2025 monitoring sample changes, while the FIVE-family direction remains near-low in both disjoint site-coordinate cohorts. The disjoint site sets are NOT random replicated wastewater treatments: some near sites in old and new cohorts share the same 2022 outfall IDs, and recording/landscape opportunity differs by recruitment cohort.

## Scientific interpretation, not a cause claim

**Supported measured finding:** The sign of a coarse WWTP-proximity association with Lycaenidae *relative standardized survey recording composition* depends on (1) what other butterfly families form the denominator, (2) which public sites were eligible based on published recording opportunities, and (3) geographic aggregation. No robust positive Lycaenidae per-survey-date reporting difference was detected in the new399 sample. The five-family relative reporting pattern near-lower is more directionally consistent across the existing two Denmark 2025 site cohorts, but cannot be assigned to wastewater biology.

**Unsupported claims:** sewage water makes Lycaenidae breed worse or better, adult protein/Na puddling differs downstream of an actual outfall, WWTP phosphorus mass changes butterfly density, or the six-region sample supports a global island/continental macroecological law. Spatial confounding by survey routes, woodland/grassland, urban setting, butterfly host plants, plant phenology, effort and EEA reporting coverage remain. A 2022 registered outfall is not measured 2025 effluent release or 2025 water chemistry.

**Next evidence threshold:** independent standardized butterfly monitoring around >=10 different VERIFIED hydrologically connected WWTP receiving habitats and physically comparable upstream/downstream/source-zero controls, butterfly/host resource survey design and receiving-water dissolved nutrient/organic substrate measurements. More recycled iNaturalist 1° shares will not close this mechanistic gap. Do NOT merge source-only exploratory drafts into unrelated chocho evidence or the gesui main branch.

## Source-locked public receipts, exact reproducibility

- [New399 SOURCE/observer opportunity with zero butterfly outcomes first](https://github.com/zuizui0223/gesui/blob/analysis/dk-ebms-2025-new399-source-coverage-20261110/results/ebms_denmark_2023_2024_WWTP_source_only/DENMARK_2025_NEW399_SITE_EVENT_ANY_OCC_WWTP20_SOURCE_GATE_BEFORE_FAMILY.json) — [run SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42767545683).
- [Original licensed 2025 family count with frozen GBIF source-taxonomy](https://github.com/zuizui0223/gesui/blob/analysis/dk-ebms-2025-new399-source-coverage-20261110/results/ebms_denmark_2023_2024_WWTP_source_only/DENMARK_2025_NEW399_SOURCE_LOCKED_ACCEPTED_BUTTERFLY_FAMILY_COUNTS.json) — [run SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42768171519).
- [2025 within-geography, 17-WWTP-component, bootstrap, remote-Pieridae leverage analysis](https://github.com/zuizui0223/gesui/blob/analysis/dk-ebms-2025-new399-source-coverage-20261110/results/ebms_denmark_2023_2024_WWTP_source_only/DENMARK_2025_NEW399_LYCAENIDAE_2022_WWTP_AGGREGATE_REGION_COMPONENT_SENSITIVITY.json) — [run SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42768830966).
- [Separate independent source-only SHA, old107 excluded and six-1° region arithmetic](https://github.com/zuizui0223/gesui/blob/analysis/dk-ebms-2025-new399-source-coverage-20261110/results/ebms_denmark_2023_2024_WWTP_source_only/DENMARK_2025_NEW399_SOURCE_BLIND_6REGION_INDEPENDENT_ARITHMETIC_QA.json) — [run SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42769157021).
- [Independent 3-versus-5-family denominator sign and exact old107 benchmark](https://github.com/zuizui0223/gesui/blob/analysis/dk-ebms-2025-new399-source-coverage-20261110/results/ebms_denmark_2023_2024_WWTP_source_only/DENMARK_2025_NEW399_THREE_VS_FIVE_FAMILY_SIGN_REVERSAL_INDEPENDENT_QC.json) — [run SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42769550575).

**Study stays Draft; no main or chocho merge; public licensed aggregate family records only.**
## Added 2026-11-10: original calendar-month × 1° region joint support (POSTHOC)

An independent original licensed eBMS v1.4 Event/Occurrence parser was used to rederive all original 2025 site-month individualCounts **only after** the frozen 399-site source and 286-site any-published-Occurrence eligibility gates. All 399 original five-family site-year counts, survey days, published-Occurrence-linked dates, original 15,606 Event rows / 21,617 Occurrence rows and original source SHA match exactly. This calendar-month analysis was conducted AFTER original 2025 source-proximity and three-vs-five-family outcomes were known. It is explicitly exploratory/posthoc, not an independent confirmatory study.

**Date/region matching rule:** retain original 2025 country×1° area × calendar-month bins with at least TWO original dates containing ANY published Occurrence in each near (<=5km) and remote (>=10km) band. Each bin's family individualCount is divided by its band-specific documented dates. The direct standardized count rate for both bands uses identical bin weights equal to the smaller of near and remote number of published dates. This weights bins by documentary overlap; it **does not** equally weight the six regions or establish true butterfly counts or recipient-water exposure.

| Matched 2025 calendar-month / 1° geography | All common months | May–August common months |
|---|---:|---:|
| Original source-matched geographic regions | 6 | 6 |
| Matched region-month strata | **32** | **24** |
| Near documented site-observation dates | 446 | 374 |
| Remote documented site-observation dates | 521 | 440 |
| Absolute Lycaenidae recorded/day ratio near/remote | **0.709** | **0.733** |
| 3-family Lycaenidae-share odds ratio near/remote | **1.819** | **2.042** |
| 5-family Lycaenidae-share odds ratio near/remote | **0.586** | **0.598** |

The sign reversal persists **within the same matched date/location support**: including only Lycaenidae, Pieridae and Papilionidae makes the relative Lycaenidae proportion positive near, but adding Nymphalidae and Hesperiidae makes it negative; the absolute per-date Lycaenidae recording remains lower near. Seasonal/geographic matching does not remove denominator dependence.

**Uncertainty and weighting:** 4,000 original 1° geographic-region-with-replacement descriptive resamples of the all-month sample produce 95% intervals for absolute rate ratio **[0.204, 1.250]**, 3-family odds ratio **[0.368, 6.234]**, five-family odds ratio **[0.263, 1.167]**. All intervals contain 1. Leave-one-original-1°-region-out results: 3-family odds ratio ranges **0.738–2.573** (changes sign), five-family **0.399–0.720**. A different original **equal-region-weight** estimate is **0.908 for 3-family**, opposite the new 1.819 overlap-day-weighted result. These are different observational estimands with six spatially heterogeneous regions, not contradictory data or significant wastewater responses. No treated-receptor chemistry, independently assigned downstream source, or causal P/N dose is measured.

**New source receipts:** [Original licensed month × 1° area full results and 4000 resamples](https://github.com/zuizui0223/gesui/blob/analysis/dk-ebms-2025-new399-source-coverage-20261110/results/ebms_denmark_2023_2024_WWTP_source_only/DENMARK_2025_NEW399_1DEG_MONTH_MATCHED_LYCAENIDAE_3VS5_FAMILY_QC.json) — [original data replay SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42772627586). [Independent source-SHA-checked arithmetic and leave-one-region reproduction](https://github.com/zuizui0223/gesui/blob/analysis/dk-ebms-2025-new399-source-coverage-20261110/results/ebms_denmark_2023_2024_WWTP_source_only/DENMARK_2025_NEW399_REGIONMONTH_3VS5_INDEPENDENT_QA.json) — [independent QA SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42772945381).

**Conclusion unchanged:** source-relative positive/negative differences in original photos and standardized recorded individuals are sensitive to denominators, geography, survey timing, reporting inclusion and weighting. They do not identify WWTP nutrient benefits or harms. Keep Draft; no main or chocho merges.
