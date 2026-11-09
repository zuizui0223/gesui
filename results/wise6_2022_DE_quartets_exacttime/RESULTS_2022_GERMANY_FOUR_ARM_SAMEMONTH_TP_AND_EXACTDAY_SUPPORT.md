# Real 2022 German river TP: four-station SAME-MONTH WWTP-interposed vs registered-zero gradients

**November 9, 2026 — original EEA2026 WISE6 numeric source fully executed and independently rechecked from original derived values.** This entire temporal check is **POSTHOC** because the original source-connected 7-quartet annual-TP result was known first. The original 2022 WWTP source and river monitoring locations were pre-frozen, but not randomly assigned. This is **not** a sewage treatment causal effect, not butterfly puddling or true Lycaenidae abundance, and not a new P mass measurement at the actual water-quality sample date.

## 1. Source-only sampling-time audit, committed before opening new within-month TP values

- Original study (Draft [PR #32](https://github.com/zuizui0223/gesui/pull/32)) frozen eight original German four-station quartets: 8 strictly WWTP-interposed two-monitor upstream/downstream river paths + 8 registered-source-zero two-monitor comparators, **32 nonreused original monitoring site IDs**, 6 original 2° blocks. The first annual TP numeric study found **7 complete quartets** (one no-overlapping-control-month).
- Original source-only group file Git blob `8bd67be25968311ca16be4026888daa783fbe69f`. Original independent 32-site numeric-QC full SHA256 `28b679184cbebf6376f1c65d7c160f7511f057ff22787e98ff3fa265b617b261`.
- Source original WISE6 2026 `DisaggregatedData` **2,024,284,753-byte** ZIP, SHA256 `2a2ef2a43a16ec068b7974a80ae5fca49e44697660a46ce3b91fb1813abe62d2`, verified **102,569,057** original individual water-quality records. Source-only first pass looked at **2022 original date and country+scheme+site ID metadata only** for the exact original 32 sites — NO numeric `resultObservedValue` or LOQ concentration.
- **Exact four-station SAME DAY** pre-outcome gate: ≥3 independently shared dates in ≥4 quartets, ≥3 original 2° blocks: **0/8 quartets qualify**. Original quartet #8 has **2** matching all-four dates, not enough. **No same-day TP contrast is reported.**
- **Same CALENDAR MONTH** pre-outcome gate: ≥3 months in all four arms in ≥5 quartets, ≥4 original 2° blocks: **7/8 quartets qualify, across 5 blocks**. The original #1 has 0 common four-site months and remains excluded without replacement. The total original all-four-arm month periods in the seven qualified quartets are **66** (7+11+11+11+10+10+6).
- [Original all-32 station sampling-date metadata/source-only code and workflow SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42601562258). Full source source file `2022_DE_8WWTP_QUARTETS_EXACT_SAMEDAY_MONTH_SAMPLE_SOURCE_ONLY_GATE.json`.

## 2. New REAL river TP numeric matched four-station monthly result

*After* the original time eligibility was committed, source-pinned EEA2026 individual 2022 original stream was re-opened and quantified station/day/year TP values computed from original `CAS_7723-14-0`, `procedureAnalysedMatrix='W'`, `resultUom='mg{P}/L'`. For original 32 source-frozen sites, **472 original sample records**, **0 below-LOQ** flags, **32/32** valid original site-year numeric quality records. All newly recomputed **annual site medians exactly reproduce** the previous original 32-site annual TP QC, confirming no new station/outcome replacement. All **472/472 original `resultObservationStatus` fields are blank**; under the Eionet reporting dictionary the `A` flag is optional confirmation for unusual results, so blank alone is not proof of invalidity, but independent EEA expert-QC status remains unverified.

Measurement grain: original 2022 TP valid same-day replicated assays → site daily median; site **month median of daily medians** → for each of the original seven quartets, only the original **four stations all having a valid month** enter the comparison. The result for each month is:

- registered-source-interposed pair: `log(TP_downstream + 0.001) − log(TP_upstream + 0.001)`;
- registered-source-zero two-site reference: identical logTP upstream/downstream gradient;
- *source-interposed minus registered-zero* observed monthly gradient; average over common months **within original quartet**, then equal weight all seven original quartets.

The original +0.001 mg P/L offset is copied exactly from the source-frozen annual TP analysis; nothing was tuned using new monthly outcomes.

| Original source-control monthly quartet | Original four-site same-month count | Mean monthly source-minus-control log TP |
| --- | ---: | ---: |
| #1 — source-only timing unsupported | 0 | **Not estimated** |
| #2 | 7 | −0.137655 |
| #3 | 11 | −0.301053 |
| #4 | 11 | +0.221896 |
| #5 | 11 | +0.094239 |
| #6 | 10 | +0.822932 |
| #7 | 10 | +0.465056 |
| #8 | 6 | +0.014727 |

**Aggregate, same original 7 quartets with equal quartet weighting:**

| Criterion | Previous annual site-median design | New four-station SAME-MONTH design |
| --- | ---: | ---: |
| Original complete source-v-registeredzero quartets | **7** | **7** |
| Average source-minus-control logTP river gradient | **+0.090448** | **+0.168592** |
| Quartets with favorable / unfavorable source contrast | **5 / 2** | **5 / 2** |
| Fixed bootstrap 95% descriptive uncertainty (quartet resample) | [−0.181608,+0.363274] | **[−0.076987,+0.442806]** |
| Exact matched-quartet 2-sided sign-flip *reference fraction* | 0.546875 | **0.296875** |
| Leave one original quartet out: mean range | [−0.00810,+0.17538] | **[+0.05954,+0.24687]** |

**Interpretation:** The same-month source-v-registeredzero comparison retains a positive *directional* average and, unlike the original annual-site medians, the leave-one-quartet-out mean stays positive. But 7 original source-controlled quartets across only **5** country×2° blocks provide wide bootstrap uncertainty **including zero** and a sign-flip fraction of **0.297**. This is **not compelling statistically robust source-specific eutrophication evidence**, much less proof of WWTP causality. One monthly quartet (#6) is relatively high (+0.823), and contrasts in #2 and #3 are negative. Positive matched contrasts need not imply actual downstream TP increased: the control reach may simply have declined more.

Crucially, **the four original sites were not sampled on exactly the same days in sufficient numbers**. Monthly matching controls seasonality more than yearly medians but cannot resolve individual WWTP effluent pulses, rain/flow or exact dates. Source input P reported annually in 2022 and the WWTP location may be on a tributary within the selected reach path; a same-month TP gradient is not a direct chemical dose. The registered-source-zero comparator means no **qualifying registered** source in a selected flow path, not no wastewater or other diffuse urban/agricultural sources.

## 3. Reproducibility, EEA quality and scope

- [Original WISE6 source 8-quartet sampling DAY/MONTH QC SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42601562258).
- [Original 36GB EEA2026 WISE6 numerical 66 quartet-month TP source/zero comparisons SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42603282019), full `2022_DE_ORIGINAL_EEA_WISE6_SAMEMONTH_REAL_TP_WWTP_SOURCE_MINUS_SOURCEZERO.json`.
- [**Independent offline source/result replay SUCCESS**](https://github.com/zuizui0223/gesui/actions/runs/42604048677), `2022_DE_SAMEMONTH_FOUR_ARM_TP_INDEPENDENT_OFFLINE_REPRO_QC.json`. Recomputed EVERY original four-arm month gradient directly from stored source-original TP values, month-bounded source-only eligibility, all seven equal-quartet contrasts, **9,999 bootstrap confidence endpoints**, exact sign-flip and leave-one-out range. All agree with committed numeric source result to numeric tolerance.
- [Official Eionet WISE6 sampling/LOQ/status data dictionary](https://dd.eionet.europa.eu/datasets/latest/WISE-SoE_WaterQualityICM/tables/DisaggregatedData). Original WISE6 site/sample IDs, source year and columns respected; extra expert-QC and source `metadata_observationStatus` flags have **not** been independently validated beyond LOQ and raw status count in this study.
- HydroRIVERS v1 **third-party converted** network SHA256 `5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f` has NOT been compared byte-for-byte or topology-equivalent to the official source shapefile. Original EEA public WWTP annual discharge masses, within-reach order, true receiving-stream flow, background agricultural and industrial P, and day-aligned original pollution dose remain unresolved.
- **No actual Lycaenidae adults/puddling photographs, Na/protein uptake, larval fitness, ant interactions or abundance were measured on these 7 river quartet sites.** The original `chocho` 22/25 Lycaenidae introduced-host **potential resource geographic range** positive, `gesui` EQR predictive positives and independent 2021–2023 iNaturalist photo-reporting correlations remain separate and untouched.

**Scientific claim at most:** temporal alignment at the same calendar month produces an observational 2022 source-interposed/control German river TP gradient that is directionally positive in five of seven quartets, but remains uncertain and noncausal. Prior four-site exact-date gate fails. No evidence for a direct wastewater nutrient benefit to Lycaenidae.

All analysis lives on a stacked Draft PR, no main merge or rewriting previously reported estimates.
