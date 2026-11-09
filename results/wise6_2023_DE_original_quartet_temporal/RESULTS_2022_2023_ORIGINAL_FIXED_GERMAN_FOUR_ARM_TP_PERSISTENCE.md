# 2022–2023 Germany: genuine observed river total-P in identical upstream/downstream WWTP vs registered-source-zero monitoring quartets

**Updated November 9, 2026. All source and numerical stages successful; exploratory/posthoc relative to original positive 2022 findings.** This is a WISE6 real river-water total phosphorus (TP) monitoring test with repeated **original four-arm river sites**, not a butterfly or mud-puddling experiment.

## 1. Exact public source and 2023 sampling-frame gate before numerical TP was opened

- Original full EEA Waterbase Water Quality ICM 2026 `Waterbase_T_WISE6_DisaggregatedData` ZIP **2,024,284,753 compressed bytes**, **36,043,700,357 uncompressed CSV bytes**, full ZIP SHA256 `2a2ef2a43a16ec068b7974a80ae5fca49e44697660a46ce3b91fb1813abe62d2`. **102,569,057** official source records streamed; 2023 Germany river-water TP `CAS_7723-14-0` with unit `mg{P}/L` **2,818 original sample rows** across all original sites.
- Freeze **all original 8 matched source/control upstream/downstream quartet site groups (32 WISE6 stations)** from the previous 2022 analysis. The 2022 site selection used original public EEA 2022 positive reported annual N/P FW outfall mass and a SHA-frozen third-party converted HydroRIVERS directed river graph, **not** original TP numerical values.
- Read **2023 source-only sample dates and station hash IDs BEFORE any new 2023 numerical TP**: 322 original 2023 TP samples across the frozen 32 sites; 26 sites have valid recorded 2023 samples and ≥4 original sample days over ≥3 months.
- Original fixed criteria: all 4 original stations sampled at least 4 days over 3 months; ≥3 exact **YYYY-MM** months with measurements at all 4 arms; ≥5 complete original quartets in ≥4 country×2° source blocks before 2023 numeric reading. **Passed: five original quartets (#1,#2,#5,#6,#7), four original source 2° blocks**. Exact same-day ≥3 days in ≥4 original quartets/3 blocks **failed: zero quartets**; NO same-day effect was estimated or invented.
- Source-only 2023 file committed ahead of numerical TP code's first run: [year 2023 source-only 32 stations / 8 quartet original date gate SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42606865044). Original 2022 8 quartet date source Git blob `c09582bee2c09ad038bd9fd16926d837644da075`.
- **Explicitly 2023 is not a separate period of verified WWTP source concentrations**: the source/interposed status is from the original **2022** EEA treatment plant report; actual WWTP 2023 month-specific P effluent delivery is unknown. A persistent 2023 TP river gradient might reflect recurring discharge, persistent agricultural/tributary sources, natural landscape structure, or study/reporter selection.

## 2. Real observed 2023 TP quantified values — five original four-arm matched quartets

Only original public 2023 WISE6 site records with TP code `CAS_7723-14-0`, original `mg{P}/L` water-matrix unit, original sampling date and exact country+site ID+scheme hash. Source value `resultObservedValue` and `resultQualityObservedValueBelowLOQ` were opened only after site-date sampling gate passed. None of the 322 original 2023 selected-site rows was below LOQ or rejected by source observation-status filter. Of the original 32 sites, **26** met four-date/three-month numeric support; all **5** prequalified four-arm quartets also satisfy numeric support.

Actual river TP site daily median → within-month median of daily measured TP → source-interposed pair downstream-minus-upstream log(TP+0.001) gradient, **minus** same-month registered-qualifying-source-zero reference downstream-minus-upstream gradient. First average months within each quartet; then average five **independent four-arm quartets equally** (months are not 12 independent ecosystems).

| Real 2023 TP source/control same-calendar-month four-arm contrast | Value |
| --- | ---: |
| Complete original four-station quartets | **5 (4 original geographical blocks)** |
| Original all-year published 2023 German RW TP samples | **2,818** |
| Samples at original frozen 32 stations | **322** |
| Average source-minus-registeredzero downstream TP log-gradient | **+0.225208** |
| Median original quartet gradient | **+0.204391** |
| Quartet directions | **3 positive / 2 negative** |
| Descriptive 9,999-quartet bootstrap 95% range | **[+0.01788,+0.43253]** |
| **Exact two-sided independent-quartet sign-flip reference fraction** | **0.25** |
| Remove-one-quartet mean range | **[+0.15449,+0.29250]** |
| Strict exact SAME-DAY four-arm test | **FAILED pre-outcome gate; no estimate** |

**Important:** the descriptive bootstrap interval excludes 0, but the more appropriate *discrete* exact 5-quartet sign-flip reference is **0.25** (no statistically persuasive evidence). Those numerical summaries are not necessarily contradictory because a 5-unit bootstrap distribution poorly represents exact small-sample sign uncertainty; no randomized wastewater treatment occurred. The 2023 positive mean direction is not a confirmed effluent-mediated river TP effect.

[Original EEA2026 2023 numeric TP source-replay and original 5-quartet same-month relative gradients SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42607880380). Exact original source/LOQ and per-quartet 2023 monthly values: [2023 original full numeric/quality result](2023_ORIGINAL_8QUARTETS_SAMEMONTH_REALTP_NUMERIC_SOURCE_ZERO.json).

For reference, the original separate **2022** same-month 7-quartet estimate was **+0.168592** with descriptive 95% [−0.07699,+0.44281] and exact sign-flip **0.296875**; these 7 and 5 quartets are not identical analysis samples and may have different seasonal windows. Never call their numerical difference a time effect.

## 3. Stronger pre-numeric-frozen seasonal year-to-year control: same quartet AND same calendar month of year

**Before seeing any 2023 numeric TP**, the original **2022/2023 source-only four-arm calendar-month intersection** fixed four eligible same-exact-site original quartets: **#2 (7 matched calendar months), #5 (8), #6 (9), #7 (8)**, spanning **three original geographical 2° blocks**. Each original site/arm of a quartet is identical in both years, and only the same calendar months in both 2022 and 2023 are included. Original source-only rule ≥4 quartets in ≥3 original 2° blocks, ≥3 same calendar-month overlaps per quartet, was **passed and committed before comparing new TP numeric**. No 2023 year-specific WWTP discharge mass was introduced.

An independent OFFLINE script reopens only the already committed old and new **per-month measured WISE6 TP medians**, checks original quarter IDs/site hashes, recomputes four-arm log-TP gradients from the *four actual concentration values*, and holds only months-of-year in both years.

| Exactly same original 4 quartets and same calendar months | Result |
| --- | ---: |
| 2022 four-arm source-minus-registeredzero mean logTP gradient | **+0.296157** |
| 2023 four-arm source-minus-registeredzero mean logTP gradient | **+0.239772** |
| **2023 minus 2022 change** | **−0.056384 log units** |
| Quartet-level change directions | **2 increases / 2 decreases** |
| 9,999-quartet bootstrap 95% interval for year change | **[−0.17100,+0.05676]** |
| Exact two-sided four-quartet year-change sign-flip reference fraction | **0.625** |

The direction of the same-month source-vs-reference TP gradient is positive in BOTH years for the same 4 original quartet population; **its estimated change is small, negative, and highly uncertain**. This is better described as **persistent environmental TP geography** than as evidence of a new 2022–2023 wastewater emission effect. A persistent WWTP source could be one mechanism, but so could persistent upstream tributary, land use, original source/monitor selection and wastewater-reporting differences. The original study does not measure daily 2022–2023 effluent delivery, and 4 independent quartets are far too few for a causal attribution.

[Independent source/date/numeric offline reproduction and exact same-month annual comparison SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42608436317) · [Full four-quartet 2022–2023 same month-year differential and every monthly original source TP value](2022_2023_SAME_ORIGINAL_QUARTETS_SAME_CALENDAR_MONTH_REAL_TP_OFFLINE_REPLAY.json).

## 4. What is positive and what remains unavailable

**Actual new empirical progress:** original 2023 individual water-chemistry measured TP data and original 2022 same-river WWTP/sourcezero four-arm sites can be linked; five sets have repeat sampled month support; four original sets have identical cross-year season/arm comparisons. The TP relative source-geometry-labelled gradient has a positive *average direction* on the same quartet group in 2022 **and** 2023, showing **temporal persistence of a measured river TP difference at these fixed sites** rather than a 2022-only artifact. This is a valid descriptive finding. It is **not** robust statistical evidence that WWTPs CAUSED the difference.

**Limitations and mandatory nulls:**
- Source/control assignment is from original EEA official **2022** WWTP outfall and **2022 annual** reported N/P mass, with strict different-reach interposition and 50km graph condition. No 2023 EEA WWTP annual P mass measured, no actual day/month effluent flux. Control sites have no qualifying *registered* intervening WWTP, not an untreated environment.
- Public HydroRIVERS v1 converted third-party GeoParquet source geometry SHA-pinned, but its topology not independently original-official shapefile verified; same-reach within-stream ordering unresolved.
- No strict same-day all-four-site comparison in 2022 or 2023 (source gate failed), so upstream/downstream monthly medians can be days apart and confounded by floods/seasonal flow.
- Source/zero river pairs were **not randomly assigned**; 2022 both-positive N/P WWTP source reporting incomplete; land cover, crop fertilizer, background rivers and tributary inputs not fully controlled at exact catchment.
- No original Lycaenidae real individual abundance, adult puddling protein/Na substrate use, host plant quality or mutualist ant behavior was measured. Original `chocho` 22/25 Lycaenidae *potential introduced-host geographic larval resource* expansion and `gesui` 2018/2022 aquatic benthic EQR prediction remain separate.
- WISE6 2023 selected 322 source rows original `resultObservationStatus` metadata was audited; blank status is not independent expert laboratory verification. No LOQ censoring in these original selected rows.

## 5. Recommendation for next true causal-physics test

Original source→receiver water TP is now a **reproducible descriptive repeated-year same-river contrast**, but inference is underpowered and lacks same-day WWTP treatment data. To assess mechanism, obtain original daily/seasonal WWTP P and in-river flow at the same monitored reaches, harmonized dates, and replicate independent river catchments or longer year series under the pre-specified same-river/sourcezero design; verify river network topology against original official HydroRIVERS. **Do not interpret historical WWTP annual P loads or neighboring Lycaenidae opportunistic photo reporting as observed insect nutrition.**

This remains an unmerged exploratory Draft. 2022 original source positive is a directional tendency, 2023 source positive likewise; both exact sign-flip fractions remain nonconfirmatory, with no evidence for a 2022→2023 change.
