# Season-controlled spatial–temporal decoupling in measured river phosphorus and benthic integrity: same 263 French river stations, 2018–2022

**13 November 2026 · Original official measured TP sampling dates, same EEA I2M2 field sites and native HydroBASINS H7×year · Source eligibility fixed before month-aligned TP coefficients · EXPLORATORY POSTHOC relative to previous spatial/temporal findings · Draft/unmerged**

## Executive scientific result

**Matching the calendar months of TP observations in 2018 and 2022 does not recover a robust association between within-river TP changes and benthic I2M2 quality changes over four years.** This survives unchanged original physical station sample (263 genuine same WQ+EEA site IDs), unchanged 73 Level7 catchment units, unchanged 24 Level6 hydrologic parent geographic blocks, fixed environmental and study-effort covariates, independent NumPy code reproduction, and an exploratory leave-one-H6-out analysis.

The earlier strong **2018 and 2022 spatial** real TP–I2M2 EQR negative correlations are retained. Neither that spatial correlation nor this null temporal relationship proves that nutrients do or do not cause benthic community change, and neither identifies specifically WWTP-derived P.

## Previously authenticated scientific baseline (PR #71)

Earlier source-first real original author Peters et al TP.2 archive (DOI `10.25400/lincolnuninz.30546371.v5`, file ID70896288, official MD5 `b596dbe9833c49426147a99f10f3fe31`, size25,003,722 bytes, original SHA256 `26f2af11f4972987ee6f7c3f497ad902820f145851f71919d4f45b589e1d809e`) and official EEA I2M2 source results identified **263 same physical WQ/EEA stations** with valid **2018 and2022** positive TP measurements on ≥4 sample dates in ≥3 months each calendar year, original country FR/method I2M2/native reach equal across years, comprising **73 H7 subbasins /24 authentic H6 hydrological parent groups**. Original same station source selected outcome-blind for 2018/22 new differences, SHA256 `2190ac905326f4f2f4f278eb0debb4a7bfda7c4c71347d4941b54cc19d3de264`. Original biological/nutrient source results in [PR #71](https://github.com/zuizui0223/gesui/pull/71).

Original unseasonally matched `Δln annual TP median` vs `ΔI2M2 EQR` across n73 H7: Spearman **ρ−0.12597**, site n263 source ρ+0.00320; original same-H7 environmental H6 5fold Ridge10 incremental EQR-change MSE improvement **−0.409%** (95% H6-parent bootstrap [−6.290,+3.782]%). 2018 year-specific TP↔EQR spatial n73 H7 ρ−0.73285 and +30.81% added predictive MSE improvement; 2022 year-specific n73 H7 ρ−0.69719 and +22.67% improvement, all same sites. TP spatial rank stability 2018 vs2022 ρ+0.954, I2M2 ρ+0.849; no reselecting sites based on TP or EQR sign.

## New SOURCE-ONLY season/month matching, frozen before new numeric effects

Original [Peters 2018/2022 official measured TP sample calendar-month support (SOURCE ONLY) run #43378522025 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43378522025); source-only original month inclusion SHA256 `a02a9f724fb1b1671a9f248ed169620bc25db8b9a5914be566f6728e25b2c4f1`.

Definition fixed before opening the new response:
- Original same **263** WQ IDs and source-verified same EEA FR/I2M2 stations in both 2018 and2022; compare only calendar month numbers **observed in BOTH years** with authentic dated, strictly positive TP measurements, without relying on station start/end metadata.
- At least **3 shared measured calendar months** in each station, plus earlier ≥4 positive samples over ≥3 months in each year; research gate ≥100 native sites, ≥30 H7 and ≥10 H6 parent groups. All 263 sites and all 73 H7/24 H6 PASS without losing field stations, so subset effects are not an artifact of dropping unsupportive original rivers.
- Original shared month support **3–12 months** per site, histogram: 3 months n1, 4 n2, 5 n10, 6 n51, 7 n32, 8 n84, 9 n10, 10 n2, 11 n17, 12 n54; **250/263 sites have ≥6 common sampled calendar months**, no original duplicate sampling dates. Original 2018/2022 *actual TP magnitudes and EQR numerical values were NOT read to select shared months*.

## Source-locked month-balanced environmental response

For each physical source station, for each month number m sampled in both years, take the genuine TP median across positive sampled dates separately for 2018 and2022. Calculate the **equal-month seasonal-matched change**
`ΔlnTP_month = mean_over_common_months [ln(median TP_2022_month) − ln(median TP_2018_month)]`.
Each shared month carries equal weight, not one observation day or one year-sampling batch. This removes observed calendar-month **composition** differences but does not match same day, precipitation, temperature or river flow conditions and does not identify a response lag.

Compare with previous `ΔlnTP_annual = ln(annual TP2022 median)−ln(annual TP2018 median)` **on the exact same 263 original stations, 73 H7, 24 H6**. All **263** original annual 2018 and2022 TP medians were independently recomputed from the original source ZIP with exact agreement to archived original 263 field-station records.

The month-adjusted and annual TP change definitions are related but not identical: n73 H7 **Spearman ρ+0.59552**, Pearson +0.49111. Thus seasonally standardized TP is substantively a different measurement summary, despite having the same underlying physical station data.

### Original aquatic source outcomes (same n73 H7)

Baseline controls held constant for both new and old TP-difference predictors: original native H7 area, source original HydroATLAS upstream `crp_pc_use` cropland, `urb_pc_use` urban, `log1p(dis_m3_pyr)` natural modeled outlet river flow, original station mean modeled discharge, n original repeated physical TP/I2M2 sites, and **2018→2022 change in TP sampling effort restricted to the same actual common months**. Original EEA method/country are constant FR I2M2. True source-frozen 24 H6-parent group5fold, train-only scaling Ridgeα10; n73 independent H7 units, never 263 pseudoindependent catchments.

| Actual measured 2018→2022 TP delta relationship | Original raw annual median on same 263 stations | NEW common calendar-month medians on same 263 stations |
|---|---:|---:|
| H7 TP change vs EEA I2M2 quality change raw Spearman | **−0.12597** | **−0.12171** |
| H6-parent 5fold EQR-change MSE gain from TP change, identical seasonal-effort baseline | **−0.2452%** | **+0.2591%** |
| 2,000 H6-parent bootstrap descriptive95 for extra TP prediction | [−6.123%,+4.051%] | **[−3.992%,+4.104%]** |
| Full 24fold leave-one-H6-parent-out extra TP prediction | **−1.1173%** | **+0.4019%**, H6 boot95 [−3.856%,+4.331%] |

New source station-level n263 original ΔTP_month versus original ΔI2M2 raw Spearman **ρ+0.01772** (site records are nested in H7), n73 source H7 Spearman −0.12171, source within-H7 station-demeaned monthTP→ΔI2M2 slope +0.01301 (unadjusted). Source 24-H6-parent-cluster descriptive adjusted coefficient of ΔI2M2 EQR per 1SD original monthTP change **−0.01700** 95 [−0.04135,+0.00735], again includes zero.

Source control containing 2018 initial measured *same-month* TP level plus ΔTP_month yields only **+0.5486%** extra heldout prediction (H6 bootstrap descriptive95 **[−1.600%,+3.397%]**). This is a secondary baseline sensitivity, not a source-specific response mechanism.

**Interpretation: The previously near-null 4-year environmental–benthic co-change is NOT explained simply by observing different calendar months in the two years.** The estimated TP change moves when month-composition is standardized, but its conditional covariation with benthic change remains weak, and the numeric sign of tiny MSE gain flips around zero within uncertainty. Do not claim “common month adjustment demonstrates positive TP effect”.

[All authentic monthly measured TP, unchanged n263 same river source and H6 geographic original outcomes SUCCESS #43378917089](https://github.com/zuizui0223/gesui/actions/runs/43378917089). The run artifact includes original n263 month-balanced TP/EQR station-level and n73 native H7 GBIF/EEA source CSV plus complete results JSON.

### Whole-native-H7 plant and butterfly source-indexed outcomes

Same original whole H7 GBIF `2017–2019→2020–2022` indexed-speciesKey `S50` plant richness change and relative `Lycaenidae` odds within FIVE accepted butterfly family occurrence reports, separate observation-opportunity-adjusted models (original Plantae, Animalia and five butterfly-family source counts, plus physical upstream environment). **These are not the same physical field station's riparian host plants or true pollination!**

| Source-index terrestrial outcome, n73 H7 with both GBIF periods | H7 raw ρ with actual ΔTP_month | Add actual ΔTP_month to 24 H6-parent group5fold prediction, signed MSE gain |
|---|---:|---:|
| Original indexed plant speciesKey rarefied richness S50 change | +0.01752 | **−1.389%**, bootstrap95 [−4.282%,+0.932%] |
| Original Lycaenidae five-family report log-odds change | −0.04859 | **−3.495%**, bootstrap95 [−6.175%,−1.000%] |

Strict 24-H6-parent leaveout TP added gains also worsen (plant **−2.091%**, Lycaenidae **−1.721%**) in the recorded source cohort. Thus no detected co-change of measured river TP and the indexed terrestrial responses; the butterfly index does NOT observe adult abundance or actual plant pollination.

## Independent numerical and EEA biological score provenance QA

A separate source-matched **manual NumPy Ridge normal equation solver** (vs prior sklearn Ridge/SVD) exactly reproduced **seven original benthic model MSE values** (baseline, monthTP, annualTP, LOHO, 2018 base) and both source original GBIF response models. Source n73 H7 raw monthTP-EQR rank and n263 site rank reproduced. [Independent numerical and EEA source-index QA SUCCESS #43379428285](https://github.com/zuizui0223/gesui/actions/runs/43379428285).

The source EEA scoring audit found both original pipelines explicitly read the **official `resultNormalisedEQRValue`**: original 2018 `scripts/gesui_2018_independent_bqe_EQR_outcomes.py` assigns it to `numeric_EQR_2018`, while original 2022 `scripts/gesui_model_2022_eea_aquatic_EQR.py` assigns the SAME original field to `normalised_EQR`. Both original source CSVs passed exact Git blob check and **263/263 identical FR/I2M2 station/year method/river IDs** were verified. This excludes a simple accidental subtraction of a raw 2018 EQR from a normalized 2022 EQR.

**Important unresolved technical point:** The original 2018 `resultEQRValue` (unnormalized) is missing/non-numeric at **263/263** of these FR/I2M2 source stations. Thus the underlying original scoring transformation and absolute interyear calibration cannot be independently reconstructed from those station raw values; a same method/field label does not itself guarantee unchanged I2M2 standardization, national reference classes or field sampling protocols. This remains a **measurement quality limitation**, not a fabricated confirmation of full calibration parity.

## Ecological synthesis and what is/isn't supported

| Same biological pathway layer | Evidential status in this actual original matching |
|---|---|
| Historic modeled registered WWTP outflow → actual receiving-river TP | **No robust source-specific prediction/counterfactual** in these original 74 French H7; actual reported WWTP P was positive at only 4 of 273 paired 2022 sites, not sufficient to fit pollutant-source attribution or causal mediation |
| **Measured river TP ↔ benthic I2M2 condition across H7** | Strong negative **spatial** association persists in 2018 (ρ−0.733) and 2022 (ρ−0.697), with positive H6-heldout out-of-region prediction using actual observed TP |
| **Measured within-site ∆TP → I2M2 ∆ over 4 years** | Weak/mixed original observational evidence even with fully matched calendar months, n263 physical same-station /n73 H7 region units, geographically heldout +0.259% [−3.992,+4.104] |
| Actual river TP change → indexed plants / relative five-family butterflies | No consistent positive predictive information in either raw annual or season-matched original H7 observation-index changes |
| Actual aquatic insect emergence, animal food, floral hosts, insect visitation, real pollination or seed set | Not measured jointly on original FR73 H7 and cannot be inferred from EQR/GBIF indices |

**Publication-safe working question:** *Why do stable spatial gradients in river nutrient concentration predict benthic integrity while contemporaneous within-river change does not?* Potential between/within variance, legacy habitats, lagged biological response, climatic discharge and standardized sampling constraints have not yet been discriminated. This is a precise boundary on **space-for-time substitution**, not proof of ecological memory or delayed effect.

Keep all original positives and nulls together, new unmerged Draft stacked on PR #71, no main/chocho modifications, no eBird or fabricated serial ecosystem-service mediation.
