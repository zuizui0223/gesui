# Lycaenidae × reported WWTP phosphorus: three-year balanced-panel DIFFERENCING audit (2026-11-09)

**Status: exploratory, strictly POSTHOC.** This test was designed and executed **after** the 2022 positive WWTP-source/photo association and the 2021/2023 retrospective analyses were already known. It is an additional falsification/interpretation audit within existing [Draft PR #27](https://github.com/zuizui0223/gesui/pull/27), not a new preregistration, treatment-control experiment, or butterfly census. No new public observations were requested.

## Immutable inputs and selection

Uses already source-hash-frozen **original** public iNaturalist family-aggregate 2021, 2022, 2023 JSONs and original 2022 EEA reported WWTP N/P/2021 WorldCover record in the same analysis branch. The [standalone Python script](../../scripts/gesui_lycaenidae_2021_2023_posthoc_photo_composition_differencing.py) reads original local JSON files, verifies original 2021/2023 source link to the 2022 land/EEA panel, validates year-specific 3-family sums, and refuses changed source cohorts.

- **42 identical country×1° cells** (DE 27, AT 15), **22 country×2°** geographic blocks; identical source-positive country×1° frame, 2021 family photo support ≥15, 2022 and 2023 family photo support ≥20; no sites are chosen by the sign of their association.
- Total original photographs of only the three available families: **2021 Lycaenidae 2,237/4,356**, **2022 2,988/5,714**, **2023 2,926/5,780**. Three-family photographic composition **is not** number of individuals, survival, reproduction, true relative species abundance, or behavioral puddling.
- Predictor in **all three** years is the **same 2022 annual EEA reported WWTP P mass** in each country×1° cell. It is neither a 2021 nor 2023 measured water concentration, nor a known treatment introduced in 2022.

## New calculation

Response for each family/year is \`log((family_photos + 0.5)/(other_two_family_photos + 0.5))\`. Unweighted cell-level OLS uses country, latitude, longitude, log(1+2021 independent Plantae photo reports), original 2021 WorldCover built fraction, cropland fraction, log terrestrial country-cell land area, plus P-source log(1+mass) standardized to one SD **within each country**. Same design is applied to 2021, 2022, 2023 and their differences on exactly the same 42 cells. No stepwise variable selection.

**There are 22 spatial groups**; independent geographic evidence must not be confused with 42 independently sampled biological populations. Descriptive confidence intervals below use **3,000 block-with-replacement resamples**, with original fixed source scaling and deterministic LCG seed 20261109; the interval is not a preregistered confirmatory test or a causal CI.

| Lycaenidae photo-share estimand | Adjusted log-odds slope / country-specific 1 SD logged 2022 P source | Descriptive 2° cluster bootstrap 95% |
|---|---:|---:|
| 2021 photo composition | **+0.0610** | [−0.2325, +0.2085] |
| 2022 photo composition | **+0.1541** | [−0.0224, +0.2828] |
| 2023 photo composition | **+0.1155** | [−0.0968, +0.2312] |
| **2022 − 2021** within-cell logit change | **+0.0931** | [−0.0919, +0.3242] |
| **2023 − 2022** within-cell logit change | **−0.0385** | [−0.2076, +0.1129] |
| **2022 − mean(2021, 2023)** photo-share anomaly | **+0.0658** | [−0.0845, +0.2410] |

The 2022 level, 2022−2021 change, 2023−2022 change, and middle-year anomaly each retain their **fitted coefficient sign** after leaving out any single 2° block (22/22), but **all their cluster-resampling intervals include zero**. Leave-one-region sign consistency is not enough for independent statistical confirmation. Annual photographic share shifts cannot be attributed to an unmeasured change in annual exposure.

## Does the positive Lycaenidae sign identify family-specific benefit?

No. The exact same compositional model on the other two butterflies in 2022 gives **Pieridae −0.0829** (descriptive 95% [−0.2547, +0.0994]) and **Papilionidae −0.2320** ([−0.4858, +0.0766]). Papilionidae has a similarly negative point estimate in 2021 and 2023. Because the three family photo shares must sum to 1, these are **compositional counterparts**, not independent ecological control outcomes. A positive Lycaenidae log-odds coefficient does not necessarily mean more Lycaenidae occur; it can arise from fewer relative Papilionidae photographic reports, different habitats or observer effort. The family fractions are also not vegetation and puddling measurements.

## Reconcile with earlier source-frozen studies

1. [2022 macro landscape PR #24](https://github.com/zuizui0223/gesui/pull/24) originally showed +0.1862 same-year P-source photo logit slope on **44 cells** without actual land-cover adjustment and +17.82% subsequent posthoc regional predictive MSE gain. This new 42-cell 3-year balanced panel changes the unit/sample and controls; **do not substitute its +0.1541 coefficient into the original 44-cell result**.
2. [WorldCover PR #25](https://github.com/zuizui0223/gesui/pull/25) retained posthoc **+13.67%** P incremental MSE prediction after land-cover adjustment on 43 cells, but 2° resampling range **[−14.18%, +33.04%]**, and the registered WWTP count benefit disappeared (−1.52%).
3. [Year check PR #27](https://github.com/zuizui0223/gesui/pull/27) found 2023 extra P prediction **−7.41%** after controlling 2022 Lycaenidae photo history; no biological annual mass/dose difference was observed.
4. [Actual river TP/WWTP source/photograph PR #28](https://github.com/zuizui0223/gesui/pull/28) found no supporting 2022 river-TP mediated landscape pathway (and not a directed point-source-to-butterfly test).

## Boundaries and decision

**Retain the narrow positive:** original reported 2022 WWTP P-source geography is associated in a positive direction with *relative photographic reporting of Lycaenidae* in the supported DE/AT 1° cell frame; the descriptive positive coefficient remains in 2021/2023 static-exposure regressions.

**Do not claim:** WWTP phosphorus benefits Lycaenidae populations, creates adult dietary protein/Na resources, improves larval host quality, causes ant mutualism, or creates a **year-over-year** increase. No directly measured puddling, nutritional chemistry, butterfly abundance, larval survival or within-reach hydrological butterfly exposure is in this panel. The absence of a robust **differenced photo-composition** result also does not prove a null sewage effect.

**Research decision:** Stop using this repeated 1° photo spatial association to argue a phosphorus-mediated adult/larval mechanism. The next meaningful test would require a standardized near-vs-upstream/downstream butterfly sampling design plus actual source/substrate chemistry (or a genuinely external independent sample with measured host/ant covariates); those source-proximity gates failed in the existing 1km/5km studies. Retain the confirmed association as bounded supplementary exploratory evidence, while the full-plants and water-quality questions remain separate.

## Reproduction

Run the [single local source-only replay script](../../scripts/gesui_lycaenidae_2021_2023_posthoc_photo_composition_differencing.py) against the exact already committed original per-cell JSONs from this analysis branch. The replay performs source-key, data-cohort and count-sum checks, the identical OLS and 3,000-block bootstraps, plus leave-one-block-out coefficients and other-family comparators. It writes a machine-readable audit JSON; the accompanying GitHub Actions workflow will preserve it as an artifact. The first source-locked numeric computation was independently evaluated against all 42 matched cell aggregates. This is an **analytical re-interpretation**, not a new biological observation.
