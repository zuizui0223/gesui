# DE/AT 2022 Lycaenidae photographic reporting vs EEA wastewater sources: new macro-scale signal and geographical falsification

**Date:** 2026-11-08. **Status:** reproducible new observational *exploratory* macroecological **photo-reporting** result, preserved on Draft PR #24. **NOT** confirmation that wastewater increases butterfly abundance, adult puddling behavior, actual N/P exposure, larval food quality or ant mutualism. The 2022 EEA reported wastewater **source** has no established spatially direct nutrient delivery to these terrestrial butterflies.

## 1. Research question and exact measurement grain

The previously completed local EEA source–butterfly photo analyses failed predeclared sample-support gates (1 of 24 original 1km pairs and 4 of 17 independent 5km 2021-botanical-effort-selected pairs). Norway standardized 2022 independent BMS Lucas sampling showed **0 near ≤2km and 2 near ≤5km** of qualified EEA outfalls. Those are still FAILED local-source support checks, not evidence of ecological absence.

The separate macroecological question deliberately shifts the unit to **country × WGS84 1° geographic cell** and response to the **fraction of public photographed, georeferenced iNaturalist butterfly observation records from Lycaenidae among three original families Lycaenidae+Pieridae+Papilionidae**, with official descendant taxon IDs. All photographic observations are reporting events, NOT independent butterfly individuals, local density or confirmed substrate feeding. Source predictors are original 2022 EEA WWTP official **reported measured annual N/P discharge masses in t/year and facility numbers**, uniquely linked to a verified FW freshwater public outfall, not pollution concentration at the observation.

The sample is explicitly **conditional on containing ≥1 qualifying reported source** and does not compare an actually wastewater-free treatment group.

## 2. Genuine pre-response source freeze

**First**, prior original source-only EEA public 2022 WWTP data: 630 German (`DE_UWWT_2022`) and 526 Austrian (`AT_UWWTD_2022`) unique reported-positive measured N/P plant sources with exactly one original FW outfall, total 1,156, original source CSV SHA-256 `926ab310904be048d5d36bbfebd5eabfe0b318279800142f2e2a86bee924f4ce`. Those were allocated to **all 53 original source-positive country×1° cells** (DE 32 / AT 21), without reading butterfly outcomes.

**Second**, 2021 iNaturalist publicly georeferenced and photographed **Plantae only** (taxon 47126) at the fixed 53 source cells, with spatial country filter Germany official iNaturalist place ID 7207, Austria 8057; research-grade, uncaptive, open location, 500m accuracy. This is an independent PRIOR-YEAR **plant photo reporting effort** proxy; it does NOT establish independent botanical abundance, exact observer effort or nutrient uptake. Threshold ≥40 Plantae photos and geographic-source diversity fixed before seeing any 2022 butterfly family data.

- **52/53** cells passed (DE 32, AT 20), **22 original country×2° regions**, all 3 source-count bands represented in each country.
- All 2021 Plantae API calls succeeded; Stage B source/effort gate satisfied and **committed to GitHub BEFORE opening 2022 butterfly outcomes**.
- QA warning computed from source-only plant data: WWTP count vs 2021 Plantae photographed observation count rank correlation **ρ +0.210** (32 German cells) and **ρ +0.898** (20 Austrian cells). Both may proxy settlement/land use and highly uneven recording effort, NOT N fertilizer effects.

[Stage B pre-outcome protocol](../../docs/LYCAENIDAE_DE_AT_2022_MACRO_REPORTING_ASSOCIATION_PREOUTCOME_PROTOCOL.md) | [actual 2021 gate success run](https://github.com/zuizui0223/gesui/actions/runs/42462421091) | [exact source-locked 52-cell 2021 cohort](2021_PLANT_ONLY_SOURCE_AND_EFFORT_GATE.json).

## 3. New 2022 butterfly-photo observations and initial within-sample positive association

Official iNaturalist public v1, observations dated **2022-04-01 through 2022-10-31**, country×1° source bbox, photographed/research-grade/noncaptive/open geoprivacy/≤500m original uncertainty, original validated descendant taxon IDs: Lycaenidae 47923, Pieridae 48508, Papilionidae 47223.

- Queries **52 source/2021-plant eligible cells × 3 families = 156** succeeded.
- Predeclared ≥20 three-family photographed records per source cell and minimum 25 cells ≥10 in each country ≥10 2° spatial blocks: **44/52** passed (**Germany 29 / Austria 15**, **22 original independent *spatial blocks*, not 22 independent whole populations**).
- Original photographed records summed across cells: **Lycaenidae 3,010**, **Pieridae 2,074**, **Papilionidae 675**, total **5,759**. Weighted pooled Lycaenidae share is 3,010/5,759 = 52.3%; **equally weighted cell mean** = 48.8%. These are observation/reporting proportions, not abundance.

Fixed exploratory OLS response `log((L+0.5)/(Pieridae+Papilionidae+0.5))` per country×1° cell, each source axis separately entered as within-country standardized `log1p(source)` plus country, latitude, longitude and `log1p(2021 Plantae photo count)`.

| Original WWTP source predictor | Within-sample coefficient per 1 within-country SD of logged source |
| --- | ---: |
| Positive-reported measured N/P unique FW **source count** | **+0.1262** |
| Sum of reported annual **P mass** | **+0.1862** |
| Sum of reported annual **N mass** | **+0.1313** |

For the source-count in-sample beta, **posthoc**, 984 valid out of 999 country×2° cluster bootstrap draws produced a descriptive 95% interval **[+0.0117,+0.2578]**. Original highest-source-count cell removal gave an in-sample beta +0.1463. These are NOT model-independent significance evidence.

Unadjusted within-country rank source-count versus butterfly photo-share: Germany ρ +0.207 (29 cells, nominal p≈.281), Austria ρ +0.213 (15 cells, nominal p≈.447). Thus the simple within-country rank trend is **weak**, and fitted multivariable/country-block results should not be overstated.

[Original 2022 source-qualified API and OLS execution success](https://github.com/zuizui0223/gesui/actions/runs/42462661003) | [source-locked raw response and entire cell-level API receipts](2022_DE_AT_LYCAENIDAE_PHOTO_SHARE_VS_WWTP_SOURCE_EXPLORATORY.json).

## 4. Posthoc regional-block heldout predictive test (not preregistered before seeing coefficients)

After seeing the initial positive coefficients, an **explicitly POSTHOC** follow-up fixed the original 44 cells, the original 22 country×2° geographical blocks, exact no-tune OLS baseline and source additions, train-only transformations, and same source-positive cell/geographic grain. Leave one whole 2° block out at each step. Baseline predictors are country, lat, lon, and *independently observed 2021 Plantae* photo reporting effort.

| Predictive addition | Baseline MSE | Source model MSE | 2° block-heldout MSE improvement | DE | AT |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2022 WWTP count | 0.188898 | 0.184925 | **+2.10%** | +0.88% | +7.60% |
| 2022 reported N mass | 0.188898 | 0.176149 | **+6.75%** | +1.05% | +32.34% |
| 2022 reported P mass | 0.188898 | 0.155243 | **+17.82%** | +10.71% | +49.70% |

These effects are cross-region predictive information **of photo reporting composition**; **not** controlled pollution experiments or validated independent field sites. Selection, country/urban source differences and spatial reporting-pattern confounders persist. No model was tuned to maximize heldout gains, but the predictive validation itself was planned only after initial positive results.

[Completed 22-block heldout test](https://github.com/zuizui0223/gesui/actions/runs/42463169543) | [full posthoc every-cell prediction and controls](2022_POSTHOC_2DEG_BLOCK_OUTOFREGION_WWTP_SOURCE_PREDICTION.json).

## 5. Honest falsification stress test on P and N mass after seeing their positive MSEs

A SECOND separately declared, **POSTHOC** stress test retained the original 44 source-qualified cells and 2022 photo responses, then removed only the maximum original source **P or N** mass cell in each country, selected **solely** from source metadata; and performed fixed random source-mass permutation controls.

| Stress, all on same underlying original 2022 macro source cohort | Reported P mass | Reported N mass |
| --- | ---: | ---: |
| Original 2° heldout MSE improvement | **+17.82%** | +6.75% |
| After removing each country's top **source-mass** 1° cell (42 cells) | **+13.49%** | +1.54% |
| Across-country-preserving source randomizations: fraction ≥ observed gain (399; with +1 adjustment) | **0.0125** | 0.0575 |
| Within-country **and 2° spatial region** restricted source randomization fraction ≥ observed (199; +1 adjustment) | **0.165** | 0.465 |

For P excluding top source-mass source-only cells, same 42-cell heldout gain DE **+8.80%**, AT **+39.13%**; for N, DE **−0.88%** and AT **+14.81%**. Thus P has more spatially robust out-of-block geographical predictive information than N **under this observational sample and fixed model**, but **P is not demonstrably the distinct causal ecological agent**. Original reported N/P masses are strongly collinear, Spearman **ρ +0.8877** across the 44 photo-eligible cells. The region-*within*-2° randomization has P tail **0.165**: the signal is **not strongly resolved to individual cells beyond 2° regional source pattern**. Country source-label permutation is a spatial negative control, *not a randomized pollution-treatment p-value*.

[Source-mass falsification & outlier test success](https://github.com/zuizui0223/gesui/actions/runs/42463399911) | [original P/N randomizations and source-mass deletion receipt](2022_POSTHOC_P_N_MASS_SHUFFLED_SOURCE_AND_MAXIMUM_REMOVAL.json).

## 6. Scientific claim boundary and interpretation

**Positive supported at the data grain:** Within this source-containing DE/AT 2022 country×1° photographic reporting cohort, registered **public WWTP source load, particularly annual reported P**, holds additional geographically transferable information about the **recorded proportion** of Lycaenidae relative to Pieridae and Papilionidae beyond a baseline of country, coordinates and 2021 Plantae photo reporting proxy. A country-permuted source-mass negative control yields a more favorable observed P prediction than most permuted candidates, though **within-2°** shuffled source can produce comparably high gains and direct landscape confounding remains unresolved.

**Not supported:** that wastewater caused Lycaenidae numbers to grow, that P mass versus N mass is a biochemically distinct mechanism, that reported mass reaches butterfly drinking sites, that photographed adults were mud-puddling, or that introduced host resources directly mediate the effect. This **cannot replace** the prior smaller-scale 1km/5km failed source-aligned photographic support audits; different units, detection and ecological source models are being tested.

**Unresolved alternative explanations:** urban/built fraction, crop cover, EEA reporting coverage and WWTP plant size, uneven one-degree cell land area, observer density and photo species identification, butterfly host plants and fragmentation, climate/latitude/phenology, adult sodium/protein availability and within-reach source transport. Prior-year Plantae photo effort only partially proxies these; in Austria it itself correlates very strongly with public WWTP count.

**Next prospective independent check** would need genuinely independently measured urban land cover, butterfly detection/observer effort or standardized butterfly transects, and independent, year-aligned nutrient chemistry/river/shore environmental data. A prior-year source 2021/2023 photo replication would test stability of *recording geography*, not demonstrate same-year wastewater dose. Do not conflate any separate `chocho` biological host geography result.

## Reproducibility

- [Original public EEA DE/AT 2022 data and 1° source grid source-only PR #23](https://github.com/zuizui0223/gesui/pull/23)
- [2022 source-positive macro exploration Draft PR #24](https://github.com/zuizui0223/gesui/pull/24)
- All scripts and original API hashes frozen in repo, no screenshot/photo/restricted data redistributed, no main merge.
- The original `chocho` result **22/25 Lycaenidae expanded potential larval host resource regions with human introduced hosts** remains a separate, unchanged positive.
