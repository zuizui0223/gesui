# 2023–24 NRSA measured nutrients, areal AFDM and chlorophyll-a: pigment per organic mass

**2026-11-14 | Draft / observational / posthoc, genuine EPA source and independent full numerical QA SUCCESS.** Do not merge, overwrite original positive/negative cohorts, or claim WWTP-specific nutrient causal effects.

## Original source gate and data availability

The authentic original EPA 2023–24 **Benthic Count**, **Benthic Metrics**, **Benthic MMI**, **Benthic Taxa List**, and their listed publisher metadata were requested using the URLs extracted from the official EPA NARS download table. In [original source-only #43453867202](https://github.com/zuizui0223/gesui/actions/runs/43453867202), the named biological data endpoints all returned **HTTP 404**. No 2023/24 benthic MMI, species count, individualCount or taxonomic EPT result was invented. The official EPA source download table still advertises the files, so this is a **publisher access problem**, not absent benthic organisms or a null biodiversity response.

The already established real, independent original [PR #73](https://github.com/zuizui0223/gesui/pull/73) 2023–24 source panel remains usable: 2,129 native UID common EPA site / real water chemistry / AFDM / periphyton-Chla measurements, strict matched **1,828** with positive concentrations, positive areal measurements and no disallowed EPA NARS flags, from **946 HUC8 watersheds /48 states**. Exact publisher original data SHA, source-UID crosswalk and strict-quality panel were revalidated, not changed after outcome signs.

## New observational hypothesis: mass vs pigment

Define `R = log10(Chla [ug/cm2] / AFDM [mg/cm2])` = `log10(Chla [ug/mg AFDM])`. This is **pigment content per standing dry organic biofilm material**, not true autotrophic biomass fraction, primary-production rate, biofilm energy flux, or pollinator ecosystem service. This ratio is a deterministic **difference of prior outcomes** (`log Chla−log AFDM`), hence the apparent new P association must NOT be presented as independent field replication or a novel independent ecological effect. It tests a distinct descriptive *contrast/estimand* on the exact same original samples.

Survey covariates: native state and eco9, latitude, longitude, Strahler stream order, actual survey year and field sampling season (sin/cos), field conductivity, and source-positive observed instream log10 total-N/total-P. No actual wastewater origin assigned. Whole HUC8 5fold heldout Ridge alpha10, train-only standardization. HUC8-cluster description 95 intervals, not causally identified uncertainty.

### Observed exact 1,828 EPA visits, 946 original HUC8

| Real response (separate same-UID models) | Adding measured river TN+TP to shared environmental baseline, heldout MSE gain | Adjusted N per +1 SD source log10TN | Adjusted P per +1 SD source log10TP |
|---|---:|---:|---:|
| **log10(Chla ug / AFDM mg)** | **+2.2359%** | +0.01725 [−0.02063,+0.05513] | **−0.09568 [−0.13406,−0.05731]** |
| log10(areal AFDM mg/cm2) | **−0.08591%** | +0.01434 [−0.02035,+0.04902] | +0.02110 [−0.00875,+0.05094] |
| log10(areal Chla ug/cm2) | **+0.81208%** | +0.03159 [−0.01211,+0.07529] | **−0.07459 [−0.11532,−0.03386]** |

- Original same-UID `logAFDM` vs `logChla` Spearman **ρ=+0.35104**.
- **Chla model additionally adjusted for AFDM**: source +1SD logTP **β−0.08914** (HUC8 cluster95 **[−0.12569,−0.05259]**), while +1SD AFDM **β+0.26408 [0.22430,0.30387]**. This *conditional observational pigment difference* is not mediation.
- A −0.09568 log10 pigment-per-mass coefficient corresponds to approx 0.802× ratio per +1SD logTP under the fitted model, but is NOT an effect of adding P, WWTP changes or nutrient manipulation. Observation-only.
- Nutrient-added relative prediction improves ratio modestly, but AFDM alone has no gain. This argues **against a simple universal nutrient→greater biofilm organic mass story in this particular geography**, not absence of eutrophication or all biomass reactions.
- The absolute source nutrient/quality coefficients retain substantial possible hydrology, shading, substrate, temperature, light, land use and detection confounding. Spatial heldout cannot fix causal source omission or method differences.

## Independent numerical validation

- [Original measured model from EPA 8 source CSVs #43454306097 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43454306097): strict original 1828/946, full AFDM/Chla/ratio model, cluster adjusted slopes, HUC8 fivefold MSE, original source hashes.
- [Independent manual NumPy algebra and Ridge/OLS #43454792597 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43454792597): reconstitutes original exact 2023 strict cohort; different matrix implementation (from initial sklearn/statsmodels) reproduces ratio baseline MSE **0.1974541102**, +TN+TP MSE **0.1930391702**, TP adjusted ratio slope **−0.0956839907**, all old component MSEs and cluster CIs within specified tolerance. Verifies **ratio slope = Chla slope−AFDM slope by algebra**. This is computational verification on the **same** original field data, not independent ecological replication.

## Consequence for biodiversity × ecosystem services

**The desired 2023–24 same-visit original stream TN/TP + AFDM + Chla + true benthic MMI/EPT is still unavailable**, because official publisher Count/Metric/MMI source URLs return 404. The separate 2018–19 actual NRSA EPT composition response (n947) and 2023–24 AFDM/Chla (n1828) remain distinct cohorts and cannot be synthetically chained as a measured pathway. No insect emergence biomass, water nutrient purification process, predator diet, host resource or actual pollination service was measured on these same visits.

The prior [PR #74](https://github.com/zuizui0223/gesui/pull/74) exact 2023/24 year-specific registered upstream WWTP→EEA DMR dual N/P gate remains **failed** (2023 N/P quantity n13 or concentration n25; 2024 n17/n38). Nothing here revises that evidence or demonstrates actual WWTP source effects.

**Next meaningful data opportunity:** recover published EPA NRSA benthic MMI/Counts from a separately authenticated official NARS API or restored original download; alternatively collect genuine water chemistry/biomass/benthic response panels from independent monitoring sources on the same source/season. Do not invent absent official data, tune retrospectively for a positive coefficient, or call simulated functions real services.
