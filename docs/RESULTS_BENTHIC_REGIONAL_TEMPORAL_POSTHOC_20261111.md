# Posthoc benthic macroecology: 2016 within-country regional signal and 2018→2022 shared-site exposure limits

**Date:** 2026-11-11 (JST)  
**Status:** Read-only, strictly **POSTHOC** robustness / source interpretation; all original 2016, 2018, and 2022 primary outputs are preserved. This is **NOT** an independent prospective test, causal estimate, or new data collection.

## What was newly checked

Using **three original files already included in the parent branch** and pinning their Git blob SHAs, independently inspected:

1. The original **2016** JRC 2014–18 **modelled domestic+industrial point-source N/P** per NUTS2 region vs real EEA benthic EQR *regional mean*, recomputed Pearson correlation with each country's region-level mean subtracted.
2. Leave-one-country-out regional association and country-level sign heterogeneity. An additional exploratory 10,000-draw country-block bootstrap, **not** a preregistered or causal confidence statement.
3. Original **2018 and 2022** annual *measured* WWTP N/P upstream-reported-load archives joined by **exact monitoring site_key** to describe EQR changes across source-reporting presence categories. Same biological assessment method restriction evaluated without redefining the reporting source as a real intervention.

## 2016 (actual biological EQR, modelled N/P sources)

- Source: 2,435 actual benthic numeric EQR sites / 12 countries, assigned to **73 NUTS2 regions**. This new audit uses **73 regional means as 73 regions**, not 2,435 independently exposed treatments.
- Equal-region 2016 EQR versus 2014–18 modelled log point N/P density: **raw Pearson** N **−0.4110**, P **−0.4440**.
- **Within-country centered regional Pearson:** N **−0.341053**, P **−0.356308**.
- **Leave one of 12 countries out:** N range **[−0.39345, −0.24053]**, P **[−0.42009, −0.25127]**; negative in **12/12** runs for both.
- **Country-block bootstrap 10,000 resamples**, descriptive percentile 95% N **[−0.54235, −0.11565]**, P **[−0.55443, −0.08267]**; posthoc dependence/exchangeability caveat applies.
- Restrict to the **8 countries with at least 3 NUTS2 regions**: N and P each **7/8** within-country correlations negative; **Sweden** is positive for both. Mean of the 8 country-specific r: N **−0.28449**, P **−0.25153**. This is not a claim that all countries agree.
- N/P point-source densities are collinear; the parent report's adjusted 2016 EQR site model estimates N-only and P-only negative slopes but **joint N+P coefficients cross zero**. Neither the regional correlation nor that model distinguishes nutrient-specific causal toxicity from other regional pressures.

**Interpretation:** the negative regional EQR association is not merely a between-country artifact or a single-country outlier. It is **still observational, coarse, method-and-habitat-confounded** and does **not** establish observed 2016 WWTP plant discharge. The published parent 2016 out-of-region forecast gains N +2.287%, P +1.402% are separate performance metrics, whose region-block bootstrap intervals include zero.

## Exact shared river monitoring sites, 2018→2022

- Source-frozen archives: **2018: 2,634 sites**; **2022: 2,475 sites**; intersection **1,167 exactly shared site keys**, not independent biological cohorts.
- Categories below use **qualifying *reported* upstream N/P loads in that year**, not an intervention, construction/closure, or independently quantified wastewater exposure in the river.

| 2018→2022 qualifying source reporting | Same sites | Mean Δ normalized EQR (2022−2018) |
|---|---:|---:|
| no→no | **1,139** | **−0.033985** |
| yes→yes | **17** | **+0.031569** |
| no→yes | **7** | **−0.034430** |
| yes→no | **4** | **+0.010495** |

- The 2018/2022 EQR classification method code is unchanged at **1,148/1,167 sites**. Restricting to these gives no→no **1,123** ΔEQR **−0.033454**, yes→yes **14** ΔEQR **+0.052648**, no→yes **7** ΔEQR **−0.034430**, and yes→no **4** ΔEQR **+0.010495**.
- Measured N-positive and P-positive exposure status flags matched for these sites within each year; this does **not** establish N vs P separability.
- Only **7 apparent report-on and 4 report-off sites** remain; their source statuses do not assert actual new/closed WWTPs. Do **not** infer a before–after treatment effect or a difference-in-differences estimate from these small, nonrandom groups.

**Interpretation:** two original annual WWTP source–EQR *predictive gains* (2018 P +0.982% N +1.190%; 2022 P +1.10% N +1.40%) are directionally reproduced in their separately designed out-of-region predictions, but cannot be converted into a **within-site temporal biological response** in this reused archive.

## Strongest defensible manuscript statement

> Measured upstream WWTP N/P reporting adds small, reproducible predictive information about European river benthic-invertebrate EQR in 2018 and 2022; a modelled broader point-source N/P regional gradient is also negatively associated with actual benthic condition in 2016 within countries. Neither result identifies the causal effect of wastewater, a specific nutrient, or a wastewater-driven within-site temporal change.

## Reproducibility and boundaries

- Read-only script: `node scripts/audit_benthic_three_year_posthoc_20261111.mjs` (repository root, Node >=18; original source Git blob SHAs enforced). The script is intentionally deterministic and does not query outside APIs.
- **2016 region JSON** Git blob `108ac711022921af8cfb2b5477d046f386d43f17`, original base [PR #56](https://github.com/zuizui0223/gesui/pull/56).
- **2018 real EQR + WWTP loads CSV** Git blob `0f30bfbc3dac13f6a182c8a910e9899c04778322`, original [PR #17](https://github.com/zuizui0223/gesui/pull/17).
- **2022 real EQR + WWTP loads CSV** Git blob `325f199e985d9862ffee5b0edb205c61a7be3833`, original [PR #14](https://github.com/zuizui0223/gesui/pull/14).
- The newly reported original-file descriptive statistics were computed from the archived aggregate JSON and actual site CSVs using an independent JS implementation; this is *not* a fresh source download or independent external lab replication.
- No changes to `gesui/main`, `chocho`, source-year PR #14/#17/#56, or any published manuscript. Keep **Draft**.

## Next real independent test (not yet executed)

Source-first select **an actually observed, contemporaneous river biological + measured municipal WWTP source discharge dataset in a previously unused country/year** with adequate upstream-exposed and unexposed rivers. Predefine region-disjoint holdout, actual routed downstream exposure, agricultural N/P and runoff/urban covariates, verified flow dilution, country×assessment-method, and an *upstream vs downstream source* negative control. Do not reuse the JRC 2014–18 model as a third measured year; do not rescue failing source gates by reweighting/choosing favorable locations.
