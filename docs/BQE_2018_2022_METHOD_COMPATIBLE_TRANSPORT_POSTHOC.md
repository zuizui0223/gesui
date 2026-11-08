# 2018/2022 biological EQR: country × assessment-method compatible transport diagnostic

**Status:** additional **POSTHOC** stress test, designed after the original positive 2018 and 2022 results and after the existing site-disjoint cross-year prediction. This document fixes the diagnostic definitions before computing its method-compatible subset results; it is **not** preregistered independent confirmation.

## Falsifiable question

Could the small positive held-out WWTP reported-load information for benthic macroinvertebrate EQR be an artifact of transporting predictions into country × biological assessment-method combinations not represented in the source-year training data?

## Fixed cohort and inputs

- Same frozen 2018 and 2022 original biological station CSVs, actual **year-specific reported** measured WWTP P/N mass, source-only EQR sites and HydroRIVERS links as Draft PR #17.
- Exact source CSV SHA-256 checks are required and original test site IDs are **excluded** from the training year.
- Fixed Ridge(alpha=10), geography (lat, lon), log river discharge and basin area, country and country × assessment-method categorical baselines, five outcomes as in original test; **no tuning**, feature changes, new data or source-selection.
- Evaluate both directions: fit 2018 / test 2022-exclusive EQR sites and fit 2022 / test 2018-exclusive sites.
- Benthic invertebrates are the target. The same protocol is applied to phytobenthos as a **specificity control**, reporting failures equally.

## Method-compatibility strata, declared before results

For each test site, define its `countryCode|procedureClassificationSystem` exact key. **Known-method** means that exact key exists among the other year's training sites. **Unseen-method** means it does not. No test outcome or load magnitude enters this assignment.

Within both directions, group and nutrient:
1. Report unstratified M0 vs M1 upstream annual P or N vs M2 physically downstream placebo MSE, reproduce committed source-locked original MSE to numerical tolerance.
2. Report n, source-positive n and MSE-based gains for **known-method** test sites and separately **unseen-method** test sites, irrespective of whether the result is positive.
3. Report the same within the **known-method × registered-source-positive** and **known-method × source-zero** groups to expose whether pooled gains reflect only one exposure stratum.
4. Provide per-country gain directions among compatible test sites; report both pooled (station weighted) and country-balanced MSE gain.
5. Use a **descriptive** 1,999-replicate country cluster bootstrap interval (fixed seed), resampling entire country blocks for the known-method subset. This CI is a sensitivity measure, **not** a confirmatory P-value (methodological choice is posthoc).

## Failure / success interpretation

- **Robustness retained:** known-method subset has M1<M0 and M1<M2 for P and N, in both temporal directions. This is a *posthoc* robustness signal, not sewage causation.
- **Compatibility-limited:** any of these conditions fails; maintain the original year-specific positive results, but narrow their cross-year transfer interpretation.
- **Cluster-limited:** a country-cluster interval spans zero or country-level improvements are uncommon. Do not claim reliable across-country generalization.

There are only two source years; report-year WWTP inventories have incomplete coverage, source-zero is **not** proven zero wastewater, and biological EQR method/country ratings may differ in calibration. A river WQ/chemical exposure signal, an intervention effect, N-vs-P separation, or changes in biological abundance/richness cannot be identified here.

**Data provenance:** [2018 source and frozen biological EQR](https://github.com/zuizui0223/gesui/pull/17); [2022 EQR original](https://github.com/zuizui0223/gesui/pull/14); [missingness and agriculture controls](https://github.com/zuizui0223/gesui/pull/16).
