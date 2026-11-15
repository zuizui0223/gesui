# Source-locked after 40-H6 opportunity: same-H6 plant–butterfly TERRESTRIAL-ONLY within-parent comparison

**2026-11-15. Exploratory posthoc stress test following prior broad positive/temporal nulls. This is not a prospectively preregistered experiment. The source-only eligibility and model are frozen before any new cross-panel fitted coefficient.**

## Fixed cohort and source bounds

- Authentic old period-paired H7 EIV/HNS source-supported cohort: 219 H7 across 86 genuine H6 from exact original 451 H7×period EIV-supported panel, source SHA256 `1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca`, original source facets SHA256 `251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4`.
- The old H6 parents with ≥2 distinct old eligible H7 are **30 H6 with 163 H7**. The 56 singleton H6 remain identified before outcome.
- Genuine source-only extra 29 H7 across 10 previously-singleton H6 (not old 232 H7), all with both years original GBIF/trait support; select **one per singleton H6** only by original prior EEA source geographic rank: 10 H7/10 H6. Verify membership and immutable source rank from covariate replay #43725920793. There are no duplicate H7 IDs or H6 parent reassignments. The 10 old originally-singleton H7 in those same H6 remain in analysis.
- Target after removing the other old singleton H6 is **183 distinct period-paired H7 in 40 genuine H6**, with an equal two-period mean of each H7. Do not count the 2 periods as independent geographic replications. A separate source-only gating script verifies this exact geography and source score support before any ecological model; fail rather than alter gates if not achieved.

## Comparable source-only environmental covariates

- Prior static modelled wastewater geography `W=log1p(max(0,upstream outlet modeled WWTP m3/day)/(UP_AREA_KM2+1))`; `area=log(local H7 area_km2)`; `crop=crp_pc_use`; `urban=urb_pc_use`; `naturalout=log1p(HydroATLAS dis_m3_pyr)`.
- Original full 232 H7/454 period rows reproduced exact to ≤1e−10 by #43725920793 using authentic 1219 W/HydroATLAS source. These 5 fields exist identically for all 29 new quality-supported H7. They are temporally STATIC and *not actual local receiving-river nitrogen*.
- New H7 often have EEA monitoring in just one original period, so former `nEEA`, original actual monitor `riverflow`, `country×original assessment method`, and original AQUATIC EQR are **not available in a comparable form**. **Do not include them** in new terrestrial-only models and **do not label** this an exact repeat of prior EEA/EQR–plant–butterfly model.
- For original and new H7 alike, reconstruct `nPlant`, `nButterfly` from ORIGINAL FIXED 6-facet GBIF raw record counts in each period, and `nAnimal` from same previous original 1219×2 GBIF index master; average `log1p` effort values between the two periods per H7. Frozen source SHA identities and missingness carried through (never zero-impute failed API).
- Original plant NUTRI = species-equal mean over exact accepted EIV numeric species for each original period, derived from same source-v2 unchanged exact 4,011 scored taxa. New H7 = same 4,011-key dictionary, no additional taxonomy imputation; only original source-supported H7. Two-period H7 means use unweighted 0.5 each period.
- Original/new butterfly HNS = species-equal mean over original fixed numeric 186 accepted speciesKey, each observed species counted once per period, averaged across two periods. Record-weighted HNS only exploratory secondary. Source-specific score coverage and scored species count are controls.

## Inference question: within the SAME parent genuine H6

**One primary response:** H7 two-period mean butterfly incidence-equal HNS. **One focal predictor:** same H7 two-period mean accepted species-equal plant NUTRI. Correct for *H6 fixed intercepts* (not a source country/method dummy) plus source W, crop, urban, naturalout, area, mean `nPlant`, mean `nButterfly`, mean `nAnimal`, mean plant EIV coverage, mean HNS score coverage, mean scored plant species count and HNS-scored species count. Standardize numeric covariates over the frozen 183-H7 cohort, do not optimize based on sign.

- Group/H6-equal WLS: each of 40 H6 contributes equal total weight (`w_i ∝ 1/number H7 in that H6`). Fixed H6 intercept; use cluster-robust SE at 40 H6 and descriptive normal95% CI.
- Source-only gate must verify full design rank, at least 183 H7/40 H6, and numeric nonmissing values from genuine source. Significant result cannot be promoted to causal effect.
- **Primary descriptive falsification:** Freedman–Lane residual-label permutation within the same H6 (9,999 draws, deterministic numpy PCG64 seed20261115), original no-PlantN weighted null fit and source explanatory variables fixed; residuals permuted only across distinct H7 within the same parent H6. Evaluate two-sided absolute W coefficient. Approximate conditional exchangeability null, *not* an experiment or phylogenetically corrected trait permutation.
- Sensitivity: original butterfly record-weighted HNS in same model and same cohort; leave-one-H6-parent-out 40 models' coefficient range (do not cherry-pick), compare sign. No multiplicity-adjusted claim for posthoc secondary.
- **Never call H6-fixed-effect fit “out-of-H6 predictive transfer”:** parent H6 dummies prevent meaningful predictions for entirely new H6. Previous H6-out result and new within-H6 inference are different diagnostics; neither implies flow of N from sewage into terrestrial host plants.
- Do not retrospectively replace which extra H7 was selected, exclude a parent H6 after sign inspection, change thresholds/controls to make a significant p, or merge prior main. Report full negative result if within-parent association vanishes.

## Limitation

GBIF source records are not detection-corrected butterfly occupancy, field plant cover, larval consumption or butterfly biomass. EIV NUTRI/HNS share underlying plant habitat-affinity dictionaries. W and most environmental covariates are static basin-wide context strongly correlated with urbanization, not measured sewage-N fraction delivered to riverbank or insects. This standalone ecological community-composition comparison cannot estimate causal mediation, trophic energy transfer or ecosystem service production. Previous within-H7 temporal delta null remains unchanged; prior original aquatic EQR ↔ HNS residual is uncertain.
