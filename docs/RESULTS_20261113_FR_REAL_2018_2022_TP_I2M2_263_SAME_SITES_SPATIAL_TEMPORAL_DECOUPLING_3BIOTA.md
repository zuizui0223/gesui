# Strong spatial but weak four-year temporal coupling of measured river phosphorus and benthic quality: genuine repeated French river stations

**13 November 2026 | Original official measured water chemistry, EEA field benthic I2M2, same physical monitoring stations | 263 original stations /73 HydroBASINS L7 catchments/24 L6 parent regions | POSTHOC exploratory source-locked test | Draft/unmerged.**

## Main scientific distinction

This is not an additional "registered WWTP causes ecosystem harm" result. Original 2018 and 2022 **measured ambient river TP** co-varies strongly and negatively with EEA benthic biological I2M2 **between** field monitoring stations and watersheds in either year, but **the change in TP from 2018 to 2022 does not predict the change in I2M2 EQR in those exact same stations/watersheds**. These are different estimands, both valid, and must not be conflated. A causal nutrient or service cascade is not identified.

The comparison builds on the existing PR #70 2022 same-station TP/I2M2 result (n273 source-matched sites, n74 French H7; raw 2022 basin rho−0.698 and +22.78% true H6 parent region-blocked predictive improvement from ambient TP). The new 4-year data were extracted from original authors' raw datewise TP archive, NOT inferred from 2022 metadata or synthetic land use.

## Source-first quality gates, real original values and source identity

**Same selected source cohort:**
- 2022 original true TP WQ site and EEA I2M2 source-only matches at native HydroRIVERS `HYRIV_ID`, distance ≤250m, 273 field-site pairs /74 FR H7, original source-only cohort SHA256 `448add2d6850cbc876c3e6b90bca278c6b711e64f9d363071e93bc0c18b84611`, source run [#43370128446](https://github.com/zuizui0223/gesui/actions/runs/43370128446).
- Original **2018 EEA I2M2** monitoring station ID and native river reach match **263/273** of those 2022 source pairs exactly. Not a nearby site substitute; **263 native EEA station IDs identical across both years**, 73 original H7, 24 original H6 parent hydrologic regions. Original EEA2018 Git blob `0f30bfbc3dac13f6a182c8a910e9899c04778322`, EEA2022 `325f199e985d9862ffee5b0edb205c61a7be3833`. [Original source report-quality+2018 biological station gate #43372659602](https://github.com/zuizui0223/gesui/actions/runs/43372659602).
- Original published Peters et al 2018–2022 TP.2 harmonized datewise ZIP file `70896288`, DOI [10.25400/lincolnuninz.30546371.v5](https://doi.org/10.25400/lincolnuninz.30546371.v5), source length **25,003,722 bytes**, publisher verified MD5 **`b596dbe9833c49426147a99f10f3fe31`**, original archive SHA256 **`26f2af11f4972987ee6f7c3f497ad902820f145851f71919d4f45b589e1d809e`**. Verified actual original date/value **positive 2018 AND positive 2022 TP measured ≥4 times in ≥3 distinct months for ALL 263 source-frozen WQ IDs**, without opening TP magnitudes or EQR scores to choose the cohort. Independently frozen 2018–22 sample eligibility SHA256 `2190ac905326f4f2f4f278eb0debb4a7bfda7c4c71347d4941b54cc19d3de264`, source [#43373061409 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43373061409).
- Original **74-H7/24-H6** spatial GroupKFold source IDs already chosen BEFORE 2022 and new 2018 biological outcome coefficients, SHA256 `294eb69b76cd4cf40f1022a5e83278529b1d2bbceac9a38931004610a23eb052`. The common **73 H7** source 5fold counts are 15/15/15/15/13. Original 24 H6 parents require prior documented source-only repair of invalid original HydroBASINS H6 geometries; do not call those original GIS files intrinsically error-free.
- Real original 2022 WQ TP medians reconstructed **exactly at all 263 original WQ IDs** from the official TP.2 datewise raw archive, compared to and matching the earlier derived 2022 source medians. No missing or incompatible original 2018 EQR/2022 EQR source magnitudes needed to be imputed or selected by sign.

**Critical WWTP effluent source gate FAILED**: Source-matched EEA 2022 annual *reported* upstream WWTP P (not measured ambient TP) positive at only **4 of 273 2022 same-site TP/EQR field stations in 2 H7 and 2 H6 parents**. Its source positivity was independently reproduced in the WQ station derivation, but both derive from the same restricted original EEA data and are not independent effluent field measurements. Source eligibility threshold predeclared 30 stations, 15 H7, 5 H6 is not met. Most zero values reflect **no positive source found in the restricted reported inventory**, NOT no WWTP/absence of pollution; therefore **NO causal wastewater→TP→benthic EQR mediation model was fit**.

## Primary new empirical macro results: *within the same watersheds over 2018→2022*

For each of the exact original 263 TP+I2M2 station pairs, original 2018 measured TP positive median and original 2022 measured TP positive median are joined to 2018 EEA `numeric_EQR_2018` and 2022 EEA `normalised_EQR` on the identical field station ID. Source method both **I2M2** and country both **France**, source river reach same. Then equal-average physical station **Δln(TP) and ΔEQR** once within each original 73 H7.

| Exact original time-change evidence | Observed original result |
|---|---:|
| **Δlog actual measured TP vs ΔI2M2 EQR**, n73 independent H7 basin averages, raw Spearman | **ρ = −0.12597**, descriptive p0.288 |
| **Δlog actual measured TP vs ΔI2M2 EQR**, 263 nested physical station pairs, raw Spearman | **ρ = +0.00320**, do not treat n263 as independent H7 ecological replicates |
| 42 H7 with ≥2 same-site pairs, 232 station pairs, pooled **within-H7 site-demeaned temporal** TP–EQR slope | **+0.03225 EQR per Δln(TP)**; H7-resample descriptive 95 **[−0.03859,+0.09607]** |
| H6-parent clustered source-adjusted slope per 1 SD **Δlog observed TP** (SD 0.26315) after original static upstream cropland/urban/mean river flow/H7 area and relative 2018/22 TP sampling count | **−0.01843 EQR**, descriptive H6 95 **[−0.05545,+0.01859]**, crosses 0 |
| True original 24 H6-parent 5fold geographic CV, same 73 H7 B0 prediction baseline MSE | **0.0143811** |
| True H6 CV B0 + actual **Δlog measured TP** EQR-change MSE | **0.0144399**; **−0.409%** predictive MSE gain (worse), parent 2000 bootstrap descriptive95 **[−6.290%,+3.782%]** |
| 24fold original H6-parent leave-one-parent-out B0 vs +ΔTP | **−1.181%** predictive gain, parent bootstrap descriptive95 **[−8.439%,+2.466%]** |
| Secondary source-adjusted B0 also containing original 2018 TP level + ΔTP | **−0.259%** EQR-change gain, bootstrap95 [−5.863%,+3.708%] |
| Mean original 73-H7 ΔEQR (2022−2018) | **−0.03338** |
| Mean original 73-H7 Δln(TP) (2022−2018) | **+0.03159** |

A null ∆ association does not prove nutrient stress has no ecological effect: yearwise medians, sampling seasonality, ecological response lags, site changes, original EEA EQR normalization and background changes were not experimentally controlled. However the original **positive spatial** TP–EQR link does **not** license claiming real 4-year within-river degradation proportional to TP changes.

[Primary source-frozen 263 repeated physical TP/I2M2 field station response 2018–2022, geographic CV and all indexed terrestrial endpoints SUCCESS #43373639391](https://github.com/zuizui0223/gesui/actions/runs/43373639391). Inspectable original 263 station-level real TP+I2M2 pairs and original 73 H7-level source/data tables plus JSON are preserved in that Actions artifact, not the entire licensed Peters datewise archive.

## Why space and time differ — identical original 263 physical stations, BOTH original years independently

After observing the temporal null, a **separate explicitly POSTHOC source-only-results numerical audit** independently reimplemented the prior original sklearn Ridge model using **NumPy manual normal equations**, exactly reproducing the original FOUR n73 H6 fivefold and leave-H6-out MSEs, original 73-basin Δ Spearman, 263-station Δ Spearman, and within-H7 demeaned ΔTP–ΔEQR slope. It then compared 2018 and 2022 physical river quality **individually at the SAME 263 field stations and 73 catchment means**, not adding new sites.

| Genuine matched-year comparison in 73 original French H7 | Raw H7 TP–I2M2 EQR Spearman | True 24-H6-parent 5fold TP-added EQR MSE gain | Posthoc 24fold H6-parent leave-one-group-out |
|---|---:|---:|---:|
| **2018 real river TP vs original 2018 same-site I2M2** | **−0.73285** | **+30.812%**, H6 bootstrap95 [+6.40,+43.54]% | **+30.668%**, H6 bootstrap95 [+6.53,+42.68]% |
| **2022 real river TP vs original 2022 same-site I2M2** | **−0.69719** | **+22.667%**, H6 bootstrap95 [+9.17,+37.65]% | **+22.556%**, H6 bootstrap95 [+8.38,+38.17]% |
| **Changes across 2018→2022 in exactly same physical river field stations** | **−0.12597** | **−0.409%**, H6 bootstrap95 [−6.29,+3.78]% | **−1.181%**, H6 bootstrap95 [−8.44,+2.47]% |

Original repeated-H7 source rank stability across years: **2018 vs 2022 logTP Spearman ρ +0.95409**, **2018 vs 2022 I2M2 EQR ρ +0.84864**. Original between-site two-year mean logTP vs mean EQR ρ **−0.74599**. This is **stable geographic stratification**, not established temporal dose-response or biological adaptation. Same sample, multiple checks POSTHOC — the 2018 and 2022 positive spatial gradients are not independent physical experiments.

[Independent NumPy exact numeric reproduction and posthoc year-separated spatial-vs-temporal comparison SUCCESS #43374058129](https://github.com/zuizui0223/gesui/actions/runs/43374058129).

## Same original H7×2 periods: GBIF plant and Lycaenidae responses do NOT establish services

The original full 454 H7×period macro GBIF/EEA study source contains BOTH 2017–19 and 2020–22 original flora speciesKey-rarefied richness S50 and original Lycaenidae five-family relative indexed occurrence ratio for all 73 H7 in the same repeated field-site 2018/22 water TP/I2M2 panel. One H7 = one ecological inference row, not 263 independent butterfly/pollinator communities.

| Real measured Δln river TP vs H7 2017–19→2020–22 indexed terrestrial change | Raw H7 rank association | +ΔTP H6 heldout 5fold MSE gain (H6 bootstrap95) | 24-parent leaveout MSE gain |
|---|---:|---:|---:|
| **Original GBIF Plantae speciesKey 50-record rarefied S50 Δ** | ρ **−0.0890** | **−2.135%** [−6.223%,+0.229%] | −2.395% |
| **Original GBIF Lycaenidae log odds among original five butterfly-family index occurrences Δ** | ρ **+0.1030** | **+0.451%** [−3.658%,+5.258%] | +0.569% |

The GBIF outcomes were indexed photo/occurrence-source diversity/record-composition proxies, **not native riparian field flora, plant hosts, flowering, adult butterfly abundance, pollinator visitation or reproductive ecosystem services**. Observing no reliable 4-year correlation with true river TP variation cannot disprove biological aquatic-to-terrestrial flows. The original +15–17% H7 **modeled WWTP-source → indexed plant S50 spatial** finding from PR #70 was in a selected n74 French subset and does **NOT** replicate as actual measured river TP or as time change. Preserve it as source/observation-geography hypothesis, not pollution-driven biodiversity gain.

## Publication-safe synthesis and scope

**Suggested narrow manuscript question**: *When do spatial nutrient–benthic quality gradients predict temporal ecological responses?* In this one country/source BQE system, two original years each show a large, negative same-site/same-watershed phosphorus–benthic integrity association with H6 geographic transfer. But the observed within-station H7 4-year changes have near-null correspondence, despite genuine repeated physical site measurements and independent numerical replay. Temporal response lag, a stronger effect of stable site context, unequal within/between variance, biological measurement noise and hydrologic seasonality are plausible and **not distinguished** here.

Evidence by ecological service chain:
1. `WWTP measured annual effluent-P → actual same-site ambient river TP`: **not identifiable** in this matched cohort (only 4 source-positive sites, reporting zeros are not untreated controls).
2. `ambient TP ↔ benthic quality`: **strong and independently numerical-reproduced cross-sectional negative association in each year**, and **weak/nondiagnostic 2018→2022 within-site change**.
3. `ambient TP/EQR → plant speciesKey S50 and butterfly five-family share`: **weak/null original indexed change associations** in original same H7 windows, not proof of no effects on terrestrial ecology.
4. `emerging aquatic insects → predators → pollination/seed set`: not jointly measured at these same original French sites/year. Do not estimate a serial SEM or claim observed services.

**Scientific priorities after this PR**: a true water-quality temporal model with *same-month sampled TP or a pre-defined moving-year average*, measured hydrologic/WWTP source completeness in France, independent 2018+22 EQR assessment harmonization review, and preferably an external region/method panel with actual field flower/host and insect emergence measurements. Choose one source/time definition before looking at additional signs. This original 263-site study alone cannot establish causes of within/between discrepancy.

No eBird; no `gesui/main` or `chocho` modifications, prior PR #64–70 remain Draft and their positive/null source results unchanged.
