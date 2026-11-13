# NRSA 2023–24: measured river TN/TP, canopy, riparian vegetation, streambed and periphyton pigment/AFDM

**2026-11-14 · exploratory POSTHOC observational ecology · source-verification and 2nd-code independent numerical audit complete. No actual WWTP source-specific effect and no measured ecosystem service response.**

## Publisher original sources and source-before-effect sample

[Exact EPA NARS official NRSA2324 physical habitat source](https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_physicalhabitat_metrics.csv), original CSV SHA256 `8a6094bb87f8a57cf45a35a1272b137e6a9327cefd158f4bee88b025eaf9b2f4`; official original metadata SHA256 `54dfd8cf0b65cec17719c6a4e5688ede73841cca1dbcc41c5afb604c68ae35c6`. **Source-only original exact same native UID join = 2,129 of 2,129** prior fixed NRSA2023–24 visit UIDs /1,041 HUC8 /48 states. 2,139 publisher PHab rows; 10 outside original source-fixed cohort; 0 duplicate PHab rows for original matched UIDs. [Original blind source/metadata QA SUCCESS #43459407538](https://github.com/zuizui0223/gesui/actions/runs/43459407538).

Before opening new numerical PHab data the source QA fixed 2 physical covariate tiers; all are genuine EPA measured original PHab metrics:
1. **Broad:** XPCAN canopy *presence fraction*, XGVEG riparian vegetated ground cover, XFC_OHV overhanging vegetation/instream cover, log10(XWIDTH) actual wetted width. Source availability for each ~2,128 original UIDs; these are field structure proxies, **not measured solar radiation or direct shade time series**.
2. **Dense-substrate:** broad + XCDENMID *mid-channel canopy density (%)* + XEMBED *streambed embedding (%)*. Original source availability ~1,304–1,305 UIDs; restricting to this tier changes the sample and survey reach/method selection. No field was chosen or removed according to effect sign. XGVEG components may overlap: not automatically clipped into an artificial 0–1 probability.

Original EPA site+measured ambient total N/P+areal chlorophyll-a+areal AFDM source n2,129 and unchanged primary numeric quality cohort n**1,828**, 946 HUC8/48 states from [PR #73](https://github.com/zuizui0223/gesui/pull/73) and [#75](https://github.com/zuizui0223/gesui/pull/75). Physical **broad n1,828/946 HUC8**, dense **n1,117/673 HUC8**, source-selected without biological-sign tuning.

## Original observed water nutrient × biofilm results after PHab geography

Response `log10(Chla[µg/cm2] / AFDM[mg/cm2])`. Same cohort and HUC8 disjoint fivefold Ridge alpha10 per tier, consistent original baseline geography/state/eco/year/season/Strahler order/conductivity; physical variables as above; ambient log TN and TP as observed from the SAME NRSA visit. Clustered partial slopes describe +1 SD logTP and are not causal N/P additions.

| Cohort / original real biofilm response | Extra TN+TP prediction after same-covariate physical environmental model | TP β per +1SD logTP, HUC8 descriptive 95% | PHab added AFTER TN+TP |
|---|---:|---|---:|
| Broad n1828 /946 HUC8, **log pigment/AFDM** | **+1.8365% MSE improvement** | **−0.08971 [−0.12928,−0.05014]** | +0.5713% |
| Broad n1828, **log areal Chla** | +0.7765% | **−0.07212 [−0.11308,−0.03115]** | +0.0715% |
| Broad n1828, **log areal AFDM** | −0.2963% (worsens) | +0.01759 [−0.01236,+0.04755] | +0.4869% |
| Dense n1117 /673 HUC8, **log pigment/AFDM** | **−0.2596% MSE improvement** (worse; uncertain) | **−0.05632 [−0.10464,−0.00801]** | **+10.3803%** |
| Dense n1117, **log areal Chla** | −0.3581% | −0.02769 [−0.07993,+0.02455] | +2.0244% |
| Dense n1117, **log areal AFDM** | −0.9453% | +0.02863 [−0.00790,+0.06517] | +2.4773% |

**Within the SAME dense n1117 subcohort, before vs after the PHab addition:** TP–pigment/AFDM adjusted slope **−0.08095 [−0.12663,−0.03527] → −0.05632 [−0.10464,−0.00801]**, and extra TN+TP source prediction gain **+1.3314% → −0.2596%**. Thus this is not merely an invalid cross-sample comparison. However physical measurements and source positive nutrient concentrations are observed simultaneously, with temporal / habitat / sampling-method confounding not resolved.

## Independent numerical and group bootstrap (same original field data)

[Original source-verified 2-tier EPA actual model SUCCESS #43459919028](https://github.com/zuizui0223/gesui/actions/runs/43459919028). [Completely independently implemented pure NumPy Ridge/OLS HUC8-sandwich + 2,000 HUC8 whole-cluster resampling SUCCESS #43460505875](https://github.com/zuizui0223/gesui/actions/runs/43460505875): reconstituted original source cohort; all 2 tiers ×3 actual responses × 4 model specification MSEs, TP/TN partial slope estimates and descriptive cluster CI match to documented numerical tolerance. *Not* new independent field observations, only code/scale verification.

**Fixed OOF prediction HUC8 cluster bootstrap descriptive95 of EXTRA nutrient gain after PHab:**
- Broad pigment/AFDM: **[+0.252%,+3.684%]** for +1.836%.
- Dense pigment/AFDM: **[−1.524%,+1.001%]** for −0.260%.
- Broad areal AFDM: [−0.983%,+0.321%] for −0.296%.
- Dense areal AFDM: [−2.037%,−0.081%] for −0.945% (model degradation).

These intervals are descriptive *resamples of fixed out-of-fold prediction errors*, not repeated retraining, confirmatory p-values, independent ecological field replicates, or causal uncertainty. Posthoc ratio is exactly log Chla−log AFDM. A positive additional predictive gain is not a direction or ecological effect size by itself.

## Conclusion bounded by data

The negative **ambient instream TP vs pigment-per-AFDM observational spatial gradient** remains under broad canopy, riparian vegetation and width controls. When actual midchannel canopy density and streambed embeddedness are also measured (smaller, selected native survey subpanel), the nutrient-only predictive improvement is NOT robust; physical environment has a much larger incremental prediction contribution. The strong decrease in extra nutrient information within the dense subpanel suggests an environmental co-gradient; not proof that sunlight causes the original TP slope or that nutrients lack any ecological effects.

**NO** EPA NRSA 2023–24 benthic Count/MMI or true EPT/raw richness endpoints were obtained (original official CSV endpoints still HTTP 404 in [PR #75](https://github.com/zuizui0223/gesui/pull/75)). **NO** contemporaneous verified wastewater source N/P (the original exact-upstream EEA/EPA DMR dual N/P gate fails in [PR #74](https://github.com/zuizui0223/gesui/pull/74)). **NO** real RSDS, productivity rate, insect emergence biomass, riparian predators, pollination services, or measured habitat manipulation. The new study is a macroscale natural gradient comparison rather than any sewage-driven multitrophic service cascade.

Follow-ups should prioritize restored official aquatic benthic original data with physical native-UID overlap, and independent measured wastewater source mass with upstream hydraulic routing, not sign-selected buffer adjustments. Never override earlier positive and null branch results, no eBird, leave DRAFT unmerged.
