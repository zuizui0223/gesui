# Posthoc four-group PHab ablation, frozen before its effects

**2026-11-14. Explicitly POSTHOC to observed PR #81 positive original six-PHab EPT external-prediction results; all four groups and contrasts fixed BEFORE inspecting group effects.** The full set of physical/human factors predicts original counted EPT fraction in three EPA survey periods, but its mechanism has NOT been identified. This diagnostic is exploratory, not a preregistered variable-selection or causal attribution study.

Use **EXACT** original source-confirmed n586/437 HUC8 2008, n1116/706 2013, n946/603 2018, canonical HUC8 and the three pre-locked original source measurements from [PR #81](https://github.com/zuizui0223/gesui/pull/81). Do not change original UID sample, field chemistry or HUC8 data; same whole-survey-year heldout target with all shared HUC8 removed from two training cycles, HUC-equal train weights and HUC-equal target MSE. Train-only category encoding and numeric scaler, Ridge α10, 2000 fixed HUC8 source-block resamples as descriptive95.

Group the exact original six available metrics by their **measurement semantics, not output sign**:
- **Substrate:** `XEMBED` measured riverbed embeddedness (%).
- **Canopy/cover:** `XCDENMID` midchannel canopy density, `XPCAN` canopy presence fraction, `XFC_OHV` overhanging vegetation.
- **Channel geometry:** `log_width` of source wetted river width.
- **Human disturbance:** `W1_HALL` proximity-weighted human stress index. This is not necessarily a physical-process causal predictor, sewage discharge, or actual nutrient dose.

For all four groups report in all three independently held-out source cycles:
1. **Standalone added value:** original geographic+ambient TN+TP baseline vs baseline+that group ONLY (NOT the full six PHab model).
2. **Conditional added value:** full six-PHab+ambient TN/TP model vs full model EXCLUDING that group (equivalently adding the group back to the other three groups). Correlations between PHab groups mean these contrasts are **not additive independent contributions** and a negative conditional added value does NOT prove ecological harm.
3. The identical source full-PHab and no-PHab model MSE must match original signed PR #81 output exactly before computing any new group contrasts; no cherry-picked group or fixed subset chosen based on signs.
4. Report group and period HUC8-equal MSE improvement percentage and whole-HUC8 fixed-prediction 2,000 bootstrap descriptive95, plus all null/worsening cases and train sample sizes. No p-value significance hunt.

**Model limitation:** this is measured predictor ablation, not `XEMBED` environmental manipulation or a cause of biodiversity change. The 2008–09 source originally lacks periphyton Chla and AFDM, and none of the cycles has verified actual municipal effluent N/P delivered at surveyed reach, insect emergence dry mass, predator food subsidy or realized ecosystem services. All data are authentic original EPA, eBird not used, maintain PR Draft unmerged.
