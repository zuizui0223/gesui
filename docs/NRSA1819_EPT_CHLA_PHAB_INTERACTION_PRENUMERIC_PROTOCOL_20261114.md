# Pre-numeric design: context-dependent EPT–periphyton chlorophyll association

**Date 2026-11-14. Exploratory/posthoc relative to known PR #77 EPT–Chla negative partial association. This protocol is fixed before inspecting any interaction coefficient.** Use no eBird. Keep unmerged existing PRs #73–77 unchanged.

## Population and physical scope
- Original EPA NRSA 2018–19 genuine same-visit source n1,221 frozen native UID from [source-success #43528570410](https://github.com/zuizui0223/gesui/actions/runs/43528570410); strict actual EPT counted fraction, periphyton areal Chla, in-stream total N/P and PHab sample n946 /603 HUC8 /47 US states from [#77](https://github.com/zuizui0223/gesui/pull/77). Reconstruct source from pinned EPA CSVs, no previous model predictions used as new field values.
- Primary observed outcome: benthic `BERW/TOTAL300` counted EPT individual fraction, **not diversity as such**, mass or adult emergence. Primary cross-biota predictor: `log10(areal chlorophyll-a µg/cm²)`, **not total AFDM or primary production rate**.
- All models include original complete PHab canopy fraction XPCAN, ground vegetation XGVEG, overhanging vegetation XFC_OHV, human disturbance W1_HALL, log10(wetted width XWIDTH), midpoint canopy density XCDENMID (%) and streambed embeddedness XEMBED (%); state, eco9, lat,lon, year, seasonal sin/cos, stream order, conductivity, in-stream TN/TP. PHab main effects mandatory even if interaction is null. All measured source columns and native units are SHA verified.
- True physical measured `XCDENMID` has one reported value >100%; include unchanged in primary and run sensitivity excluding >100 only, never truncate silently. Ground vegetation XGVEG can exceed one; it is not a probability. Original repeated SITE_ID may have more than one native UID; whole HUC8 CV and unique-site sensitivity handle it.

## Fixed interaction falsification (posthoc, not confirmatory)
Same **n946** observations and 603 heldout HUC8 groups for every comparison. Before looking at sign, test one moderator at a time:
1. `XCDENMID`: midchannel canopy density (physical shade proxy, not measured solar radiation);
2. `XEMBED`: streambed embeddedness;
3. `XPCAN`: canopy presence fraction;
4. `XGVEG`: riparian vegetated ground index;
5. `log_width`: log wetted stream width.

All log Chla and moderators standardized mean 0, SD1 on full analysis cohort for interpretable coefficients; **train-only** centering and scaling in heldout Ridge model (interaction transformed from raw columns with fold-specific fit, not leaked heldout outcome). For inferential contrasts, covariates include linear main effects + both measured river TN/TP and `logChla²` to reduce spurious interaction from a curved Chla main trend. No sign-based outcome transformation or moderator selection.

- Model B: original PR77 EPT fraction ~ spatial/region/river PHab + measured TN/TP + logChla.
- Model Q: B + logChla².
- Models I1–I5: Q + logChla×each of the 5 fixed PHab moderators separately. Model Iall: Q + all five interactions.
- 5-fold **whole HUC8 holdout** Ridge10 MSE improvement over Q and fixed 2,000-bootstrap on heldout HUC8 group errors; no leakage, cluster resampling describes error stability not independent dataset.
- Separate OLS HUC8 clustered **descriptive** β and confidence interval for each moderator interaction with main effects and squared term, plus Holm familywise multiple-test adjustment (5 prespecified tests); sample-size guards n946, HUC8 >=500, states >=40 and each CV fold >=100 visits.
- Effect scale: change in EPT counted fraction per 1SD log Chla at the 25th vs 75th *source PHab modifier*, holding model covariates constant; do not call it treatment effect.
- Controls/sensitivities: first unique physical SITE_ID visit, >=250 benthic counted subsample, excluding the single raw XCDENMID >100 recorded anomaly, and year-specific panels if n>=300; preserve all outcomes including nulls and weak signs. Evaluate only strongest *source-defined* overall interaction model without search-optimizing knots or quantile thresholds. Stop if evidence insufficient; avoid false confirmation.
- The raw cross-biota Spearman ρ≈−0.049 from previous PR77 is low, and standalone extra Chla predictive MSE +1.09% has 95% cluster range spanning zero. Any interaction is exploratory and should not upgrade this to causality.

## Scientific interpretation boundaries
The question is whether the **conditional cross-biota association varies with measured stream habitat**. Neither modifier is randomized; Chla and EPT may respond jointly to unmeasured flow, algal composition, season, contaminants or both, and input nutrient sources are not proven wastewater-specific. N/P effluent original 2018 source support failed as documented in [#74](https://github.com/zuizui0223/gesui/pull/74). No true insect emergence biomass, pollination, ecosystem service or sewage causal chain inferred.

Official [EPA National Aquatic Resource Surveys data page](https://www.epa.gov/national-aquatic-resource-surveys/data-national-aquatic-resource-surveys), original publisher SHA and date/site ID source #77.
