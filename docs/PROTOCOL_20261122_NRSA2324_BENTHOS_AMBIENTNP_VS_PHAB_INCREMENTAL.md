# EPANRSA2023–24 source-locked posthoc additive drivers: original nutrients versus physical river habitat
2026-11-22, after examining the 4 ambiguous canopy×ambient TP interaction results in this draft branch. **This is explicitly posthoc, not a confirmatory replication.**

## Scientific question
Does the available measured benthic community generalize better across whole HUC8 basins from the two measured chemistry variables (original stream TN/TP) or from original habitat width and bank/midchannel canopy variables? Unlike PR139 these are not measured sewage-delivered exposures. The new question uses the exact four original strata/outcomes already audited in the interaction run with source SHA and all rows preserved, no result-sign selection.

## Complete-case cohort locked before examining additive gains
Reuse `scripts/gesui_nrsa2324_tp_canopy_benthos_interaction_20261122.py` source validators and `cohort()`. Four original groups: BERW MMI n1165/693 HUC8, BERW EPT n1165/693, BETB MMI n725/476, BETB EPT n725/476. Same original five-table EPA SHA and date/UID gate, all valid TP/TN/cond, canopy, wetted width, original site coordinates, state and ecoregion; no refilling missing or pooling methods. Must re-match original main-effects MSE precisely before new conclusions.

## Fixed model families on IDENTICAL observations and folds
- B = geography/stream order/year/season/stream conductivity/state and eco9.
- B+N = B plus observed log10(TN) and log10(TP).
- B+H = B plus original log10 wetted width and canopy percent, separately XCDENMID BERW and XCDENBK BETB.
- B+N+H = all main effects without interaction; this must equal earlier interaction study's reference model.
No interactions, nonlinear additions, outcome transformations, or class-specific variable searches in this follow-up. Fit train-only numeric standardization and category one-hot, Ridge alpha10; identical deterministically assigned five HUC8-disjoint folds. Report group-equal and visit-weighted heldout MSE for all four models. Key effects: nutrient gain AFTER habitat and habitat gain AFTER nutrients; also each after B. 2000 HUC8 fixed-OOF bootstrap percentile95 of each percentage gain. Report all four biological outputs, no p-hunting or cherry-pick. This is incremental predictive value, not cause or mediation.

## Interpretation
Do not claim observed in-stream phosphorus proves sewer/WWTP source; do not infer shade mitigation or aquatic-to-terrestrial flux from main effect gains. Compare qualitatively with preexisting 2008/2013/2018 cross-era EPT findings but do not call this fourth *independent experiment*, since it reuses NRSA and model architecture. A tiny HUC8 OOF gain with CI crossing zero is uncertain.
