# Further test of 2023 HUC8 source-novelty: substratum and hydrological domain
2026-11-22, exploratory posthoc follow-up to PR #142 (its ecological outcome signs and prior 214 matches already observed).

**Question:** Does measured *substrate composition* (silt/sand/coarse particle availability) or river hydrological size, which may not be represented by XEMBED alone, explain the difference in externally trained EPT predictions between previously sampled and source-novel 2023 HUC8? This is a possible explanation, not a causal mechanism.

## Step 1. Unblinded context, but source-only new variables
- EPA original signed 2023 SiteInfo / PhysicalHabitat / WaterChem / Periphyton files; reconstruct the exact prior 1,165 UIDs /693 HUC8, and the 214 source-only matched pairs from PR142 with unchanged HUC8-pair identity SHA256.
- Before examining any additional matched biological residuals, inspect explicit EPA physical field metrics: **PCT_FN**, **PCT_SAFN**, **LSUB_DMM**, **LRBS**, **XDEPTH**, **XFC_NAT**; site hydrology **MEAN_FLOW**, **WSAREASQKM**; site **URBN_NRS23** and riparian disturbance **W1_HAG**. Report actual presence, non-NA fraction, ranges and invalid/sentinel counts, including method. Do not silently interpret source missing as zero; no EPA outcome values read.
- Fixed two scientifically distinguishable blocks, provided enough genuine source coverage:
  1. **Granulometry**: PCT_FN (percent fines), PCT_SAFN (percent sand+fines), LSUB_DMM (log geometric mean bed particle diameter in the original EPA source definition). XEMBED is a related but **non-equivalent** physical habitat metric; its field collection protocol mechanically sets embeddedness to 100% for sand/fines, so XEMBED may conflate burial with particle mixture.
  2. **Hydro-landscape**: log10(MEAN_FLOW [modelled cfs]), log10(WSAREASQKM [km²]), and binary urban category URBN_NRS23. The flow is modeled mean annual flow, not 2023 sample-day discharge. URBN_NRS23 is categorical and must not be converted to numeric.
- Source GO for each *block* if all required fields are present, meaningful, and observed for >=80% of the **original 214+214 source-only matched HUC8** on each arm. For genuine quantitative rematching, source-GO additionally requires >=70 match pairs />=5 eco9 under fixed Euclidean caliper **3.5** after adding the block, and overlap in ecoregion between groups. Report STOP without relaxing thresholds.

## Step 2. Source covariate divergence and matched falsification
- For source-go blocks, compare unadjusted 2023 HUC8-equal source group environment differences and differences in the original 214 pairs (same pre-outcome matched SHA). Compute paired shifts using fixed metrics and report source missingness per arm.
- Re-match on the original 7 source-matching features (lat,lon,stream order, river logTN/logTP, XEMBED, log width) **PLUS** exactly one chosen GO block; never choose a model based on the outcome sign.
- HUC8 equally averaged source features, numeric variables standardized using only the original 551 older TRAINING HUC8. For the new substrate block, old training original source may **lack the same measured variable**, so if identical training-period source variables are not available, source-only standardize these supplemental covariates using the **2023 original combined 693-HUC** marginal distribution, without outcome values. Report which standardization was applied. Caliper 3.5 and minimum-pair GO remain fixed.
- Once a block's source-only matching is frozen in an independently signed run, evaluate **unchanged** original past-era-trained EPT predictions, with 2023 train-target HUC8 zero overlap, and differences in paired HUC8-equal MSE reduction for (a) all six PHab after in-stream N/P, and (b) XEMBED alone after N/P and after the other five PHab.
- 2000 fixed-prediction bootstrap matched-pair 95% for a difference of gains and, separately, *absolute* squared error reduction (no pooling of normalized MSE effects with independent data). Negative or inconclusive results stay reported.

## Interpretation and hard limits
- Four EPA survey cycles are within one national survey family, not global datasets. The 214 new HUC8 are *new only to the previously eligible original analytical source panels*, not never-sampled rivers. EPA source sampling/quality flags may differ by time or target group. Additional matching does not randomize past source inclusion.
- Sediment particle fraction and XEMBED are not experimental treatments; any loss/gain of residual difference cannot prove sediment mediates historical sample-selection effects or sewage causality.
- No source of delivered POTW phosphorus has been measured; this must remain separate from PR #139's original sewage exposure STOP.
