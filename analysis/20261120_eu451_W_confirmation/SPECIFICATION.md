# EU451 wastewater–butterfly confirmation: frozen specification

2026-11-20. Written and committed before executing this focal comparison. This is a posthoc confirmation on previously analysed observations, not a preregistration on unseen data and not independent replication.

## Correction to the premise

The existing `analysis/20261115_eu454_plant_nutri/h6_model_if_gate.py` already includes `crop`, `urban`, `riverflow`, `naturalout`, and `area` in both P_GEOG and B_GEOG. The 20261115 butterfly family decomposition also includes these in COMMON. Earlier analyses were not adjusted for urban cover alone. The 451 observations are H7 × period rows, representing 232 H7 catchments and 90 H6 parents, not 451 distinct catchments.

The existing 451-row plant–butterfly analysis focuses on adding plant NUTRI after W and geography. This confirmation instead explicitly tests adding W after the same non-plant-response controls. No claim that it introduces cropland adjustment for the first time.

## One primary comparison

- Cohort: exactly the original 451 source-quality-eligible H7 × period rows, original plant gate and original butterfly panel; no new selection, countries, species, or family-specific subgroup.
- Outcome: `primary_incidence_mean`, the species-equal known larval-host nutrient-habitat affinity of the five butterfly families combined. No switching to record-weighted, within-family, or Lycaenidae outcomes after results.
- Focal exposure: original static W, standardized to one SD of this same cohort. It is modelled wastewater-associated source geography, not measured delivered nitrogen.
- Numeric controls: `crop`, `urban`, `riverflow`, `naturalout`, `area`, `nEEA`, `nAnimal`, `nPlant`, `nButterfly`, `primary_record_coverage_raw`, `primary_host_completeness`, `NUTRI_coverage_all_Plantae_speciesKeys`, `log1p_eiv_scored_taxa`.
- Categorical controls: original `country_method` and `period`.
- Do not add Plant_N_incidence or EQR to this primary model: the question is total conditional W information beyond landscape/observation conditions, not information after a proposed biological mediator. Do not add family proportions: the outcome is whole-community composition. Climate and elevation are not in this frozen adjustment set; do not claim all confounding removed.
- Estimator: reuse the original h6_model_if_gate.py H6-equal WLS, H6 cluster-robust standard error with small-sample correction, normal-approximation 95% interval. Report W coefficient in original HNS units per W SD.
- Predictive comparison: baseline is all listed controls; augmented model adds only W. Reuse the original deterministic 5-fold H6 hash split, training-only scaling/categorical encoding, H6-equal Ridge alpha=10. Paired resampling of 90 fixed H6 mean squared errors, 2,000 draws with original seed 20261115, yields the descriptive interval for percentage MSE improvement. It does not include full training-model uncertainty.

## Fixed decision rule

1. If W coefficient has a positive 95% lower bound AND H6-heldout MSE improvement has a positive 95% lower bound: support positive landscape-adjusted association plus additional H6-heldout predictive information, within this cohort.
2. If coefficient is supported but predictive improvement is not: support association only; generalizing extra information remains uncertain.
3. If predictive improvement is supported but coefficient is not: predictive information only, no stable positive linear association established.
4. If neither meets its criterion: this comparison does not support a positive independent W signal. Do not infer zero wastewater effect, cropland causation, or a proven common nutrient source.
5. Preserve a supported negative result if present. Never translate failure to support positive into a positive landscape-only causal conclusion.

No additional model/outcome search in this confirmation. Source gate/hash/rank/leakage/numerical checks are validation, not new hypotheses. If eligibility or design checks fail, stop the focal fit and report that reason without relaxing criteria.

## Outputs and scope

Save the frozen-specification commit, source/code hashes, exact covariates, coefficient and predictive result in JSON, then a Japanese manuscript note and PR update. Virginia biosolids data remain candidate material for a separate field study, excluded from this paper's ecological analyses. French sludge/geography downloads already in flight may finish caching, but their spatial analysis is deferred while this requested confirmation is completed.
