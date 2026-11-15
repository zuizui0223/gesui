# Posthoc adversarial covariate-block and suppressor diagnostic for source-authentic 183 H7/40 H6

**Locked after seeing the main adjusted positive β=+0.1181 and high source conditional VIF_N≈16.4. New science diagnostic is intentionally POSTHOC, not a preregistered robustness confirmation. All variants and unfavorable signs must be retained.**

Source: exact unmodified `GENUINE_TERRESTRIAL_ONLY_183_H7_40_H6_PAIRED_H7_MEANS.csv` SHA256 `c0840b1c5c6e85342fea41d57dc40fa5ec295fabf6bcd6908a4072557db9a736`. Genuine 183 different H7 nested in 40 official H6, each H7 mean over 2017–2019 and 2020–2022, 173 original +10 source-prior rank selected. No new observational site or sample/trait score modification.

**Question:** is positive species-equal butterfly HNS versus species-equal plant NUTRI only induced by controls / suppressor terms? Disaggregate fixed13 covariates into **three source-meaningful blocks** and test all eight nonempty/full combinations with exact same H7 rows, H6 group-equal WLS/H6 fixed intercepts.

1. Environment block (`W,crop,urban,naturalout,area`) — static modelled wastewater geography and HydroATLAS, not actual dose.
2. Recording opportunity block (`nPlant,nButterfly,nAnimal`) — log1p source indexes, not individual abundances.
3. Trait score support block (`plant_coverage,hns_coverage,hns_host_completeness,plant_scored_species,hns_scored_species`) — source availability / possible response detection selection, not a direct causal mediator.

Report the 8 combinations: H6-FE alone (no covariate block), each of 3 blocks alone, each of 3 pairwise combinations, all 3 blocks. Use original full-cohort N SD for all β and original H6-equal weights for all 183; H6-cluster descriptive 95, weighted plant-N VIF/condition number, rank. **No model will be chosen based on sign or p.** Not an independent multiple-hypothesis confirmatory test. Primary prior fully adjusted result retained.

Report unadjusted within-H6 **demeaned** plant-N/HNS Pearson covariance weighted across H6 parents, exact same 183 sites. Contrast with the full adjusted β. Reconfirm beta_Full equals +0.11807615189599471 to numerical tolerance, otherwise STOP. If direction changes when recording effort/trait score support blocks are added, highlight as a susceptibility to conditional-selection/collider effects. If stable, still no independent nutrient mechanism.

No new ecological response or source year measurement; wastewater dose, leaf chemistry, local species abundance, larval feeding, energy flow or causal SEM are not observed. Previous within-H7 temporal cochange null and no within-H6 unique W signal remain.
