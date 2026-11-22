# Butterfly within-genus nutrient-affinity recording test

Status: new exploratory analysis authorized by the user's 2026-11-22 「進め」. The biological coefficient has not been fitted when this specification is committed. Existing whole-community results were already known; this is not independent preregistration or a new validation dataset.

## Question and source correction

Within the same true H6 parent catchment and period, do the relative recorded identities of congeneric butterfly species align with local plant nutrient habitat affinity? The historical p=0.6864 permuted **plant** NUTRI within **plant** genera while butterfly HNS was fixed. It did not test within-butterfly-genus sorting.

## Fixed inputs

- Base repository commit: 6a0ebac1274f23cd202e2f515176b6d77d264b54.
- Original 183 H7 / 40 H6 geography; periods 2017-2019 and 2020-2022; no new locations or taxonomic imputation.
- Original old facets SHA256: 251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4.
- Original additional facets SHA256: e076f6c8e2d93282be19aea3c75ae8431302318f068b9d75e278e2a21c20aafb.
- Original paired-H7 panel SHA256: c0840b1c5c6e85342fea41d57dc40fa5ec295fabf6bcd6908a4072557db9a736.
- Original 303-key HNS table SHA256: 00c1a36a9237a3d1a6d77334aa8baa5dc8b329a49925c2e31f033fa2d302ecdb.
- ACCEPTED butterfly records, numeric HNS, no actual-used-alias conflict; preserve current frozen genus names. The known full eligible roster is 186 species / 88 genera. No genus score filling.
- Original EIV crosswalk 11,511 keys / 4,011 numeric accepted plant species. Rebuild same-period species-equal plant NUTRI from raw facets; compare its two-period H7 mean, and all reconstructed original butterfly CWMs, to the original panel within 1e-10.
- Original 13 controls: W, crop, urban, naturalout, area, nPlant, nButterfly, nAnimal, plant_coverage, hns_coverage, hns_host_completeness, plant_scored_species, hns_scored_species. Use the already frozen **two-period H7 mean values for all controls**, repeated in both periods. This choice preserves the published control set and does not invent missing period-specific animal effort. The focal plant NUTRI is same-period. Standardize X across all 366 source cells; standardize controls across all 183 original H7. No outcome-driven control selection.

## Approaches and chosen estimand

1. A decomposition of community HNS into between/within butterfly-genus terms is easy but its within term can still change merely because locally biased genera change frequency.
2. A full conditional-subset model can condition on a site-genus total but needs numerous species intercepts and empirical regional pools.
3. **Chosen:** pair-stratified conditional logistic composite likelihood. It conditions out arbitrary species-pair baselines separately for each H6 and period and directly answers the within-genus question. It changes the estimand to relative recording among discordant congeneric pairs; it does not estimate occupancy or abundance.

For each unordered scored butterfly pair within a genus, fix A/B by lexicographic speciesKey (never by HNS). For every H6 × period, keep H7 cells where exactly one of A and B was recorded. Define y=1 when A alone is recorded, y=0 when B alone is recorded. Condition on the number of A-only cells within each H6 × period × pair stratum. A stratum is informative only when both outcomes occur. Both recorded / neither recorded cells provide no discordant-pair information and are not called absences. Retain originally equal-HNS pairs within variable-HNS genera so permutation eligibility is not changed by the trait labels.

Linear predictor, after the stratum intercept is conditioned out:

`eta_i = (HNS_A - HNS_B) * (beta * z_plant_NUTRI_i + sum(gamma_k * z_control_ik))`.

The conditional likelihood sums over all subsets with the stratum's observed success total. Its product across reused pairs is a **composite** likelihood, so ordinary independent-observation standard errors are not used. Each species is coded recorded/not-recorded regardless of record count. Give each unordered pair equal total **stratum weight**: weight=1/(number of informative H6-period strata for that pair), with a common normalization. This does not imply equal information per stratum or individual.

Primary effect: beta and exp(beta), the conditional A-vs-B recording odds ratio for +1 original-source SD of plant NUTRI when A has HNS one point greater than B, holding the original controls fixed. Primary scope: this fixed observed species roster, empirical geographic overlap, and two periods.

## Source and numerical gates

- Input hashes, species-family identity, unique H7-period keys, all366 original cells and original CWM replay must pass.
- At least20 informative H6 and10 informative butterfly genera for the main inference. These are pragmatic sufficiency cutoffs, not proof of power. All source eligibility counts are reported before fitting.
- Center each design column within its conditional stratum to diagnose rank. Remove only invariant/exactly redundant nuisance columns by deterministic QR; report every removal. The focal column must add rank beyond nuisance columns, or report identification STOP. No ridge penalty or new covariate is substituted after seeing a coefficient.
- A fit must converge, have finite coefficients and a positive identified Hessian. Report gradient and conditioning. Do not silently regularize or discard a failed full model.

## Uncertainty and fixed additional checks

1. Primary95 interval: H6-clustered sandwich of the weighted composite score, CR1 G/(G-1), t quantile with G-1 degrees of freedom. Both periods, H7 and all pairs remain together in each H6 cluster. This conditions on observed species and does not establish transfer to all butterfly genera.
2. H6 nonparametric bootstrap399 draws, seed20261122; refit with entire-H6 multinomial weights and frozen source scaling/pair weights. Report every failed replicate. Percentile interval is reported only if at least95% converge; otherwise bootstrap uncertainty is unresolved.
3. Descriptive trait-label reference999 draws, seed20261123: shuffle HNS once globally per species, within each butterfly genus among species in any informative pair. Keep cells, pair orientation, eligibility, controls and weights fixed; refit the full same model. Include formerly tied pairs. Reference p=(1+number of |beta_perm|>=|beta_obs|)/(1000). This is **not** an exact beta-only null in the presence of arbitrary trait-control relations. Any failed reference fits are reported; do not calculate a complete-reference p by silently dropping them.
4. Leave one butterfly genus out, refit fixed model; report all coefficients, including failures. This is influence analysis, not out-of-region prediction.
5. Fixed sensitivities: (a) retain only species with documented-host numeric coverage>=0.8; reapply identical support gates and report STOP if too small; (b) remove the four butterfly-derived controls nButterfly, hns_coverage, hns_host_completeness, hns_scored_species, retaining the primary eligible pairs and scaling. These do not replace the primary result.

No post-result search over widths, families, model forms, thresholds or directions. Host breadth/moderation and individual Lycaenidae species effects remain later tasks; do not choose examples based on the new coefficient.

## Interpretation and deliverables

Report coefficient, odds ratio, uncertainty, exact effective sample and all contrary checks whether positive or negative. A positive result is spatial alignment of butterfly relative **recording** within genus, not population increase, evolutionary constraint, causal nutrient transfer, larval feeding or wastewater effect. Failure to reject is not equivalence. Do not compare the new p with historical plant-genus p as a direct effect-size contrast.

Save new code, tests, source receipts, aggregate result tables and a Japanese results report on this isolated analysis branch; keep raw third-party bulk files in the reproducible cache, not in Git. Create a reviewable draft PR; do not merge or alter the current manuscript.
