# Same-H7 butterfly HNS vs static wastewater geography — 2026-11-15

**Exploratory observational posthoc macro analysis, not a BACI/downstream-upstream wastewater test, not a causal mediation/SEM, not a meta-analysis.** Same frozen EU catchments as prior PR84-86: 454 H7×period source records, 232 H7, 90 authentic H6, 2017–19 and 2020–22, 222 H7 repeated. Each period is NOT independent from same H7.

Source HNS scores: original 303 butterfly taxonomy keys retained; accepted+numeric subset 186; 117 not numeric/unaccepted are not zeros. Five-family GBIF source index records n4,400,147, speciesKey-known n4,318,269; HNS-assigned original indexed records 91.064%. The 63 host-coverage≥80% butterfly species capture only 13.676%, so their response cannot replace primary cohort.

Mean community HNS per basin-period is the original indexed butterfly-species record-weighted mean of larval-host habitat N indicator; distinct-species incidence-weighted alternate is evaluated. Neither is actual butterfly abundance, leaf N, adult Na uptake, local larval feeding, or wastewater N/P chemistry.

## Important values; W beta is adjusted difference in HNS per +1 SD modelled source index
| Model/source rows | W beta [H6 cluster95] | H6 out-of-parent heldout MSE gain % [fixed-error bootstrap95] | Country-equal gain % |
|---|---|---|---:|
| Main record-weighted n454 | +0.13006 [0.06261,0.19752] | +7.6208 [-6.2997,19.2449] | +14.7299 |
| Geographic terms + record-weighted n454 (posthoc) | +0.10802 [0.03888,0.17717] | +5.3740 [-5.9011,14.7795] | +6.4153 |
| Geographic + species-incidence equal n454 | +0.03418 [-0.00168,0.07003] | +1.4323 [-5.2622,6.4544] | -0.0151 |
| Geographic + all-six-group taxon-ID source coverage≥80% n390 | +0.12588 [0.05543,0.19633] | +10.9956 [3.1712,17.1017] | +1.3320 |
| Temporal Δ record-weighted for repeated H7 (n222) | +0.03532 [-0.01622,0.08686] | +5.4120 [0.9631,8.8577] | +5.6319 |
| Temporal Δ equal-incidence (n222) | +0.02843 [-0.00725,0.06412] | -0.2962 [-2.5954,1.2482] | -2.0047 |

Same-H7 temporal sign count: 112 increase, 110 decrease; raw mean ΔHNS -0.01725, median +0.00355. W is static in both periods: a temporal prediction does NOT show wastewater treatment change caused biotic turnover.

Country- and era-transfer heterogeneity: FR-only beta +0.01754 [-0.03668,0.07175] and H6 gain -7.286%; non-FR beta +0.26369 [0.09972,0.42766], H6 gain +10.183% but CI [-11.399,+26.834]; geographic + source80 country-equal +1.332% [-11.798,7.833]. Generalization is not uniform.

Known missing trait scores materially limit inference: adversarial species-specific unknown-host N assignments permit W beta from **-0.14956 to +0.32447**. Setting every unknown host to 1 vs every unknown host to 9 yielded mechanically identical conditional slopes because completeness was a linear covariate: DO NOT call these two distinct robustness confirmations.

Models: H6 group disjoint validation, train-only scaling, Ridge alpha 10, H6 equal fitting/error units, H6 clustered coefficient; country holdout removes shared H6, 3000 bootstrap fixed out-of-fold group-error resamples. Coefficients and prediction intervals are different diagnostics, not causal effects. Urban covariate and W strongly correlated (source rho ~0.7806).

Independent local source QA: 4,400,147 total indexed family records, and 37,785 positive species source rows; all 454 aggregation checks largest deviation 1.78e-15; alternative Ridge max prediction deviation 5.40e-13; 20 model specifications numerical-QA checked. 14/14 new HNS-macro unit tests pass. Verified local HNS macro ZIP 7,868,436 bytes, SHA256 3986e74781d2ea0934cef5b08a2b8fb86d034aff53105b41d379f6580d160c21.

Same-scale is basin period, NOT riverbank plot or same EEA monitoring site for every butterfly occurrence. See PATH_ANALYSIS_IDENTIFIABILITY_20261115.md for stage-gated path tests.
