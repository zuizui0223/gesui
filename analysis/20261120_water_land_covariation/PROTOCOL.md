# Water–land covariation: analysis specification

Written before fitting the new models, 2026-11-20. Exploratory follow-up to previously inspected data; not preregistration.

Question: do aquatic and terrestrial indicators covary across matched catchments, or change together between two observation windows, and does static wastewater-related geography explain part of that covariance?

Use the frozen FR73 physical-station-matched 2018/2022 water source, original 454-row butterfly panel and verified species facets, and eligible plant summaries. Join by H7 and period; retain original H6 parent clusters. No new occurrence collection or source substitution.

Primary cohort: all 73 catchments, five-family species-equal host nutrient score. Lycaenidae: separate sensitivity restricted to both-period original support (10 scored records, 2 species, 70% records scored). No family specificity claim. Spatial analysis averages the two periods per catchment; temporal analysis subtracts early from late. Water years are 2018/2022, while terrestrial windows are 2017–2019/2020–2022: no instantaneous or annual synchrony claim.

Responses: logTP, minus EQR (reversed for display only, not a nutrient-specific metric), plant incidence EIV NUTRI, all-butterfly incidence HNS, and supported Lycaenidae incidence HNS. Do not pool them into one ecological scale. Primary cross-realm pairs: minus EQR with plant and with all-butterfly HNS. TP comparisons are secondary; within-realm pairs and Lycaenidae are context/exploration.

Fit parallel weighted least-squares responses, with equal total weight per H6. Report residual Pearson correlations and standardized residual covariances under (i) intercept only; (ii) observation controls; (iii) observation controls + W; (iv) observation controls + landscape; (v) observation controls + landscape + W. Observation controls: plant and butterfly record opportunity, plant score coverage, butterfly score coverage, water site support (spatial) or change in water sampling opportunity (temporal). Landscape: crop, urban, area, riverflow and naturalout; use period means. Temporal observation controls are differences. W is period mean, not change in exposure. Climate covariates are unavailable in these frozen inputs; do not claim adjustment for climate.

For nested W models report the paired change in residual covariance and correlation with an H6 cluster-bootstrap interval. A decrease is an association-based decomposition, not causal mediation; an increase or interval spanning zero is retained. No percent explained when baseline covariance is near zero. Also report response-specific W coefficients and shared fitted covariance. Refit all models in each of 3,000 H6 resamples, retaining paired resamples across specifications. Use fixed original outcome scales; drop only exact constant controls, flag rank failures. Report marginal 95% intervals and simultaneous intervals over the four primary EQR cross-realm correlations (2 pairs × 2 spatial/temporal modes) per specification. No independent-study meta-analysis.

Validation: frozen hashes, unique keys and cohort sizes; correct period matching; paired difference identities for water measurements; full model rank; fitted/residual covariance decomposition; a synthetic shared-driver check; deterministic seed; source hashes and runtime versions in receipt. Inspect results and plot before integrating the report. Keep all null and contrary results. Existing PR branch only; no merge.

Implementation: one Python analysis script, aggregated CSV/JSON results and figure, and a Japanese results note linked from the manuscript index. This is the execution plan for the user-authorized comparison; no unrelated analyses.
