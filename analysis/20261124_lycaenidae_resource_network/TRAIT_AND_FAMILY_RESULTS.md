# Source-verified trait partitions and stricter family composition falsifier (2026-11-24)

Exploratory reuse of the original GBIF H7 panel and static catalogue; no new field observations. Both analyses use H6-disjoint training-only 2017–19 geographical null, exact within-stratum expected graph, early EQR added after original covariates and early network excess, Ridge10, H6-equal fixed-loss bootstrap. Plant-side symmetric term includes the average of early and late consumer compositions; a separate early-fixed consumer sensitivity excludes late consumer information. All model fits are observational, not wastewater effects.

## Original remote trait partition

[Successful run 45381603761](https://github.com/zuizui0223/gesui/actions/runs/45381603761), artifact 14912974425, ZIP SHA256 `4cd03196f4fee45fdbc8ebc7e5c839555c961a8cdac29e4c87472578719d9fe8`. The five output files were extracted byte-for-byte from this digest-verified ZIP. Source gate 45381165806/14915055021. Trait projections operate **within plant family** on matched hosts, using orthogonal projections of the matching leverage, not actual management exposure or a causal mediation. Life-form 537 exactly scored keys, 403 eligible keys in 16 families; disturbance/management 378 scored, 320 eligible in 20 families. Both periods >=50% matched eligible resource inclusion mass; quality n173/77 H6; all n222/89 for life form and 221/89 for other domains. The stricter 80% coverage has only 12, 0, 0 quality H7 respectively; no strict 80% result. Trait-identity gain, quality: life form +0.37% [-4.00,4.41]; management +0.30% [-1.50,2.41]; disturbance +2.86% [-0.70,6.42]. All: -1.12%, +1.58% [0.08,3.80], +2.87% [-1.14,6.67]. No trait-versus-residual contrast survives Bonferroni-3. Failure to identify these mechanisms does not prove ecological traits irrelevant. No measured nitrogen fixation, agricultural use, biomass or management events.

## New family/resource partition, local source replay

Reconstructed all 454 original H7-periods, original 303 butterfly and 990 host keys, 1565 edges, 222 paired H7/89 H6 and quality 173/77. Verified the six ZIP digests and original 499-permutation network replay. The stricter null has exact within-family row-sum conservation (also retains narrower training-only strata); tested each of 71 plant families. The symmetric decomposition of plant-side change into amount and identity uses the mean of endpoint leverage and standardized host inclusion masses, with empty group composition carried from the other period so zero mass has zero contribution. All decompositions reproduced to 1e-12. An all-host mass term is expected number of catalogue host taxa represented in 50 GBIF records, not biomass or physical abundance. For 13 endpoints x two source cohorts x two consumer definitions, 52 fits with identical fixed H6 folds and two model implementations; 4 paired fixed-loss contrasts per cohort/mode. Bootstrap reuses sorted H6 losses so small CI differences relative to previous report are Monte Carlo ordering, not changed point estimates.

| Symmetric plant-side endpoint | quality 173 gain [95%] | all 222 gain [95%] |
|---|---:|---:|
| total | 7.328 [0.840,13.492] | 8.228 [2.614,14.210] |
| all-host mass | -2.595 [-4.440,-0.229] | -0.440 [-2.434,0.263] |
| all-host identity | 9.309 [1.989,16.568] | 8.477 [2.489,15.394] |
| family mass | 2.691 [-1.100,5.676] | 4.096 [0.236,8.692] |
| within-family identity | 3.866 [-2.226,9.977] | 3.785 [0.076,8.171] |
| Fabaceae delta | 1.387 [-1.093,4.598] | 1.612 [-0.286,3.917] |
| non-Fabaceae delta | 3.410 [-2.923,11.170] | 6.017 [-0.113,13.716] |
| actual catalogue graph delta | 5.038 [1.039,9.176] | 5.447 [0.868,10.566] |
| null graph delta | 0.902 [-0.512,2.419] | 1.802 [-1.038,4.811] |

Symmetric all-identity minus all-mass gain: +11.904 pp (Bonferroni-4 [0.666,22.159]) quality and +8.917 pp [1.167,18.811] all. Fabaceae minus non-Fabaceae and within-family versus family mass contrasts not supported. Exact graph versus null contrast not robust after correction. Early-fixed all-host identity all 222: +11.297% [2.649,20.496], quality 173 +7.488% [-2.076,15.121]; corrected paired identity-minus-mass contrast for all narrowly crosses zero. The coefficient of EQR for symmetric total and identity is negative in every fold. No claim that plant family differences alone explain the association or that legumes mediate it. Fold-specific outcomes and fixed-error intervals, not independent field replication or full uncertainty. Source association remains observational; direct measured TP and WWTP-specific chain remain unsupported.

Local numerical check: original source and identities verified, independent sklearn predicted values agree within 1e-8. See script and full losses, contrasts and results. Earlier results retained. Main unchanged.
