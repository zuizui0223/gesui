# Temporal resolution of adult bait identity — 2026-11-25

Frozen protocol 0b46e007afebcd213291c2c2aff888f1ee3982e4, following PR176 binary carrion/banana positive. Salazar-Carrión et al. (2026), Ecuador 2017, DOI 10.1017/S0266467426100510; same previously opened field experiment, not independent replication. Exact original signed source artifact15011284503 SHA256 5adaf08823a3d56e01d82fa2c15f6b2275a976fb07fde0488dde4a288d05c780 and all three native workbook hashes verified by byte-identical prior audit d5f1823345e627c9f6a6216a7828fc5c8719352ac39292c3c3dadc25f73ad7df. Same 27 taxa, unique verified first captures, same 14 event blocks. All 27 retained; 24 have >=3 animal-bait captures per half. Animal capture counts: early443 (freshwater145, marine136, shrimp162); late329 (freshwater97, marine119, shrimp113). Source gate passes, no new outcome-based selection. No new larval host, ant or adult substrate availability observations. Preserve source ID and stratum flags conflicts, single-site, unknown sex and deployment-opportunity limits.

## Conditional 3-category animal-bait prediction
Training-only smoothed family×forest-stratum offsets, symmetric penalized multinomial species intercept, precision2; native genus baseline secondary. Forward primary, reverse descriptive, 2000 paired fixed-event-loss resamples, 4999 conditional training label permutations preserving three-category margins in family/genus×stratum×event groups.

| Direction | Baseline | Species increment % [conditional95] | Permutation p | Improved events |
|---|---|---|---|---|
| early->late | family×stratum | -3.940990 [-5.772257,-2.215276] | .8132 | 0/7 |
| early->late | plus genus | -3.498123 [-4.259216,-2.682770] | .4214 | 0/7 |
| late->early | family×stratum | -4.491124 [-7.506547,-1.244102] | descriptive | 1/7 |
| late->early | plus genus | -3.175673 [-5.082447,-1.139047] | descriptive | 1/7 |

Primary fails. Forward family baseline mean logloss 1.099772515 -> 1.143114444 (gain in nats -.043341929); species-equal gain -4.503393%; all leave-one-event gains -4.391682 to -3.184820%; all leave-one-species gains -4.470425 to -3.318468%. Negative gains are predictive loss worsening for the frozen model, not absence of physiological discrimination. Observed independent direct-record optimizer max error 1.268e-8; permutation gradient residual max1.36e-14. Independent arithmetic replay of all 4 gains, intervals and 2 p-values max absolute error2.22e-14.

## Post-result common-denominator descriptive comparison
Exact chain rule for four-class log loss: L4 = Lbinary + I(carrion)*Lconditional. Combine prior saved binary and new conditional heldout event losses with exact same original event capture counts and conditional fractions. No new fits. Same 835 verified first captures; event-equal denominator. Prior source archive SHA256 4a1642d78d0f7c9de6e8d18d11344ff252622dffaba76e5f21e78c22402df0e8 verified.

| Model with species identity at which resolution | Forward four-class logloss / gain % | Reverse four-class logloss / gain % |
|---|---|---|
| neither (family/stratum both levels) | 1.293973 / 0 | 1.292301 / 0 |
| broad only | 1.261899 / +2.478715 | 1.253913 / +2.970578 |
| fine only | 1.334603 / -3.139961 | 1.336773 / -3.441272 |
| both levels | 1.302529 / -.661245 | 1.298384 / -.470693 |

This decomposition is an exact descriptive architecture comparison made after outcomes, not a preregistered contrast, new p-value, equal-effect test or proof of a universal hierarchy. Earlier binary +11.5983/+13.5943% gain and conditional three-class gains use different denominators; the common four-class summary solves that descriptive denominator difference. Candidate interpretation: identity contains temporal information at broad carrion/fruit resolution but does not improve prediction of individual animal bait types in this cohort. Only one lycaenid, mostly Nymphalidae; no ant integration, actual ingestion, physiological nutrition, sewage effect or independent regional generality. Stop screening this opened dataset; next distinct evidence needs external adult behavior data with larval host/ant context. Main unchanged.

## Reproduction
Complete local derived archive gesui_carrion_resolution_20261125.zip SHA256 e51cbc0a1c501b7c03bb2c708f2c404b14cc73e22b3735c9a5920fb1aec782f0, 29 members. audit.py SHA256 dcb88bb98a871c50796130b5ffea25e79ebaef42b59fed32859cdd59daf68837; results.json 3aa36536c8976cf432c83b9f39ad79804613f11971e10635d226f6eee4c6d7bf; hierarchical_loss.py 8016ff01311745523bc0a42b3ec69d1ad916e339b56a368313bb46f653d97b87; hierarchical_summary.csv 5ba089559d90f6ac1eba17915ff480f7a472dc24625b81a886f25c3cd600285d; independent_loss_replay.csv 064e71b16f1ecf15ce6fb4c06ca6fec8a038fc43b192e5f7552d8d09f81f9a66. Complete local code and numerical tables in archive; GitHub protocol and report only; remote numerical replay not performed. No raw specimen or image source redistribution.
