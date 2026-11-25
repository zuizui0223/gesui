# TREAM repeated monitoring: nested community persistence and limited environmental transfer (2026-11-25)

Frozen protocol b344d20fe4e0e5efc442a32e23b3d470110565e6; code 088710e870a46f95e5e761e2608b602b266aeb6e; [Actions 45542432633](https://github.com/zuizui0223/gesui/actions/runs/45542432633) SUCCESS. Signed derived artifact 14989139284 ZIP SHA256 46a0f4c1f31bd2f8ee3724d0a1c04dfb886529e829a8f9441d3ef4fbce46f581. Publisher [TREAM DOI](https://doi.org/10.5063/F1NG4P4R), CC BY 4.0. Raw CSVs SHA256 and publisher checksums checked; exact four source-only derived tables byte-identical to signed scoped source artifact 14989507257. Strict global source STOP preserved, scoped 563-site gate passed. No original files altered.

## Positive observation support
| Universe | Sites / site-years | Countries / projects | Family universe (detected) | Retained positive source records / all order-label records in used years |
|---|---|---|---|---|
| Ten fixed insect orders | 546 / 6058 | 9 / 12 | 86 (77) | 62573 / 128757 |
| EPT sensitivity | 474 / 5239 | 9 / 12 | 29 (26) | 31142 / 59394 |

Same exact label-support rosters and SHA256 identities reproduced. No invalid or zero abundance entries among clean retained records. Each of four fixed windows has >=2 positive-support observed years per retained site. Both universes pass minimum site/country/project support. Retained rows represent conditional recorded family incidence, not abundance, biomass, true richness or whole river community. Global exclusion of 35 source-ambiguous family labels removes substantial taxonomic coverage; excluded major families may influence results. Do not generalize to all invertebrates. Known sampling-count/duration inconsistencies in some source metadata not repaired; observed annual keys used instead. Eleven of twelve projects have >=2 sites and are informative for within-project variation; uneven country support. Prior UK/France excluded; physical overlap with other European monitoring panels not exhaustively ruled out.

## Primary comparison, equal project and equal site within project
Family binary incidence conditional on clean retained positive records, normalized to sum one per observed year; averaged across matching years in each window. Orthogonal family-coordinate projection replaces coordinates with within-order mean (between-order component), remainder is within-order family allocation deviation. Same family universe and site-year support for both, not comparing different indices. Half reproducibility = 2 times weighted cross-half spatial covariance divided by sum of early/late spatial variances, not ICC, not ordinary Pearson correlation. Remove each project's spatial mean separately at each window. Positive exponential project multipliers, 2000, seed 20261125; nominal descriptive 95% intervals.

| Universe | Total | Between order | Within order | Between minus within |
|---|---|---|---|---|
| 546 insects, within project | .581801 [.516796,.642148] | .633956 [.550024,.711094] | .574540 [.510974,.633921] | +.059415 [.021964,.108075] |
| 474 EPT, within project | .656643 [.521861,.750618] | .622023 [.444376,.757180] | .659416 [.527909,.750773] | -.037393 [-.105213,.035118] |
| insects, all space | .721340 [.661679,.768359] | .815697 [.738143,.869804] | .703292 [.647781,.746129] | +.112405 [.077134,.147321] |
| EPT, all space | .761903 [.671865,.818144] | .822254 [.698350,.891804] | .753589 [.666534,.808942] | +.068665 [.012987,.115395] |

Within-project leave-one-informative-project ranges for fine component insects .554738–.590920, EPT .605290–.712285. Paired coarse-fine contrast insects .045123–.073607, EPT -.060092 to -.015845. All 2000 draws estimable. Mean cross-window covariance divided by mean window spatial variance remains positive: within-project insects B .517208 [.450210,.597364], W .440780 [.375230,.500005]; EPT B .562490 [.432362,.674178], W .565334 [.458594,.647607]. Orthogonal projection assigns .883787 [.865108,.904361] and .927191 [.914143,.938392] of within-project family-coordinate spatial variance to W for insects and EPT, respectively. This partition depends on the Euclidean family weighting and number of coordinates, NOT a universal share of biological diversity or causal processes.

Primary finer-level persistence supported for both universes. Coarse-over-fine persistence criterion passes for ten insect orders but not within-project EPT sensitivity. Do not claim universal taxonomic ordering. Sustained fine spatial correspondence in repeated field monitoring is consistent with, but does not independently replicate, prior species-specific host-graph correspondence in butterflies.

## Secondary matched environment comparison (state and change)
Four source-native upstream urban/crop cover and prior-12-month climate variables, exact same selected annual keys. Train first two windows, leave test country entirely out; project-center heldout late response and predictors, so this is heldout centered correspondence, NOT deployable absolute-state or prospective change forecast. Seven eligible test countries for insect, six for EPT. Fixed Ridge10, equal project weights, train-only scaling, fixed-country loss bootstrap 2000. MSE gains relative to zero centered predictor; point estimate [nominal fixed-loss 95%], percentage points for the paired state-minus-change gain contrast.

| Universe/component | State gain % | Change gain % | State minus change pp |
|---|---|---|---|
| insects total | -0.986 [-3.764,.198] | -10.358 [-34.989,-.033] | +9.372 [-2.534,34.194] |
| insects B | +1.014 [-.897,3.200] | -6.789 [-20.745,.402] | +7.803 [.286,20.913] |
| insects W | -1.267 [-4.499,.028] | -10.733 [-36.775,-.061] | +9.466 [-3.156,35.723] |
| EPT total | -2.072 [-9.103,.633] | -32.958 [-99.479,-.260] | +30.886 [-7.154,97.181] |
| EPT B | -.424 [-7.961,3.711] | -45.328 [-122.514,-.471] | +44.904 [-3.894,121.292] |
| EPT W | -2.213 [-9.754,.499] | -31.975 [-97.657,-.218] | +29.762 [-7.885,96.045] |

Neither state nor change shows a robust positive increment for any component. Only secondary insect B paired gain contrast is positive, but positive state prediction itself not established, and no cross-universe replication. No model tuning or outcome-selected predictor additions. This is not absence of environmental filtering, universal temporal decoupling or evidence for measured nutrient effects. These are four non-exhaustive environmental summaries; no measured sewage or nutrient dose. Persistent geographic composition can coexist with limited transfer of these particular environment models. Strongest manuscript role is external data-type sensitivity of persistent fine-scale community composition, not a wastewater-specific causal mechanism or extension of GBIF forecast accuracy.

## Independent local numerical QA
Both signed ZIP digests verified. Recomputed all 72 persistence summary rows (point, multiplier intervals and deletion ranges) using normalized indicator projection and independently formed half-mean Gram matrices; direct weighted vector covariance point checks. Independently refit all 78 country/mode/component losses with augmented weighted least squares instead of original normal equations; regenerated all 18 environmental summary rows. Maximum absolute discrepancy 3.56e-13. No independent raw abundance reconstruction locally; original-source computation in Actions. Full derived matrices and remote code in PR165 artifact/repo, independent local audit code and CSVs in conversation ZIP. No original data reposted in Git. Main unchanged.