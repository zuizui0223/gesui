# Nonoverlapping time-window spatial covariance and change diagnostic (2026-11-25)

Protocol fixed at 4fb255120b52ee5ede9dabe35e1f3dad53a226cc before this calculation, after all preceding outcomes were inspected. Actual local run on 20 hash-verified ZIPs; same PR160 quality and inclusive common-46 rosters. Previous source reconstruction rechecked against the sealed 2020–22 butterfly levels and original France263 site means. Computed 1,294 diagnostic rows and 134 summary metric rows. For every covariance calculation, original H6-equal sufficient-moment result independently reconciled to directly centered weighted arrays (absolute tolerance 2e-13, relative 1e-10). No new field observations, no prediction fitting or threshold tuning.

Four windows split into disjoint early 2014–19 and late 2020–25 means. For each fixed species, cross-half covariance divided by average within-window spatial variance is a signed descriptive persistence fraction, not an ICC. E–L Pearson spatial correlation separately calculated. H6 equal, H7 equal within H6; 2,000 H6-block bootstrap draws seed 20261125 with fixed estimable taxa. Medians over species, not independent species replication.

| Cohort | Index | Cross-half fraction [nominal 95%] | E–L spatial r [nominal 95%] |
|---|---|---|---|
| 103 H7/53 H6 quality | Butterfly | .5938 [.5156,.6767] | .7771 [.7274,.8443] |
| same | Ferrer-Paris host | .5418 [.4641,.6294] | .7299 [.6686,.7917] |
| same | Clarke host | .5593 [.4757,.6296] | .7188 [.6587,.7837] |
| 176 H7/75 H6 all | Butterfly | .5824 [.5145,.6588] | .7736 [.7227,.8382] |
| same | Ferrer-Paris host | .5773 [.4862,.6231] | .7218 [.6775,.7796] |
| same | Clarke host | .5860 [.5164,.6364] | .7451 [.6940,.7866] |

## Association across disjoint time windows versus association of changes
Symmetric cross-half association = (Cov(X_E,Y_L)+Cov(X_L,Y_E))/2 divided by sqrt(mean spatial variance of X half means × mean spatial variance of Y half means). Half-change association is Pearson correlation of L-E differences. These coefficients have different normalizations, so contrasts are descriptive comparisons of standardized statistics, not tests of equal causal effects. The median of paired taxon contrasts is not the difference of two marginal medians. No inference of direct interaction or temporal transmission.

| Cohort | Pair | Cross-half association [95%] | Half-change r [95%] | Paired contrast median [95%] |
|---|---|---|---|---|
| quality | old host/butterfly | .1997 [.1126,.2443] | .0510 [-.0005,.0943] | .1246 [.0372,.1848] |
| quality | Clarke host/butterfly | .2367 [.1811,.3284] | .0050 [-.0304,.0534] | .2000 [.1277,.3106] |
| all | old host/butterfly | .2015 [.1039,.2346] | .0784 [.0216,.0997] | .0946 [.0244,.1545] |
| all | Clarke host/butterfly | .2452 [.1822,.3027] | .0373 [-.0017,.0761] | .1576 [.1037,.2562] |

Within-H6 geographic demeaning (quality 17 informative multi-H7 parents, 45 estimable taxa): old cross .1490 [.1187,.2446], change .0128 [-.0485,.0701], paired .1788 [.0776,.2700]; Clarke cross .2563 [.1923,.3298], change -.0199 [-.1184,.0449], paired .2573 [.1402,.3882]. All sample 24 informative parents: old paired .0385 [-.0160,.1443], Clarke .1399 [.0501,.2531]; old inclusive within-parent contrast uncertain. Not an independently validated fixed-effect model and does not adjust finer geography, climate, observation practice or host degree. Singleton parents retained as zeros.

## Same France river/catchment frame (two source periods only)
61 H7/23 H6 exact inclusive intersection, 42 estimable butterfly taxa, 2018 and 2022 river observation years and matching 2017–19/2020–22 terrestrial windows.

| Pair | Cross-period association [95%] | Change r [95%] | Paired contrast [95%] |
|---|---|---|---|
| measured river logTP / benthic EQR | -.6861 [-.7739,-.5835] | -.1338 [-.3622,.0667] | -.5524 [-.7739,-.3071] |
| old host / butterfly | .1703 [.0791,.2699] | .0250 [-.0267,.0669] | .1635 [.0163,.2511] |
| Clarke host / butterfly | .2418 [.1429,.3419] | .0000 [-.0616,.0484] | .2047 [.0779,.3234] |

Full original FR73/24 H6: logTP/EQR cross -.6817 [-.7722,-.5709], change -.0868 [-.3003,.1457]. Within-H6 France61 (13 informative parents) old/Clarke paired contrasts .1058 [.0055,.2249] and .1846 [.0271,.2815]; TP/EQR -.5553 [-.9769,-.2012]. Strict France33/18 H6 intersection is smaller: old host/butterfly paired contrast .1364 [-.0051,.2619] uncertain, Clarke .2026 [.0548,.3355]; TP/EQR -.6477 [-1.0151,-.1749]. No comparison of effect sizes across compartments or equal-window two-versus-four biological stability.

## Interpretation
Persistent spatial correspondence survives temporally disjoint observations, not merely averaging the same windows. Weak or uncertain change correspondence contrasts with those spatial associations; in inclusive old catalogue the half-change correlation itself is small positive, so avoid declaring no temporal association. These are source-panel retrospective descriptions, conditional H6 resampling, not causal tests, independent prospective replication, measured biomass, wastewater attribution, or equivalence tests. Existing stronger controlled forecast and EQR-network nonreplication stand. New evidence complements manuscript Fig 2/4 by providing matched descriptive statistics, not proving the proposed ecological mechanism. Exact river station and whole-catchment GBIF observations are different grains.

## Archive
Local ZIP gesui_disjoint_time_diagnostic_20261125.zip SHA256 6cb0d8a6657cca3aea8a596bfdda00f3c30ba37e170350689ec6e9d501c2ceca. audit.py 297debe14a2890d83270eb3eeb6381d091129ff3dd9cfd6630417cabcefa211a; summary.csv f620654af7833111587b64b107f50b237a67153b9b341d280dab2a876ce37a7e; species_diagnostic.csv 74a4bb9337703adcd5a238c3f6e631337197ce14a697e98b9966b5a494e41d24. Code and full tables in conversation ZIP, report and protocol committed; no remote numerical replay, no main merge.
