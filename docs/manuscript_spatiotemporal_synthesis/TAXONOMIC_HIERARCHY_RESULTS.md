# Nested taxonomic identity audit (2026-11-25)

Protocol fixed at 39c2064fc217e052e2705d0a0b703e7180fa7a8f before this calculation, after the preceding outcomes were opened. Actual local computation on the same 20 hash-verified source ZIPs and one original flora taxonomy artifact (inner crosswalk SHA256 ab45f36c0dad047442bb958bcef3135b45014bafec8a70811bcafb6109a02679, 11,511 original species keys). 1,176 candidate plant keys, 472 source-defined genera; 990 old names from verified GBIF match genus, 186 extra accepted names from the original exact canonical crosswalk. All 520 overlapping accepted old species genus names agree. Common 46 accepted butterfly keys, 32 source-defined genera. No live taxonomy rerouting or inferred phylogeny. Original 263/519 edges preserved. Reproduced all eight earlier cohort/graph/scope paired actual-minus-family-null median contrasts to 1e-12. No previous model refitting.

## Plant taxonomic level
For each H6 separately construct two exact nested expected host adjacency matrices, preserving the number of host links in family × early observed H7-window incidence strata (F), or additionally plant genus (G). Actual A, null F and G use identical early mean plant inclusion; Y is later mean butterfly inclusion. Exact nested strata and row-degree preservation verified. H6 equal, H7 equal within parent, same fixed estimable species; 2,000 parent bootstrap draws. Primary comparison within H6 is the species median of paired Pearson r differences (not the difference of marginal medians).

| Cohort | Catalogue | Actual r | Family null r | Genus null r | Genus minus family [95%] | Actual minus genus [95%] |
|---|---|---:|---:|---:|---|---|
| quality 103 H7/53 H6, 45 taxa | Ferrer-Paris | .1530 | .0982 | .1221 | .01621 [.00000,.04939] | .01673 [.00000150,.03197] |
| same | Clarke | .2241 | .1426 | .2029 | .02795 [.00000,.06362] | .00000 [-.00013,.00159] |
| all 176 H7/75 H6, 46 taxa | Ferrer-Paris | .1487 | .0764 | .1176 | .01752 [-.00251,.04131] | .00547 [-.00002,.01613] |
| same | Clarke | .2106 | .1272 | .1983 | .03489 [.00640,.06592] | .00049 [-.000001,.00165] |

Neither component meets the pre-fixed both-catalogue positive lower-bound rule in the primary quality within-parent panel. The all-spatial sensitivity supports genus-minus-family in both quality and inclusive cohorts, but is not the within-parent primary. On a different, additive common-denominator covariance decomposition (not Pearson increment), the between-plant-genus component in quality is positive in both catalogues: .04276 [.00583,.07008] and .04937 [.01260,.07970]. Within-plant-genus covariance .04727 [.00724,.06183] old and .00375 [.00000,.00893] Clarke. This supports a descriptive contribution of genus-level plant structure, not a causal taxonomic filter or replacement of the primary decision. Component medians need not add to the median total; equality checked at each taxon and bootstrap draw.

Exchangeability materially shrinks from family null to genus null: quality movable available links old 87.04% -> 54.69%; Clarke 81.26% -> 36.47%; all old 85.13% ->51.08%, Clarke 79.42% ->34.28%. Nonmoving links and sparse information can limit discrimination; no absence-of-species-level-effect claim. Fixed-taxon resamples with every primary r estimable: 1,175/2,000 quality, 784/2,000 all within H6. All 2,000 produce a finite median among estimable taxa; conditional intervals do not include support uncertainty.

## Butterfly taxonomic level (secondary)
Use only 23 species in nine multi-species genera, excluding singletons, standardized separately for each species and variable. Subtract H6 means; decompose average matched-species correlation exactly into covariance of the congeneric mean standardized fields plus covariance of the congeneric residual fields. Family-excess index is actual host index minus expected family × early-incidence index, then standardized; not the same as the host index Pearson increment. H6 equal within species and species equal; denominator for within-genus residual correlation is pooled residual SD. Independent direct weighted-array decomposition reconciled to 1e-12.

| Cohort | Catalogue | Predictor | Total mean r | Between butterfly genera | Within butterfly genera [95%] | Within-genus residual r [95%] |
|---|---|---|---:|---:|---|---|
| quality 103/53 | Ferrer-Paris | actual | .2647 | .1596 | .1051 [.0690,.1359] | .2135 [.1430,.2797] |
| quality | Clarke | actual | .2545 | .1758 | .0787 [.0390,.1112] | .1815 [.0959,.2500] |
| quality | Ferrer-Paris | family excess | .2234 | .1324 | .0910 [.0409,.1221] | .1697 [.0777,.2313] |
| quality | Clarke | family excess | .1755 | .1302 | .0452 [.0090,.0802] | .0945 [.0195,.1647] |
| all 176/75 | Ferrer-Paris | actual | .2145 | .1358 | .0787 [.0451,.1136] | .1566 [.0893,.2218] |
| all | Clarke | actual | .2322 | .1556 | .0767 [.0301,.1130] | .1700 [.0677,.2431] |
| all | Ferrer-Paris | family excess | .1718 | .0996 | .0722 [.0390,.1095] | .1373 [.0743,.2049] |
| all | Clarke | family excess | .1878 | .1265 | .0613 [.0206,.0975] | .1225 [.0417,.1882] |

These descriptive associations among congeneric species persist for both catalogues and both cohorts, even after subtracting the regional family-null resource index. Only 17 quality and 24 inclusive H6 contain multiple H7; only 1,228/2,000 quality and 1,259/2,000 inclusive within-parent bootstrap draws have all required variances nonzero; intervals are explicitly conditional on those valid draws. No multiple-comparison correction or prospective replication. The association is not independently adjusted for earlier butterfly composition, and prior regional-pool early-butterfly partial association was uncertain in both quality catalogues. No demonstration of demographic response, within-genus adaptation, phylogenetic causality, or a dynamic host effect.

## Interpretation
The previous PR105 genus-null result concerns the plant NUTRI–butterfly HNS trait association, not this species-specific actual host-composition index and not a prohibition on consumer species-level spatial structure. PR134 W–plant NUTRI within-plant-genus score sorting is also a different endpoint and predictor. This new result rules out neither those associations nor their limitations. Do not claim a universally genus-only consumer response; it is possible to retain species-specific geographic host correspondence within congeneric butterflies even when within-genus plant nutrient-affinity sorting does not transmit to the butterfly HNS endpoint. Earlier EQR–network cross-catalogue nonreplication, dynamic change/coherence failures and population confounding remain. This is supplemental descriptive taxonomic sensitivity, not a new identified sewage pathway or causal selection mechanism.

## Reproduction
368 species/graph/cohort/scope rows; 72 plant summary rows; 256 regional support rows; 64 butterfly genus summaries and 368 genus-membership rows. 20 source archive hashes verified; additional original crosswalk inner SHA256 verified, outer ZIP hash recorded in receipt. Local full derived archive gesui_taxonomic_hierarchy_20261125.zip SHA256 258d9e3325f00e1d65ac67a0fd71b423a3d3879c89738a620d4d8ff8baad612b. audit.py SHA256 4b5c96e7a2f74a47df12d1d4bb3dc647c82c509a9140375457270c26762e487b; summary.csv 2c9c2f41449e1e3264ee2e7e479a76ce24bbd2fcaca48128a9955166a33123bc; butterfly_genus_summary.csv 3baa0d8c97af561b8b336599951920371ce3cb9d2f4d6f5b056f692a4d05a61e. Full code/tables in conversation ZIP, report and protocol committed. No raw archive redistribution, no remote numerical replay, no main merge.
