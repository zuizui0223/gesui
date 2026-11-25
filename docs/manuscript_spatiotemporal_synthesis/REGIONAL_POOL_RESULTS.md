# Regional candidate-pool-conditioned host identity diagnostic (2026-11-25)

Protocol frozen before this calculation at 11755c99e61ca2197729787b2a4ff454de6236b6, after prior analyses had opened the same outcomes. Actual local computation: 20 SHA256-verified input ZIPs, exact common 46 Lycaenidae and 1176 candidate plant keys (990 original +186 extra), 263 Ferrer-Paris and 519 Clarke edges. Same common PR160 roster. For each cohort and H6 separately, stratify plant labels by family and exact number of occupied H7-window cells in 2014–16 and 2017–19; exact expected adjacency under within-stratum exchange preserves each butterfly's degree and stratum counts (checked for all). Both real and null degree-normalized host indices use only early plant occurrence inclusion. Response: 2020–22/2023–25 mean rarefied butterfly inclusion. No later plant data enters the null, no later butterflies enter its construction. This is a regional **observed candidate-pool surrogate**, not a complete independent regional species inventory, and not a test of matched physiological niches.

## Primary result: within-parent association with late butterfly inclusion
H6-equal and H7-equal within H6, demean variables within each parent. Conditional 2000 H6-block bootstrap seed 20261125, fixed estimable taxa. Paired contrast median across taxa (not subtraction of marginal medians). Both quality graphs pass the frozen positive lower-bound rule.

| Cohort | Graph | Taxa | Actual r [95%] | Expected null r [95%] | Paired difference [95%] |
|---|---|---:|---|---|---|
| quality 103 H7/53 H6 | Ferrer-Paris | 45 | .1530 [.1122,.2317] | .0982 [.0478,.1520] | .0449 [.0082,.0889] |
| quality 103/53 | Clarke | 45 | .2241 [.1394,.3216] | .1426 [.1007,.2224] | .0674 [.0126,.0900] |
| all 176/75 | Ferrer-Paris | 46 | .1487 [.0822,.1945] | .0764 [.0410,.1338] | .0485 [.0022,.0702] |
| all 176/75 | Clarke | 46 | .2106 [.1421,.2778] | .1272 [.0722,.1959] | .0538 [.0218,.0779] |

Only 17/53 and 24/75 H6 have multiple H7 and contribute nonzero within-parent information; singleton H6 retained as zeros. Quality all-space paired differences old .0440 [.0099,.0726], Clarke .0309 [.0079,.0609]; inclusive .0231 [.0133,.0695] and .0297 [.0138,.0557]. Local regional null mean pool size 328.47 candidate plants per H6 quality and 311.91 inclusive; 87.04%/81.26% of early available links actually movable within strata in quality, 85.13%/79.42% inclusive (old/Clarke). Null not degenerate.

## Critical sensitivity: condition on early butterfly composition
Host excess = actual−expected null index. Within-H6 partial correlations with late butterfly inclusion, residualizing against the null index, then also early butterfly inclusion. These are same-panel descriptive partial associations, not causal mediation or predictive model validation.

| Cohort | Graph | Partial r, null index only [95%] | Partial r, also early butterfly [95%] |
|---|---|---|---|
| quality | Ferrer-Paris | .1045 [.0736,.1971] | .0244 [-.0271,.0903] |
| quality | Clarke | .1496 [.0954,.2698] | .0490 [-.0206,.1242] |
| all | Ferrer-Paris | .1014 [.0344,.1414] | .0178 [-.0135,.0715] |
| all | Clarke | .1814 [.1039,.2063] | .0611 [.0231,.1223] |

Quality within-H6 partial estimates use 45 old and 44 Clarke taxa; inclusive 46/45. After early-butterfly adjustment neither quality graph shows a clear association. All-space Clarke adjusted r .0324 [.0049,.0928] quality, .0412 [.0066,.0895] inclusive; old uncertain. Do not rescue the joint result by selective choice of catalogue or cohort.

## Interpretation and manuscript role
Actual catalogue host identities retain a modest positive spatial correspondence with later butterfly composition beyond the specified local family/availability-matched null, even within H6. This strengthens the existing descriptive persistent host-resource correspondence and is not merely coarse regional pool overlap under this particular null. The additional association attenuates after conditioning on prior butterfly state. It does not demonstrate independent dynamic driving, a causal environmental filter, realized feeding, biomass, nutrient transfer, or wastewater effects. Shared fine-scale habitats, species pools and observation practices remain possible. The adjusted result does not establish mediation or prove absence of a causal host effect. Distinct catalogue graphs share the same observations; no independent field replication or formal permutation significance test. Previous four-window change/coherence failures and cross-catalogue river-EQR nonreplication stand. Include as supplementary host-specificity sensitivity alongside the main manuscript temporal result.

## Reproduction
All covariance results reconciled to directly centered weighted arrays, partial correlations independently checked by weighted least squares. 368 taxon rows, 48 summaries, 4 support summaries, 256 regional graph support rows. No original forecasts refitted. Local full archive gesui_regional_pool_diagnostic_20261125.zip SHA256 74063e78029910b27e0bd1452493997e9192e2a60c7e48be9b34d55c2118a824; audit.py 390a8b6e452f1286256dff079746e85c35e67e80733a30aa524770ff7555dad0; summary.csv 04f762e3c10504a5ddc8747731cb8ef18d01f816f4b789db99efabb43f9564d3; species_diagnostic.csv db5e961efb23577fce271fc14560c129207057b314641cded9126a8ce02a196d. Full code and numerical tables in conversation ZIP, report and protocol in Git. No remote replay, no main merge.