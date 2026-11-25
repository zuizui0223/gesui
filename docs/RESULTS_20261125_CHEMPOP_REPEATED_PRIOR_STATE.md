# Repeated English river monitoring: prior-site state, nutrients and modeled effluent — 2026-11-25

Protocol fixed at dc7453e3ae03e49a272e559281afeab90e4d7df6 before this analysis. Workflow [45546744470](https://github.com/zuizui0223/gesui/actions/runs/45546744470) on 9ab5bca31659be138b9d0b8cbc85e04779bfd956 SUCCESS, signed artifact 14991936109 SHA256 350182684a24a9ee095fcf5020a2ffa755da4d8349c7ea635ed9468db21f265a. Original 35 source-member hashes, published score dictionary and exact 5000 analyses / 2624 site-years /591 sites/7 regions reconstructed against earlier signed receipts. Exact original 85 observed family scores and order grouping preserved.

## Fixed repeated-site source gate
Four windows 2011–12, 2013–14, 2015–16, 2017–18. Require at least one eligible observed year in every window, static EDF, and unambiguous native site identifier. 139 candidates in seven regions; exclude the eight southern-region sites under minimum ten sites per region. Main passes: **131 sites, six regions, 944 site-years, 1866 original biological analyses**, window site-years 249/247/252/196. Exact roster SHA256 5dd31f4839292160872bc40f3aaa55fdcd496924463deac803d3514261ae0fd3. Stable recorded-abundance-method cohort fails its preregistered source support: 91 candidates, 83 sites/four regions after region filter, 600 site-years/1188 original biological analyses. No outcome fit for this cohort. Recording method includes log abundance estimates, not all actual counts. No claim of gear/effort equivalence.

## Historical-only prediction results
Train 2011–12 -> 2013–14 on other regions, predict 2015–16 -> 2017–18 in entirely held-out region. P=other-site regional score mean, PS=historically learned convex mix of P and own-site score. B=source background/recording controls, E=static log1p EDF, C=three measured nutrient-associated water chemistry window summaries at predictor time. Fixed Ridge10 residual correction, equal regions and sites, identical total score weight and linear design for additive order-between/within-order score components. Target-window outcomes or chemistry not used in predictor fit. Target-region prior observations are available predictors. Learned own-site weights .77855–.89095 (six folds), not a causal memory coefficient. Ordinary paired region bootstrap 2000, fixed forecast losses, descriptive 95% intervals. Only six geographic units in one country. BMWP equal-family score components, not family-vector composition, abundance, population change or material flux.

| Contrast | Total | Between order | Within order |
|---|---|---|---|
| PS vs P | +71.465 [64.747,75.395] | +72.793 [64.916,77.092] | +41.436 [7.507,56.797] |
| PSBCE vs PBCE | +51.180 [32.753,60.315] | +49.241 [35.499,59.513] | +41.982 [5.761,57.003] |
| PBCE vs PBE (chemistry after EDF, without prior own state) | +12.448 [7.206,18.221] | +13.488 [5.529,20.376] | +1.908 [-.044,5.709] |
| PSBCE vs PSBE (chemistry after EDF and prior own state) | +.568 [-.026,1.167] | +.060 [-1.594,1.865] | +.160 [-1.965,2.259] |
| PSBCE vs PSBC (EDF after chemistry and prior own state) | -.387 [-.979,.054] | -.220 [-.629,.039] | +.638 [-.131,1.801] |

All values relative MSE improvement percentages on same response/support, not percent changes in biology. Primary within-order EDF contrast FAILS fixed positive lower-bound criterion. Improves five of six regions, leave-one-region range +.114 to +.948%; this does not rescue primary. Prior state after full measured covariates has positive nominal intervals for all three components; 5/6, 6/6, 4/6 regions respectively. PS vs P positive all six regions for total and between, five for within. No new tolerance-label search; PR163 p=.3548 remains. No stable-method outcome sensitivity because source gate failed.

## Post-result same-denominator arithmetic
Using stored regional losses, compare the absolute MSE reduction attributable to adding chemistry before and after adding own-site state, divide BOTH reductions by PBE MSE, and apply the same paired 2000 region bootstrap. This removes differing relative-gain denominators but is an exploratory arithmetic contrast, NOT causal mediation or independent test.
| Response | Chemistry increment before own state, pp | After own state, pp | Before-minus-after, pp [95%] |
|---|---:|---:|---|
| total | 12.4482 | .2441 | 12.2041 [7.0872,17.7311] |
| between order | 13.4884 | .0265 | 13.4619 [6.2769,19.7420] |
| within order | 1.9076 | .0913 | 1.8162 [-.8718,5.3820] |

For EDF after chemistry, with common PBC MSE denominator, within-order increment changes -.1839 to +.3731 pp; difference before-minus-after -.5570 [-.8240,-.1435]. This descriptive contrast does not establish a positive conditional EDF increment, and is not a substitute for the failed primary criterion or the earlier score-specificity null. Total and between-order EDF contrasts uncertain. Full 18 common-denominator rows preserved.

## Interpretation and limitations
Past own-site condition retains substantial predictive information after the particular measured nutrient summaries and static modeled effluent fraction. Chemistry's predictive increment for total and between-order scores is substantially smaller when earlier local scores are available. Prior score can summarize stable habitat, legacy environmental effects, demographic processes and monitoring practice, and may lie downstream of actual exposure. No identified causal ecological memory or proof of mediation, no absence-of-effect claim, no wastewater dose. Compared with TREAM own-site composition forecasts (PR168/169), different response and sampling grain: do not pool gains or compare biological effect sizes. This repeated-panel subset represents 131/591 original sites, not the whole ChemPop network. Earlier spatial results remain valid within their original estimands, not independent confirmation; stable-method sensitivity was unavailable.

## Local QA
Rehashed all three signed result archives, verified 35-member receipts and original score/order identities, independently reconstructed main and stable-method rosters and 3930 site-loss rows to 18 region/response summaries. All 45 model contrasts, 2000 paired bootstrap intervals and deletion ranges reproduced (maximum absolute numerical discrepancy 1.14e-13). New common-denominator table and script in audit archive. Original source-level model execution in Actions; local QA is saved-loss arithmetic, not raw-source refit. Local archive gesui_chempop_prior_audit_20261125.zip SHA256 33142fceeac0e807b4eae29cae713aa601fb8950d3c1af18803ad91db30b398b; audit.py b6bf46ef7f364bc02f7e4af4d09819260eb3b6cf1f09fd1134b26b38b50e5df4; common_denominator.csv ab66b1a18573b575fde489b0ae78afced74b242e3194d7ace9bb5dc90a2eedbb. Source/model code and report on draft PR170, no main merge.