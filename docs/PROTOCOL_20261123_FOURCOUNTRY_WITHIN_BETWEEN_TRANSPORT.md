# External 4-country EPT × fine sediment: within-country signal vs country-mean ecological composition
Date 2026-11-23 JST | **POSTHOC numerical audit**, chosen after the original four positive heldout MSE gains in PR #146. No significance claims; this is not preregistered confirmation.

## Question and independent source
McKenzie et al. *Global Change Biology* 30:e17084 (2024) [doi:10.1111/gcb.17084](https://doi.org/10.1111/gcb.17084) **already established** region-dependent macroinvertebrate responses and the relevance of EPT indicators to fine sediment. The genuine Loughborough original workbook doi:10.17028/rd.lboro.24615570.v1 (CC BY-NC 4.0; SHA256 `a02573bef0e11ea19c9465a78ebecc130fb9c39d76d35a79863c2efd5debcd48`) has 6,491 `ANALYSIS_ID`-paired original visual fine sediment cover observations and *relative taxon abundances*, not raw counted individuals or a physical site identifier. Original E/P/T taxon assignments were fixed by an external taxonomy-source gate before opening the ecological relationship (original 302 header names, upper bound for ambiguous `Insecta` source label). Existing [PR #146](https://github.com/zuizui0223/gesui/pull/146) held out each country entirely and observed tiny source-conditioned MSE improvements AU+1.814%, BR+2.396%, NZ+1.448%, UK+1.785% using fixed Ridge α=10 with one train country one vote.

**Question:** Does this small cross-country positive prediction mostly arise from *within-country* ecological co-variation of relative EPT and fines, or from comparing countries with different mean sediment cover and different mean EPT relative abundance?

## Exact source and prior numerical replay gate
- Source must be original publisher Figshare XLSX 4,640,542 bytes; MD5 `1e4d6e37b73c663406e5567dd79d9452`; SHA256 as above. 4 named country sheets, original row totals AU3099 BR232 NZ342 UK2818, original ANALYSIS_ID one-to-one, no original physical SITE_ID, 302 original taxon labels.
- Import unchanged original source-gated EPT-to-order relative proportions and original `fine/100` as in `scripts/gesui_external_fourcountry_fixed_country_holdout_EPT_20261122.py`. Download the independently signed source-only original UNKNOWN-mass gate and hierarchy-corrected source taxonomy, and **the original signed PR146 four-country heldout numerical artifact** from Actions #45032490453. Recompute and verify all **four original full-gradient gains, baseline MSE, full MSE and train slope**, on original lower and upper taxonomic bounds, at numeric tolerance 1e-10 before performing anything new. STOP if not reproduced. Never impute unknown taxon mass or site IDs.
- All country-heldout train sets comprise exactly the *other three countries*; country-equal train weight (each country's individual sample weights sum to 1), original sediment fraction in [0,1], response original confirmed EPT relative abundance.

## Freeze decomposition and fixed falsifiers
For a target held-out country T, train-country means of sediment `s` and relative EPT `y` are `μs, μy`, with each training country weighted 1/3. Exact training source covariance identity:
`Cov_bal(s,y) = mean_j Cov_j(s,y) + Cov_j(μs_j, μy_j)`,
likewise `Var_bal(s) = mean_j Var_j(s) + Var_j(μs_j)`. Report the within-country and between-country components, their signs and each source country's raw population mean/covariance. Use population moments (not independence claims).

**Fixed models at original α=10:**
- M0 = original train-country-equal grand mean of `y` (identical to signed PR146 baseline).
- M_pool = original alpha-10 ridge penalized slope from *total pooled* country-equal covariance, unpenalized grand mean intercept; must numerically reproduce signed original PR146.
- M_within = ridge penalized slope from **within-country demeaned covariance/variance only**, same original train-country-equal grand mean intercept. Train-only means and centering, no heldout target labels/means used for fitting.
- For M_pool and M_within report 4 heldout country MSE changes vs M0. Never cherry-pick model variants or write within-only as the previously published M_pool result.

**Exact gain identity**, for either trained β:
`MSE0 - MSE1 = 2β Cov_T(s,y) + 2β(μs_T-μs)(μy_T-μy) - β²[Var_T(s)+(μs_T-μs)^2]`.
Report each component in **heldout baseline-MSE percentage points** and independently assert component sum reproduces the original gain to tolerance 1e-10. First component = within *target* association; second = country-mean mismatch correction; last = nonnegative predictor penalty. `μy_T` appears **only in after-the-fact error attribution**, NOT model training or features. This decomposition identifies a mathematical driver of MSE gains, not a causal ecological mechanism.

**Sensitivity grid fixed now, not tuned:** alpha `[0, 1, 10, 100]` for both M_pool and M_within, same country-equal training, in both original entire fine 0–100% source domain and original 0–50% common support domain. Additionally one visit-equal rather than country-equal training at original alpha10 as a transparent weighting sensitivity; this is not a new target model selection procedure. All countries with valid source support reported. For taxonomy upper bound `y_upper=y+unresolved Insecta mass`, repeat alpha10 decomposition and 4 heldout MSE gains to confirm the original source-uncertainty bound.

**Predeclared descriptive decision rules, not statistical tests:** a sign mismatch between train-within slope and target within-country covariance in any country challenges an *every-country* transferable response direction; if in >=3/4 countries the absolute mean-offset term exceeds the absolute within-target term, the original country-holdout gains cannot be presented as direct evidence of a mainly within-country ecological signal. If M_within α10 fails to improve heldout MSE in >=2/4, avoid claiming robust within-country predictive transfer. No p-values or CIs without true physical site IDs.

## Boundaries
These 6,491 ANALYSIS_IDs are not verified independent river sites and represent pre-normalized relative assemblage counts, not abundance or biomasses. Country-level mean differences may reflect taxonomic sampling, environment, or source domain shifts; decomposition cannot identify the cause. Four countries are only four independent country domains. This is **not** an independent replication of EPA HUC8 original-source-membership × XEMBED interaction; deposited <2mm cover and EPA embeddedness are different environmental constructs. No wastewater effluent dose, land transfer or ecosystem services measured. Maintain sewage-exposure STOP in [PR #139](https://github.com/zuizui0223/gesui/pull/139).
