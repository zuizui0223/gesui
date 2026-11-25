# Matched chemistry/order resolution — verified outcome (2026-11-25)

Actual run [45530295378](https://github.com/zuizui0223/gesui/actions/runs/45530295378), job 138923627828, succeeded on original code 8ddee67333ceefdbbb16c309bf14291634a44238. All 35 original data-member SHA256 values checked against previous signed source receipts, original BMWP score reference checked, same 5000 analyses/2624 site-years/591 sites/7 regions. 85 observed scored families, 79 exchangeable in 12 nonsingleton order groups, Phryganeidae and Taeniopterygidae source-ambiguous and retained as singletons. Original total WLS and four regional prediction estimates reproduced by workflow. Local audit of signed output ZIP 14985171063 SHA256 9c04699d4fab6ba27b3e5ba99e6c17670b0c69bb912584dbe9a697391185edfc, original chemistry 14599583276 and order 14598879538 verified hashes and all receipt identities, reconstructed all 12 gain/interval rows from stored 7-region losses, all 4999 null draws and p-value, and coefficient additivity. No local re-download of raw EIDC archive or independent row-level predictions.

## Primary chemistry hypothesis: within-order NOT supported
| Response | Joint three chemistry variables after EDF gain (%) [fixed-region 95%] | Improved regions |
|---|---|---:|
| total BMWP incidence CWM | 10.558940 [5.062079,14.289563] | 6/7 |
| between-order mean score component | 11.840724 [6.935631,15.981726] | 6/7 |
| within-order family score deviation component | -0.410884 [-1.163157,0.385056] | 2/7 |

Observed within-order partial R2 .00154750 (direct independent solve .00154750), label-null median .00826841, 95% null quantiles [.00077362,.03389785], upper-tail p=.936. Both prespecified conditions fail. Total standardized orthophosphate/nitrate/ammoniacal-N beta -0.066419/-0.111957/-0.112870; between-order -0.057730/-0.116808/-0.113360; within-order -0.008690/+0.004851/+0.000489, all three within-order nominal intervals include zero. Preserve original chemistry whole-community positive, no absence-of-effect claim.

## Secondary finding: static EDF conditional on three chemistry summaries
| Response | EDF after nutrients gain (%) [fixed-region 95%] | Improved regions |
|---|---|---:|
| total | -0.518565 [-0.761173,-0.253323] | 0/7 |
| between-order | -0.568302 [-1.264757,0.055442] | 3/7 |
| within-order | +1.718449 [1.283008,2.258632] | 7/7 |

Within-order EDF beta -0.019920 [-0.038826,-0.001014]; total +0.007091 [-0.045447,+0.059630]. EDF after geography only within-order gain +2.252681 [1.295434,3.356617]. All comparisons other than original primary joint chemistry rule are secondary; no multiplicity correction. Leave one region from fixed predictions (no refit): within-order EDF after nutrients +1.528342 to +1.861373%. Seven region gains .9077/1.8041/2.9227/1.9799/2.4142/1.2747/1.0829%.

Posthoc stored-error arithmetic, express each component's absolute mean MSE reduction divided by same total baseline MSE: EDF after nutrients total -0.518565 percentage points = between -0.481682 + within +0.289843 + cross remainder -0.326725. For nutrient block after EDF total +10.558940 = between +10.186936 + within -0.060357 + remainder +0.432360. This is not summing component-specific relative MSE gains. Identical linear Ridge designs imply cross-error covariance remainder, but stored outputs alone do not independently verify row-level additivity. Thus aggregate BMWP score can conceal a predictive association with finer taxonomic score deviations. Do not call this a demonstrated community turnover or a causal sewage effect.

EDF is static modeled effluent fraction, the three chemistry measurements are biological-date window summaries, not wastewater-specific chemical doses or all possible water-quality/physical habitat variables. Residual geography, unmeasured compounds, habitat, method and observation differences remain viable. Prior within-order EDF score-label null without nutrient adjustment p=.206 is a different model. The new within-order EDF conditional result has not yet been tested against an order-preserving score-label null. This is a targeted secondary question, not permission to optimize endpoints until positive. No cross-country generality or actual abundance inferred. No claim that plant genus and benthic order are directly equivalent. Main unchanged.

Local audit ZIP gesui_chempop_resolution_audit_20261125.zip SHA256 cd29292f0aefd1e65233fd860f8a217503464208948857d706ca6ec9fe2da639; includes audit.py, receipt, 12 gain rows, 84 region loss comparisons, 16 decomposition rows and 84 fixed-loss deletion rows. Source code and original workflow already in PR162, derived numerical audit code in local ZIP only.
