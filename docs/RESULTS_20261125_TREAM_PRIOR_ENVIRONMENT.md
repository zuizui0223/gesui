# TREAM prior-site composition and additional environmental information — 2026-11-25

Fixed protocol 2ffad834bfbb09e0e450aface9f014d66428e4e4 before local computation. Exact previous insect 358 sites/7 projects/5 countries and EPT 320/6/4; all five input archive SHA256s verified. Reproduced all previous P and PS predictions, historical weights and stored site losses exactly (max abs discrepancy 0). Four windows, fit 2007–09 -> 2010–12 on all other projects, score 2013–15 -> 2016–18 in heldout project. Only T0/T2 predictor-time environmental fields used, not target-period environment. P=other-site project mean; PS=historically estimated own-site mix with P. C=two climate variables; L=upstream crop and urban cover; CL=joint four variables. Project-centered predictors, weighted RMS training scaling, fixed Ridge10 residual correction; same total fitted prediction projected to order-between and order-within subspaces. Project-equal, site-equal within project. Paired positive project multiplier intervals, 2000 fixed-loss draws, not frequentist confidence guarantees. Posthoc on previously inspected observations; no independent field replication.

| Universe, total composition | Gain % [nominal95] | Improved projects |
|---|---:|---:|
| insect PS vs P | +23.00145 [16.24586,31.93155] | 7/7 |
| EPT PS vs P | +33.70713 [19.49443,46.25836] | 6/6 |
| insect PSCL vs PCL | +23.06582 [16.34656,31.78875] | 7/7 |
| EPT PSCL vs PCL | +33.94292 [19.59939,46.09389] | 6/6 |
| insect PSCL vs PS | -.27906 [-1.03812,.17835] | 4/7 |
| EPT PSCL vs PS | -.08039 [-.57779,.34840] | 4/6 |
| insect PSCL vs PSC | -.28925 [-1.00093,.04484] | 2/7 |
| EPT PSCL vs PSC | +.01526 [-.30528,.23176] | 3/6 |
| insect PCL vs P | -.36297 [-1.83524,.52417] | 5/7 |
| EPT PCL vs P | -.43763 [-2.76791,.94471] | 4/6 |

Both joint primary criteria FAIL. Total prior-site incremental gain remains positive after adding the measured landscape/climate predictors. There was no clear gain from these four environment variables with or without own state; common-denominator after-minus-before environmental increment is not a mediation test. No model tuning. Linear environmental composition scores sum one, but some coordinates are negative: PSCL insect 4984/30788 coordinates (all358 sites); EPT 1243/9280 (315/320 sites), minimum -.038964/-.116589. The pre-fixed Euclidean simplex projection sensitivity for PSCL vs PS +.01362 [-.39768,.29425] and +.17151 [-.02912,.52053]; PSCL vs PSC -.06696 [-.42152,.11108] and +.20081 [.00566,.43281]. Neither meets joint criterion; EPT-only simplex land increment cannot rescue primary. All 132 contrasts retained. No statement that prior composition causally mediates environment, absence of environmental filtering, or memory mechanism. TREAM climate/land cover are not measured nutrients or wastewater dose. Need separate matched chemistry/effluent study for direct extension. EPT is overlapping sensitivity.

QA: 132 summaries, 78 project loss rows, 32544 site loss rows, 78 fold diagnostics, 16 score-validity rows. Independent augmented least squares matches normal-equation predicted corrections to 2.22e-16; exact orthogonal loss additivity. Original source global STOP preserved, scoped acceptance passed. Local full archive gesui_tream_prior_environment_20261125.zip SHA256 586dd83e592c9fb3b974343c9e7c3d686a4987ec5a8f83570098d85bcaac0df9; audit.py 782a5b5a4004a133bdba7166a23c78a591f3ebb2d70fe1e1684c38b13ccb8e72; summary.csv e9d323d5001a186712e8d159e348f08d43825aed85186554aa9f141a70122a6c. Report and protocol in Git; full code and tables in conversation ZIP; no remote numerical replay or main merge.