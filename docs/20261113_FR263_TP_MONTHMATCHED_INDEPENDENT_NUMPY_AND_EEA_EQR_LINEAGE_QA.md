# Independent original EEA EQR lineage and measured TP monthly model numerical replay (2026-11-13)

**POSTHOC quality/replication audit AFTER prior month-matched outcome already observed as near-null (73 H7 rho−0.1217; H6-heldout +0.259% MSE gain, 95 [−3.99,+4.10%])**. It must not be portrayed as independent field-site replication.

## Original official biology response provenance

Source script `scripts/gesui_2018_independent_bqe_EQR_outcomes.py` reads official original 2018 EEA `resultNormalisedEQRValue` into `numeric_EQR_2018`; its `raw_EQR_2018` is the *different* EEA `resultEQRValue`. Source script `scripts/gesui_model_2022_eea_aquatic_EQR.py` reads official original 2022 `resultNormalisedEQRValue` into `normalised_EQR` and raw `resultEQRValue` into a separate raw response. Both archived true original EEA site/year source CSVs are SHA-pinned:
2018 Git blob `0f30bfbc3dac13f6a182c8a910e9899c04778322` and 2022 `325f199e985d9862ffee5b0edb205c61a7be3833`.
This proves both **analysis source columns represent EEA's normalized EQR** and original 2018 raw EQR is not incorrectly being subtracted from normalized 2022. It does NOT prove unchanged sampling protocol, reference class, scoring calibration, or inter-year lab processing just because both labels are I2M2. Count original 2018 FR I2M2 raw-vs-normalized differences as a consistency audit, no sign-based filtering. Source original 263 physical EEA station IDs match both years and original 2022 method I2M2, country France, original native `HYRIV_ID`.

## Independent month-matched response arithmetic

Source common-month eligibility original SHA256 `a02a9f724fb1b1671a9f248ed169620bc25db8b9a5914be566f6728e25b2c4f1`, source gate [#43378522025](https://github.com/zuizui0223/gesui/actions/runs/43378522025): all exact 263 field sites, 73 native H7, 24 H6, ≥3 shared actual sampled calendar months per site. Original source-only TP.2 official ZIP SHA `26f2af11f4972987ee6f7c3f497ad902820f145851f71919d4f45b589e1d809e`, all original measured TP annual medians of 263 sites reproduced.

Use only original source-frozen month-matched RESULTS/CVS published by [#43378917089](https://github.com/zuizui0223/gesui/actions/runs/43378917089), no source/year/site re-selection. Independently recompute all seven exact original baseline + month-matched TP, baseline + raw annual TP, and sensitivity baseline + 2018 monthly level TP Ridge10 MSE values on original same 73 H7 /24 H6 hydrological group fivefold and 24fold leaveparentout with NumPy normal equations rather than original sklearn SVD. All raw TP rank correlations and original n263 station-repeated 2018/22 EQR difference must match previous source to 1e−9 or STOP.

No reinterpreting the original small month-matched H6 CV +0.259% [−3.99,+4.10] as evidence of nutrient causation or environmental service response. If independent code disagrees, report the discrepancy and withhold original effect until resolved.

No merges, no modifications to `gesui/main`, `chocho` or previous original PRs.