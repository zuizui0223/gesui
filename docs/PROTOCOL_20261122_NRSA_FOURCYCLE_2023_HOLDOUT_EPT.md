# Four-cycle EPA NRSA EPT: 2023–24 external-era transfer protocol

**2026-11-22; exploratory, posthoc to known 2008–19 and 2023–24 associations.** This is not an independent preregistered experiment or a wastewater-source/terrestrial transfer result.

## Target and source opportunity before response fitting
Reuse the original three-era source-harmonized EPA data and *exactly* the previous 2,648 UID strict survey cohorts for 2008–09 BERWW (586 visits/437 HUC8), 2013–14 BERW (1,116/706) and 2018–19 BERW (946/603), source SHA and taxonomic quality assertions from [old PR #81](https://github.com/zuizui0223/gesui/pull/81). This study **does not** change old datasets or their previously reported 3-cycle results.

For **entirely held-out 2023–24**, acquire genuine EPA SiteInfo, WaterChem, PhysicalHabitat, BenthicCount and separate PeriphytonChlorophyll files with SHA256 pinned as in prior 2023 EPA source QA. Use original same UID/SITE_ID/VISIT_NO/DATE_COL joining; only BERW + WADEABLE and original counted `TOTAL300` (the same response definition as 2013 and 2018); do **not** pool 2023 BETB or 2008/2013/2018 boatable. Define source candidate based on actual valid chemistry TN/TP original blank NARS quality flags, positive mg/L-converted stream TN/TP and field conductivity, six measured original PHab fields, original HUC8, state and eco9. EPT response is only considered after independently auditing source join and method coverage.

Six authentic PHab features: original `XEMBED`, `XCDENMID`, `XPCAN`, `XFC_OHV`, `W1_HALL`, log10 field `XWIDTH` (retain native field units/definitions, do not rescale physical percentage if ambiguous). Other controls: lat/lon, stream order, season sin/cos, log conductivity, native state/eco9, same-source observed river log TN/TP. Require at least 500 strict target visits, 300 independent HUC8, >=30 states and >20 original source HUC8 not seen in old eras before numeric fitting.

**Two target cohorts, before looking at new model directions:** (i) all original source valid 2023–24 wadeable, (ii) source intersection with 2023–24 same UID original periphyton chlorophyll record (original eligible quality flag), to approximate old 2013/2018 periphyton-selected source. Retain both, no tuning to observed result.

## Core independent 2023 target test
Train on the 2,648 original previous cycles **excluding every HUC8 shared with the target 2023 cohort**. Inspect and report remaining training n and HUC8; train GO if >=500 visits and >=300 HUC8. Without any target tuning, use original PR #81 exact train-only continuous StandardScaler + state/eco9 OneHotEncoder(handle_unknown='ignore'), Ridge(alpha=10) and **train HUC8-equal sample weights**; no cycle year as an extrapolation predictor. Target is **unseen survey years 2023–24 and unseen HUC8** (actual target HUC8 may have appeared in raw training data but those HUC8 are removed completely from model fitting).

Models held fixed:
- B: geography+river order+season+conductivity+state+eco9.
- B+NP: B+measured river log10(TN/TP).
- B+PHab: B+all six original field habitat measurements.
- B+PHab+NP.
- Field width block **alone after B+NP**.
- Field width **conditional after B+NP+five other PHab**.
- Bed embeddedness **alone after B+NP** and conditional after all other PHab.
- Canopy/cover block XCDENMID + XPCAN + XFC_OHV after B+NP.
- Human disturbance W1_HALL after B+NP.

Primary metric: **equal independent target HUC8 heldout MSE gain**; secondary UID-weighted MSE. Report both increment NP after PHab and PHab after NP, plus width/bed/canopy block support, not cherry-pick strongest. Compute 2000 fixed prediction whole-HUC8 bootstrap descriptive95; this is conditional on fitted model, not full inferential uncertainty. For novel target-only original HUC8 that never appeared in any old source, show separately if >=100 HUC8 and >=150 visits. Repeat target predictions with visit-equal rather than HUC8-equal training as sensitivity. Report observed target MSE relative to own HUC8 mean outcome variance and category novelty.

## Controls and stop conditions
Original 2023 site HUC8 normalization strictly 8 digits; old historic harmonizer already handles 2008 7-digit and 2013 'H' prefixes. Reject unresolved geography, duplicate source UID, date mismatch or missing official source. Report sample attrition separately for 2023 full and 2023 periphyton-intersected subsets. Never turn 2023 absence of effluent DMR into no sewage. The old full-six PHab model's positive transfer in three earlier survey periods remains prior evidence, not a new test repeated as if independent. If 2023 source support fails, STOP: no numeric coefficients, no false effect. Interpret feature signs or prediction gains strictly as observational and potentially protocol-/time-confounded.

**Only after source and model results:** decide whether the 2023 independent-era target supports the old pattern that physical field habitat adds information after measured river nutrients. Do not claim a new ecological law or cross-ecosystem sewage mechanism from a single national survey series.
