# Source-to-receptor reality check: directed upstream 2022 WWTP P does not improve held-out measured river TP prediction

**9 November 2026 — completed on a new, SHA-locked original 2022 source/chemistry panel.** This analysis directly addresses the missing *hydrologic source connectivity* of the earlier **country×1°** WWTP-P versus Lycaenidae *photographic-reporting* results, but is **NOT** a test of butterfly nutrition, larval fitness or causal wastewater treatment. All model construction after prior Lycaenidae-positive discovery is explicitly POSTHOC. Positive and negative checks are archived equally.

## 1. Actual independently observed 2022 water chemistry with physical source direction

- EEA original Waterbase WISE6 2026 `DisaggregatedData` archive, original public ZIP SHA256 `2a2ef2a43a16ec068b7974a80ae5fca49e44697660a46ce3b91fb1813abe62d2`. WISE6 **2022 actual river Total phosphorus**, CAS_7723-14-0, water matrix `W`, unit **mg{P}/L**, original per-site day median then 2022 station-year median, no substitution for below-LOQ, and ≥4 independently quantified days in ≥3 distinct 2022 months. 187 original station IDs were source-frozen *before* numeric TP values were read; **184** passed original numerical QC.
- EEA UWWTD `TreatmentPlants_2022` original **5,492** both-positive **measured-reported** annual N/P facilities joined by the full reporter+facility ID to `DischargePoints_2022` original 24,416 points; original **4,057** unique one-outfall **freshwater** sources inside the Europe source bounds. Original **2,958** discharge points could be unambiguously snapped within 500m to the independently distributed river network after checking competing reaches within 2 metres; 1,088 were farther away and 11 ambiguous.
- HydroRIVERS Europe v1-derived **THIRD-PARTY GeoParquet** is SHA-frozen at `5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f`. Its underlying source geometry/topology has **NOT** been independently verified against the official HydroRIVERS shapefile. This limitation is explicitly retained.
- Original WISE6 station coordinates are from official 2026 public SpatialObject source joined on `country|site identifier|site identifier scheme`. Strict nearest river segment ≤500m; reject ambiguous competing reach lines ≤2m of the closest; no within-reach source/order guessing. Original discharge-to-monitoring path follows explicit HydroRIVERS `NEXT_DOWN` within **50km** of original directed river length, originating source reaches excluded. Distinct **downstream 50km** WWTP sources are retained as predeclared transport-negative control.

**Source-only original 187 stations:** 173 uniquely snapped, 113 downstream of ≥1 qualifying upstream registered original WWTP and 60 with zero *qualifying registered* upstream source, 14 off-river/ambiguous. Positive source original reaches 111 and 20 country×2° blocks. Source connectivity is not sewage chemistry measured at the receiving water station.

**Numerical quality/site selection independently committed without using TP magnitudes:** 170 of 187 sites remain after strict quantified-date/LOQ and graph/flow source QC, 113 with ≥1 upstream registered qualifying source, 57 original registered-source-zero. DE **55 exposed / 40 zero**; AT **58 / 17**. Two Austrian original river reaches each had 2 eligible WISE6 water stations; these were collapsed *by median of station-year medians* to the **168 independent reach** primary panel, 111 source-positive and 57 source-zero, 21 country×2° heldout blocks. Every reach had original modeled mean discharge and drainage area and independent ESA WorldCover2021 country-source-cell built/crop/land area. All original source-support criteria PASSED.

Sources:
- [187-site directed source-only original outcome-blind result](https://github.com/zuizui0223/gesui/blob/analysis/wise6-2022-TP-directed-WWTP-reaches-20261109/results/lycaenidae_2022_measured_river_P_bridge/2022_WISE6_187_TP_STATIONS_DIRECTED_ORIGINAL_EEA_WWTP_SOURCE_ONLY_QC.json).
- [Original numeric-QC-only 170-site source-panel SHA-locked BEFORE testing TP magnitudes](https://github.com/zuizui0223/gesui/blob/analysis/wise6-2022-TP-directed-WWTP-reaches-20261109/results/wise6_2022_riverTP_directed_WWTP_source_only/2022_NUMERIC_ELIGIBLE_TP170_REGISTERED_UPSTREAM_WWTP_SOURCE_PANEL_FREEZE.json).
- [Source-only 170-site gate workflow SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42586255469).

## 2. Real TP environmental model with hydrologically upstream and downstream registered-source controls

**Response:** log(1 + 2022 actual measured TP mg{P}/L), 168 original independently identified receiving HydroRIVERS reaches, equal reach weight. Same original country×2° region left out (21 blocks) in turn; all source centering and scaling fit within the training data for each test region. Baseline is **country, lat, lon, log1p(original modeled mean river discharge), log1p(original upstream basin area), ESA WorldCover2021 measured built/crop fraction and country land area**. Each official 2022 WWTP annual source axis entered **separately**. The annual WWTP source **P tonnes/year** is an indicator of reported origin mass, NOT in-river delivered P or sample-date exposure.

| Independently fixed source covariate added to baseline | Held-out baseline MSE | Extended MSE | Relative heldout MSE improvement (positive better) | 2° block error bootstrap descriptive 95% |
| --- | ---: | ---: | ---: | ---: |
| **Strictly upstream ≤50km registered WWTP P mass (t/year)** | 0.00206085 | 0.00209967 | **−1.88%** | **[−5.35%, +0.52%]** |
| Strictly upstream registered WWTP N mass (t/year) | 0.00206085 | 0.00211451 | −2.60% | [−6.68%, +0.06%] |
| Number of strictly upstream original qualifying registered WWTPs | 0.00206085 | 0.00213542 | −3.62% | [−8.43%, +0.30%] |
| **STRICTLY downstream registered WWTP P only — negative-control placebo** | 0.00206085 | 0.00211094 | **−2.43%** | [−6.92%, +1.40%] |
| Strictly upstream annual reported P mass ÷ original modeled river flow — potential mg/L equivalent, NOT dose | 0.00206085 | 0.00208134 | −0.99% | [−10.75%, +3.73%] |

**Within-country heldout upstream reported P MSE improvement:** DE **−3.59%** across 95 independent receiving reaches; AT **+0.15%** across 73. Under 199 original-source-mass randomizations *within countries*, **91.5%** of randomized registered P placements produced equal or better out-of-block prediction improvement than the actual original P geography; within-country×2° region randomizations this proportion is **82.5%**. These are **source-label negative controls** and are NOT randomized WWTP operations or a treatment-effect p-value.

**Descriptive real river TP without attributing cause:**
- DE original reach-wise median measured TP: upstream-source-positive **0.107 mg P/L** (n=55) versus no registered qualifying upstream source **0.11575 mg P/L** (n=40).
- AT: upstream-source-positive **0.03415 mg P/L** (n=56 independent reaches) versus registered-source-zero **0.0175 mg P/L** (n=17).
- Directions differ between countries; neither is a valid causal effect estimate because land use, upstream industry/agriculture, discharge routing, sampling, and reported source coverage differ.

Original [regional-block environmental validation SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42586833456). [Full source/reach-level outcome, every heldout prediction and both randomization receipts](https://github.com/zuizui0223/gesui/blob/analysis/wise6-2022-TP-directed-WWTP-reaches-20261109/results/wise6_2022_riverTP_directed_WWTP_source_only/2022_REAL_TP_168_REACH_UPSTREAM_WWTP_P_VS_DOWNSTREAM_PLACEBO_BLOCKCV.json). The primary code executed after source-only sample freeze with the exact original input blob checks.

## 3. Ecological implication and limitations

**Correct positive conclusion:** The original measured 2022 water-quality source and original public EEA WWTP freshwater outlets can be mapped into a **strictly directed river network**: **168 independently source-classified reaches**, both DE and AT, with enough original measured TP site-year QC and negative-control support for a cross-region observational test. This is a materially better environmental measurement than co-membership of 1° photography/source cells.

**Correct environmental result:** In this preselected 2022 source/photo geographic cohort, the *annual reported positive WWTP P load* does **not** produce a measurable additional improvement for out-of-region predicting real river TP after original river mean flow, catchment area, geography, and 2021 measured land-cover controls. The source-negative P-shuffle and downstream-source placebo do not uniquely support positive physical source delivery. A particular upstream WWTP might nonetheless affect river chemistry locally; heterogeneous original sources, unreported outlets, within-reach ordering, time averaging, flow dilution, natural & agricultural phosphorus, and transformed third-party river topology limit inference.

**Not established:** sewage → higher real TP at all sampled rivers; real TP → increased Lycaenidae populations; actual adult puddling sodium/protein N use; larval host nutrient/fitness effect, or a mediating causal chain. Original 2022–2023 iNaturalist *photographic reporting composition* remains a different and spatially persistent detection process. Standardized butterfly 2022 WWTP-outlet sampling lacked positive support in prior studies. The original `chocho` **22/25 Lycaenidae** introduced-host **potential-resource geography** finding and `gesui` 2018/2022 benthic macroinvertebrate EQR incremental prediction gains remain independent and unaltered.

**Publication blockers not silently cleared:** independently verify third-party converted HydroRIVERS reach topology and source coordinate positions against the original official HydroRIVERS v1 shapefile; review original EEA WISE6 `resultObservationStatus` plus other expert-QC flags beyond LOQ; validate official WWTP reported measured annual P mass units against reporter dictionary; consider proper date-aligned WWTP effluent chemistry, station flow and nutrient retention. This is a **posthoc environmentally negative predictive result**, not an experiment disproving wastewater eutrophication.

## 4. Why the next national river-pair analysis is a separate study

An independently source-frozen extension to all DE/AT 2022 date-eligible river TP stations (~301 public source-geolocated sites) and strictly same directed river upstream/downstream monitoring pairs exists in a separate Draft PR chain, avoiding the original selection of 42 country×1° *photo-opportunity* regions. It must be evaluated as **a new national population with fresh, independently frozen site/year physical source rules**, not retroactively substituted for these 168 reaches because this photo-region sample yielded weak WWTP-P environmental prediction.

[National original source-only pair gate (separate branch, not part of this analysis)](https://github.com/zuizui0223/gesui/tree/analysis/wise6-2022-national-DEAT-river-TP-pairs-20261109/results/wise6_2022_national_DEAT_TP_sourceonly).

**This Draft report explicitly distinguishes the field data support success from the lack of a transferable reported-P environmental signal.** No `main` or `chocho` manuscript alterations.
