# EPA NRSA 2018–19 / 2023–24: nutrient–periphyton mass–benthic composition at native UID

**Date:** 2026-11-13. **Status: Draft exploratory; all ecological coefficients computed from original EPA CSVs; independently audited numerically.**
**This is NOT WWTP effluent exposure attribution, nor a proven ecosystem-service chain.** All findings are observational spatial associations.

## Key independently audited results

Same original EPA native river site/visit units (UID); pre-numerical source gate and pre-existing 2018 source-frozen Benthic-UID panel; region, state, coordinates, stream order, sampling year/season and conductivity baseline. HUC8 whole-catchment heldout fivefold Ridge10. Primary excludes original NARS quality flags rather than treating them as zero. Response-specific MSE baseline is separately fit, **a predictive gain is not an effect size or proof of causation**.

| Period, response | Strict joint visits | HUC8 groups | Adding measured TN+TP to background, heldout MSE improvement | TN adjusted beta per +1 SD logTN [HUC8 cluster descriptive95] | TP adjusted beta per +1 SD logTP [95] |
|---|---:|---:|---:|---|---|
| 2018–19 benthic EPT **counted fraction** | 947 | 603 | **+1.7211%** | **−0.03109 [−0.05645,−0.00572]** | −0.02086 [−0.04296,+0.00124] |
| 2018–19 periphyton **log10 areal chlorophyll-a ug/cm²** | 947 | 603 | **−0.9651%** | +0.02842 [−0.03280,+0.08964] | −0.02302 [−0.08050,+0.03447] |
| 2023–24 periphyton **log10 areal AFDM mg/cm²** | 1,828 | 946 | **−0.08591%** | +0.01434 [−0.02035,+0.04902] | +0.02110 [−0.00875,+0.05094] |
| 2023–24 periphyton **log10 areal chlorophyll-a ug/cm²** | 1,828 | 946 | **+0.8121%** | +0.03159 [−0.01211,+0.07529] | **−0.07459 [−0.11532,−0.03386]** |

Interpretation: within the measured **river nutrient geography**, EPT relative composition has a negative TN signal in the older cohort, 2023–24 areal chlorophyll has a negative TP signal, and 2023–24 areal AFDM has **no robust additional nutrient prediction**. Therefore no universal measured **nutrient-driven increase of biomass with decreased biodiversity** has been demonstrated. No measured WWTP-derived N or P was connected to these NRSA stations.

Sensitivity accepting advisory quality flags but excluding Q/NR/ND, same fixed model:
- 2018–19 1,162 UID /696 HUC8: EPT TN+TP **+1.7053%**; areal chlorophyll **−1.0074%**.
- 2023–24 2,089 UID /1,029 HUC8: areal AFDM **−0.0539%**; areal chlorophyll **+1.1859%**.
- 2018–19 EPT additionally adding **measured areal chlorophyll** beyond TN+TP in same cohort changed baseline gain from **+1.7211% to +1.8657%**; adjusted Chla EPT coefficient **−0.00951** [−0.02401,+0.00498] zero-spanning. This is NOT observed mediation or transfer biomass.
- The 2018 and 2023 sites are NOT assumed identical physical sites. The chlorophyll chemistry result is not a measured before–after at the same river.
- EPT standardized subsampled *fraction* is not absolute insect abundance or adult emergence export.

## Critical source authenticity and measurement incompatibilities

1. Official EPA files from [NARS source site](https://www.epa.gov/national-aquatic-resource-surveys/data-national-aquatic-resource-surveys). For 2018/19, original pre-outcome fixed 1,221 NRSA benthic native UIDs are retained from [PR #58](https://github.com/zuizui0223/gesui/pull/58) run #43034107486. Full source SHA256 rechecked for eight original EPA CSVs, source IDs not reselected after coefficient signs.
2. 2018/19 chlorophyll-a original **UID** did not align: only 1 UID directly matched 2,112 original site UIDs; the source-only **SITE_ID + VISIT_NO + DATE_COL** one-to-one crosswalk resolved 2,109 additional records and 1 original exact UID (2 remain unavailable). This resolution preceded testing numeric effects, and source-level n1037 joint 2018 biomass/chem/chla IDs is an *availability* count, not the benthic-complete-case count.
3. **2018/19 AFDM** is mostly reported in volume-normalized **mg/L** (`RESULT_VOL`, 1000/1038) or bulk grams (`RESULT`,38/1038). **2023/24 AFDM** is original **areal mg/cm²** (`RESULT`,2104/2130). No unapproved conversion; 2018 AFDM was **not** modeled as comparable biomass. 
4. Areal chlorophyll-a reported in **ug/cm²** (2018 2017/2112 source rows, 2023 2128/2130). Chla is pigment standing stock, not measured gross primary production.
5. 2023/24 official webpage lists a benthic-count CSV but exact publisher link `nrsa2324_benthiccount.csv` and its metadata returned **HTTP 404** when authentically downloaded on GitHub Actions. Missing does NOT mean zero EPT. Thus no 2023–24 same-visit complete triad of AFDM + EPT + nutrient, and no numeric 2023 EPT effect should be claimed.
6. Same-UID numerical quality support: 2018 primary n947 visits in 603 catchments /47 states /881 physical sites; 2023 primary n1828 visits in 946 catchments /48 states /1,660 physical sites. No missing values or absent N/P sources imputed as zero.

## Reproduction, source validation, limits

- Original source-only multi-era schema, SHA and UID overlap: [GitHub Actions SUCCESS #43417456903](https://github.com/zuizui0223/gesui/actions/runs/43417456903).
- Strict source-only composite mapping, original categorical NARS flags/units and bridge: [SUCCESS #43417864266](https://github.com/zuizui0223/gesui/actions/runs/43417864266).
- Measured EPA original-sourced outcome fit: [SUCCESS #43418686337](https://github.com/zuizui0223/gesui/actions/runs/43418686337).
- **Separate independent sklearn Ridge + statsmodels** using temporary same-source panel (not a new field-data validation), exactly matched both eras × 2 quality rules × 2 outcomes × every MSE and per-SD slope within numerical tolerance: [SUCCESS #43419198293](https://github.com/zuizui0223/gesui/actions/runs/43419198293).
- Original full JSON fits and QA are Actions artifacts with source SHA/UID provenance. Original ecological source rows were not committed; ephemeral audit panels are not artifact-uploaded.
- Remaining confounders: upstream land use and runoff, light/shading and substratum, hydraulic disturbance, water quality contaminants, nutrient co-variation, nonlinear responses and measurement dates. HUC8 holdout / cluster descriptive intervals do not identify mechanisms or causal source effects.
- The 2023–24 negative TP↔Chla association could reflect any of these, not necessarily nutrient-induced suppression. Do not claim temporal trends or ecological dysfunction from it.

## Science decision and next source

**Decoupling of endpoints is suggested as an observational difference, not mechanistically demonstrated.** To prove the proposed WWTP→nutrients→biomass→EPT→terrestrial consumer-service sequence, find a genuine source-matched **same river/year** panel with (i) verified upstream measured plant effluent N/P, (ii) real site TN/TP and receiving-water flow/dilution, (iii) areal AFDM + benthic counts, and (iv) emergence biomass and consumer response. All stepwise mediation statements remain unsupported until actual same-unit measurements exist.

Preserve original European EEA benthic EQR results and the 2018/2022 French TP-I2M2 negative spatial gradient from other PRs; they are different source populations. No eBird, no site selection on positive effect, no main merge.
