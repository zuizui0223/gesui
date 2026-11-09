# Directed freshwater WWTP source→actual 2022 river TP, Germany and Austria

**2026-11-09 · Completed source-only and station-level source-versus-receptor prediction. Observational and POSTHOC, not a randomized wastewater phosphorus test.** This update proceeds beyond the previous country×1° WWTP source mass / 2022 river TP / Lycaenidae *public photo* association by using **river-reach routing**. The source-only directed exposure was Git-committed **without opening TP magnitudes**, before this model joined the already verified numeric WISE6 station-year medians.

## New independent geographic/hydrologic evidence

**Original reporting and river source:** Official 2022 EEA UWWTD original 5,492 facilities with both-positive **measured-reported** annual N and P discharge loads; 24,416 original discharge points; 4,057 direct-FW outfalls with exactly one unique full reporter+facility key within the pinned Europe study bounds. The exact original 31-page EEA REST response digest sequence was verified to match the earlier independently preserved `gesui` source-only 2022 UWWTD PR #9.

**River network:** Independently hosted HydroRIVERS EU v1 converted GeoParquet original SHA256 `5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f`; 199,545 original directed reaches in study geometry. **Not independently matched byte-by-byte or topologically against official HydroRIVERS shapefile**. Original directed `NEXT_DOWN` and `LENGTH_KM` define **strictly upstream, different-reach only, ≤50 km**, source snap ≤500m and TP receptor snap ≤500m. **Same-reach sources are excluded** because outlet versus TP site longitudinal order is unresolved. New WISE6 site coordinates matched public original site country+ID+scheme and the original 187 pre-numeric TP source registry; **none of the 2022 numerical TP concentrations were used to choose source or monitoring points**.

| New official source and directed snap | Number |
|---|---:|
| Original unique direct freshwater both-positive reported measured N/P WWTP outfalls | **4,057** |
| Original outfalls unambiguously snapped ≤500m to pinned rivers | **2,967** |
| Original independent WISE6 2022 TP source-qualified stations | **187** |
| WISE6 original monitoring stations uniquely snapped ≤500m | **174** |
| WISE6 stations with both strict numeric TP QC and unique directed reach | **171** |
| Numeric QC + routed **registered strictly upstream source-positive** stations | **114** |
| Numeric QC + routed **no qualifying registered upstream source** stations | **57** |
| Unique routed monitoring reaches among the 171 | **169** |
| Unique source-positive river reaches among 114 source-positive stations | **112** |
| Original country×2° geographic source-positive blocks | **20** |

Among the 171 eligible actual 2022 WISE6 TP station-years, **Germany 95** sites (55 registered strict upstream source-positive / 40 registered source-zero); **Austria 76** (59 source-positive / 17 source-zero). The sample passes the pre-source thresholds ≥20 each pooled, ≥5 positive each country, ≥5 independent country×2° blocks and ≥15 distinct exposed reaches. Recorded stations with no qualifying registered upstream source are **not guaranteed untreated**, because EEA reporter coverage, multiple-outlet ambiguity, unmeasured/zero N/P, other phosphorus inputs, and even unregistered wastewater remain.

[Strict directed 187-site source-only full original QA success](https://github.com/zuizui0223/gesui/actions/runs/42578347508). The original 187-site source-only reach/flow/WWTP P/N and downstream placebo SHA records are in [2022 directed source QC](2022_WISE6_187_TP_STATIONS_DIRECTED_ORIGINAL_EEA_WWTP_SOURCE_ONLY_QC.json).

## Same 171 actual numeric TP stations: directed environment source prediction, with reverse-direction negative control

Source-dependent TP model was specified before rejoining the TP outcome to the directed source registry. Response: `log1p(actual 2022 WISE6 TP site median in mg{P}/L)`, using the already independently processed WISE6 original 2022 sample-date/LOQ quality rules (no below-LOQ imputation). Baseline: original country, station latitude/longitude, **HydroRIVERS modeled mean flow** and upstream drainage area, independent 2021 ESA WorldCover country×1° built-up/cropland/country-terrestrial land area. Predictors centered/scaled **on training only**, source log1p standardized **within each country and training only**. Spatial holdout: leave one original **country×2° block out**, **21** groups. Unweighted independent monitoring sites (169 distinct original reaches; two reaches host multiple sites), original 187 source sites not resampled by TP concentration. Original 2022 mass reported from WWTPs strictly upstream only; downstream WWTP source annual P mass is the **directional negative-control**.

| Additional candidate source to fixed geographic + hydrologic baseline | MSE baseline | MSE + candidate | 2° block-heldout MSE improvement |
|---|---:|---:|---:|
| **Strict upstream original registered WWTP P mass** | 0.00202926 | 0.00206659 | **−1.84%** |
| Strict downstream registered WWTP **P-mass placebo** | 0.00202926 | 0.00207445 | **−2.23%** |
| Upstream WWTP P mass / modeled HydroRIVERS Q | 0.00202926 | 0.00204608 | **−0.83%** |
| Strict upstream WWTP **N mass** | 0.00202926 | 0.00208087 | **−2.54%** |
| Strict upstream qualified WWTP **count** | 0.00202926 | 0.00210414 | **−3.69%** |

**Negative heldout MSE improvement means the source predictor increased prediction error.** None of the original source P/N/count predictors yields a positive pooled region-heldout improvement in actual measured river TP after the geographical, modeled flow/catchment and country land-cover controls.

- Strict upstream P original source: DE **−3.60%** MSE gain across 95 sites, AT **+0.24%** across 76; original country×2° prediction-error block-bootstrap descriptive 95% interval **[−5.35%, +0.62%]**. This does not establish a robust directionally source-specific TP effect.
- Strict downstream P placebo: DE **−2.25%**, AT **−2.20%**, pooled **−2.23%**, 95% descriptive block interval **[−6.48%, +1.50%]**. The upstream result is not much stronger than the reverse-direction source-control. Apparent original P/flow index also fails pooled positive gain (−0.83%).
- Raw, unadjusted source-positive versus registered-source-zero **site median TP**: DE **0.107 vs 0.11575 mg P/L** (55 vs40 original stations), AT **0.0323 vs 0.0175 mg P/L** (59 vs17). These unadjusted differences are in **opposite directions across countries**, and are NOT source causal differences or comparable randomized catchments.
- Independently recomputed the exact per-site original 2°-block predicted errors for all **five** candidate source models; reported MSE gains match arithmetic to floating point tolerance across **171 unique station hashes**, **169 river reaches**, **21 blocks**.

[Actual 2022 numeric river TP + original directed upstream source / downstream negative-control out-of-region predictive test SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42579054225). Full original 171 site source and 2° holdout predictions are in [source-specific region CV QC](2022_REAL_WISE6_TP_DIRECTED_UPSTREAM_WWTP_P_VS_DOWNSTREAM_PLACEBO_CV.json).

## Scope and interpretation

This specifically tests **year-aligned reported WWTP source-to-recipient water chemistry geography** with directed river reach topology, but it is STILL **observational**: source annual reported P tonnes at an original outlet do not equal 2022 in-river TP mg{P}/L at individual dates, and actual 2022 river discharge, WWTP removal efficiency, retention, diffuse agricultural sources, hydrological barriers, post-source in-stream processes, and monitoring site selection are unresolved. HydroRIVERS mean Q is modeled climatology, not local same-day observed streamflow. Source-positive stations may be urban, and source-zero sites may receive unregistered or nonqualifying wastewater. The third-party river conversion's topology has not been directly reconciled with official original vectors. Multiple site monitors along one river corridor are not fully independent organisms.

**What is supported:** We can now (a) link 2022 official WISE6 **184** quantified TP site-years to precise public WISE positions, (b) source-match **171** stations to actual directed rivers, and (c) find **114** with one or more qualifying registered WWTPs strictly upstream 50km. This overturns the old *data-source unavailable* barrier and establishes a real source-connected environmental comparison.

**What is NOT supported:** On this source-qualified station panel, official WWTP P annual reported loads do **not** improve heldout prediction of measured TP beyond geographic/hydrologic/WorldCover land variables. This does **not** prove WWTPs cannot raise TP in a specific river; spatial blocked predictive gain is not the same estimand as causal source→river increment, and annual source mass, river mixing and catchment confounding remain. The prior coarse 38-country×1° river-TP/photo connection and the original 2021–2023 Lycaenidae **public photograph-recording** geography still cannot be rebranded as a wastewater-driven butterfly benefit or P ecological mediation.

## Additional publication gate

The EEA WISE6 original 2022 `resultObservationStatus`, `metadata_statusCode` and `metadata_observationStatus` codes at the exact 187 preselected TP monitoring sites are being independently audited. Before publication, their values/semantics and expert QC flags must be reviewed. Initial numeric TP QC has already excluded 53 of 2,471 original at-site samples below quantification, without zero/half-LOQ substitution, and 3 of the 187 original stations with insufficient valid measured dates or censoring. No additional flag-specific exclusion has yet been decided on unknown status meanings.

## Next actually discriminating ecological mechanism test

Rather than shifting photographs or regression thresholds again, a follow-on source-only screen should identify **upstream and downstream actual WISE monitoring-station pairs along the same directed river**, with at least one qualifying reported WWTP outfall **strictly between** the station reaches, and distinguish those from within-river monitoring pairs with **no registered intervening source**. Fix station pairs and source support before inspecting their TP concentration difference. Such a paired test is more mechanistically local than the current 171-station unrelated-river source comparison. If insufficient station pairs, STOP rather than tune distance thresholds posthoc.

Original user-related research positives remain separate: `chocho` **22/25 Lycaenidae** potential *introduced larval-host resource-geography* expansion, `gesui` independently positive 2018/2022 aquatic macroinvertebrate EQR prediction association, and the earlier exploratory country×1° 2022 WWTP-P Lycaenidae **photographic reporting** geographic association (attenuated by 2021–2023 historical butterfly photos and not explained by actual TP). No adult protein/sodium puddling or larval survival was measured.

## Source references

- [Official EEA UWWTD 2022 source](https://www.eea.europa.eu/en/datahub/datahubitem-view/6244937d-1c2c-47f5-bdf1-33ca01ff1715)
- [Original EEA Waterbase Water Quality ICM 2026](https://www.eea.europa.eu/en/datahub/datahubitem-view/fbf3717c-cd7b-4785-933a-d0cf510542e1)
- [Original source-only reported WWTP river routing, gesui PR #9](https://github.com/zuizui0223/gesui/pull/9)
- [New real 2022 WISE6 measured river TP and originally coarsened Lycaenidae photo trial, PR #28](https://github.com/zuizui0223/gesui/pull/28)
- [Official HydroRIVERS site](https://www.hydrosheds.org/products/hydrorivers); river topology source used above is a source-pinned **third-party** converted mirror, not this official file.


## Completed original 2022 TP sample status QA and worst-case U-flag sensitivity

Additional source-only original 2026 WISE6 status audit read the **same complete 102,569,057 source lines** a second time, but only raw source status codes (not numeric TP values) at the pre-frozen **187** WISE6 real 2022 TP sites (**2,471 original river Total-P sample rows**):

| Original reporting/quality status | Original sample rows |
| --- | ---: |
| Original metadata_statusCode = stable | **1,043** |
| Original metadata_statusCode = valid | **1,428** |
| Original metadata_observationStatus = A | **2,462** |
| Original metadata_observationStatus = U | **9** |
| Original resultObservationStatus = A | **1,043** |
| Original resultObservationStatus blank | **1,428** |
| Original below-LOQ resultQualityObservedValueBelowLOQ = 1 | **53** |

All nine metadata **U** records occurred in **nine distinct original DE stations**, just one per station among otherwise normal reported observations. These are original code values, never silently reinterpreted as zeros/valid or propagated across countries. Source uses original station and unique source ZIP SHA; [original independent status audit SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42579391197), [full original per-site original code histogram](2022_WISE6_TP_ORIGINAL_187_SITE_RESULT_AND_METADATA_STATUS_ONLY_QC.json). The Eionet WISE data dictionary confirms below-LOQ observations must be flagged and values should not be used as observed concentration; our 53 censored rows were excluded with no half-LOQ substitution: https://dd.eionet.europa.eu/datasets/latest/WISE-SoE_WaterQualityICM/tables/DisaggregatedData . The **metadata U-code semantics remain to be confirmed specifically from the current official 2026 dictionary** before publication, rather than blindly assuming they are independently validated low-accuracy samples.

A deliberately conservative, explicitly POSTHOC sensitivity removes **every station with even one original U row**, not just those nine individual measured dates. This yields **162** independently TP-qualified 2022 WISE6 stations and actual public directed reaches: **112** source-positive, **50** registered-source-zero; AT59 positive and DE53 positive, **110** positive distinct reaches, **20** positive original country×2° regional blocks. Same exact original 2° regional holdouts, country, true site locations, modeled discharge and catchment, and ESA WorldCover2021 built/crop/land controls:

| Original directed 2022 WWTP source variable | Initial 171-site relative MSE gain | Strict no-any-U 162-site relative MSE gain |
| --- | ---: | ---: |
| **Strict upstream original EEA P reported load** | **−1.84%** | **−1.68%** |
| Strict downstream original EEA P **directional placebo** | −2.23% | −2.62% |
| Strict upstream P divided by modeled discharge | −0.83% | **+0.15%** |
| Strict upstream original N reported load | −2.54% | −2.27% |
| Strict upstream qualified outfall count | −3.69% | −3.38% |

The original 171-site model was independently re-executed from all frozen sources before filtering, and every originally reported MSE gain matched exactly. On the stricter 162-site panel, **P total source does not improve prediction**, and the small +0.15% increase in P/modelQ source proxy lies within a highly uncertain descriptive block-bootstrap interval **[−9.31%, +5.21%]**. No evidence of a positive robust annual source-P transport effect is established.

[Full status-exclusion audit Action SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42580204403) · [strict original nine-flag station exclusion and heldout results JSON](2022_WISE6_DIRECTED_WWTP_TP_EXCLUDE_9_ORIGINAL_META_U_SITES_POSTHOC.json).

**Publication qualification:** this excludes original U-marked sites as a conservative stress test, not a formal proof all remaining original source rows pass every additional expert EEA QC criterion. Main original 171-site and strict 162-site analyses must be reported together; source-positive/negative labels remain nonrandom and do not identify causal effects.
