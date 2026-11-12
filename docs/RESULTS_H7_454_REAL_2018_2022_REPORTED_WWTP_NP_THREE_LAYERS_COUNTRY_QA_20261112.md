# Real 2018/2022 EEA annual reported upstream N/P source versus EQR, plant S50 and Lycaenidae, SAME H7×period

**2026-11-12 | Genuine original-data source-reconciliation and ecological regression | EXPLORATORY POSTHOC | Draft, no merge.**

## Main source upgrade beyond static HydroWASTE

Prior PR #65 found a spatial **negative** association between the static mixed-quality HydroWASTE outlet-pressure model and actual EEA river benthic EQR (+7.33% marginal H7-hash spatial holdout gain across 454 same H7×year rows), while GBIF-indexed plant and Lycaenidae responses were weak. This source cannot identify annual effluent mass or actual load at river station.

Here, original *reported* EEA 2018 / 2022 upstream WWTP **annual N and P tonnes/year** are rejoined to original **same EEA monitoring station** through the prior HydroRIVERS upstream ≤50km source matching, scaled by the **original modeled mean river discharge DIS_AV_CMS**. The parent source-matching study had many monitoring stations without a matched positive reported WWTP source; those should NOT be assumed genuinely pristine or WWTP-free.

Conversion from upstream reported t/year divided by modeled m³/s to **input mg/L-equivalent** (hypothetical full delivery, no river retention), anchored to 365.25 days/year:

`mg L^-1-equivalent = (reported source t year^-1 / modeled river m3 s^-1) * 1e6 / (365.25 * 86400)`.

This is NOT measured chemical TN or TP concentration at the monitoring station, and it cannot distinguish agricultural loads, flow dilution and operational treatment effects. True P and N source reports are strongly correlated.

## Source-only gate before opening this new model's outcome

- Original 454 fixed original H7 basin×period rows, 232 unique H7 basins and 15 countries, 227 early /227 late. These are the same original 2018/2022 EEA EQR, genuine GBIF plant 50-rarefied speciesKey S50 and five-family Lycaenidae log-odds response rows from prior PR #65.
- Original 3,525 monitored source station-year contributions in those 454 rows, after original strict station match and valid modeled flow. **132 original station-year source reports positive for upstream N or P**.
- **56/454 H7 basin-periods** have strictly positive matched reported upstream WWTP P or N; **38 unique H7 basins**, 28 positive-period rows early and 28 late. Positive reports span **13 original countries**, with France only **4** of the 56 (but France 299 of 454 total source period rows). The other **398 source-periods have no matched positive reported load**, NOT a source-confirmed absence of sewage.
- Real EEA source original Git blob 2018 `0f30bfbc3dac13f6a182c8a910e9899c04778322`, 2022 `325f199e985d9862ffee5b0edb205c61a7be3833`. Original H7 source SHA `a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89`. PRE-new-effect signed source-only matching manifest SHA256 **`012ad79feabc8c11e0d74f803730b5b7913b96a64a6af26fb0c97263d9aa7c26`**.
- Original source-only eligibility threshold ≥30 positive H7-year, ≥20 unique H7 basin positive, ≥10 positive in each original era, and ≥3 original countries passed **before reading newly joined biological effects**. [Original source-only annual N/P gate SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43222864099).

## Same-source three response models (source first, but known original static sign)

A fixed original baseline includes previous year, country FE, original benthic assessment-method FE, area, modeled river flow, matched EEA monitoring-site counts, original GBIF Animalia/Plantae/five-butterfly reporting counts, and original HydroWASTE **static** registered/modelled source pressure. Exact 454 original sample and original baseline 5fold Ridge alpha10 H7 disjoint folds were **reproduced to numerical precision before adding new annual N/P covariates**.

New source predictors were predeclared separately: reported-positive binary, `log1p(site mean P µg/L-equivalent)`, `log1p(site mean N mg/L-equivalent)`, and joint P+N. Source-positive zero here means no positive *reported matched source*, not true zero effluent. 2000 H7 basin-bootstrap descriptive uncertainty intervals. H7 hash CV can leak nearby watershed geography, so **NOT** a fully spatially blocked test.

| Same 454 original H7×period rows, separate outcomes | Additional reported N 5fold heldout MSE gain | Additional reported P | Additional P+N |
|---|---:|---:|---:|
| **Actual EEA benthic EQR** | **+0.9276%** | **+0.3756%** | +0.6267% |
| GBIF plant source-index rarefied S50 | −0.3230% | −0.3920% | −1.8699% |
| GBIF original five-family Lycaenidae log odds | −0.8435% | +0.0143% | −1.5916% |

**EQR N** original 2,000 H7 cluster-resampling 95% predictive gain **[+0.0668%, +1.9457%]**; original N partial adjusted H7 cluster coefficient **−0.01296 original EQR units per original 1-SD indexed N intensity**, descriptive 95 **[−0.01997, −0.00596]**. Original equal-country-total weight EQR/N coefficient **−0.00889** [−0.01762,−0.00015]. P-alone adjusted EQR coefficient −0.00823 [−0.0208,+0.00434], MSE gain bootstrap [−0.416%,+1.201%]. Reported P/N indexed source log correlation **r=0.936**; cannot assign effect to N alone vs P.

**The simple binary "any matched upstream report" variable is worse** as an EQR predictor (−0.220% extra OOS MSE gain). Source intensity rather than registered presence appears predictive only in a narrow subset of source-positive stations.

Main fully original **reported annual EEA N/P into the same 3-layer H7 ecological responses SUCCESS** [GitHub Actions #43223533195](https://github.com/zuizui0223/gesui/actions/runs/43223533195). Stored source-locked original result JSON and original macro cohort are downloadable as that run artifact. No data were simulated.

## Strict posthoc independent numeric/code robustness: the positive N signal is STRONGLY country heterogeneous

A second separately implemented original source-locked Python script reproduces the exact original 3-outcome H7 CV MSE for **the baseline and N-added models**, independent of first fitting code. **All subsets and top-intensity caps below were inspected AFTER seeing original positive N source result**:

| EQR original source support condition | H7-hash heldout MSE gain from adding original reported N |
|---|---:|
| Whole original 454 H7-period source | **+0.928%** |
| 2017–2019 only (n227, 28 source-positive) | +0.921% |
| 2020–2022 only (n227, 28 source-positive) | +0.319% |
| **France-only (n299, ONLY 4 source-positive)** | **−2.493%** |
| **France-excluded (n155, 82 H7 basins, 52 source-positive)** | **+3.618%** |
| Country/method mixed river-basin-years excluded n442, n48 source-positive | +0.692% |
| Top positive N intensity 95th percentile winsor | +0.978% |
| Top positive N intensity 90th percentile winsor | +0.807% |

The signs are positive except within France, and the French source/reporting coverage is almost absent. **Do NOT call this a pan-European or France-replicated N biological result.** The "overall positive N" is mainly the non-France reporting-positive geography (N source: 52/56 positive basin-periods non-French), even though 299/454 sample rows are French. The large France-only regression learns essentially no reported point source and suffers added predictor variance.

Source-positive-only (n56/38 H7 basins) original raw (UNADJUSTED) Spearman rank diagnostics: source N vs EQR **−0.298**, P vs EQR **−0.267**; source N vs plant S50 **+0.417**, P vs S50 +0.371; source N vs Lycaenidae family log odds **−0.306**, P vs butterfly −0.315. This is a very small and **selected** positive-source set across varying countries and years. The plant-positive rank association **does not persist as a stable new prediction** in the full modeled/recording-adjusted n454 panel. No claim of sewage positively affecting plants or butterflies is justified.

[Original 3-outcome/5fold numeric independent exact replication and period/country/source-leverage diagnostics SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43224180601). These are diagnostics on the same original field data, NOT new observational replication.

## Macroecological implication, without causal overstatement

- Replacing static globally modeled outlet pressure with actual European time-specific **reported annual nutrient SOURCE mass** retains a small negative benthic EQR predictor sign **where source reports are recorded**, particularly outside France. This is coherent with prior reported-load EEA 2018/2022 small predictive gains and US NRSA river measured TP vs EPT composition — but different observation grain and original datasets.
- The same original reported source N/P environmental predictor **does not consistently improve both GBIF-indexed terrestrial biota responses** after recording opportunity controls. This is an intriguing **differential predictive sensitivity across measurement layers**, NOT a demonstration that terrestrial biology is insulated from wastewater or that ecosystem service cascades fail.
- **Source reporting failure dominates the next inference problem**: positive source N/P in only 56/454 H7-periods and 132/3525 station-year records; absence of report cannot establish zero emissions or control sites. Source regimes concentrated outside France, country/method and land use confounded, and P/N highly correlated. More aggressive regression on sparse data will not solve this.
- Remaining *critical physical design*: true station-specific upstream all-WWTP active source exposure (not basin outlet), original HydroATLAS local/upstream agriculture, artificial land cover, human population, and receiving water modeled flow, plus ideally contemporaneous measured river TN/TP. HydroATLAS source catalog itself is globally available but multi-GB download and has not been joined in this PR.
- Measured field emergence, predator assimilation, plant host, butterfly pollination or seed set are still unavailable on these matching basin-periods. Never infer causal benefit/damage from this index.

**All previous confirmed positives, nulls and already-frozen source studies preserved. Draft / unmerged / no eBird / no main or chocho changes.**
