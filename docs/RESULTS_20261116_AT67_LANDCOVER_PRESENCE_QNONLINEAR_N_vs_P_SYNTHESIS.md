# 2026-11-16 final source-authentic Austria67 synthesis: WorldCover land use, WWTP source reporting, nonlinear flow and N/P identity

**State:** All FOUR new same-source AWS/ESA/GitHub Actions analytical checks SUCCESS and signed. They concern the SAME 67 original 2022 Austrian water field sites, and are exploratory/posthoc, not four independent ecological replications. Their evidence is quantitative **spatial predictive associations** and mechanism ambiguity, not wastewater treatment causality, aquatic-terrestrial material transfer, terrestrial biomass or butterflies.

## 1. Real measured source and same-site environmental units

- Exactly original 2022 Austrian genuine BMLUK FW1 national annual unfiltered river total-P and EEA2022 benthic normalized EQR at 68 exact physical stations, same original river Q/area, 5 original 3-degree geographic blocks, original EEA UWWTD annual reported 2022 P and N on river model source graph within 50km UPSTREAM, plus downstream noncausal negative controls. GZÜV official PKINT identity crosswalk revealed one country-vs-EEA geolocation discrepancy of 954.614m at ATFW54110117, excluded BEFORE terrestrial source data, leaving **67** (original geographic regions 15_3:9,15_4:14,15_5:21,16_4:13,16_5:10), 52/67 eligible N/P upstream report-positive and 15 with no matching eligible reported source (not chemically clean).
- Original BMLUK 2022 genuine FW1 unfiltered river TP XLSX SHA256 17de22544e9d043f83a928154847fb3fcd938894ce7007cfd06c000444b6eeb2; original EEA2022 2475 site BQE river data Git blob 325f199e985d9862ffee5b0edb205c61a7be3833. ESA official WorldCover2021 v200 10m public class COGs, CC BY 4.0, original fixed 1km radius around the actual 2022 monitored WATER site, class fraction crop40/built50/tree10/grass30, and 250m source context. All **67/67** official source landcover circle quality gates pass with >=20,000 real 10m valid landcover pixels and >=1,000 near pixels; no sites selected on EQR/TP/P sign.
- These WorldCover2021 class fractions are LAND USE COVARIATES for water-site landscape, not proof of a 2022 CHANGE in natural riparian vegetation. The separate genuine S2 L2A 2022 near-250m and outer250–500m same-natural-class optical NDVI source smoke STOP remains 5/10 < preregistered 6/10. **No NDVI–water biological response coefficient was estimated.**
- In all new 67-site fits, evaluation is complete out-of-geographical-region leave-one-out (5 blocks), train-only numeric scaling and Ridge10; ecological conditions in different river reaches may have spatial autocorrelation. Five original blocks are too few for reliable conventional cluster-robust confirmatory statistical inference. All fixed-errors region resampling intervals are descriptive only.

## 2. Authentic ESA2021 landscape confounding test: river chemistry and benthic ecology

[Original signed full67 landscape source and predictive model #43888461430 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43888461430). Compare original B = log1p river annual Q, log1p catchment upstream area, lat, lon, then B+4 true local WorldCover classes, then original reported annual upstream P/Q index.

| SAME67 out-of-region response and incremental predictor | Added heldout MSE gain | Regions improving |
|---|---:|---:|
| Actual 2022 river lnTP: 2021 near-field WorldCover landcover → B | **+30.7749%** | 5/5 |
| River lnTP: upstream reported 2022 P/Q AFTER B+landcover | +4.5822% | 3/5 |
| River lnTP: downstream reported P/Q placebo AFTER B+landcover | −66.8876% | 2/5 |
| 2022 benthic EQR: near-field WorldCover landcover → B | **+26.5822%** | 5/5 |
| Benthic EQR: actual river lnTP AFTER B+landcover | **+23.5179%** | 5/5 |
| Benthic EQR: upstream annual eligible P/Q AFTER B+landcover | **+15.0307%** | 5/5 |
| **Benthic EQR: upstream P/Q AFTER B+landcover+actual river lnTP** | **+8.3883%** | **5/5** |

Source reported mass dilution has mathematical equivalent 0.03170979 × official annual P tonnes/year / annual mean river Q (m³/s), **not** measured WWTP phosphorus actually reaching the receiving-station water. 2021 local cover does NOT capture full upstream 50km farm/urban modifications, sediment, multiple sources or longitudinal landscape change.

Original EQR source EEA normalized quality index is NOT invertebrate biomass, individual counts or even necessarily species richness. River 2022 TP is annual mean from all nutrient sources, not same-day benthic sample and NOT source-specific sewage contribution.

## 3. Binary WWTP-report presence vs magnitude (posthoc negative control)

[Original same67 source #43889331093 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43889331093), preregistered as a NEW post hoc stress only after observing the landscape result, no source population change. Benthic EQR baseline always includes original water lnTP and 2021 landcover.

| EQR extra predictor beyond real river TP+landcover | Out-of-block MSE improvement |
|---|---:|
| eligible UPSTREAM N/P positive-source report presence indicator | **−3.4181%** |
| original upstream annual reported P mass t/year | **−1.2651%** |
| original upstream annual reported P/Q potential concentration | **+8.3883%** |
| **P/Q AFTER report presence** | **+10.4405%** (5/5 regions) |
| report presence AFTER P/Q | −1.1015% |
| actual measured river lnTP: P/Q after landcover+report presence (separate outcome) | +1.9698% (3/5 regions) |

This excludes the simplest interpretation that source P/Q prediction merely codes the existence of an eligible reported WWTP. It does NOT exclude nonlinear river hydrology, source detection bias, treatment clustering, land use, river physical habitat, persistent long-term nutrient histories or other pollutant stresses. The P/Q after presence is not a measured mg P/L impact.

## 4. Q-only nonlinear hydrologic falsification

[Original same67 Q-only curves source #43890241923 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43890241923). Pre-specified FIVE Q-only covariate shapes: ln(1+Q)^2, ln(1+1/Q), both, ln(1+0.03170979/Q) unit pseudo source, and Q-inverse source-only 99th-percentile-winsor transform. All source regions and original background unchanged.

| Added to EQR model already adjusted for real lnTP + near-2021 landcover | Mean out-of-region MSE gain |
|---|---:|
| Q-only square | +1.75% |
| Q-only inverse | +1.21% |
| Q-only square+inverse | +2.04% |
| **UPSTREAM reported P/Q AFTER Q-only square+inverse** | **+13.153%**, 5/5 regions |
| UPSTREAM P/Q after inverse Q alone | +12.351%, 5/5 regions |
| Same downstream Q-normalized P source after inverse Q | −1.989% |

Thus an especially simple nonlinear Q-only covariate does not recover the observed reported P/Q indicator's out-of-geography EQR information. However all fits use the SAME source 67 stations, and the Q-only bases are only 5 arbitrary smooth alternatives, not a complete control for river habitat, Q regime, discharge variability or causal identification. For actual measured river lnTP, P/Q extra gain after Q controls ranges about +5–10% but has opposite-sign heldout region effects (only 3/5 regions improve), hence not a stable source→chemistry effect after geographic covariates.

## 5. Actual published reported 2022 source N vs P (nutrient specificity falsification)

[Original same-source N vs P 2022 run #43890958636 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43890958636). The official 2022 UWWTD source pool reported N and P from the SAME eligible facilities, positive-source presence status overlaps for **all 67/67** real sites. Source upstream reported annual N vs P t/year Spearman **ρ=+0.970698**; N/Q vs P/Q Spearman **ρ=+0.945272** and Pearson +0.899428.

| Benthic EQR: same source67 background + actual water lnTP + local 2021 landcover | Out-of-region added MSE gain |
|---|---:|
| Reported 2022 upstream P/Q only | +8.3883% |
| **Reported upstream N/Q only** | **+12.1991%** |
| P/Q AFTER N/Q | **−1.3166%** |
| N/Q AFTER P/Q | **+2.8979%** |
| Raw upstream P annual mass | −1.2651% |
| Raw upstream N annual mass | +0.9431% |
| DOWNSTREAM P/Q placebo | −13.2084% |
| DOWNSTREAM N/Q placebo | −2.5515% |

In stream lnTP itself, upstream N/Q adds +6.3701%, P/Q +4.5822%, adding P/Q after N/Q −3.3951%, and adding N/Q after P/Q −1.4577%. Therefore there is NO empirical identification of a specifically phosphorus-driven biological mechanism. The remaining signal is a **composite source reporting/dilution and potentially multiple environmental-stress indicator**. It could reflect N, P, co-transported organic chemicals, habitat/morphology, urban infrastructure, land-use selection, and river long-term histories. No actual river TN was obtained in FW1 sheet and these source N/Q figures are not a genuine measured receiving-water nitrogen concentration.

## 6. What this means for the user's one-paper question

**Strong empirical AQUATIC part (source-unit consistent):** genuine official 2022 upstream WWTP reported eligible N/P mass and real hydrological network within 50km, plus actual measured 2022 river TP and EEA benthic normalized EQR, same real water-station geography. In independent France73 2018/2022 same-site river TP/I2M2, out-of-parent river TP→EQR predictive signal was also strong. Those are separate populations and different methods, not treatment replicates.

**The aquatic mechanism is still only predictive:** upstream P/Q and N/Q source indicators predict EQR after local landcover and actual river TP, but the specific chemical mediation path is not identified and the same source N and P are highly correlated; no source-specific delivered chemical fraction, contemporaneous in-stream nutrients from the outlet plume, pre/post facility treatment, same flow event, independent reference river, full upstream landcover, river sediment/temperature, or real WWTP/source intervention.

**The terrestrial biological branch is NOT supported as a same-source response:** original Austrian 67 STAC 2022 scenes exist, but natural-cover matched same-source S2 NDVI smoke FAILED 5/10 against fixed 6/10; WorldCover2021 used here is only historical environmental covariate. EU H7 basin-level plant nutrient EIV NUTRI and butterfly host nutrient HNS covary spatially but genuine 2022 stream TP→those terrestrial scores showed no stable region transfer, and both are trait-affinity reference scores, not actual plant P uptake, butterfly food intake, biomass or ecological service. A manuscript that claims a completed wastewater→river→plants→butterflies food-web subsidy is not defensible.

**Publication decision:** a narrow **aquatic ecological source/condition predictive study** investigating conditional source load/dilution/habitat proxies, true measured nutrient and benthic EQR is now supported as a focused observational manuscript. A single cross-compartment wastewater ecological subsidy or causal plant–butterfly response paper is **NO-GO** until new same-source terrestrial data and experimental/quasi-experimental source control exist. A possible shared theoretical framing is that anthropogenic point sources provide **multi-component environmental footprints**, but do not convert the positive water predictors into measured water–land nutrient transfer.

## Reproducibility and review state

Original signed source/analysis runs:
- WorldCover full67: https://github.com/zuizui0223/gesui/actions/runs/43888461430
- source-report presence vs mass: https://github.com/zuizui0223/gesui/actions/runs/43889331093
- nonlinear source Q-only placebo: https://github.com/zuizui0223/gesui/actions/runs/43890241923
- source N vs P chemical-specificity stress: https://github.com/zuizui0223/gesui/actions/runs/43890958636
- natural-cover class-matched S2 2022 riparian SMOKE fixed source STOP: https://github.com/zuizui0223/gesui/actions/runs/43885473411

Source code: analysis/20261116_AT67_riparian_satellite/WorldCover2021_67_WWTP_PQ_confounded_streamTP_benthic.py, WorldCover_AT67_WWTP_report_presence_vs_Pdose.py, WorldCover_AT67_upPdivQ_vs_Qnonlinear_stress.py, WorldCover_AT67_NQ_vs_PQ_reported_WWTP.py. Positive and negative, fold-by-fold source results in original signed Actions JSON; original raster COGs and third-party published government XLSX/GeoPackage not rehosted. Draft PR118 is stacked onto draft PR117. main and chocho remain unchanged. No actual ecological or academic manuscript submitted.
