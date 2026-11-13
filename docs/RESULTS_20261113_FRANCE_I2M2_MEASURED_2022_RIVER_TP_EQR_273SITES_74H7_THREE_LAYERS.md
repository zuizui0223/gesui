# Same-site measured 2022 river phosphorus and benthic quality: a physically matched, catchment-scale three-layer macroecology audit

**13 November 2026 | New original field water-quality–benthic result, source geometry frozen before the TP–EQR test | All analyses exploratory/posthoc relative to prior wastewater–EQR knowledge | Draft/unmerged.**

## What changed from the older "wastewater geography predicts benthic EQR" finding?

Original PR #68/#69 found that **registered/modelled wastewater flow to a H7 basin outlet** provided **+11.1%** source-conditional prediction of the measured EU benthic quality index across 454 H7×period rows, but the predicted signal shrank to **+1.4%** when only confirmed upstream, distinct river-reach outfalls were counted, and **+1.87%** when 352 confidently ordered same-reach upstream WWTP×EEA field-site pairs were added. Downstream facility *geography* was similarly predictive despite not being a simple direct upstream water source. Thus strong modeled source geography was not a valid proxy for *receiving river nutrient concentration*.

Here we used **genuinely measured river TP** at a physical station on the same river reach and calendar year as genuine **2022 benthic quality EQR**, then compared corresponding catchments to opportunistic GBIF botanical and Lycaenidae observations. This is a qualitatively better measurement of the environmental water-quality mediator, while still **not identifying sewage-specific P** (agricultural/urban runoff and natural sources are major confounds).

## Provenance and exact source-frozen sampling, before new numeric TP/EQR effects

- Peters et al. original source real measured 2022 river TP station medians from DOI `10.25400/lincolnuninz.30546371.v5`, n1,291 stations meeting ≥4 positive 2022 measurements across ≥3 months, derived original source Git blob SHA1 `a75a4cf140de64fd6038079a3fac6756a461f3f8`. These are measured 2022 receiving-stream TP concentration **mg/L**, NOT treatment plant 2022 reported P tonnes/year.
- Original actual EEA 2022 benthic macroinvertebrate biological EQR n2,475 with authentic native `HYRIV_ID`, latitude/longitude, `procedureClassificationSystem`, source Git blob SHA1 `325f199e985d9862ffee5b0edb205c61a7be3833`. Real benthic quality index, not EPT biomass or richness.
- A true source-only stage matched the **identical original HydroRIVERS reach `HYRIV_ID`**, original station points at a **predeclared ≤250m** spatial separation, and one-to-one nearest source station pair by independent geographic distance/lexical ties. **TP and EQR numerical values were NOT accessed to pick source pairs**. Across ALL EEA2022 original sites, 508 WQ TP stations share a `HYRIV_ID` with EEA and 473 unique WQ–EEA pairs remain within 250m (471 FR, 1 NO, 1 AT).
- Original prior real multiclass macro panel fixed **232 H7 native EU/US basins ×454 basin-period rows** from earlier PRs; **227** of these H7 had 2022 EEA EQR. Restrict to original, not selected by new TP–EQR sign: **273 physically close real TP/EQR station pairs in 74 independent H7 native subbasins, all FR, all with SAME biological assessment method `I2M2`**; 239 source site pairs have literally identical coordinates.
- In same source-fixed 2022 H7 population, actual **TN** measurement had only 1 site/1 H7 overlap, so a pan-European or even French TN–EQR analysis is NOT feasible. Do not pool with TP to achieve a positive nutrient effect.
- Original geographical source gate before opening TP/EQR numbers [`SOURCE ONLY SUCCESS #43370128446`](https://github.com/zuizui0223/gesui/actions/runs/43370128446), exact 273 original matched-source manifest SHA256 **`448add2d6850cbc876c3e6b90bca278c6b711e64f9d363071e93bc0c18b84611`**.
- Original 74 H7 correspond to **24 genuine H6 parent hydrological basin groups** (official source parents required explicit previously documented GEOS repair on invalid original HydroBASINS L6 geometry in the broader original 232 study). The 5 source-frozen independent GroupKFold hydrologic block sample sizes are **15/15/15/15/14 H7**. Actual 74-fold IDs saved before numeric coefficients, SHA256 **`294eb69b76cd4cf40f1022a5e83278529b1d2bbceac9a38931004610a23eb052`**. [Outcome-blind group feasibility/fold freeze SUCCESS #43370462196](https://github.com/zuizui0223/gesui/actions/runs/43370462196).

**Scientific inference grain:** one 2022 **H7 watershed** contributes one original site-equal mean normalized 2022 EEA I2M2 EQR and one mean of `log(actual TP mg/L)` at exactly the same **273 source-matched physical monitoring stations**; individual stations are not 273 independent hydrological ecosystems. Same 74 H7×2020–22 mapped onto existing original source whole-H7 GBIF plant S50 and Lycaenidae five-family reporting composition. Those terrestrial observations are **not field-sampled at the TP/EQR station**.

## Primary genuinely observed river-quality/benthic result

| Original physically compatible TP–I2M2 EQR comparison | Observed result |
|---|---:|
| TP 2022 at exactly source-matched river monitoring sites | **273 source station pairs** |
| Same original 74 H7 watershed Spearman, `mean log TP` vs mean EQR | **ρ = −0.69833** |
| Source 273 true physical station Spearman TP vs EQR, descriptive only | **ρ = −0.54587** |
| Source `within H7` demeaned original station TP vs EQR regression slope, in basins with at least 2 sampled stations | **−0.11634 EQR per unit ln(TP mg/L)**; 2,000 H7 cluster bootstrap descriptive95 **[−0.16043,−0.07206]**; 42 H7 /241 original site pairs contribute |
| Adjusted EQR vs actual 2022 TP, per **1 SD mean log TP** (SD 0.76366) | **β = −0.10076 EQR**, 24 H6-parent cluster descriptive95 **[−0.13610,−0.06541]** |
| Baseline B0 after original upstream cropland `crp_pc_use`, upstream urban `urb_pc_use`, modeled natural outlet flow, site-modeled flow, H7 area and original EQR/TP station count | **H6-blocked 5fold MSE = 0.029539** |
| Add actually observed **2022 in-stream TP** to the fixed baseline | **MSE = 0.022809**, incremental **+22.7826%** predictive gain |
| 2k H6-group source-bootstrap descriptive95 5fold gain | **[+9.1618%,+37.2604%]** |

All basins are French and all quality methods are native **I2M2**, avoiding cross-national EQR method pooling. The TP 2022 original observed same site range **0.0075–0.72 mg/L** and EQR range **0.0119–0.9860**. The same-year median TP sampling is spread over months, not necessarily benthic index sampling date.

[Original 2022 station co-location, measured TP→I2M2 EQR and native 74-H7 three-layer source models SUCCESS #43370847415](https://github.com/zuizui0223/gesui/actions/runs/43370847415). Original machine-readable H7 74 and station-paired 273 source CSV plus report JSON are in run artifacts.

### Independent numerical validation and stricter transfer

A second independently implemented **NumPy normal-equation Ridge10** exactly reproduced all **three outcomes ×3 core model specifications = nine** prior sklearn H6-blocked MSE values, all 273 original site rank correlations, and original same-H7 within-basin slope.

**Leave ONE of 24 original H6 hydrological parent groups out at a time (24fold LOHO, posthoc):**
- Actual river TP→EEA 2022 I2M2 EQR MSE incremental improvement **+22.271%**, H6-parent bootstrap descriptive95 **[+7.592%,+38.296%]**;
- Original modeled HydroWASTE **H7 basin-outlet** flow index→EEA EQR improvement **−2.722%** (worsens);
- Adding modeled outlet flow **after actually measured TP** worsens EQR prediction **−3.270%** (H6 95 [−8.685%,+1.337%]); conversely adding **real TP after original WWTP index** improves EQR **+21.857%** (95 [+7.470%,+37.364%]).
- **Physical source bridge check:** Original HydroWASTE modeled H7 outlet flow index vs actual measured 2022 TP basin means raw Spearman **ρ = +0.3571**; but after upstream cropland/urban/outlet flow and river flow controls, modeled W predicts measured TP only **+2.713%** under original 5fold (descriptive H6 bootstrap [−5.21%,+9.33%]) and **+0.033%** under 24-H6-parent LOHO (descriptive95 [−8.10%,+6.17%]). **Not strong evidence of registered wastewater causing observed river TP**.

[Independent 9-model NumPy arithmetic, 24 H6 parent leaveout, and wastewater→real measured TP source-specificity SUCCESS #43371398070](https://github.com/zuizui0223/gesui/actions/runs/43371398070). This is independent code, **not independent field data**.

### Positional uncertainty sensitivity — stricter coordinates AFTER observing the primary TP–EQR sign

Retained only original identical `HYRIV_ID` river station coordinate pairs within **≤10m**, not entire 250m support; reaggregated EQR/TP from the same previous source pairs, no new observations or changed H6 grouping method. This **posthoc** geometry sensitivity retained **255 station pairs, 72 H7, all 24 H6 parents**, and found H7 TP–EQR **ρ = −0.67492**, within-H7 TP–EQR station slope **−0.12517**, and environmental baseline + true TP H6 5fold MSE gain **+18.946%** (2k H6-bootstrap descriptive95 **[+4.162%,+34.601%]**). The strong negative TP–EQR relationship is not generated solely by joining stations that are 100–250m apart.

[Original strict ≤10m sensitivity SUCCESS #43371662688](https://github.com/zuizui0223/gesui/actions/runs/43371662688).

## How the terrestrial layers behave: not a serial ecosystem service chain

Original genuinely native same-H7×2020–22 GBIF response fields were present in original 74 H7:

| Different response observation, same 74-H7 temporal/spatial unit | Add measured river TP (H6-group CV MSE gain) | Add modeled registered HydroWASTE H7 outlet flow | 24 H6 leave-parent-out sensitivity |
|---|---:|---:|---|
| **True measured I2M2 benthic EQR** | **+22.783%** | −0.890% | TP +22.271%; W −2.722% |
| **Indexed rarefied plant speciesKey S50 (not field flora)** | −3.761% | **+17.493%**, descriptive H6 bootstrap [+5.688%,+24.925%] | TP −3.938%; modeled W **+15.438%** but 95 [−0.129%,+23.838%] includes zero |
| **Lycaenidae share of five original GBIF butterfly family recorded occurrences (not real butterfly abundance)** | +0.229% | −6.648% | TP −0.472%; W −3.605% |

Plant S50's modeled W/outlet source association in the original selected **74-FR-H7 environmental monitoring sample** is worth retaining as a reproducible **posthoc observational predictor**. It has no direct correspondence to measured river TP or actual same-site flowers. It may measure different human habitat/vegetation/GBIF observation geography despite original recording-opportunity adjustments. The broader PR #65/#67 all-EU n454 showed plant S50 modeled W near-zero/worsened CV, so the n74 selective French panel result cannot be claimed globally robust or sewage nutrient-driven ecological benefit.

Lycaenidae record composition did not show stable improvement with the actual river TP proxy. Butterfly index is not adult visit rate, larval food availability or realized pollination.

## Honest science synthesis

**Measured ambient river TP ↔ benthic biological condition is the strongest direct environmental association in this entire `gesui` macro chain so far**: physically co-located original 2022 water chemistry + benthic quality, same country/method, fixed H7 catchment means, same-H7 spatially structured H6 group CV and independent numeric rechecks. It is **NOT** novel proof that WWTP effluent specifically caused that TP, that TP caused the benthic decline, or that butterfly/pollination ecosystem services were propagated through the same units.

**Original chain edge measurement levels and conclusions:**
- `WWTP actual source-reported P/N → actual in-stream TP/TN`: **NOT corroborated broadly on shared monitoring H7 units**. The registered HydroWASTE static H7 outlet flow index does NOT substantially improve actual measured 2022 TP under strict 24-H6-parent leaveout (+0.03%, highly uncertain). Original reported annual true EEA N/P source was positive for just 56/454 older H7×periods, country-selective, and no full treatment status verification.
- `Observed in-stream TP → river benthic EQR`: **strong negative observational co-response** in the French 74-H7, 273 same-station source frame, +22.27% strict LOHO predicted EQR gain, source H6 cluster negative beta.
- `TP/Benthic condition → H7 plant index and Lycaenidae family index`: measured TP is not a strong out-of-region predictor of observed H7 GBIF terrestrial layers in this matched sample, and direct realized resource/flower/host/pollinator information is absent.
- `Aquatic emergence/consumer subsidy → pollination/seed set`: no same-H7×2022 measured functional data. Existing Oregon 17-stream adult emergence chronology and Canadian Dryad riparian studies are independent locations, not joinable source rows.

Source co-location ≠ wastewater-specific causation. No historical source/reporting imputation, no synthetic multilevel mediation, no "ecosystem service positive chain" declaration. Future ecological identification should source actual treatment-plant effluent 2022 P/N plus receiver discharge chemistry and H7 source-network operations at the SAME French exact monitoring stations, and independently measured host/pollinator endpoints.

Retain all original PR #64–69 positive and negative outcomes. No merging to `gesui/main`, `chocho`, or replacing original source datasets; no eBird.
