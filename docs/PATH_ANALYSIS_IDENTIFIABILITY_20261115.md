# Structural path analysis: what source data actually identify

This document is a verification gate, not a fitted causal SEM. Avoid forcing unrelated study units into one ecological services mediation chain.

## Target DAG (hypothesis, not measured causal arrows)

Measured effluent fraction Q_eff/Q_down (volume/volume) and mass loading N/P (kg/day) are distinct and should not be confused; annual nutrient mass divided by discharge volume is concentration (kg/m3), not effluent fraction.

Exposure -> water chemistry (nutrients, Na, temperature, DO, pharmaceuticals) -> periphyton and aquatic benthos (index and measured productivity separately) -> measured **emergent adult insect biomass** -> observed riparian predator diet / abundance.

Exposure -> reach groundwater/overbank soil N or hydrology -> field host-plant cover and leaf chemistry -> observed larval abundance and host use -> butterfly community host N-affinity and reproductive success.

Exposure -> local wet mineral/organic patch -> verified adult butterfly proboscis contact and Na/N/DON intake -> adult aggregation (separate from larval occupancy).

Predator abundance/diet may influence butterfly mortality but prey switching makes sign ambiguous. Do not draw predator -> butterfly HNS as an observed mediation without evidence.

## Present measured-grain gate
| Path edge | Existing available data | Status |
|---|---|---|
| W (static H7 outlet geography) ↔ EQR | EEA H7 basin-period aggregated EQR, historic EEA source N/P/outfalls partly | Eligible for noncausal co-response association only, confounding and methods remain |
| W (static H7 outlet geography) ↔ butterfly HNS | 454 original common H7 periods, 232 H7, 90 H6 | Observational HNS association verified in attached ledger, not dose-response at adults |
| W ↔ plant source S50 and candidate hosts | Original GBIF plant species keys in same H7×period plus historic host catalogue | Potential resource opportunity; source recording not floristic plot/field host usage |
| Ambient chemistry ↔ benthic EPT role/tolerance | US EPA NRSA across 3 survey eras (PR81-83) | Predictive associations; no verified wastewater-origin dose |
| Aquatic larvae -> emergent insect dry mass | Not available at same 454 H7×period source | **STOP causal mediation** |
| EPT emergence -> riparian predator feeding | Some separate published source/isotope studies; not same H7×period | **STOP shared SEM**; evidence synthesis only |
| Soil/water N -> host plants -> actual larvae | No reach-matched root-zone/leaf chemistry and larval occupancy | **STOP causal mediation** |
| Na/organic patch -> adult butterfly puddling | No observed species×sex×substrate behavior at original H7 | **STOP direct effect** |

## Recommended modelling sequence
1. Same H7×period multi-response OOS diagnostics of EQR, plant S50 and butterfly HNS, shared covariates and H6 disjoint splits; label association not mediation.
2. Separate repeated H7 Δ diagnostics for butterfly HNS, plant S50 and EQR, preserve EEA physical site/method eligibility and 222/169 distinct cohorts. Do not infer changing source from static W.
3. Only fit actual SEM once nested plots/reaches/outfalls and source chemistry, emergent mass, predators, butterfly behavior are observed at compatible space/time with repeated independent discharge systems; use multilevel model with outfall random effects and spatial/season controls.

Reject pseudo-replication: 454 H7×period rows are not 454 independent basins; 232 unique H7 nested in 90 H6; 2 periods share sites. EQR is method-dependent ecological quality ratio, not benthic species biomass. GBIF records are not individual butterflies. Morphological, functional and reproductive outcomes should not be fabricated from these proxies.
