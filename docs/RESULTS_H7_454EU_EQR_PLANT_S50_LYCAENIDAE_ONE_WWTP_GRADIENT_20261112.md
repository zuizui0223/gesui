# Same-watershed macro comparison: registered wastewater geography versus aquatic benthic condition, observed flora diversity and butterfly composition

12 Nov 2026 | EU/NA source provenance, EUROPEAN SOURCE PANEL ONLY | Exploratory POSTHOC original-data results | unmerged Draft

## Why this matters

A single original HydroBASINS Level7 catchment × 2017–19 or 2020–22 period row is compared across THREE distinct response observations. The environmental predictor is a single static original modelled HydroWASTE cumulative outflow pressure to the basin OUTLET. The core question is not whether wastewater *causes a chain*, but whether a common anthropogenic source geography predicts aquatic and terrestrial biodiversity indexes similarly or differently.

This spatial analysis was designed AFTER the source-first 169-basin temporal co-response gave nearly zero correspondence; therefore results cannot be described as prospectively confirmatory.

## Original identical geography and source quality

- 234 original European HydroBASINS L7 river subbasins with real 2018 and 2022 EEA benthic original EQR monitoring, 468 potential basin-period records. Real GBIF original native polygons with plant indexed speciesKey 50-draw expected rarefied richness S50, and the real five accepted butterfly-family indexed counts with Lycaenidae log odds.
- Source-only covariate requirements yield **454 analyzable basin-periods in 232 distinct native L7 basins**, original 227 2017–19 +227 2020–22; 222 basins have both periods. Original 15 country codes; 299 of 454 rows are FR (65.9%). Original country and EQR method source labels retained, including 11 border-country-mixed and 12 method-mixed basin-periods.
- Original cross-source geographic/source payload SHA256 07fa9744cd4e9574587a7b25c4ed4fa7f3183e377cac61c28723a4f24a86c875; original source crosswalk SHA256 a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89. Original EEA Git blob for 2018 0f30bfbc3dac13f6a182c8a910e9899c04778322, 2022 325f199e985d9862ffee5b0edb205c61a7be3833.
- Original static HydroWASTE publisher genome = 58,502 registered plants, about 90.65% facility status NOT REPORTED. Outflow locations aggregated and routed along original Level7 NEXT_DOWN to the basin OUTLET. 2018 vs2022 true WWTP annual loads and discharge actually reaching interior original river station are NOT measured here.
- Original same source baseline for all 3 outcomes: country, original modal river biological assessment method, observed period, original local basin area, original monitoring station n, original Animalia/Plantae/five-butterfly-family GBIF indexed record opportunity counts, and original river DIS_AV_CMS (fully available at observed stations). Additional one source variable: log1p(original registered-modelled WWTP outflow / (source upstream area+1)). Entire source is observational and does not adequately control agricultural/urban land-use confounders.

Five deterministic H7 basin-hash (both period observations from same H7 remain in same fold), training-only predictor standardization, Ridge alpha 10, original equal basin-period MSE and 2000 H7-basin cluster resamples. This is NOT a geography/contiguous-region or country-heldout prospective evaluation.

## NEW original same-gradient 3-layer results

| Original measured or indexed response | Additional basin-ID-heldout MSE gain from static WWTP spatial gradient | Standardized 1SD source gradient effect, outcome original units (H7-cluster descriptive 95%) |
|---|---:|---|
| **actual aquatic EEA benthic EQR** | **+7.328%** | **−0.07806 EQR** [−0.10614, −0.04999] |
| Plant original GBIF 50-record rarefied indexed speciesKey S50 | **+0.127%** | −0.2255 S50 [−0.4609,+0.00994] |
| Lycaenidae original index log odds among five butterfly families | **−0.498%** | +0.01502 log odds [−0.03830,+0.06834] |

Corresponding standardized outcome effect: EQR **−0.418 SD/SD of original source**; plant S50 −0.115 SD/SD; Lycaenidae +0.042 SD/SD. But causal interpretation is NOT valid.

H7-basin bootstrap **descriptive** 95% gain: EQR **[−2.73%,+15.77%]**, plant [−2.07%,+2.63%], Lycaenidae [−1.52%,+0.49%], each crosses zero. Equal-country contribution weighting weakens the EQR gradient effect to **−0.0270 EQR per SD** [−0.06674,+0.01269], crossing zero. Plant corresponding country-balanced +0.2818 [−0.2038,+0.7673] and Lycaenidae −0.02661 [−0.1521,+0.0989]. Thus the source-pressure/EQR negative association is **predictively suggestive**, not strong universally replicated ecological impact.

Source-only original real three response fitting and results, original 454 native basin-period rows: [GitHub Actions #43215754868 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43215754868).

## Independent Python reimplementation and subgroup diagnostics (all POSTHOC)

A completely separate second original source-file computation reproduced **all 3×2 original baseline vs added-WWTP Ridge5fold MSE numbers** to numerical tolerance. No new biological survey or independent country validation was performed, only reproducibility.

| Fixed after-source subset | Original EQR MSE improvement from static WWTP | Original plant S50 improvement | Lycaenidae improvement |
|---|---:|---:|---:|
| Whole 454 observed basin-periods | **+7.33%** | +0.13% | −0.50% |
| 2017–19, n227 | **+14.27%** | −0.68% | −0.47% |
| 2020–22, n227 | **+3.47%** | +2.72% | −0.71% |
| France-only n299, 152 H7 | **+12.72%** | +1.06% | −0.14% |
| Other countries n155, 82 H7 | **+3.17%** | **−7.81%** | **−4.15%** |
| Exclude both country/method-mixed basin-periods, n442, 230 H7 | **+7.05%** | −0.13% | −0.48% |

**Independent original-saved numeric model reproduction and source-domain validation SUCCESS** [#43216359735](https://github.com/zuizui0223/gesui/actions/runs/43216359735).

For EQR the signed *incremental predictive improvement* was positive in every reported posthoc subset. However its magnitude varied markedly with time and country composition, and the subgroup results were inspected after the original effect was known. Plants and Lycaenidae lack comparable stability. These diagnostics are not 5 new independent experiments, since the same study sites and original observations are reused.

## Combined interpretation with the n169 genuine repeated-station analysis

The **temporal co-response** study [#43214039321](https://github.com/zuizui0223/gesui/actions/runs/43214039321), using ONLY the exact 1,023 same-station same-method EQR pairs across 169 source-eligible H7 basins (148 FR), found near-zero unadjusted EQR–plant Delta S50 rho −0.054, EQR–butterfly Delta log-odds rho+0.035, plant–butterfly rho+0.013. Adding plant and/or butterfly temporal change to adjusted EQR change models **worsened** H7-heldout MSE (−0.28%,−0.97%,joint−1.20%). Independent second script verified all ranks and CV values [#43216686503](https://github.com/zuizui0223/gesui/actions/runs/43216686503).

**Consistent cross-study evidence pattern:** The broad coarse source geography of registered wastewater plants may carry information about *spatially distributed river benthic quality*, while opportunistic terrestrial flora and butterfly index observations do not provide a stable same-basin temporal covariance with EQR over these fixed two windows. A detected predictive difference across observed response layers is **not** evidence of actual ecosystem service decoupling, because no emergence flux, actual plant hosts, butterfly pollination, or shared observed WWTP dose was measured.

## Main threats and required next macro source

- Same-country/method and GBIF observer sources concentrate in France; the EEA normalized EQR original 2018 vs 2022 protocol may have harmonization artifacts. Raw original site EQR n by method/year should be audited for stability before any publication claim.
- Static modeled outlet pressure is not actual nutrient mass at the original biological station; estimated treatment operation status mostly missing and confounders agriculture, urban cover, climate, urban human recording concentration remain. A strong EQR association may reflect those other pressures.
- GBIF rarefied richness is computed from indexed speciesKey records with original 50-observation draw, not field abundance; the five-family butterfly composition is not standardized butterfly visits/nectar supply/pollination. Great catchment coverage does NOT eliminate GBIF detection bias.
- Source response layer comparison does not create a serial SEM and cannot identify mediation.
- Priority next quantitative work: original Level7 **within-basin cropland and built-up fractions / streamflow**, original directed true station-specific effluent dose using HydroRIVERS if possible, country-heldout test or alternative continental river biomonitoring, and independent source vertebrate/insect functional outcomes. Keep repeated basin×period labels and signed nulls.
- Preserve original global WPI coastal plant results separately; no eBird; do NOT merge or overwrite old PRs or gesui/main/chocho. Remain Draft.

Working title: *Wastewater source geography predicts river benthic integrity more consistently than indexed terrestrial biodiversity across European watersheds*. This remains a descriptive **hypothesis and exploratory pattern**, not a verified causal ecosystem-services paper.
