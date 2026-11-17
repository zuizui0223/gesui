# 2026-11-17 | REAL same-2018 EEA benthos × Copernicus Riparian Zones: source PASSES, geography association and cross-boundary nulls

**Verified source-based ecological result, not a causal wastewater→terrestrial mechanism.** This new common-site comparison was designed after original 2018 PR #17 water EQR results existed, so it is exploratory/posthoc, NOT an independent preregistered test. The source-only sample and full quality rules were frozen before attaching RZ land cover to actual 2018 EQR for this comparison.

## Authentic sources, final grain and signed runs

1. **Original 2018 2018 true EEA benthic station × reported-measured wastewater N/P hydrologic source**: 2,639 pre-EQR river station candidates / 15 countries / 204 upstream positive reported sources, original 2018 source-year frozen original site CSV SHA256 `5a074f5182711c4ab74ef9bdacfa3918e9fcb434dd0b75c67072a0f06459c3e1`; authentic later 2018 normalized benthic EQR 2,634 sites, source bytes SHA256 `2a1a72d43656aca588a832dfae13860054c01a651e76bf25f984c5a374709b02` (PR #17).
2. **Independent genuine terrestrial source**: Copernicus/EEA *Riparian Zones LCLU 2018*, 56 publisher Level4 land-use classes, original polygon REST queries, reprojected EU LAEA EPSG3035, fixed 0–250m circle and 250–500m annulus from each exact original river monitoring site, actual mapped polygon union area. Not direct current biomass or actual riparian groundwater nutrient uptake. Official product https://land.copernicus.eu/en/products/riparian-zones/rz-land-cover-land-use-2018 (DOI 10.2909/2afca4ec-76e2-4155-b4e6-7460f9f6ae01).
3. **Original all-site land source [Actions #44004357868 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/44004357868)**, source full 2,639-site derived artifact **14293433193**, original ZIP SHA256 **00b1b72b3812d9c541740d3542da82cf4717265ed014df3a6852508865ce3f6b**. Frozen >=70% classified polygon coverage BOTH zones: **986/2,639 technical map-eligible sites**, **14 countries,49 geographic 3° blocks**, **175 positive eligible reported-source upstream sites**, far above full source GO thresholds 200/6countries/10blocks/20positive.
4. Additional prior frozen >=10% actually mapped *terrestrial area* in BOTH zones (exclude wholly aquatic/mostly water): 8/986 fail -> **978 source-eligible**. Exact original 2018 EQR join leaves **976 distinct river monitoring sites, 13 countries, 49 3° geographic blocks, 169 positive reported-upstream source sites**. No many-to-many joins, no missing numeric 2018 EQR imputed. 5 original pre-outcome geographic folds reused (heldout n=220/178/121/242/215).
5. **Genuine ecological prediction [GitHub Actions #44005855019 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/44005855019)**, derived artifact **14293429762**, original ZIP SHA256 **a898195e9136a5d801b2ad1e29aee06c15a613ea8638d15aa9c5ea0d2a12cb32**. Original actual 2018 landcover/source 2018 benthic EQR and source reports on identical physical station; no water-to-terrestrial resource flux is implied.

## PRIMARY: original year, real EQR, fixed five geographical heldout folds; Ridge alpha=10

The terrestrial predictor labels refer to mapped **woodland/forest + grassland + heath/scrub + wetland** class groups from Level1 3/4/5/7; **forest plantations and managed grasslands included**, so NOT a measurement of pristine vegetation or production. Category 8 (open water) removed from terrestrial fraction denominator. Baseline B = lat/lon/log river Q/log upstream area + country + country×original BQE method. LandContext L adds natural-like vegetated, built, cropland fractions **250–500m** from actual water site. Near−outer `D` = comparable noncrop vegetated fractions at **0–250m minus 250–500m**, not plant foliar N or abundance. All fold scaling training-only.

| 2018 benthic EQR: additional predictor vs previously frozen baseline | heldout site-weighted MSE gain | 49-region equal MSE gain | descriptive95 fixed observed region errors |
|---|---:|---:|---:|
| L = broader 250–500m landscape composition **after B** | **+3.423%** | **+9.434%** | **[+1.876,+17.005]** |
| D = near–outer local vegetation class fraction **after L** | **−0.248%** | **−0.349%** | **[−1.257,+0.493]** |
| 2018 original eligible **UPSTREAM P t/year after L** | **+1.363%** | **+1.241%** | **[−0.703,+3.158]** |
| 2018 original eligible **UPSTREAM N t/year after L** | **+1.424%** | **+2.272%** | **[−0.295,+4.486]** |
| D added AFTER L + upstream P | **−0.106%** | **−0.264%** | **[−1.205,+0.613]** |
| UPSTREAM P added AFTER L + D | **+1.503%** | **+1.325%** | **[−0.685,+3.302]** |
| DOWNSTREAM P placebos after L | **+0.114%** | **+0.117%** | **[−1.685,+1.853]** |
| UPSTREAM P × D interaction after L+P+D | **+0.250%** | **+0.096%** | **[−0.381,+0.715]** |

**Positive result:** genuine landscape composition at the *same river observation point* improves 2018 aquatic macroinvertebrate EQR prediction beyond original country, assessment method, hydrogeography. This may be due to shared habitat or catchment pressures, not demonstrated nutrient transfer, and is not terrestrial biomass production.

**Not supported:** source-specific wastewater terrestrial vegetation response or a directional water→land material/biotic ecological pathway. Neither the closer-river vegetation contrast, its source P interaction nor P-specific source prediction after broad landcover has an original 49-geoblock fixed-loss uncertainty interval clearly above zero. P and N source reports covary and do not isolate P or N causality. Downstream-placebo after L is uncertain, not a conclusive negative control.

## Explicit cross-boundary terrestrial-as-response falsifier

**Outcome** original mapped `D=near−outer` noncrop vegetated fraction at 2018 stations, predicting its spatial variation with broad land cover and hydrology plus reported upstream source: **+0.018% heldout site-weighted gain**, region-equal **−0.430%**, descriptive95 **[−1.201,+0.236]**. Downstream source similarly adds −0.092% site, region-equal −0.081% [−0.358,+0.162]. No evidence that a 2018 WWTP-report load index improves 2018 *mapped river-proximate terrestrial composition* prediction at this spatial grain. This is not a biological null for actual riverbank groundwater/foliar N or emergent insects, which are unmeasured.

## HIGH priority selection bias: source mapping opportunity is non-random

**Do not interpret +3.423% in 976 stations as the original 2,634-station population without qualification.** The official Riparian Zones polygon product covers selected river corridors, not all EEA monitoring points. The original RZ source gate passes **175/204 (85.784%)** sites with eligible reported upstream P **vs 811/2,435 (33.306%)** reporting-negative/uncaptured original sites. Source map-eligibility odds ratio ≈**12.084**. After the distinct frozen terrestrial-land-area 10% requirement, **171/204 vs 807/2,435** remain pre-EQR source-eligible; 2 further reported-source-positive sites lack valid 2018 EQR, leaving 169 reported-source-positive in the actual 976 biological sample. The physical river-network selection/reporting process and national procedures differ sharply; no-report is NOT clean water.

Spatial model results are conditional on this selected study frame. The 49-region equal confidence limits are **descriptive resampling of fixed out-of-fold errors**, not confirmatory independent-study CIs, and do not account for site-selection propensity or model refitting. Good landscape EQR prediction does not mean ecosystem services have increased.

## Manuscript relevance and what remains

The same-year, same-physical-site observational bridge is a concrete macro-ecology advance beyond earlier EU447 trait-habitat geography and France73 river chemistry cohorts, but **the new terrestrial endpoint is LCLU class composition**, which could reflect watershed urbanization, agriculture, morphology and observer/inventory differences. An aquatic→terrestrial sewage subsidy/food chain is still **not identified**. The current narrow paper can discuss *multicomponent anthropogenic source pressure embedded in landscape mosaics* and distinct aquatic biological status, preserving negative tests and selection bias.

Newly posthoc country-heldout and source-bias falsification (not used to choose this primary model) is being separately implemented as `analysis/20261117_RZ2018/eea2018_RZ2018_country_purge_source_selection_QA.py`, and will be reported only after successful original Action. The original 2022 Austria S2 5/10 smoke STOP, 2022 HRVPP official-public STAC STOP and 2018 plant/host trait scores from a *different EU H7 source population* are not merged into the same causal graph.

Reproduction: `docs/PROTOCOL_20261117_EEA2018_SAME_RIVER_RIPARIAN_LCLU_SOURCE_GATE.md`, `docs/PROTOCOL_20261117_2018_RZ_LCLU_BENTHIC_EQR_SPATIAL_VALIDATION_IF_SOURCE_GO.md`, full-source `analysis/20261117_RZ2018/eea2018_original_riparian_lclu_source_pilot.py`, model `analysis/20261117_RZ2018/eea2018_RZ2018_benthic_EQR_geoblock_model_if_fullsource.py`. No original EEA/CLMS raw polygon or unlicensed occurrence source replicated in repository; signed derived output artifacts remain on GitHub.


## 2026-11-17 posthoc full-country transfer: source-verified caution about continent-wide generalization

[Original-source Actions **#44006778058 SUCCESS**](https://github.com/zuizui0223/gesui/actions/runs/44006778058), derived artifact **14295235790**, ZIP SHA256 **d3af4869597f04c10b5f27a8189171e5a3dc5b4989964fc77575b54c06d5472c**. Independently reconstructs genuine 2,639 RZ source rows and 2,634 genuine numeric 2018 EQR sites; exactly 976 site-year matches, 13 source country groups; 2018 original original source/biology SHA locked. Whole-country disjoint test (train on all other original countries, test unseen full country) uses original Ridge10, same source variables and train-only standardization. This was designed AFTER 5 geographical fold results; strictly *posthoc*, not new independent confirmation.

| Adding to benthic EQR predictive model | All 13 country-purged test rows site-weighted MSE gain | Country-equal mean-MSE gain | Descriptive95 of resampled FIXED 13 country errors | Positive countries |
|---|---:|---:|---|---:|
| Original mapped 250–500m landcover after original B hydrogeo controls | **+2.399%** | **+1.861%** | [−3.270,+10.119] | **6/13** |
| 2018 eligible-reported upstream P after landcover | +3.239% | +2.123% | [−5.976,+7.181] | 8/13 |
| 2018 eligible-reported upstream N after landcover | +3.518% | +2.101% | [−5.964,+6.769] | 8/13 |
| 0–250m vs 250–500m mapped vegetated contrast after broader landcover | +0.055% | −0.169% | [−1.441,+1.142] | 8/13 |

**Science downgrade:** the positive 49-geoblock fixed-error interval for broad landcover **does not translate to a zero-excluding whole-country interval**. The 13 national groups are extremely unbalanced (e.g., FR 455 sites, SE 2, SI 7, HR 8), with original reporting effort/methods/geography varying. Broad landscape composition predicts aquatic status on familiar source countries across geographic blocks, but external country-transfer universality is **NOT established**. The 2018 report-source incremental gains and direct near−outer contrast remain country-heterogeneous and uncertain. Do not claim wastewater caused any terrestrial reorganization from these source proxies.

**Independent original source-selection QA validated** 2018 RZ polygon eligibility **175/204 = 85.784% upstream eligible-report-positive** against **811/2,435 = 33.306% no-eligible-report** (odds ratio **12.0838**), so the eligible subset is non-exchangeable even before modeling ecological responses. No regression approach can make RZ missingness fully ignorable without verifying external sampling/river hydrology. 2018 reference-year bank LCLU remains a landcover class composition (including plantations, managed grasslands), not field vegetation biomass, leaf N or aquatic insect emergence. There is no actual physical wastewater-to-terrestrial-trophic flux measured.

The signed source code is `analysis/20261117_RZ2018/eea2018_RZ2018_country_purge_source_selection_QA.py`, source-staging task `.github/workflows/gesui-20261117-RZ2018-country-purge-selection.yml`. The original fivefold geographic result is preserved rather than overwritten by this new validation.
