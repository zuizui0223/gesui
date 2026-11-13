# Spatial transfer does not equal country transfer: European freshwater EQR, terrestrial plant observations and butterflies along registered wastewater gradients

**13 November 2026 | Original source-pinned geographic validation | Exploratory, POSTHOC | Draft / unmerged.** This study reuses the genuine same 232 HydroBASINS Level7 river basins × 2017–19/2020–22 = 454 original basin-period rows from PR #65 and original HydroATLAS + EEA reported NP in PR #66. No new ecological sample is introduced. New evidence is **stronger geographic holdout validation**, not causal ecosystem service measurement.

## Question and pre-existing known observation

Prior source-verified result [HydroATLAS adjusted original 454 three-layer model #43227058486](https://github.com/zuizui0223/gesui/actions/runs/43227058486) used original H7 basin-ID hash heldout and reported:
- Adding **modeled HydroWASTE static outfall flow to H7 basin outlet** after original **upstream** cropland, urban cover, natural modeled flow, country, EEA BQE method, river-site mean discharge and original GBIF record-opportunity covariates improved EEA benthic EQR holdout MSE **+8.08%**; after LOCAL crop/urban controls +1.35%, uncertain.
- Adding original **2018/2022 report-source upstream annual N mass scaled by modeled receiver river flow** improved EQR holdout MSE another +0.80% (95% descriptive H7 bootstrap [+0.040%,+1.731%]); however actual reporting-positive source in only **56/454** original H7×periods and **38 H7** (France just 4 positive of 299 rows).
- The GBIF-original botanical S50 is **rarefied 50 indexed speciesKey observations** (NOT field flora), and Lycaenidae is a five-family indexed butterfly RECORDING proportion (NOT abundance/flower visitation/pollination). Both lacked stable extra static-WWTP prediction.

H7 hash heldout may overstate spatial prediction because neighboring basins can appear in train and test; original method/country distribution highly uneven. Thus this next test checks **source-native Level6 H6 hydrological parent group** and **leave-one-country-out**. No response source/year/covariate thresholds changed.

## Source-first geography and source failure: important to retain

Official HydroBASINS v1c Europe Level6 original archive SHA256 `6a4d9bcbadacbf6e32936970f84ec77c21809948085d79d1d4388964d607fa`, and official original Europe Level7 original ZIP SHA256 `b9b810f4db271f0fa6b39fecf01784041b948944b3c3f296a393766388e86b44`. Original source geometry **has 207 INVALID Level6 polygons**. Without any repair, **216/232 target H7** could be uniquely mapped to a containing H6, leaving 16 unassigned, so the ORIGINAL strict unchanged-sample source gate **FAILED**. This failure was preserved [Original strict no-repair SOURCE gate #43338063245](https://github.com/zuizui0223/gesui/actions/runs/43338063245). **Do not call the subsequent source-sensitivity fit prospectively prevalidated.**

Before inspecting the new H6 blocked numerical effects, a separate, transparently post-failure source-only technical sensitivity used **GEOS make_valid on only the 207 originally invalid H6 features**, never repaired or changed H7 biology/geometry or already valid H6 features. All invalid features repaired, restoring the 16 original missing H7 assignments. Every original H7 polygon assigned ≥95% area uniquely into a true source H6 parent. Preserved full **232/232 original H7 /90 unique original H6 parent IDs**, with parent-heldout fold groups 17/15/20/20/18 parents (60/90/123/106/75 original basin-period rows). Stage0 geo SHA256 `372b800b85abe2f395241987cea2000ff901cf5aff31bf1d1fb825304f2495fc`. Original source-only explicit correction [SUCCESS #43338311684](https://github.com/zuizui0223/gesui/actions/runs/43338311684).

## Actual 90 H6 parent heldout numerical results — identical 454 original H7×period biological rows

Train-only feature standardization, exact Ridge10, SHA256(H6 true parent ID) %5, every original H7 and both periods held out with its entire true native H6 parent. Source environmental "upstream" means **HYDROATLAS original upstream** cropland/urban fractions plus `log1p` natural modeled basin outlet discharge, controls for the same original country, BQE method, observational period, actual EEA river station DIS_AV_CMS, H7 local area and GBIF recording opportunities. Original static W is cumulative modeled HydroWASTE flow pressure to BASIN OUTLET; actual source N is original reported annual upstream mass divided by station river modeled discharge converted to hypothetical mg/L-equivalent; no actual measured TN/TP.

| Original response on identical H7×periods | Add static modeled WWTP after upstream landuse, **Level6 holdout MSE gain** | Then add source-reported annual upstream N | After **local** landuse, add static modeled WWTP |
|---|---:|---:|---:|
| **Measured EEA benthic EQR** | **+11.129%** | **+0.687%** | +2.726% |
| GBIF botanical speciesKey-rarefied S50 | −0.196% | −0.338% | −1.204% |
| GBIF Lycaenidae five-family indexed relative composition | −0.762% | +1.441% | −0.848% |

H6 parent 2000 resample descriptive95 for additional static W EQR gain **[+3.944%,+18.267%]**; EQR reported N [−0.123%,+1.637%] crosses zero; local crop/urban-source EQR static W [−3.419%,+8.051%] crosses zero. The Lycaenidae source N 1.441% gain H6 bootstrap95 [−0.070%,+3.313%] crosses zero, and with the local landuse +1.403% [0.030%,+2.951%], **but country-heldout local-landuse N worsens predictive performance (−1.271%)**. Thus it is NOT a stable universal N-butterfly benefit and certainly not pollination.

These are **observed out-of-parent predictive differences** in a source-dependent model, NOT "effects of sewage" or "measured ecosystem services" and the 90 H6 parent groups come from an after-failed original geometric repair sensitivity.

[Original 90-parent real H6 watershed group CV full EQR/plant/butterfly 454 same-source analysis SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43338763077).

## Strongest limitation: country generalization fails for modeled W

Country heldout groups were created with one stable majority country assigned per genuine original H7 so both periods of a basin stay held out together; this yields **11 effective native H7 country groups** (15 distinct period-level original source country labels, 14 source-period labels reassigned to stable H7 group to prevent leakage). Models for country-heldout OMIT country fixed effects, since a previously unseen country's effect is non-estimable, but retain EEA BQE-method dummies (unseen methods remain unsupported) and original environment/GBIF/river sources. Country transfer is therefore a different validation/model specification from H6 blocked within-country CV.

| Original response | Modeled static WWTP after upstream landuse **row-weighted country-heldout MSE gain** | Annual reported N added afterwards, row-weighted country-heldout |
|---|---:|---:|
| EEA EQR | **+0.953%** | +3.404% |
| Plant S50 | **−6.764%** | −1.254% |
| Lycaenidae five-family log odds | **−0.345%** | +1.367% |

Across **six** source country groups with ≥10 original native periods and ≥5 H7 basins (FR, HR, IT, LT, NL, NO), an **equal total country weight** gives **EQR static W −7.989%** MSE gain (i.e. worse), and annual reported N +9.119% (the latter very heavily influenced by sparse positive source reports, not independent crosscountry mechanism).

Actual original EQR per-country static W heldout gains are highly variable: FR n299 **+10.35%**, IT n48 **−5.06%**, LT n30 **−7.00%**, NO n22 **−17.26%**, NL n12 **−4.24%**, HR n16 **−5.83%**. Thus row-weighted gain is dominated by FR. Reported N similarly varies: FR n299 and only 4 N-positive **−1.27%**, IT **+11.89%**, LT **+12.68%**, NL **+33.70%** (11 of 12 NL source-periods N-positive), NO **+1.86%**, HR **−9.15%**. Do not interpret large % from tiny country groups as strong biological effects; baseline MSE varies and original monitoring method differs across countries.

**Bottom line:** Some modeled WWTP source geography predicts European river benthic status within related hydrologic parent groups after environmental controls, but the exact spatial predictive function is **not robustly transportable to a new country**. Plant indexed richness and butterfly reporting show no comparable stable cross-layer pattern. A regional source-reporting/data-method interaction is plausible; a pan-European biologically causal pollutant sensitivity is **NOT established**.

## Independent numerical audit using different implementation

A separate original-data reader with scikit-learn `Ridge(alpha=10,fit_intercept=True,solver='svd')` instead of the first study's manual NumPy Ridge normal equations exactly reproduced **all 3 biological responses ×6 shared predictor sets ×2 heldout strategies** source MSE values. The same 90 original Level6 parents and 11 stable source-country groups were used, no independent ecological field sample.

[Independent sklearn numerical and source identity QA SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43339114428).

## Publication-safe interpretation and next identification step

**Working question**: *When does reported and modeled wastewater-source geography contain transferable information about benthic biological status rather than simply regional differences in river hydrology, human land use and observation/reporting methods?*

- Observed H7-level repeated temporal EQR×plant×Lycaenidae **co-change** was near zero in the much smaller n169 exactly same physical monitoring station / assessment method panel. Do not conceal this negative result.
- Mapped HydroWASTE static flow to the **H7 outlet** does not identify a real upstream delivered exposure at the original interior EEA field station; source plant operating status mostly missing.
- Real report-source annual 2018/2022 EEA P/N supports only **56/454 H7 periods, 38 basin IDs**, with France registering 4 positive source-periods; 398 source-index zeros are not untreated controls.
- HydroATLAS crops and artificial built cover are **static historical covariates**, not 2017–22 measured land use transitions; unmeasured agriculture run-off or human observation effort can still confound.
- National EEA EQR methods differ, river networks are not causally independent and GBIF plant rarefied S50 and Lycaenidae indexed family mix are not ecosystem services, and no adult emergence, nectar, seed set or pollinator observation exists on the same 454 catchment-period rows.
- A genuinely independent future macro evaluation requires (i) monitored receptor-level hydrologic effluent routing and measured delivered nutrient concentrations, (ii) national method-linked consistent EQR/BQE standardization, (iii) another set of countries or decades held out before fitting, and (iv) original river-watershed-aligned field plant/pollinator consumer observation data rather than opportunistic GBIF alone.

No eBird, no merge to gesui/main or chocho, old positive EQR/N and old null plant/butterfly results preserved. All analyses explicitly observational exploratory/posthoc after previous results were observed.
