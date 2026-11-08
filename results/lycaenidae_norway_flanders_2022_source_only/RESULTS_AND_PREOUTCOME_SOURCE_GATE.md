# 2022 Norway/Flanders × wastewater reported measured N/P × butterfly monitoring: independent source-only verdict

**Date:** 2026-11-08. **Scientific status:** actual public publisher sources, SHA-locked Event-only audit, then original EEA 2022 discharge-source reconciliation. **No butterfly species-level outcome, abundance, puddling, or aquatic EQR numeric values read. No wastewater effect estimated.**

## 1. Independent positive source-quality results

Publisher-original Darwin Core Archives were downloaded, verified, and their Event core was inspected without opening any `occurrence.txt`, species, or individualCount extension.

| National monitoring source | 2022 original events | Source geolocated 2022 | Observed coordinate uncertainty | Source-independent biological sample support |
| --- | ---: | ---: | ---: | --- |
| **Norway NINA** CC BY 4.0 | **4,411** | **4,199** | **5m** | **70 independent Lucas grids**, 20 source sections each |
| **Flanders INBO Meetnetten** CC0 | **3,573** | **3,573** | **7,071m** all records | 0 sites pass ≤500m source-position gate |
| **Flanders historical** CC BY-NC | **0 in 2022** | 0 | N/A | dataset stops 2020 |

Norway source correction is crucial. All 1,400 originally extracted section locations are not independent! Reconstructing the original grid identity from `stateProvince|municipality|locality` gave **exactly 70 groups × 20 distinct original section locationIDs**. All 70 groups have ≥2 source-dated 2022 repeat visits and compact ≤2.5km source spatial footprints; **68/70 have 3 distinct dates**. Their 2022 source segment positions have published accuracy **5 m**. Each grid, not section or visit, is the independent monitoring unit. Thus the study's source-side statistical unit is defensible before reading any response. [Source-only independent unit code and SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42455790667).

The withheld Flanders original accurate coordinates are never reconstructed or reverse-engineered. Source publicly reports 7,071m uncertainty, so the designed ≤500m local freshwater outfall test **cannot** use these 2022 sensitive observations, even though they are rich repeated observations. The dataset's public broad UTM grid is not exposure to a specific outfall.

Archive original outer ZIP SHA256 receipts are in `PUBLIC_NORWAY_FLANDERS_2022_EVENT_SITE_SOURCE_QC.json` and `NORWAY_LUCAS_2022_INDEPENDENT_GRID_SITE_QC.json`. Source-country licences retained individually.

## 2. 2022 EEA measured effluent N/P source data explicitly reconciled

Official EEA Waterbase UWWTD 2022 public source REST reports **5,492** WWTP entries with both positive annual `uwwNDischargeMeasured` and `uwwPDischargeMeasured` mass; original public discharge inventory **24,416** points. Join full original `repCode|uwwCode` uniquely to exactly **one** direct freshwater `dcpWaterBodyType=FW` outlet and quality-check coordinates. Reporter codes were identified from the original source-side country and geography before any butterfly outcome: **Norway `NO_UWWTD_2022` and Flanders/Belgium `BE-2022`**.

| National reporting source | Unique eligible original 2022 FW outfalls |
| --- | ---: |
| Norway (same 2022 reporting cycle) | **43** |
| Belgium (same reporting cycle) | **15** |

Report-year WWTP annual mass is not a chemical concentration/dose measured at a butterfly puddling patch, and individual laboratory measurement dates have not been externally verified. Outfalls with absent N/P values or multiple outlets were excluded; registered-source-zero is not a river free from sewage.

## 3. Actual Norway BMS Lucas grid × EEA FW outlet proximity, prespecified gate

For each frozen **70** independently located Lucas grid (all 20 source 50m sections), compute nearest source outfall distance to 43 Norwegian `NO_UWWTD_2022` qualifying FW outlet coordinates. **Conservative rule:** near if grid centroid distance + maximal original section radius + published GPS uncertainty ≤2km (secondary ≤5km); remote if centroid distance − section radius − GPS uncertainty ≥10km.

| Source-only sample | Independent 2022 Lucas grids | Counties |
| --- | ---: | ---: |
| Strict near registered both-positive N/P FW outfall ≤2km | **0** | 0 |
| Secondary near ≤5km | **2** | 1 (Viken) |
| Conservative remote ≥10km | **63** | 7 |
| Intermediate / ambiguous 5–10km | 5 | — |

**Predeclared 2022 source support gate**: ≥10 independent near grids **and** ≥10 remote grids, with near/far groups both present in ≥2 counties. **Fails** for both 2km and 5km bands (near 0 or 2, only one county at 5km). Note: proximity is a terrestrial landscape proxy, not directed water delivery or adult puddling.

**Decision:** do **not** open Norwegian original species/count occurrence data or fit a Lycaenidae wastewater association. The measured-source spatial sampling overlap is inadequate under the fixed 2022 criteria, so the outcome is unidentifiable here. Do **not** call it a negative biological association.

[2022 Norway/Flanders official source-only BMS–WWTP audit GitHub Actions SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42455925114).
Full `2022_NO_BE_BMS_GRID_VS_EEA_WWTP_SOURCE_ONLY_QC.json` and all 70 exact original-grid distance records versioned and SHA-locked.

## 4. Existing positive results preserved and next defensible source

- Original frozen `chocho`: **22/25 Lycaenidae (88%)** have increased potential larval-host resource geography from human-introduced hosts; not actual larval reproduction or sewage-caused abundance.
- Original `gesui`: independent 2018/2022 EEA benthic macroinvertebrate EQR prediction gains from upstream measured-reported N/P remain supported, but separate from butterflies; 2020 third-year source coverage insufficient.
- Previous 2022 butterfly photo 1km source-pair support **1/24**; preregistered prior-year plant-effort 5km cohort butterfly support **4/17**. Both source-support failures retained and not reoptimized.
- Original Norway 70-grid reconstruction is independently usable for **other predeclared, same-study 2022 ecological questions**, subject to species-access/effort/rights, but cannot currently establish a spatial wastewater relation.
- A defensible future wastewater × Lycaenidae **standardized population** comparison needs transect event data from countries with many 2022 registered measured effluent outlets (e.g., EEA reporters **DE, AT, HU**) and valid fine public coordinates, or a coordinator-authorized eBMS dataset. Central eBMS original individual count/coordinates have data access constraints and cannot be silently acquired. This is a source/data-access condition, not a promise of ecological result.
- A separately verified controlled photo/field study can test adult organic protein/Na substrate use directly, versus larval host quality and ant specificity. WWTP total N/P mass is not protein or Na.

**Only source eligibility was studied; all ecological effects remain untested.** Keep branch Draft and do not merge into `main` or alter `chocho` GEB manuscript.

## Provenance
- [Original Norway monitoring public source](https://www.gbif.org/dataset/aea17af8-5578-4b04-b5d3-7adf0c5a1e60)
- [Original Flanders Meetnetten public source](https://www.gbif.org/dataset/63938753-1fec-4c08-ae39-e9f8a6576521)
- [Source-only Event core QA SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42455632989)
- [Source-only independent Lucas-grid QC SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42455790667)
- [Original 2022 EEA WWTP × 2022 monitoring-source QC SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42455925114)
