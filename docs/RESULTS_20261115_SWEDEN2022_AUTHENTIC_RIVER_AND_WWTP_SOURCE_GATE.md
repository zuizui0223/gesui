# 2022 Sweden SeBMS point sites × genuine EEA wastewater outfalls × authentic HydroRIVERS — source-only audit

**2026-11-15 JST / source-only / noncausal / frozen primary HNS NOT RUN.**  
This is continuation of existing Draft PR #95 and its exact frozen v1.0 riparian source/reference protocol. The biological study was **not** resumed; no HNS species counts, butterfly traits, leaf chemistry, adult puddling, or 8-arm effects were read or fitted. All 264 Swedish point-site results derive exclusively from 2022 public Event-core locations and 2022 EEA registered source geography, with original HydroRIVERS v1 river centerlines.

## New first-party provenance and independent recheck

- Sweden source points: exact [SeBMS 2022 point-only 264-source site registry](../results/20261115_riparian_source_only/SEBMS_2022_POINT_SITE_LOCATIONS_SOURCE_ONLY.csv), original Git blob **a56d143a3474027092a438399e6c2e59066d911f**. Original SeBMS publisher DWCA v1.13 full ZIP SHA256 **384a2e5145df29268fe2de79fde9e196fc113e68adfc922a81acd58069437348**; 2022 May–September **1,635 point visits /264 distinct sites**, **127 with >=5 distinct visit dates**, all reported geolocation uncertainty 25 m. No occurrence extension opened.
- EEA original report: official 2022 [UWWTD public REST TreatmentPlants_2022](https://water.discomap.eea.europa.eu/arcgis/rest/services/WISE_UWWTD/TreatmentPlants_2022/MapServer/0) and paired DischargePoints_2022 layer, **exact Swedish reported code SE369**. Server-side filter returned **25 reporter plants with both positive reported N and P** and **430 national discharge records** before one-to-one join. Only **14 matched single-FW source outlets** meet the *previously fixed source definition*. **Do not label the excluded 416 points as zero sewage.** Neither annual reported N/P mass nor source coordinates demonstrate measured 2022 nutrient concentration at the butterfly points.
- **Official original** HydroRIVERS v1 EU SHP ZIP from https://data.hydrosheds.org/file/HydroRIVERS/HydroRIVERS_v10_eu_shp.zip, publisher ZIP SHA256 **500da7d36ceee0aa4c82dd625f5c53d80ba2f7ac14ad022cd3a6601bf55068e5**; **79,439** native reaches in extended Sweden/border spatial extent, retaining native HYRIV_ID/NEXT_DOWN/LENGTH_KM and original river-line geometry. Official source validated; the original copyrighted shape network is **not redistributed** here.
- Prior independently source-audited [PR #21](https://github.com/zuizui0223/gesui/pull/21) 2022 site-to-outfall proximity CSV original Git blob **bbc1c423dd95eeddda224d8326caba13e61af06f**. Among the new 264 centers, **225** matched earlier 2-visit source geographies; all matched nearest-distance values to the same 14 official outfalls re-calculated **exactly** (maximum difference 0 km to saved original precision). **39** newer centers are not covered by the old 2-visit/geographic selection; never infer zero distance from those missing earlier rows.

### Source-only audit GitHub Actions SUCCESS

[Verified original source-only rerun #43694358102](https://github.com/zuizui0223/gesui/actions/runs/43694358102). Original numerical JSON and detailed geometry output **signed** workflow artifact ID **14153822488** was independently SHA-checked before replaying into Draft branch: [archive replay SUCCESS #43694715664](https://github.com/zuizui0223/gesui/actions/runs/43694715664).

- [Exact original-run 264-site river distance and source report table](../results/20261115_riparian_source_only/SW2022_ORIGINAL_RIVER_SITE_DISTANCE_SOURCE_ONLY.csv), original artifact file SHA256 **cbf62975ccdb4d8f52ba5056799c51a1b54fb6da1cb821310f876a82ddf8eff7**.
- [Six original-site/outfall candidate relations, no species](../results/20261115_riparian_source_only/SW2022_ORIGINAL_RIVER_DIRECTION_SOURCE_RELATIONS.csv), original artifact SHA256 **bbf608164bd93a7b2d3efe6a5ba0782f2af3a551cc6a27e2529ca6243853054f**.
- [Source-only machine-readable frozen STOP result](../results/20261115_riparian_source_only/SW2022_ORIGINAL_RIVER_SOURCE_FEASIBILITY.json), original artifact SHA256 **0cf39b1a07294d6a0498a2df5376a881be65a29ba8d5061c86d93014a441ffb8**.
- [Executable original river and source auditor](../analysis/20261115_riparian_source_only/audit_sw2022_original_river_and_outfalls.py). Frozen protocol Git blob **b9ce415c36d398eee77d8318a20aef4873e0bf44** unchanged.

## Newly measured SOURCE geography only

Among the **127 source point centers with >=5 May–September 2022 visit dates**, and requiring the *entire original 25 m public positional uncertainty disk* to fit the frozen distance bands for the **nearest any mapped HydroRIVERS river**, source-only counts are:

| Geographic source-only class | 2022 repeat-5 centers |
|---|---:|
| **Entire 25m circle** within 0–200m of nearest mapped river | **13** |
| Entire 25m circle within 500–1500m from nearest river | **38** |
| Boundary/transition (not contained in either locked band) | **34** |
| No original mapped river within 1700m | **42** |
| **Total** | **127** |

Distances to *nearest river anywhere* alone do not show that it is the **receiving river for an EEA source** or that the sampled host plants or individual butterflies contact that river. A survey center is not an individual butterfly's observed position. The 25m point-count area is a source-only conservative uncertainty footprint; its biological sampling footprint needs separate confirmation.

**14 original qualifying EEA freshwater reported source outfalls** were source-georeferenced, and **12** snapped to original HydroRIVERS within the previously used 500m rule. For all >=5-visit SeBMS centers, we screened all these snapped sources for necessary Euclidean proximity within **7.5 km** (broad source-only superset of a 1–5km network path plus <=1.5km terrestrial offset and <=0.5km source snap). Result:

- **6 distinct repeat-5 butterfly survey site centers** paired with **4 distinct official reported outfalls**, source-to-site Euclidean distances **2.981–6.771km**.
- Among these 6: **2 river-near**, **2 inland**, **2 boundary or intermediate band**. Both prior <=5km outfall-near sites are **inland or boundary**, not source-verified river-near.
- Native NEXT_DOWN graph **full-reach route-length** screening found **zero** source→site downstream **1–5km approximate** candidate relations, and zero strict upstream qualifying relations; the six are classified **ROUTE_OUTSIDE_1TO5_OR_UNRESOLVED**.
- **Critical caveat:** this is a *conservative reach-level approximate source-only route*, not a point-along-the-line shortest path. Whole-reach length may overstate a source-to-point path, and endpoint source/survey nearest-reach assignment, small missing tributaries and the expanded-country clipping can obscure connection. Zero qualifying approximate paths is NOT proof that no real source reaches the terrestrial survey locations. It is a hard stop for claims using the **current verified** source data.

## Pre-existing protocol decision

**Frozen biological 2022 8-arm study remains STOP_DATA_OPPORTUNITY_NOT_YET_DEMONSTRATED.**

The exact, immutably frozen primary requires **24 independent matched reported-source/reference reaches × 8 site-position cells, >=12 genuine H6 parent basins, >=3 countries**, with upstream/downstream real river routing, 0–200m near and 500–1500m inland, 2022 date-aligned >=5 observer-date opportunities, >=20 butterfly five-family records and >=2 **scored** Lycaenidae per cell and >=70% HNS record-source coverage. The Swedish 127 repeat-5 sites do not meet this independently source-verified 8-arm opportunity. Sweden is **only one country**, and no new 2022 fine-precision source/reference sample was acquired in another two countries. None of the species-based conditions has been tested; no fake zero was imputed for sites without source confirmation.

The previous [DK 2025 directional gate](https://github.com/zuizui0223/gesui/pull/42) (only two independent downstream sources) and [DE 2022 nationwide source/photo gate](https://github.com/zuizui0223/gesui/pull/36) (five old 5km four-arm quartets, different geometry) remain separate, prior *unsatisfied* studies. They must not be combined with these Swedish results as if they were 2022 frozen 8-cell replicates.

## Citable acquisition implications, no new inferential test

Independent two-country *standardized 2022 event/site* data with precise original section coordinates is now the main bottleneck, not additional H7 ridge/SEM fits. Public [Flanders older all-species BMS](https://ipt.inbo.be/resource?r=bms-occurrences) has temporal coverage ending in 2020 on GBIF and publishes generalized 5km positions; newer [Meetnetten priority butterflies](https://www.gbif.org/dataset/63938753-1fec-4c08-ae39-e9f8a6576521) generalizes sensitive sites into 1/5/10km cells. [eBMS data access](https://butterfly-monitoring.net/index.php/ebms-data-summary) offers research-access pathways; precise original site/section/event data can only be considered with applicable permissions. Neither published coarse centroid nor requested precise data is an eligible sample without the unchanged 2022 source-only gate and matching.

**No 2022 Lycaenidae numerical HNS from individual sites was read, no confidence interval/p-value for the main outcome was estimated.** Do not interpret this source gate as evidence that wastewater does not transfer to terrestrial insects. Keep PR #95 Draft; no main/chocho merge, no new branches, no eBird/TinyFish.
