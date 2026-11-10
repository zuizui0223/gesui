# Denmark eBMS 2025: registered wastewater outfall proximity rarely verifies downstream river-path exposure

**2026-11-10. Original source/geometry-only Stage 0 COMPLETED. STOP: directional biological comparison unsupported.** This is retrospective source feasibility after the original 2025 NEW399 Lycaenidae, 3-versus-5 family inversion and 2021 WorldCover adjustment were already known. **NO butterfly family/taxon counts were accessed in this new physical source gate**. Absence of verified river connection is not biological evidence of wastewater absence, and near-source site proximity is not a WWTP nutrient effect.

## Original fixed sources and trust boundaries

- Denmark **2025** eBMS licensed CC-BY-NC 4.0 official v1.4 Event-only QA, original published public sites; same exact original 399 new-site frame (previous analyzed 107 site coordinates excluded), **286 source-event-supported site coordinates with >=2 original any-published-Occurrence dates** (70 originally <=5km from a registered qualifying WWTP, 92 >=10km, 124 intermediate). The 2025 taxonomy/individualCount file was **not opened** in this source stage.
- Official EEA 2022 registered single-freshwater-outfall DK2022 both measured-N/P-positive source records, **203** distinct original outfall SHA IDs, from the originally archived pre-2025 butterfly Event-only inventory. Actual 2025 WWTP release, N/P/organic matter, stream concentration and seasonal flow are NOT measured.
- Reused original **third-party** converted European HydroRIVERS v1 GeoParquet from gesui PR #9, exact original SHA256 `5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f`, downloaded and hashed in a separate GitHub Action. **This conversion has not been independently verified against authoritative HydroRIVERS official source bytes**; all connectivity estimates are *provisional model-based directed river graph results*, not official measured exposure.
- Parsed the original 15,045 third-party HydroRIVERS directed reaches within expanded Denmark/border bounding box (53–59 N, 6–17 E). Snapped site point centroids and registered outfall point geometries independently within **500m of stream segment** in EPSG:3035; tied reaches sorted deterministically by distance and reach ID. `NEXT_DOWN` directed reach IDs determine approximate river-path direction, with source/site on the SAME river reach **explicitly excluded due to unknown order within reach**. A 500m site-centroid snap does not guarantee monitoring transects follow a river.

## New actual original-source findings

| Original source eligibility measure | Numerically verified result |
|---|---:|
| All selected 2025 eBMS site centroids (including middle/remote) | **286** |
| eBMS sites snapping to ANY HydroRIVERS centerline <=500m | **60/286 (21.0%)** |
| EEA FW registered outfalls snapping <=500m | **125/203 (61.6%)** |
| Original <=5km outfall-near eBMS site centroids | **70** |
| Of original 70 near eBMS sites snapping <=500m | **20/70 (28.6%)** |
| All source/site straight-distance <=5km candidate pairs (multiple outfalls per site possible) | **92** |
| Source/site pairs without a site-to-river snap <=500m | 69 |
| Pairs with source outfall not snapped <=500m | 2 |
| Both snapped but no verified directed route <=50 network-km | 3 |
| Same original reach with within-reach ordering UNKNOWN | 4 |
| Strictly **DOWNSTREAM connected within 10 network-km** | **10 source×site pairs** |
| Strictly **UPSTREAM connected within 10 network-km** | **4 source×site pairs** |

Crucial pseudoreplication: the **10 downstream pairs are exactly 10 unique eBMS survey-point centroids but only TWO registered WWTP outfall IDs AND TWO downstream HydroRIVERS reach IDs**. One outfall contributes **7** sites to one destination reach, the other **3** sites to a second. They cover just **two original 1° geographic cells**. The four upstream placebo sites likewise arise from **two different registered outfalls and two river reach IDs** across two 1° cells. So ten butterfly sites are **not ten independent wastewater exposure systems**.

Only **7/70** original near sites are <=100m from any HydroRIVERS line, and **16/70** <=250m. Modelled <=500m snapping is a permissive physical proximity gate; the point does not prove actual butterfly access to downstream streambank vegetation or nutrient-rich water. The dataset contains no independent riparian butterfly transect segment geometry, larval host plants, ants or adult mineral/protein puddling observation.

## Original fixed causal-inference stop conditions (established BEFORE this graph had been evaluated)

| Geometry/source support threshold (precommitted) | Required | Observed |
|---|---:|---:|
| Strict downstream site coordinates within 10 river-km | >=12 | **10** |
| Distinct registered downstream source outfalls | >=8 | **2** |
| Independent 5km spatial outfall components | >=5 | **2** |
| Original 1° regions with downstream support | >=3 | **2** |
| Strict upstream site coordinates | >=8 | **4** |
| Distinct original upstream outfalls | >=4 | **2** |
| Original 1° regions with upstream support | >=3 | **2** |

**Outcome Stage-1 is NOT PERMITTED.** No registered-source/downstream-vs-upstream effect coefficient, p-value, butterfly family trend, effluent chemistry or adult water/substrate benefit was computed. Changing 500m snap, 5km source proximity, 10km directed reach-path limit, or including same-reach as downstream would be a **new study** and cannot retroactively rescue this source-frozen gate.

## Implications for the previous Lycaenidae/WWTP project

The prior broad 2025 NEW399 study showed a **pooled three-family composition OR 2.270** near registered outfalls, but **absolute Lycaenidae documented count rate ratio 0.833**, and a five-family relative composition OR **0.772**. Six-region/actual measured 2021 landcover-adjusted site-equal 3-family OR was **0.511**, while 5-family OR was **0.377**. Those results remain **faithful observational reporting relationships**, not proven sewage effects. The new river graph shows why the point-source-to-butterfly exposure mechanism is particularly weak: only two potentially directed source–monitor reach linkages can be distinguished with existing Event locations, even after expanding to 399 NEW sites. That does NOT mean actual freshwater effluent is unavailable to all other terrestrial butterfly habitats; smaller streams are omitted in HydroRIVERS and point centroids may not represent survey paths.

**Scientific decision:** Stop reusing this dataset to fit new causal wastewater-butterfly models. Preserve the cross-denominator/geographic/vegetation Simpson-like observation sensitivity as a legitimate bounded result. A source-to-receptor ecological mechanism needs a newly designed or independently obtained riparian transect/source pair population with exact streambank survey geometry, >=10 genuine independent receivers and control rivers, repeated season-specific field water/substrate chemistry, larval-host vegetation and observation efforts. Treat an outfall/river as a physical water transport path only, NOT a butterfly dispersal route.

## Reproducible original source receipts (all Draft, no main/chocho merge)

- [Original source-only 286 eBMS × 203 EEA × 15,045 river reach directed classifications](https://github.com/zuizui0223/gesui/blob/analysis/dk-ebms2025-hydrorivers-outfall-exposure-gate-20261110/results/ebms_denmark_2023_2024_WWTP_source_only/DK2025_NEW399_286SITE_203WWTP_SHA_LOCKED_HYDRORIVERS_DIRECTED_SOURCE_FEASIBILITY.json), [original verifiable river source run SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42779582374).
- [Independent original site/outfall/reach SHA integrity, mutually exclusive direction classes, 10 sites collapsing into 2 reach IDs](https://github.com/zuizui0223/gesui/blob/analysis/dk-ebms2025-hydrorivers-outfall-exposure-gate-20261110/results/ebms_denmark_2023_2024_WWTP_source_only/DK2025_NEW399_286SITE_203WWTP_HYDRORIVERS_SOURCE_DIRECTION_OFFLINE_AUDIT.json), [offline source QA run](https://github.com/zuizui0223/gesui/actions/runs/42780093058).
- [Fixed design, thresholds, interpretation restrictions](https://github.com/zuizui0223/gesui/blob/analysis/dk-ebms2025-hydrorivers-outfall-exposure-gate-20261110/docs/DK2025_ORIGINAL286_EEA2022_HYDRORIVERS_DIRECTED_SOURCE_ONLY_GATE_20261110.md).

**No old source/outcome file overwritten, no original licensed eBMS individual source data, no original full GeoParquet committed, no biological outcome accessed in this new source-only procedure.**