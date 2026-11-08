# Wastewater × Lycaenidae 2022: new 5km independent-effort cohort result (2026-11-08)

**A new source-backed observational FEASIBILITY result, not a biological wastewater effect.** Original prior 24-pair 1km study remains unchanged: one photo-sufficient pair. This new 5km cohort was independently defined from source data **before** 2021 effort and 2022 butterfly photo outcomes were queried.

## Reproducible gates and data

### A. Freeze freshwater-registered source vs no qualifying registered source

- Original EEA 2022 BQE river source-only 3,282 benthic site records; **394** source-positive under the qualifying 50-km upstream reported-measured N/P network, **2,888** no qualifying registered source.
- Source-only matched country and 3° region with log river discharge/upstream basin area and site separation ≥10 km, without using any butterfly or botanical observation outcome. Pre-locked **78 pairs / 156 sites / 10 countries**, capped at 16 pairs per country. Original source-only SHA and CSV committed before observational effort.
- The last group is NOT a verified sewage-free control and the 5km buffer is not a direct wet-substrate exposure measure.

### B. 2021 non-focal public Plantae photo effort gate

- Official iNaturalist public observations with taxon Plantae, georeferenced photo, research grade, uncaptive, ≤500 m location error, 2021, 5km radius for **all 156** previously chosen sites.
- Prior locked criterion: at least **8** plant photos on each side and smaller/larger photo count ratio ≥**0.20**, then rank by smaller 2021 botanical count, cap 8 pairs/country and 32 total; ≥12 pairs / ≥3 countries required before any 2022 butterfly query.
- Source-backed result: **17/78** effort-eligible, 17 selected, **6 countries** (AT3, BE1, IT8, LU1, NL3, PL1). All 156 queries succeeded. Median lower-side 2021 plant count **31**, median paired ratio **0.60**. **Predeclared effort gate PASSED**. Source records committed before querying butterflies.
- Plant reporting effort is not butterfly reporting effort; conditioning on prior-year plants can itself select access/land-use conditions and change target population.

### C. Separate 2022 butterfly photo opportunity audit

- Independent official descendant taxa: Lycaenidae `47923`, Pieridae `48508`, Papilionidae `47223`. Same 5km radius, year 2022, research-grade georeferenced photo, uncaptive, public GPS accuracy ≤500m.
- All **34 site × 3 family = 102** API queries succeeded.
- Raw pooled source-associated photographic recording counts (NOT individuals or puddling):

| 2022 butterfly photo records | Registered upstream qualifying source | No registered qualifying upstream source |
|---|---:|---:|
| Lycaenidae | 50 | 73 |
| Pieridae | 52 | 40 |
| Papilionidae | 20 | 11 |
| Total | **122** | **124** |

- The previous 1km cohort had 36 three-family photo observations across 48 sites; the new 5km effort-selected cohort has **246** across 34 sites. These are very different target areas and source-effort designs and **cannot** be called a tenfold biological increase.
- **PREDECLARED within-pair support gate**: ≥5 photos from any of the three families on *both* stations in ≥12 pairs and ≥3 countries. **Only 4/17 pairs** met photo counts, distributed IT2 / LU1 / NL1. **GATE FAILED.**
- The four pairs' count ratios and photo shares exist as source diagnostics, but the gate disallows pooled inferential share effects or significance. No wastewater-caused Lycaenidae increase/decrease is inferred; the pooled 50 vs 73 counts are **not** an effect estimate.

**Important distinction:** The much larger pan-European 2022 Lycaenidae photo pool (40,726 in 44–64° N, -12–26° E) is NOT a 40,726-individual census and is spatially disconnected from many eligible EEA reported-outfall receptor stations. Spatial and behavioral matching remains the bottleneck, not global photo availability.

## Scientific interpretation and next source that can break the observational bottleneck

**Switch validation source type rather than posthoc moving the radius or significance threshold on the existing photos.** European Butterfly Monitoring Scheme (eBMS) uses fixed-visit transect counts; open national sampling-event datasets are registered in GBIF (e.g. Sweden public CC0 and Denmark public CC-BY-NC). Raw eBMS centrally merged count/transect data are governed by coordinator approval; do not scrape private data or assume full 5km coordinates can be downloaded. A new independent official-source eligibility audit can read only freely licensed *event/visit site metadata*, freeze overlap with 2022 EEA WWTP outfalls and require actual visit denominators before opening species-specific Lycaenidae abundance counts.

- [eBMS access policy](https://butterfly-monitoring.net/ebms-data%20access)
- [GBIF SeBMS sampling-event dataset](https://www.gbif.org/dataset/be77e203-486c-4651-91b9-8347968b728c)
- [GBIF eBMS Denmark sampling-event dataset](https://www.gbif.org/dataset/caa38f94-ddeb-43e7-9a27-1ac46d7ddc5c)

## Artifacts, immutable source provenance

- [New 5km river source-only 78-pair freeze](https://github.com/zuizui0223/gesui/actions/runs/42453851208)
- [Pre-2022 2021-Plantae effort 17-pair freeze](https://github.com/zuizui0223/gesui/actions/runs/42453941319)
- [Same-year 2022 butterfly family photo counts and support-gated analysis SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42454226280)
- `2022_NEW_5KM_PREPHOTO_SOURCE_ONLY_PAIR_QC.json`
- `2021_PLANT_EFFORT_SELECTED_PAIRS_FROZEN.json` and same-stem CSV
- `2022_BUTTERFLY_PHOTO_SHARE_AFTER_PREYEAR_PLANT_EFFORT.json`

Previously verified `chocho` potential larval-host geography **22/25 Lycaenidae expanding under introduced host resources** is retained but is **not** a sewage or adult behavior effect. Also retain 2018/2022 `gesui` EQR held-out predictive gains. All new source work remains Draft and unmerged; do not imply an ecological effect from source-support failure.
