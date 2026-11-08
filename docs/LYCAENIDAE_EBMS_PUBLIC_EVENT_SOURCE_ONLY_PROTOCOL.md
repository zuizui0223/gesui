# Source-only eBMS alternative to photographic sampling: GBIF national event archives (2026-11-08)

**Objective:** determine whether accessible **standardized butterfly transect VISIT locations and dates** can support a source-blind, 2022 EEA wastewater comparison before any Lycaenidae species/count outcomes are accessed. This is an independent source feasibility branch, NOT a reinterpretation of failed 1km or 5km photo gates.

## Why necessary
Previously frozen 2022 EEA BQE × iNaturalist 1km photo comparisons had just 1/24 adequately observed pairs; independent 5km 2021 botanical-effort screen retained 17 pairs, but only 4/17 had ≥5 focal-family 2022 photos on both sides (threshold 12). Those are source observation-coverage failures, not evidence for/against a wastewater biological effect.

## Open, attributable candidate archives

- Swedish Butterfly Monitoring Scheme (SeBMS), [GBIF dataset be77e203-486c-4651-91b9-8347968b728c](https://www.gbif.org/dataset/be77e203-486c-4651-91b9-8347968b728c). Publisher Department of Biology, Lund University, **CC0**. Source IPT [publisher](https://www.gbif.se/ipt/resource?r=lu_sebms); original DwC-A candidate `https://www.gbif.se/ipt/archive.do?r=lu_sebms`.
- EBMS Denmark Transects, [GBIF dataset caa38f94-ddeb-43e7-9a27-1ac46d7ddc5c](https://www.gbif.org/dataset/caa38f94-ddeb-43e7-9a27-1ac46d7ddc5c). Publisher Biological Records Centre, **CC BY-NC 4.0** (attribution and noncommercial limits); source IPT [publisher](https://ipt.brc.ac.uk/resource?r=ebms_denmark_transects) and DwC-A `https://ipt.brc.ac.uk/archive.do?r=ebms_denmark_transects`.
- Full shared central eBMS transect data are **not** openly downloadable by default; [eBMS scientific data access requires a coordinator-approved request](https://butterfly-monitoring.net/ebms-data%20access). Do not scrape restricted central records or assume event-level point coordinates exist in all national archives.

## Prospectively required stage A; strictly NO species-level outcomes

1. Download only official public IPT published DwC-A archives, hash outer ZIP and `meta.xml`, inspect `rowType`, column names and source license. Do not open any `occurrence.txt` species/count data.
2. If `Event` core exists, read only original **eventID, parentEventID, locationID, country, eventDate, decimalLatitude/Longitude, samplingProtocol/Effort, survey date/time and any coordinate uncertainty**. Reject blurred/invalid locations. If core is Occurrence-only, report source infeasibility rather than reading species.
3. Count distinct year-2022 visits and independent site IDs per country and public coordinate resolution, preventing visit duplication. Do not infer absent butterflies from missing species records.
4. Before inspecting Lycaenidae outcomes, combine only event-location metadata with 2022 EEA *public source-only* outfall locations and hydrologic network, ensuring upstream/downstream/outside-reach logic makes biological sense for butterflies (their adult movements are not river-water transport).
5. Require a source-only minimum of **≥20 spatially independent suitable transects exposed and ≥20 suitable control transects, at least two countries or independent regions, with repeated visit denominators**. If not met, STOP before butterfly abundance/outcome inference. This threshold is exploratory and should not overwrite any earlier failed photo gate.
6. If met, freeze geographical matching (country/biogeographic zone, urban/land cover, month/year, sampling effort), source SHA and site IDs before opening species-level **Lycaenidae count** outcomes. For true tests use per-visit counts and offset for effort, not GBIF occurrence records alone.

## Scope limitation
Even an event-based Lycaenidae difference near WWTPs would be observational and cannot prove puddling, sodium/protein uptake, host plant quality, or nitrogen causal effects. WWTP reported annual N/P mass is not environmental exposure concentration. Preserve legacy chocho 22/25 Lycaenidae introduced-host resource expansion separately.


## Fixed prereceptor source-only proximity bands before SeBMS species/count values

- Use original SeBMS `event.txt` only, accepting 2022 events of `Fixed route Pollard walk transect` as the **primary standardized protocol**, with `Point site count` as a separately reported survey-type secondary group; group by original `locationID`, not by photo upload.
- Require unique original open coordinates and at least **2 distinct 2022 survey eventDates** per site; reject missing or uncertainty exceeding **500m**. Sites in 44–64°N and −12–26°E only.
- Source is the EEA official **2022** both-positive *reported measured* N/P annual load source (5,492 qualifying WWTP reports before outfall linkage), joined on the full `repCode|uwwCode` to exactly one direct freshwater `FW` discharge point with original coordinates. Report missing, duplicate and ambiguous source points, and exact underlying API query metadata.
- The fixed primary **source proximity** is ≤**2 km** nearest qualified FW outfall, and registered-source-remote comparison ≥**10 km**. A separately **prespecified descriptive sensitivity** is 2–5 km; no bands are retuned after species counts are examined. These circles are not measured downstream chemistry or verified butterfly exposure; no "effluent-feeding" claim is permitted.
- Require at least **20 independent coordinate site IDs** in near ≤2km and 20 in far ≥10km with repeat 2022 survey dates for primary support; alternative ≤5km may be described as a distinct *secondary* feasibility class if at least 20 near ≤5km and 20 far ≥10km. Require same protocol and at least two counties in each group for any comparative model.
- Matching by original source-only Swedish county and locale, survey protocol, source-site latitude/longitude and 2022 visit frequency must be frozen before opening species/count outcomes. Report if county overlap is impossible.
- Distances must be calculated from original coordinates in a metric Swedish CRS or geographic Haversine and checked against coordinate uncertainty. Do not assume that a nearby outfall irrigates a nearby butterfly transect, or that offshore/estuary outlets are freshwater sources. Do not infer "absence of sewage" from a missing EEA registered WWTP.

**Physical-source gate ≠ butterfly ecology result.** If the support fails, keep the published SeBMS source-only dataset inventory as a useful metadata result and do not open Lycaenidae species count outcomes.
