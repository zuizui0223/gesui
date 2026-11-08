# 2022 Lycaenidae × registered WWTP: source-matched public photo opportunity audit

**2026-11-08 — Source-backed data-quality / feasibility outcome, NOT a sewage butterfly effect and NOT a confirmed puddling analysis.** This is a new independent `gesui` observational feasibility pilot, separate from the successful historical `chocho` resource geography and the 2018/2022 EQR models.

## 1. Source-blind pairing completed before photo queries

Original 2022 EEA river biology site metadata was taken from [2022 source-only registry](https://github.com/zuizui0223/gesui/blob/625af95fc12e01e9c796230ca4709be72225ab50/results/eea2022_biology_bqe_20261108/SOURCE_ONLY_BQE_SITE_REGISTRY.csv), with Git blob `16f6b531e29ed09f070f9fb86f7d8feb899108ed`. **Numeric EQR scores were never loaded.** Among 3,282 source-eligible benthic river station candidates, 394 had a qualifying upstream official 2022 report of both positive measured N/P mass from a single freshwater outlet within the directed 50-km river network and 2,888 had **no registered qualifying** upstream source.

Stations were matched *without any butterfly observations* on country and 3° geographic block, log river discharge, drainage area and distance, with 3km minimum separation between any selected source/control station. A fixed hash ordering selected **24 pairs (48 locations)**: Poland **13 pairs**, Luxembourg **4**, Austria **3**, Italy **2**, Netherlands **2**. This is a *source-only* sampling frame; no-registered-source is not no-sewage. Site location derives from nonrandom river monitoring locations.

[Source-only 24-pair definition and SHA verification](https://github.com/zuizui0223/gesui/actions/runs/42447192871).

## 2. Same-year public iNaturalist photos — complete 24/24 pair queries

For each frozen location, query the official iNaturalist v1 API with taxon descendants `Lycaenidae=47923`, `Pieridae=48508`, `Papilionidae=47223`, year **2022**, radius **1.0 km**, photo present, georeferenced, public/open geoprivacy, research-grade, uncaptive, location accuracy at most 500m. These filters are photo-*opportunity* proxies, NOT a sample of butterflies, confirmed drinking, or puddling.

All **144** API family × site queries returned source totals, covering **24/24 pairs**.

| Source group | River stations | Stations with at least one photo in three taxa | Lycaenidae photo observations | Pieridae | Papilionidae | All three |
|---|---:|---:|---:|---:|---:|---:|
| Upstream official registered WWTP source | 24 | **4** | **10** | 9 | 1 | 20 |
| No registered qualifying upstream source | 24 | **4** | **2** | 13 | 1 | 16 |

The 10 versus 2 Lycaenidae photo counts **must not** be described as a favorable sewage association. The 36 total butterfly photo observations are clustered in **only eight station buffers**; one country, Luxembourg, accounts for **32/36** three-family records across the entire 48-site panel. Because no observer effort, vegetation/host resources or urban conditions were equivalently measured, the source-present group is not causally identified.

**Prospectively fixed data-support gate:** at least **12 pairs** with ≥5 photo records among all three families at **both** matched sites, before any quantitative share comparison. Exactly **1/24 pairs** met this criterion, all in Luxembourg. Therefore **do not fit a pooled WWTP–Lycaenidae proportion model, effect estimate, significance test or puddling outcome**. The four 2022 European Lycaenidae observations with `q=puddling` remain a separate, not photo-verified candidate list and were not included as a puddling outcome here.

[Completed 144-query public-photo count and matching audit](https://github.com/zuizui0223/gesui/actions/runs/42447303358) and [full per-station results](https://github.com/zuizui0223/gesui/blob/analysis/lycaenidae-wastewater-stage-contrast-20261108/results/lycaenidae_2022_source_matched_river_pairs/2022_BUTTERFLY_PUBLIC_PHOTO_OPPORTUNITY_PAIR_QC.json).

## 3. Existing positive conclusions preserved

These source-data support failures do **not** refute original `chocho` **22/25 (88%)** Lycaenidae gaining potential introduced-host resource geographic coverage; their scale is separate. Nor do they alter the `gesui` independently positive 2018/2022 benthic EQR predictive gains from upstream reported N/P mass. Neither established wastewater–butterfly biology.

**Scientific conclusion actually supported:** Source- and behavior-aligned public butterfly photo evidence in a close (1km), historical (2022), hydrologically defined sample is **too sparse** for this particular test, despite the much larger full-European 2022 Lycaenidae photographic opportunity of **40,726** observations. Scale and detection process are limiting steps.

## 4. Next reproducible approach, WITHOUT posthoc rescuing results

A separate, explicitly **new exploratory protocol** may select photo-rich regions *based on prespecified target-group observation effort and land cover* (blind to Lycaenidae share/puddling labels), with 2022 EEA source exposure independently restored, and require matched land-use/urbanization, photo timing, open GPS accuracy and actual wetted substrate observations. This changes the estimand to effort-rich accessible river sites and must be presented as a **new cohort**, not a sensitivity that rescues the 1km 24-pair failed gate.

For adult behavior, actual photographic image-based annotation with blind independent reviewers (wet surface contact/proboscis vs flower feeding vs resting/aggregation) is mandatory; no positive wastewater effect is inferred from the iNaturalist mixed-purpose observation field 17428 or free-text `puddling` hits.

**Frozen data contracts:** `2022_PREPHOTO_SOURCE_MATCH_QC.json` and `2022_PREPHOTO_EEA_BQE_SOURCE_MATCHED_RIVER_SITE_PAIRS.csv`, then separately `2022_BUTTERFLY_PUBLIC_PHOTO_OPPORTUNITY_PAIR_QC.json`. Keep all failed and successful stages with exact original source hashes, and PR Draft unmerged to `main` and without altering `chocho`.
