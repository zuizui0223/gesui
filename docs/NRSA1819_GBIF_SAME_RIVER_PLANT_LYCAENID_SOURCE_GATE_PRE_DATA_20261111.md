# GBIF terrestrial bridge: original NRSA river cohort, pre-GBIF-recording SOURCE gate (2026-11-11)

**Date**: 2026-11-11 JST. **Strictly source-first** posthoc research question following the *known* EPA NRSA river TN/TP→EPT result, not preregistered independently of it. Keep open Draft PR #58. No `main` changes.

Aim: establish whether authentic public **Plantae**, **Lycaenidae** and **Animalia observation-opportunity** records have adequate same-place / same-time coverage at **fixed original US NRSA 2018–19 river sites** for a later aquatic–terrestrial macro-analysis. The word ecosystem services requires direct processes and will NOT be inferred from GBIF count associations.

Input: frozen original **1,221 EPA UID** SOURCE-only site manifest, SHA256 `fbf53b8379d4e6564dbb51b3705b02fdd03e4a2caacad06a2314d643c8412e91` from SUCCESS Actions #43034107486. Study *site* is distinct physical `SITE_ID`, not 1,221 independent locations. The later nutrient/EPT numerical outputs are not used to select or filter sites.

**Freeze a deterministic, source-only 360 unique NRSA physical SITE_IDs** prior to fetching GBIF counts, ordered using original native `AG_ECO9` and SHA256 salted site hash. Specifically, sort within eco by digest of `NRSA2018_19_GBIF_v1|SITE_ID`, then round-robin across nine ecoregions, take first 360. Preserve chosen original site IDs, HUC8, eco9, state and precise 2018–19 point coordinates only in short-lived Actions artifact. No resampling after counts to achieve desired positive effects.

Source: genuine official `api.gbif.org/v1/occurrence/search?limit=0&geometry=...` aggregate count only. Each radius5km WKT geodesic approximation polygon (24-sided) around original river point; year **2017–2021 inclusive**, `hasCoordinate=true`, `occurrenceStatus=PRESENT`. Three separate taxa with the SAME filter:
1. `taxonKey=6` Plantae;
2. `taxonKey=47923` Lycaenidae;
3. `taxonKey=1` Animalia record counts, only a broad recording-process comparator, **NOT a negative-control unaffected taxon**.

Exact request URL and full HTTP JSON SHA256 captured; HTTP error = **missing**, never count=0. GBIF species occurrence != flower visits, host use, insect flux, or standardized butterfly counts. The butterfly and plant ecology *must* be controlled for record effort in any later model. The study is exploratory and cannot attribute a nutrient response to actual wastewater without measured WWTP source and directed exposure.

**Locked source-only support gate BEFORE looking at these new GBIF counts:** 250/360 sites all three queries succeeded; >=100 unique sites with >=50 original Plantae records AND >=20 Lycaenidae records; >=15 US states and >=5 AG_ECO9 categories *among those well-supported sites*. If any fail -> stop, no aquatic→terrestrial regressions or mediation. Do **not** relax distance, dates, taxa, source filters, or response denominators to recover sign. Include nested sensitivity at 10km only after 5km gate pass and with a new explicit protocol; never claim expansion changes service levels.

The initial outcome to report is **data coverage only**. A later joint model (if pass) needs effort-controlled within-HUC8 and heldout-region tests comparing `EPR larval fractions ~ water TN/TP` and `GBIF Lycaenidae photo abundance/record share ~ EPT/potential plant hosts`, with a proper background observation model; no mechanistic mediation unless emergence/diet/pollination directly measured.
