# SeBMS v1.13 2022 event/source-only acquisition: point sites, date coverage and position precision

**2026-11-15 / source-only QA / biological HNS data NOT accessed / primary STOP unchanged.**
Existing frozen one-shot riparian HNS protocol SHA (Git blob): b9ce415c36d398eee77d8318a20aef4873e0bf44. No amendments to core year, 8-arm source/reference geometry, distance bands, HNS endpoint, 24-pair source gate or inferential plan.

## Actual original source obtained and pinned

- Publisher: [Swedish Butterfly Monitoring Scheme (SeBMS), Lund University](https://www.gbif.se/ipt/resource?r=lu_sebms&v=1.13), CC0, published **2026-10-16** as **v1.13**. Current publisher reports **39,466 Event-core, 201,386 Occurrence-extension, 238,737 eMoF-extension** records for the entire archive. The GBIF landing page may still cache v1.12 from 2025; **the newer verbatim publisher source is controlling**.
- Reproducible exact source download URL: https://www.gbif.se/ipt/archive.do?r=lu_sebms&v=1.13
- Full source **DWCA ZIP SHA256 384a2e5145df29268fe2de79fde9e196fc113e68adfc922a81acd58069437348**. Exact Event-core member event.txt SHA256 **422b06fdb13c77d25ec96471aa348b7c357692622d87ad4af7215a975578c6bc**; metadata parsed through publisher meta.xml.
- Only the **39,466 Event-core** rows were opened and parsed. The **201,386 Occurrence** rows and the **238,737 eMoF** rows were NOT opened. **No butterfly species/taxa, larval host HNS, individualCount, or model response were read or calculated**. Original source provenance can be independently checked by running the script with immutable protocol and archive version pin.
- Original [GitHub Actions 2022 event QA #43692977024 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43692977024) and expanded [point-site source QA #43693414744 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43693414744). Source-only output [JSON](../results/20261115_riparian_source_only/SEBMS_2026_V113_EVENT_SOURCE_ONLY.json); public coordinate-only [264 point sites CSV](../results/20261115_riparian_source_only/SEBMS_2022_POINT_SITE_LOCATIONS_SOURCE_ONLY.csv); [source-only Python code](../analysis/20261115_riparian_source_only/audit_sebms_2026_113.py).

## Source evidence: observed 2022 summer visit opportunities

All counts below are from **real publisher Event-core rows**, not model predictions, not biological effects, and not verified WWTP-linked cells.

| Original 2022 May 1–Sep 30 event source | Point counts (25m published uncertainty) | Fixed-route transects (1000m published uncertainty) | Total |
|---|---:|---:|---:|
| Event-core rows | **1,635** | 805 | 2,440 |
| Distinct 2022 site locationIDs with stable exact coordinates | **264** | 178 | 442 |
| Sites with >=5 distinct dates from May–September | **127** | 68 | 195 |
| Meet frozen original maximum 100m coordinate-uncertainty *at the event center* | **1,635 / 1,635** | 0 / 805 | 1,635 |
| Original 2022 full-year event-core visits | | | 2,677 |

The published uncertainty is **25m** for all included 2022 summer "Point site count" events, and **1000m** for all "Fixed route Pollard walk transect" events. The SeBMS publisher defines point counts as a **25m radius observed for 15 minutes**; a 25m sampling footprint around a center is **not a precise observed location for an individual butterfly**, and a 0.5–3 km transect centroid is not a 100m butterfly coordinate.

Original 2022 event ID/site-date identity and a source-only site-year CSV were regenerated from the exact Event-core using **original locationID + exact latitude + longitude** as source grouping. The 264-row CSV uses SHA256-truncated hashes instead of raw locationIDs, and contains only public point-center latitude/longitude, original coordinateUncertaintyInMeters, original samplingProtocol, May–September event and distinct visit-day counts, and status flags. No taxon identities, observer details, owner IDs or butterfly outcomes are included. **All 264 rows, 127 repeat-site counts and 1,635 event totals were independently checked against the Event-core reader's full-sample aggregate.**

## What this DOES and DOES NOT accomplish for the frozen eight-cell test

**Source-only positive finding:** an authentic 2022, 25m-radius point-count scheme with repeat visits exists. **127 publicly georeferenced point sites** have at least 5 distinct summer survey days under their original reported site positions, and no lack of precision is demonstrated at this first site-source gate.

**Crucial untested constraints:**
1. No official reported **2022 municipal WWTP** source has been hydrologically connected to these 264 point centers, and no true H7/H6 basin assignment or source-vs-reference river reach pairing has been verified. Zero verified pairs means **not yet evidenced**, not geographically absent.
2. No entire **25m footprint uncertainty disk** has been tested against the frozen 0–200m bank and 500–1500m inland bands, and 1–5 network-km downstream/upstream positioning is unknown. If the center is inside a 200m river bank band but within 25m of its boundary, the whole sampled footprint may straddle the band.
3. The original 2022 butterfly speciesKey number, >=20 five-family source records per cell, >=2 scored Lycaenidae per arm, and >=70% Lycaenidae HNS-source coverage are not evaluated. **Do not use 1,635 event visits as 1,635 butterfly records, independent treated streams, or 1,635 matched local cells.**
4. Sweden alone supplies only **one of three required countries**. Even if every one of 264 point sites were river-linked, the frozen three-country and 24 independently matched source/reference pairs criteria would still need independent Germany/Belgium/other countries and event/site/biological availability. Thus **primary HNS STOP persists**.
5. Some current publicly listed schemes are unsuitable for the frozen 2022 precise-point requirement: [Denmark EBMS 15-minute counts](https://www.gbif.org/dataset/99f9976a-dfdd-4359-bf4a-5538f29458cf) started in **2023**; [Flanders transects](https://www.gbif.org/dataset/63938753-1fec-4c08-ae39-e9f8a6576521) generalize original locations to **1/5/10km UTM grids**; the public [UKBMS occurrence dataset](https://www.gbif.org/dataset/1e266c3d-92ef-4d5a-8e4a-c04742c772c3) currently gives a temporal coverage ending in **2020**. Treat these as publisher metadata access constraints, not proof that no exact 2022 data exist under authorized research access.
6. [European Butterfly Monitoring Scheme data access](https://butterfly-monitoring.net/index.php/ebms-data-summary) describes standardized multi-country monitoring data (30+ schemes). Restricted site/section coordinates and sampling-event participation could potentially meet the primary 2022 source gate if formally accessible, but **no such cross-country source set has yet been acquired or approved**. No fabricated sites, no unlicensed uses.

## Decision

**SEBMS_2022_EVENT_METADATA_ELIGIBLE_AS_CANDIDATE; PRIMARY_HNS_STOP_NOT_YET_8CELL_SOURCE_QUALIFIED.**

The **sole next allowed source-only operation** is directed river/outfall matching of 264 published point centers, plus independent 2022 event-location availability checks in at least two additional countries. Without demonstrated ≥24 non-overlapping eight-cell source/reference matched pairs, ≥12 genuine H6 and ≥3 countries under immutable frozen filters, no primary HNS calculation or ecological sewage-specific interpretation. This is a source opportunity finding only.

Keep PR #95 Draft, no new branch, no main or chocho merges, no eBird/TinyFish.
