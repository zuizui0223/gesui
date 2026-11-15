# 2026-11-15: Denmark and Belgium 2022 public butterfly sampling Event-core originals

**Status: ORIGINAL SOURCE-ONLY STOP. No butterfly species/HNS counts accessed.** The single frozen riparian 2022 source/reference biological protocol (Git blob b9ce415c36d398eee77d8318a20aef4873e0bf44) remains unchanged. No opportunistic relaxation of year, eight-arm support, study area, distance bands, sample count, HNS definition or source control.

[Reproducible original-source GitHub Actions #43696277771 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43696277771). Original Event core member event.txt and its meta.xml metadata ONLY; original Occurrence and MeasurementOrFact extensions deliberately not opened. Full source [numeric QA JSON](../results/20261115_riparian_source_only/DK_BE_2022_ORIGINAL_EVENT_ONLY_GATE.json) and [exact event-only reader](../analysis/20261115_riparian_source_only/audit_dk_be_2022_original_event_only.py).

| Actual official original source | Denmark eBMS v1.4 | Belgium Flanders Meetnetten priority butterflies v1.40 |
|---|---|---|
| GBIF dataset | [Denmark transects](https://www.gbif.org/dataset/caa38f94-ddeb-43e7-9a27-1ac46d7ddc5c) | [Flanders priority monitoring](https://www.gbif.org/dataset/63938753-1fec-4c08-ae39-e9f8a6576521) |
| Exact published original archive SHA256 | 2cc64ec352754223db06a92eef5b01803968eac94bae34057c17f1abdb6da79b | e0563615e6fa90dbfcdaeed5400fe87a11a6f30fb6ad9365b3c4a7e7ea803a34 |
| Native Event-core member SHA256 | 06b582ac3e640e4e18b81a45423713703d7524788f75b6dcfae70a74e2ad8da7 | 425cd46663d139895c202bdd73c17dfb410b1d316443e500080d27fff1653996 |
| All original Event-core rows | 15,606 | 40,233 |
| 2022 full calendar year | **125** | **3,573** |
| **2022 May–September** | **125** | **3,491** |
| Approximate published coordinate groups | **22** | **59** |
| Coordinate groups with five or more distinct summer event dates | **14** | **16** |
| Events with published uncertainty <=100m | **125 / 125** | **0 / 3,491** |
| Published coordinateUncertaintyInMeters | 50m | 707m (2,038); 3536m (917); 7071m (536) |
| Published samplingProtocol | **Missing all 125** | butterfly transects (3,491) |
| Event ID | Present all 125 | Present all 3,491 |
| Published dataGeneralizations and informationWithheld | Neither available as an original Event-core field | **Both populated 3,491 / 3,491** |
| samplingEffort and samplingSizeValue as original Event fields | Missing | Missing |
| Frozen 2022 0–200m river-near eligibility | **NOT VERIFIED** | **PUBLIC POSITION PRECISION FAIL** |

## Two real metadata and design findings

1. The Danish eBMS public description states the project began in June 2023, but the original publisher v1.4 Event-core **really has 125 events dated to 2022**, plus older records (2011 and 2018–2021). The existence of these events is now proven by original source bytes, **not evidence they followed a common standardized 2022 butterfly survey protocol**. Original 2022 samplingProtocol is completely missing; coordinate uncertainty 50m can describe an event/transect centroid rather than exact butterfly or section position. The 14 coordinate groups with five dates therefore cannot be promoted to valid near-river 200m spatial study sites. Earlier blanket statements that no 2022 Danish event rows exist were inaccurate.
2. Flanders priority GBIF registry temporal coverage ends 2023, whereas the exact 2026 source Event-core contains 2024 (3,863 events) and 2025 (3,638 events) as well. Registry dates lag actual archive rows. In 2022 specifically, all 3,491 summer event records have public UTM **generalized** positions from 707m to 7,071m uncertainty. They cannot resolve the 0–200m versus 500–1500m primary bands with <=100m source position accuracy. Also, this monitoring frame intentionally targets priority species; other butterfly species are incidental, so taxon-composition sampling opportunity needs independent review before any comparison with Sweden or Germany.

## Wider registry and 2022 access implications

A separate [11-dataset GBIF registry metadata-only audit #43695886459 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43695886459) verified source types, licence, current publisher dates and archive endpoints. [Original registry JSON](../results/20261115_riparian_source_only/GBIF_2022_OTHER_COUNTRY_DATASET_REGISTRY_SOURCE_ONLY.json) has the receipts. The Finnish agricultural butterfly dataset updated 2026-10-24 is an **Occurrence**, not Event-core dataset (event coverage requires separate proof). The Dutch nationwide NDFF butterfly monitoring entry is **Metadata** only; the 2022 scheme existed but exact event locations and visits require the provider. UKBMS public GBIF is an Occurrence dataset with registered temporal coverage to 2020. Older Flemish all-butterfly transect events end 2020, different from the Flanders priority source. The older Spain Sierra Nevada Event archive ends 2020. Publication year must never be taken as an observation year.

## Eligibility decision

**STOP_SOURCE_OPPORTUNITY_NOT_VERIFIED**. Neither DK nor BE contributes even one newly verified eight-cell source/reference 2022 pair. No 2022 freshwater outfall routing, exact same-H7 reference matching, genuine H6 grouping, source-side larval host scores, adult sodium/protein drinking, or ecological effects were calculated. A centroid plus an event ID is not a point observation near a wastewater-receiving river.

The frozen inferential study requires >=24 independent matched source/reference pairs, >=12 genuine H6 parent basins, >=3 countries, eight locations per pair, at least five valid survey dates and 20 five-family records per location, two scored Lycaenidae species, 70% trait score coverage and exact spatial/directional source quality. The unmatched and generalized published sources are a **data-availability limitation**, not evidence of no terrestrial transfer or negative butterfly response.

The next legitimate operation is **custodian permission to original precise 2022 point / 50m-section river position data and survey protocol/effort**. See [consent-first eBMS / German TMD / Dutch NDFF request prepared but NOT sent](20261115_EBMS_2022_SOURCEONLY_DATA_ACCESS_REQUEST_NOT_SENT.md). Only after legally authorized point-level source/effort eligibility could the one-shot HNS outcome run. All results remain Draft PR #95, no new PR, no merge to main or chocho, no eBird/TinyFish.
