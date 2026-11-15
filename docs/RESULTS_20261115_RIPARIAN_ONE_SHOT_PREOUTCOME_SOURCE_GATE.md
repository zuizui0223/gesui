# Source-only feasibility decision after frozen 2026-11-15 one-shot Lycaenidae HNS protocol

**Source-only status: STOP_DATA_OPPORTUNITY_NOT_YET_DEMONSTRATED.** The biological primary estimator was **not run**. This is NOT an ecological negative result.

## Successfully verified original-source files

- Frozen protocol Git blob: b9ce415c36d398eee77d8318a20aef4873e0bf44. It was verified byte-for-byte unchanged.
- DK source-only river/site network JSON original Git blob: 16a89e29f6bfe4aa6c39542b28f06e85c3880823, exact head of prior PR #42.
- DE original 73-pair opportunity audit Git blob: 1c1af1b5be5d061e8d7d2ca901f09156f7cc40e1, and pre-butterfly STOP JSON Git blob ec76f9a33f4021650c0a1af4cb0b8ddf971b8acf, exact head of prior PR #36.
- Original 454 H7-period HNS panel SHA256: 1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca. **Only field names and H7/H6 grouping columns were extracted**; HNS values were not calculated or inspected for the local contrast.

The complete re-execution of this audit succeeded via [GitHub Actions #43691928617](https://github.com/zuizui0223/gesui/actions/runs/43691928617). The source-only machine-readable QA is saved as [FROZEN_GATE_SOURCE_ONLY_RESULTS.json](../results/20261115_riparian_source_only/FROZEN_GATE_SOURCE_ONLY_RESULTS.json). The executable checker is [audit_frozen_gate.py](../analysis/20261115_riparian_source_only/audit_frozen_gate.py). Historical source-only files are fetched at immutable commits and verified against their original Git blob SHA-1. No HNS mediation or outcome model is imported or fitted.

## Actual original-source feasibility findings

| Source | Independent source support/QA | Why it does not qualify for 2022 frozen 8-cell primary |
|---|---|---|
| Danish standardized EBMS 2025, original 286 monitoring centroids and 203 EEA2022 outfalls | 92 <=5 km source/site relation candidates; 69 sites unsnapped; 4 same-reach order unknown. Strict downstream 10 sites **from only two outfalls and two river reaches**; strict upstream 4 sites from two outfalls/reaches | Butterfly observation year **2025**, not locked **2022**; 500 m centerline snap does not show 0–200 m *occurrence-level* distance or the 500–1500 m inland control |
| Germany original 2022 source+photo opportunity prior gate | 73 source directional river pairs, 100 sites. Prior 5 km Plantae photo effort: 35 pairs supported, but only **5 independently assembled four-arm geo-quartets**, below the former fixed support requirement; **zero 2022 butterfly outcome API requests** were made | Prior 2/5 km photographs and four-arm station pairs are **not the new 0–200 m/500–1500 m eight-cell protocol**. Exact observer-day support, river-distance band occupancy and matched reference 8-cell sets are **unknown**, not demonstrated |
| Euro H7-period GBIF HNS from PR #95 | 454 H7-period rows /232 H7 /90 H6; source genus/family counts available | 169 fields but no occurrenceID, exact occurrence latitude/longitude, eventID, recorder/date, coordinate uncertainty or reach-band assignment. Cannot recover locked 2022 May–September near/inland arms |

**Source-only STOP means only that no complete matched, verified 2022 eight-cell pairs are documented in these audited legacy source files. It does not mean that zero such pairs exist across GBIF or monitoring programmes.** Therefore the 24 independent matched source/reference pair gate has not been demonstrated, and the study remains in the acquisition stage.

The fixed 24-pair × 8-arm design requires **192 qualified cells**. Its independent necessary minima are **3,840 five-family georeferenced observation records** (20/cell), **960 cell-specific observer-date opportunities** (5/cell), and at least **384 scored Lycaenidae species-in-cell incidences** (2/cell). These are arithmetical **necessary, not sufficient**, and do not imply 960 globally independent observer-date identities. The design further requires ≥12 genuine H6 parents, ≥3 countries, river topology, nonoverlapping habitats, matching and score coverage. The counts above are not proof that any of these requirements can be achieved with opportunistic records.

## One forward acquisition lead, not a protocol change

Official [GBIF Swedish Butterfly Monitoring Scheme](https://www.gbif.org/dataset/be77e203-486c-4651-91b9-8347968b728c) is a public **sampling-event** dataset with records from 2009 onward (including 2022), point-count and transect events, dates, and sampling effort; dataset licence CC0. Its original records may have site centroids, not precise transect sections; the centroid cannot automatically satisfy a 100 m location-uncertainty gate or the frozen 0–200/500–1500 m riparian assignment. It is therefore a *candidate for source-only coordinate/visit coverage checks*, **not** an approved eligible unit or a valid third-country replicate.

Likewise, [Flanders event monitoring records](https://www.gbif.org/dataset/63938753-1fec-4c08-ae39-e9f8a6576521) may generalize sensitive site coordinates to **1/5/10 km cells**, explicitly incompatible with the frozen <=100 m per-record accuracy where so generalized. The [eBMS data access policy](https://butterfly-monitoring.net/ebms-data%20access) provides a scientific research licensing route for precise event/section information where the custodians permit it. Those data must still be traceable to eligible 2022 point-georeferenced GBIF records to be used in the **already frozen** endpoint; otherwise using them as a new endpoint would be a separate, explicitly revised design, not an unnoticed substitution.

**Do not broaden radius, swap year, relax pair count or observer-day criteria, or calculate 2022 Lycaenidae-HNS by source treatment group before a source-only independent access and eligibility manifest proves that the freeze gate passed.** Do not claim adverse or beneficial terrestrial resource transfer from this source-only outcome.

No new branch/PR, no merge to main or chocho, no eBird or TinyFish.
