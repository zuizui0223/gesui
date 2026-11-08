# Source-locked legacy chocho Lycaenidae findings (NOT wastewater results)

This audit uses only 239-species frozen resource geography and its independent family summary at exact chocho commit e8b8f86c7a889b6a3224a4ef69ffe5e0c7738660.
The two original Git blobs are checked byte-for-byte using Git SHA-1. No new GBIF or EEA ecology outcome is queried.

| Summary | Existing evidence |
|---|---:|
| Lycaenidae with introduced-host expanded potential resources | **22/25 (88.0%)** |
| Family median additional resource regions | **73** |
| All butterfly species with expanded potential resource opportunity | **206/239** |

| Illustrative species (posthoc) | Native-only regions | Contemporary including introduced hosts | Added | Original host-taxonomy quality |
|---|---:|---:|---:|---|
| *Lycaena phlaeas* | 140 | 224 | +84 | adequate |
| *Polyommatus icarus* | 128 | 207 | +79 | adequate |
| *Brephidium exilis* | 242 | 267 | +25 | adequate |
| *Pseudozizeeria maha* | 16 | 275 | +259 | NOT adequate |

**Interpretation:** These are potential plant-host resource regions, not actual butterfly population change or wastewater responses. The extreme *Pseudozizeeria maha* ratio (16→275) has a non-adequate host taxonomy support flag and should **not** anchor the main claim.

Unresolved: whether local nutrient availability, adult organic nitrogen puddling, ant mutualism, host quality and effluent contamination interact. Every such connection is an **unmeasured hypothesis** at this stage.

Original sources:

- https://github.com/zuizui0223/chocho/blob/e8b8f86c7a889b6a3224a4ef69ffe5e0c7738660/data/frozen/figure_sources/anthropogenic_species_metrics.csv
- https://github.com/zuizui0223/chocho/blob/e8b8f86c7a889b6a3224a4ef69ffe5e0c7738660/provenance/reviewer_defenses/results/butterfly_taxonomic_family_expansion_v0.1.json
