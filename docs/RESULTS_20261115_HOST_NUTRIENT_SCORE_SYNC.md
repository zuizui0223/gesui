# Butterfly host nutrient-affinity scores — source-verified local result sync

Computed 2026-11-15 from the frozen Ferrer-Paris butterfly larval host source, gesui candidate links, chocho original230 host catalogue and Tichý et al. European Ellenberg-type N (EIV Zenodo v2, DOI https://doi.org/10.5281/zenodo.7427088, CC BY 4.0). Full source retrieval GitHub run https://github.com/zuizui0223/gesui/actions/runs/43663192859.

**Estimand:** HNS_s = mean of numerical EIV N over distinct source-known larval host plants. Equal weight per known host; not field host abundance, leaf nitrogen, adult puddling, measured wastewater exposure, response sign, or evolutionary adaptation.

| Catalogue | All roster rows retained | Numeric HNS available | ≥80% documented-host N coverage |
|---|---:|---:|---:|
| gesui original GBIF speciesKey | 303 | 190 | 63 accepted species keys |
| chocho original frozen names | 230 | 147 | 25 names |

63 acceptable gesui species keys with >80% host N coverage do not represent the 303-species original population; in H7 source data they account for only 13.676% of all raw butterfly source records. HNS coverage must be retained and missing taxa never assigned zero or no effect.

| Example accepted gesui taxon | HNS | Numerically matched / documented host concepts |
|---|---:|---|
| Aglais urticae | 8.20 | 2 / 2 |
| Araschnia levana | 8.20 | 2 / 4 |
| Lycaena dispar | 7.15 | 2 / 3 |
| Pieris napi | 6.3773 | 22 / 33 |
| Pieris rapae | 6.0111 | 36 / 56 |
| Lycaena phlaeas | 4.525 | 4 / 6 |
| Polyommatus icarus | 3.1714 | 7 / 11 |
| Plebejus argus | 2.30 | 2 / 2 |

Interpretation: different species within Lycaenidae have very different host nutrient-affinity. This does NOT assign positive/negative/no-effect biological sewage response classes. Adult mineral/protein/Na feeding is a separate axis; **301/303 adult records were not curated** and missing behavioral evidence is unknown, not negative.

Original Tichý source Excel MD5 3519d41639647dde8d80042fc98b2ef5; SHA256 b3fa7c3a6e9d356e52a45472549f807967bd0c3189e58e7e8a22e0389c29cc4e. Original external indexed 624 numerical NUTRI matches agreed exactly (max absolute difference 0). Source included 8,908 taxon rows and 7,195 numeric N values; source publication summary 7,193 discrepancy is retained, not silently corrected.

Critical source limitation: Aglais io has one extracted host Humulus lupulus, despite known Urtica dioica use in other experimental literature. A 100% within-catalogue trait coverage is NOT completeness of biological host knowledge. Global chocho host-plant source also contains non-European hosts; do not impute European local hosts.

Local reproducibility ZIP SHA256 6f3927afa5bd0f489b8e2ed2f437dee6812a29ef4edda3e39e389a89a6feabac (17,076,785 bytes); numeric gesui CSV SHA256 00c1a36a9237a3d1a6d77334aa8baa5dc8b329a49925c2e31f033fa2d302ecdb (256,411 bytes). See the accompanying SHA ledger for full filenames and original source endpoints.

Local fresh 2026-11-15 test verification: 13/13 host score tests pass. Numeric score examples above were inspected from original accepted-key CSV; not a new global field study.
