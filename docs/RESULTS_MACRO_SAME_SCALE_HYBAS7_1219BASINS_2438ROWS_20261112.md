# Macro source-unified HydroBASINS level-7 two-period results (2026-11-12)

Source-first, no cross-taxon ecological effect estimates yet. All older confirmed numerical results remain unchanged.

## Actual original data linked to one hydrological geography

Official HydroBASINS L07 Europe and North America archives obtained from data.hydrosheds.org and SHA256 audited. Original study coordinates and no numerical biological responses were used.

| Original measured source | Original rows | Valid HydroBASINS L07 assignments | Unique L07 subbasins |
|---|---:|---:|---:|
| EEA Europe 2018 numeric benthic EQR / reported WWTP N/P | 2634 | 2422 | 406 |
| EEA Europe 2022 numeric benthic EQR / reported WWTP N/P | 2475 | 2224 | 329 |
| EPA NRSA USA 2018-2019 native UID actual TN/TP/benthic sites | 1221 | 1173 | 718 |

Europe 2018/2022: 234 common subbasins with monitoring in both periods. Of 1167 exact repeated original river station keys, 1064 map to the SAME valid HydroBASINS L07 polygon; 0 mapped to different valid polygons; 103 lack at least one valid assignment.

Source geometry warnings: 336 EU and 367 North America genuine source polygon features were initially invalid and skipped WITHOUT silent repair. As a result, 212 EU2018, 251 EU2022 and 48 US source records remain unmatched. These must remain missing for any future effect analysis unless independently source-corrected before outcomes; do not treat them as EQR/biological zeros.

Official EU file zip SHA256 b9b810f4db271f0fa6b39fecf01784041b948944b3c3f296a393766388e86b44. Official NA file zip SHA256 4e5150f0fe3a2c3d17b21938ec5c6ca1959d136497f572a0203af6b1842e501d.

[Original EEA + EPA + true river basin source QA](https://github.com/zuizui0223/gesui/actions/runs/43165326244)

## Same analytical unit generated: actual HYBAS_ID_L07 × period

Fixed 3-year bins: 2017–2019 (observed EU2018 + USA2018-19), 2020–2022 (observed EU2022; US not observed).

Original 1219 distinct mapped subbasins × 2 periods = 2438 generated master table rows, including explicit missing combinations. 1453 basin-period rows have a real river biology source: 735 Europe (406 early /329 late), 718 US (early only). 234 Europe subbasins are observed twice.

This is a physical geographic SOURCE table, NOT an ecological outcome table. All current plant, butterfly, adult emergence, pollination and unmapped future wastewater exposure response fields are NULL rather than absent species or measured zeros. We have not created actual fully measured basin-period records containing bottom fauna, plant taxa, butterflies and source-specific WWTP loads together.

[Reproducible two-period native HydroBASINS macro source master CSV and provenance](https://github.com/zuizui0223/gesui/actions/runs/43165738190). Master CSV SHA256 a9a8eaa5c328026e46118810f9fbde29cc24451a594ad9437ba2222201bf5e0a.

## Critical original global WPI port 2015 vs 2018-2022 inland aquatic geography result

Archived worldwide 1055 port point locations, original EU/USA river stations, and the SHA-pinned US original UID source were all assigned to the *same H3 r5 cover* (roughly ~250 square km H3 hexes, actual area varied). EU2018 1730 hexes; EU2022 1588; US NRSA 1066; port WPI 974; EEA2018/2022 overlap 944 hexes.

Coastal WPI port centers share ONLY 6 r5 hexes with EEA2018 river records, 7 with EEA2022, and 0 with the US NRSA panel. Co-location is NOT downstream water delivery and port original GBIF 2011-2024 observation facets are not period-specific or basin-specific.

**Conclusion:** original 1055-port 2015 WWp plant-photo analysis CANNOT be used as the same-site inland macro river-plant-butterfly panel. Preserve it as a distinct coastal study; for a genuine common-basin 2018-2022 macro study need fresh plant/butterfly records in each actual river catchment and correct time bin.

[Original equal-spatial-support source-only H3 audit](https://github.com/zuizui0223/gesui/actions/runs/43164858721).

## Fixed next analytical operations — no coefficient seeking

1. Source-audit HydroBASINS invalid polygons and geographic coverage by method/state. Try make_valid as a documented source-only sensitivity, never overwrite the original no-repair result.
2. Worldwide sewage exposure: georeference 58502 registered global HydroWASTE WWTP estimated outfalls (static mixed measured/imputed flow) to the *same* HydroBASINS L7 IDs and upstream HydroRIVERS reaches. Keep annual real 2018/22 Europe reported nutrient discharge separate from modeled flow. A nearby plant is NOT necessarily an upstream source.
3. Native GBIF plant and butterfly occurrence source queries MUST use each same HYBAS_ID_L07 polygon and original fixed periods 2017–2019 / 2020–2022. Original port 5 km GBIF facets cannot be cloned as new basin data. Count independent accepted species and use observation-effort corrected/rarefied plant richness, potential host species, and Lycaenidae / Papilionoidea standardized relative reporting. Do not use GBIF zeros for failed requests; preserve source quality and provenance.
4. Fit a hierarchical multivariate macro model only on actual shared basin-period support: environmental WWTP intensity or measured TP/N, within-method EU EQR, separate within-method US EPT, plant candidate resource richness and butterfly presence/detection with different observation likelihoods. National assessment methods, HUC8 river-order and observation effort must be controlled. EQR and EPT fractions are NOT biologically interchangeable even if both numerically between 0 and 1.
5. Predeclare HydroBASINS L7 primary, L6 spatial scale sensitivity, and holdout by independent basin/country, not row. Evaluate whether directional signs and source opportunities persist when spatial scales differ. Never convert statistical association into serial ecosystem-service causation without measured adult emergence/predator diets or seed set on shared catchment/time units.

No new biological effects have been calculated in this branch. No merging into gesui/main, chocho, previous positive/negative PRs. No eBird.
