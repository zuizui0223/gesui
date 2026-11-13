# HydroRIVERS receiver-specific WWTP source exposure: fixed source-first design (13 November 2026)

**New posthoc physical routing source audit BEFORE any new receiver-level ecological numerical outcome model.** Preserves existing PR #67 232 genuine EU HydroBASINS L7 basins / 454 basin-periods and all previous EQR, plant S50, Lycaenidae, HydroATLAS results.

## Original source evidence and exact geography

- Actual original EEA biological **2018** `benthic_invertebrates_2018_WWTP_source_real_EQR.csv` and **2022** `BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv` both include **original monitoring-reach `HYRIV_ID`**. Source original Git blob SHA1 2018 `0f30bfbc3dac13f6a182c8a910e9899c04778322`; 2022 `325f199e985d9862ffee5b0edb205c61a7be3833`. Monitoring station IDs are mapped to the same fixed authentic H7 basin/year cohort through original source-only [run #43165326244](https://github.com/zuizui0223/gesui/actions/runs/43165326244), SHA256 `a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89`.
- Genuine HydroWASTE original v1 `HydroWASTE_v10.zip` official author/WWF linked Figshare source [source audit #43198783733](https://github.com/zuizui0223/gesui/actions/runs/43198783733), original ZIP SHA256 `6e19d8f31d95d533b2b29056b823eeb915d933fa22ea090f234dc38f6427d044`, original CSV SHA256 `af040d9b4fdc93bab2eb8d1b52985f5bacc6b0db9f16795c29e31c10fb615dd4`. There are 58,502 registered facilities, 58,123 nonmissing `HYRIV_ID`, `WASTE_DIS` **m³/day** estimated/observed with `QUAL_WASTE` provenance; operational status not reported for most facilities. Static, **not** actual 2018/2022 annually measured effluent N/P or active operation.
- New original official [HydroRIVERS v1 Europe/Middle-East SHP](https://data.hydrosheds.org/file/HydroRIVERS/HydroRIVERS_v10_eu_shp.zip) source, ~68 MB ZIP. Required authoritative native `HYRIV_ID`, `NEXT_DOWN`, `LENGTH_KM`, `DIS_AV_CMS` field names verified at source. No replacement with basin center or 2015 port point. Never share raw standalone HydroRIVERS source ZIP or its DBF file publicly; only derived receptor exposure and source hashes (HydroSHEDS license redistribution constraint).
- Original identical **454** 2018/2022 H7×period source study period rows (232 native H7) from [real three-response source #43215754868](https://github.com/zuizui0223/gesui/actions/runs/43215754868), originally 3,525 EEA monitored station-year source entries. Do not choose a 234/218/170 river cohort by looking at exposure signs.

## Source-only stage 0: exact HYRIV reach topology, no ecological outcomes

Download official HydroRIVERS EU v1 shapefile ZIP (not global and not RiverATLAS). Read native DBF attributes, native reach IDs, `NEXT_DOWN` directed downstream chain; validate unique `HYRIV_ID` with no cycles, source flow/length, original EEA 2018/22 station-reach IDs and original HydroWASTE outfall `HYRIV_ID`. Source receiver reach ID must exist in this exact network. HydroWASTE global facilities not in EU network must be marked **outside_source_region**, not dropped without counts. Do not inspect or use any EEA EQR, Plant S50 or Lycaenidae magnitude/sign in source stage.

For every unique monitored receiver reach:
1. Compute **strict upstream** cumulative original HydroWASTE source counts and `WASTE_DIS` at *other* river reaches that flow to the receiver's reach using `NEXT_DOWN` chain, **excluding same-reach facilities** because position within the reach is unresolved.
2. Compute **strict upstream <=50 km full-reach length** count/discharge by traversing reverse directed reach edges, counting the complete original length of each upstream source reach; this is a conservative *source opportunity* measure, not travel time, dilution, or guaranteed delivery.
3. Compute source registered **operational-only** (source STATUS explicitly Operational) and original `QUAL_WASTE=1` (reported actual effluent) separate pressure, without silently excluding 90% status Not Reported.
4. As a directional **negative control**, count outfalls on distinct reaches downstream within 50 km of the monitoring reach. They cannot be a simple upstream source under a correctly directed network; however unequal hydrology, river size, urban concentration and source reporting make this a geographic comparator, not strict randomized negative control.
5. Candidate discharge fraction = strict upstream sum `WASTE_DIS` m³/day divided by original monitoring `DIS_AV_CMS × 86400` m³/day (both long-term model estimates). Exceedance >1 must be retained and QC-flagged, not artificially capped at 1. No nutrient measured from the `WASTE_DIS` field.
6. Aggregate site-equal averages across original distinct physical source stations inside **the SAME original 454 H7×period rows**, keeping original annual EEA 2018 and 2022 station inventories, source status missing distinct from zero, same source IDs in repeated years yield static exposures. Never call this a dynamic 2018→2022 discharge intervention.
7. Store native source station-reach id, original H7/time, source site counts and source reported quality counts, derived per-station georeferenced exposure only (no raw HydroRIVERS/HydroWASTE redistribution), source original ZIP/CSV digests and exact model dataset SHA.

**Predeclared physical source gate**: ≥90% of the original 454 H7×periods with ≥1 valid genuine HydroRIVERS receiver reach, ≥95% of all original 2018/22 station-year HYRIV_ID in the fixed 454 panel present in original HydroRIVERS network, at least 150 H7×periods with ≥1 **strict upstream distinct-reach** registered original source, and ≥75 distinct original H7 with upstream source support. If fail, **STOP WITHOUT ECOLOGICAL REGRESSION**; keep original results. Any row not passing routing remains missing, never no-WWTP/zero. No chosen 50km radius revision after outcome signs.

## Possible stage 1, only after source0 PASS

On **exact original 454 rows** and original same outcome values, test incremental receiver-specific source variables compared with **prior original Level6 source parent heldout** baseline after HydroATLAS upstream crop/urban/natural flow, EEA method and country and GBIF record effort:
- original all-upstream static registered outfall discharge fraction; original strict 50km fraction;
- documented operational-only and original direct-reported-flow quality sensitivity;
- downstream 50km direction-control source intensity.
- EEA actual benthic EQR separately from GBIF source-indexed rarefied plant S50 and Lycaenidae relative five-family indexed records. EEA 2018 vs 2022 annual reported N/P should stay a separate source-specific layer; no double counting or substituting HydroWASTE estimated flow for TN/P.

**Critical publication interpretation:** source geography alone (river topology, WASTE_DIS modeled flow and source status) cannot identify source nutrient delivery, contaminants, biological causality, spatially transferable response or pollination/consumer benefits. Country-heldout transfer failure of original PR #67 must be explicitly retained. New original country checks are POSTHOC, never prospective confirmations.

No eBird; leave all original gesui/main, chocho and prior source-positive/null PRs unchanged. New work remains unmerged Draft.
