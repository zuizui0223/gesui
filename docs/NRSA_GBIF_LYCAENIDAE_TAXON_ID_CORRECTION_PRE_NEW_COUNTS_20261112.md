# 2026-11-12 GBIF Lycaenidae taxonomy-source correction — before any revised butterfly counts

**Identified error:** the previous NRSA-GBIF360 workflow used `taxonKey=47923` for Lycaenidae, but 47923 is **iNaturalist**'s taxon ID, not GBIF's. The legacy GBIF Backbone taxon ID for Lycaenidae is **5473**; current GBIF COL XR has alphanumeric identifiers and is a distinct taxonomy. Verified via [Wikidata taxon cross-reference](https://www.wikidata.org/wiki/Q158717) and [GBIF's July 2026 migration guidance](https://data-blog.gbif.org/post/catalogue-of-life-taxonomic-backbone/). The initial **342/342 numeric zeros** for "Lycaenidae" cannot be treated as biological absences or true source coverage.

**Unchanged source geography and a priori eligibility:**
- Exact original *360 NRSA US river sites*, originally selected without the nutrient/EPT effects or new GBIF counts.
- Original 2017–2021 and 5-km source polygons, existing exact original individual request geometry, `limit=0`, `hasCoordinate=true`, `occurrenceStatus=PRESENT`.
- Preserve original immutable Plantae and Animalia GBIF count responses, exact URLs and SHA: no new source request for these taxa and no re-selection.
- Only re-query the **original fixed 360 Lycaenidae source URLs**, replacing the invalid `taxonKey=47923` with verified legacy GBIF **`taxonKey=5473`**, with explicit GBIF Backbone `checklistKey=d7dddbf4-2cf0-4f39-9b2a-bb099caae36c` to prevent unannounced COL XR mixing. Verify original GBIF species taxonomy name and rank BEFORE any new response.
- Missing HTTP requests remain missing, not zero. No nearby radius, year or region selection rescue.

**Preserved initial gate:** at least 250/360 original sites all 3 taxon queries succeed; >=100 of those have >=50 original Plantae indexed counts **and** >=20 corrected Lycaenidae indexed counts; those >=100 sites span >=15 US states and >=5 of the original nine eco regions. If fail, STOP: no aquatic/terrestrial response regression. If pass, freeze the exact eligible site IDs and source receipt hashes before first cross-taxon regression; comparison with the already known EPT result will be explicitly posthoc, not a new independent experiment.

**Do not use `iNaturalist taxon 47923` as a GBIF taxon key again.** Distinguish GBIF source record counts from larval host availability, adult butterfly density, pollination, wastewater effects or measured services. Note 57 original query errors (not zeros): 18 Lycaenidae, 15 Plantae, 24 Animalia; failure audit source [Actions 43047272367](https://github.com/zuizui0223/gesui/actions/runs/43047272367). Initial incorrect ID source [Actions 43036787547](https://github.com/zuizui0223/gesui/actions/runs/43036787547).

No change to source-frozen real EPA EPT/TN/TP or earlier EU EQR. Keep Draft, never merge automatically.
