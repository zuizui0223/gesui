# After official Part2 SOURCE STOP: separately frozen Part1 single-sample EEA TNTP source

2026-11-16, post-result extension. Original Part2-only source gate failed: 12,139 WQ stations in original 232 H7, 23,717 candidate station×year windows; strict Part2 TNTP valid only 54 station×year (of which 32 country-matched), giving primary 5 original H7×period/3 H7/2 H6/1 country. Existing frozen Part2-only STOP remains. Part1 is a DIFFERENT original source table and will not retrospectively relabel Part2 as comprehensive.

## Pre-Part1 original source and immutable source geography

- Preserve byte-for-byte complete ORIGINAL_2026_EEA_WQ_STATIONS_IN_ORIGINAL_232_H7_SOURCE_ONLY.csv from successful Actions #43816554985, SHA256 c5cde7bd8134b254c5e95b06849292fa260fd6f01832e6b6d2e1994a903c5910. Do not remap geographies or adjust geographic thresholds after viewing source chemistry.
- Official EEA 2026 Waterbase ICM Part1 DisaggregatedData 1,639,231,716 bytes, original published ZIP SHA256 ff779521749b40b952cce35a4671240c3b6c47d87f32a8fef9240a4ac503a082, ZIP DEFLATE64 method9 (stream with pre-installed 7z x -so). Original full CSV 96,597,294 source records, 28,458,399,518 bytes uncompressed (verified earlier PR #107). No full raw source rehost.
- Target original WQ station-year full key (country, monitoringSiteIdentifier, monitoringSiteIdentifierScheme, year 2018/2022) from the frozen geography. phenomenonTimeSamplingDate is the genuine Part1 sample date; its year must be exactly 2018 or 2022.
- River RW/water W; actual total TN EEA_31615-01-7 and TP CAS_7723-14-0; mg{N}/L, mg{P}/L. All source values numeric >=0 and finite, belowLOQ must not be true. Require >=3 distinct original water-sample records per nutrient in each WQ station and year. Use original sampleIdentifier for de-duplication; where absent, conservatively use original exact date+value. If same original sample key carries conflicting values, disqualify that nutrient×station×year rather than cherry pick. No nitrate or phosphate replacement.
- Compute simple original annual arithmetic mean per TN and TP from that same station's source valid samples. Do not merge separate TN-only and TP-only stations, do not set missing to 0.
- H7×period: equal station weight, median of official annual values across valid original TNTP-paired water stations in the same genuine H7 and year. Primary station country must match original GBIF H7 majority country, any-country only secondary source opportunity count.
- Do NOT condition source station/year choice on WWTP registry, W, EQR, plant NUTRI, butterfly HNS, previously published sign of slopes or actual model prediction.

## Source-only hard support gate

Pre-fixed exactly the same independent primary minimum as Part2: >=30 distinct H7×period, >=20 distinct original H7, >=20 true H6, >=3 WQ countries with same-station/year strict TN and TP. If insufficient, STOP_PART1_SOURCE_CHEMISTRY_SUPPORT, NO ecological effects. If passed, record source opportunity as such; any ecological response/causal WWTP claim requires separate pre-specified analysis.

## Combining with Part2

A Part1+2 source union is not two independent measurements: for each SAME station×year and each nutrient, prefer a strictly valid Part1 annual mean; only fill from Part2 when original Part1 missing or invalid and Part2 valid. Never average source table values, blindly combine station-year counts or switch preference based on matching biological outcome. Part2 rerun to derive exact original same station/year Part2 values is permitted only when reproducibly SHA-locked; the original Part2 site-year H7 period medians alone are NOT enough for deduplicated union. Report Part1 support and Part1+Part2 support separately and preserve the old failed Part2-only outcome.

## Interpretation

These data are ambient river TN/TP at official Water Quality stations inside the same native HydroBASINS L7 polygons as original GBIF species lists, NOT the EEA benthic EQR monitoring stations, not actual WWTP delivered dose, not soil/root-zone/leaf N and not observed feeding. WQ year 2018 is compared to three-year GBIF 2017–2019 and 2022 to 2020–2022; not simultaneous organism observation. An H7 polygon shares geography, not direct water→plant energy transfer. Any cross-compartment correlation must be treated as posthoc observational and evaluated with genuine H6 heldout groups and country sensitivity.
