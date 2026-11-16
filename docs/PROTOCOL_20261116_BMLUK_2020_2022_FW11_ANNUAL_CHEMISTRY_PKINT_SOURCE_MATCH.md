# 2026-11-16 original BMLUK 2020–2022 yearbook FW1 2022 chemical station source recovery (before opening means)

**Independent national publisher source, specified before reading workbook cells or target station chemical values.** The official [Wassergüte Jahresbericht 2020–2022](https://www.bmluk.gv.at/themen/wasser/wasserqualitaet/wassergute-jahresbericht-2020-2022.html) publishes a 587 kB ZIP with 2020–2022 results. Publisher report section 9 documents:
- **FW 1.1: Jahresmittelwerte und Perzentile für ausgewählte Parameter für das Jahr 2022 (Wiederholungsbeobachtungsjahr)** — actual 2022 annual chemistry summaries at some river stations. May contain phosphorus/nitrogen, but do not infer total TN/TP until field audit.
- **FW 2: GZÜV operative Fließgewässermessstellen 2019–2021** — NOT a verified 2022 biological EQR table; do not treat prior-year statuses as 2022.
- 2022 sediment/fish chemical results are not macroinvertebrate EQR site-year samples.

This original publisher source is distinct from EEA Waterbase ICM2026 Part1 and Part2: original 96,597,294 individual sample rows and 6,450,192 annual source rows both contained **ZERO 2022 exact original IDs** for 225 targeted GZÜV PKINTs near 117 positive 2022 annual-reported UWWTD N/P outfalls. EEA official WQ spatial archive matched only 10/225 national GZÜV points at <=10m, all SAME PKINT IDs; 18 at <=100m, not certified same station. Do not merge 2 different source schemes by geographic closeness.

## Source identification and eligibility, prior to numeric values

1. Recover source ZIP ONLY from BMLUK's official original 2020–22 page link. Verify actual URL, host, ZIP SHA256 and original member list. Stage relevant FW1 workbook temporarily for spreadsheet-source QA; no raw publisher workbook committed in Git, no original complete GZÜV GPKG/EEA files rehosted.
2. Inspect exact workbook sheet names, original column headers, station-ID column, date/year and units, whether each row is an annual mean or threshold. First source support and missingness before interpreting measured nutrient response.
3. Join by exact national GZÜV PK to original 225 selected PKINTs where PKINT = 'AT' + PK, verified crosswalk SHA256 040d610b6172d165787775406ea50c0f73f34db1db3fcb6c8576898736146dde. The 225 stations were preselected at radius <=2km of 117 2022 original N/P-report-positive outfalls; DO NOT infer same downstream receiving reach from distance.
4. Primary source availability gate: if workbook contains real 2022 annual mean TOTAL N (Gesamtstickstoff/N ges) in mg{N}/L and TOTAL P (Gesamtphosphor/P ges) in mg{P}/L at same original physical station, report source availability counts; do not treat nitrate/ammonium/phosphate-P as equivalent to totals. If totals are absent, report which distinct species exist and year/grain.
5. Preserve source row counts outside the original 225 (publication coverage), all 225 denominator and all 117 outfall denominator. Annual means are NOT proof of >=3 original independent sampling dates, not a source-specific delivered effluent dose, and not 2022 macrozoobenthos at the same site.
6. Do not fit any EQR, NDVI, plant NUTRI, butterfly HNS or ecosystem service effects from an annual report without first verifying same-year source benthos and directed outfall receiving river plus clean matched reference.

## Next scientific decision

If >=20 independent original 2022 outfall systems have actual measured 2022 river chemistry at valid stations, next verify genuine source NEXT_DOWN river delivery and actual 2022 macrozoobenthos dates with independent untreated river controls. If source rarity persists, use authorized public H2O/WISA five-step export or official helpdesk request for 225 original PKINTs. Missing EEA/annual report rows do NOT establish lack of national samples.

Preserve existing EEA Part1/Part2 original 225-ID STOP, 18/225 national FW3 2022–24 monitoring-list presence (only 1 chem+MZB/PHB both marked), and 68/68 exact GZÜV-EEA2022 benthic ID matches. National FW3 2022–24 and FW1 2020–22 have distinct grains; cannot be pooled into causal sewage linkage.
