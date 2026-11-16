# Source-only EEA 2026 WISE6 2022 water TN/TP for original 225 Austrian GZÜV PKINT stations (freeze)

Date: 2026-11-16. **Prespecified BEFORE opening any 2022 water nutrient results at the target 225 locations.** This is a source-availability test, not a biological response model or wastewater dose assignment.

## Exact prior authentic population

Official Austrian BMLUK GZÜV Fließgewässer original *FWMST.gpkg* file SHA256 4b87a495cd71bf2c413494b80158e8cfc91bb0de3e6e2559c186298388986fcf, 3,788 government river station points, official unique station ID `PK` and EU-interoperable `PKINT`. Exactly original **117** 2022 measured reported annual N+P EEA UWWTD positive outlet facility IDs have **225 unique** GZÜV station PKINT positions within 2,000 m. Source-only derived CSV from signed [Actions #43865663130](https://github.com/zuizui0223/gesui/actions/runs/43865663130): filename `BMLUK_GZUEV_2022_AT145_2KM_SOURCE_STATION_PKINT_MATCHES.csv`, SHA256 **040d610b6172d165787775406ea50c0f73f34db1db3fcb6c8576898736146dde**, 225 rows /225 unique original ATFW-prefixed site IDs /117 outlets. Every row has original matching EEA2022 WWTP measured N/P annual loads, and source distance; those annual facility loads must **NOT be treated as actual river TN/P concentration** or downstream delivery. The original 117 outlets themselves are a subset of genuine AT145 2022 annual N+P facility sources with authentic directed river-line geometry, not selected for outcome signs.

Original 2022 EEA benthic BQE in current main has **68 official AT sample-site IDs**, each matches GZÜV PKINT, 67/68 original location pairs within 5 m but one station `ATFW54110117` differs by 954.614 m; flag it and do not assert exact source location from shared ID. National WISA Benthosdata could contain more 2022 samples, but that raw source not yet verified.

## New independent actual WATER chemistry source, not WISA private database

Official EEA Waterbase WQ ICM July 2026 original public source used in previously source-verified PR107/111:
- **Part1** 2026 official original WISE6_DisaggregatedData-csv.zip 1,639,231,716 bytes SHA256 **ff779521749b40b952cce35a4671240c3b6c47d87f32a8fef9240a4ac503a082**, full true CSV **96,597,294** original individual rows, 28.46 GB expanded, archived in public EEA source Nextcloud WebDAV. Zip member DEFLATE64 method 9; genuine file streaming requires preinstalled 7z (stdout) after complete source SHA check. Original PR107 source parser and #43817344055 official Part1 source read completed previously.
- Year **2022** (exact original `phenomenonTimeSamplingDate` calendar year), true country `AT`, station `monitoringSiteIdentifier` = exact GZÜV `PKINT`, scheme **euMonitoringSiteCode** as the primary source (also report any other original schemes as uncertain, don't combine same named IDs across schemes), waterbody RW and matrix W. Keep *only* true total nitrogen original determinand **EEA_31615-01-7** in **mg{N}/L** and true total phosphorus **CAS_7723-14-0** in **mg{P}/L**; nitrates/phosphates are not substitutes.
- No lower-than-LOQ observations, negative/nonfinite/missing values, invalid units. Each nutrient×PKINT×2022 must have **at least 3 valid distinct original sampling DATES** and **at least 3 valid unique sample identities** (sampleIdentifier + exact date if available, else exact date+value). Preserve exact distinct sample dates, source scheme, invalid observation count, conflicts, official source digest and raw full-original row count. No user/observer personal information.
- **Do NOT save chemical concentrations in results** during this source-only phase: evaluate value validity and sample-count eligibility only. No EQR, NDVI, NUTRI, butterfly HNS, individual biomass/abundance coefficient. All missing values are UNKNOWN, never impute zero.
- Additional Part2 annual site-year report data may be checked separately, but Part1-only gate must be preserved as first independent source finding. Do not treat Part1 and Part2 summaries as independent sample replicates.

## Pre-locked scientific source support and explicit limits

Report per original source 225 PKINT sites:
1. any matched official 2022 Part1 record, and actual 2022 original TN/TP analyte sample-record counts, distinct sampled days and valid dates;
2. true 2022 both TN+TP ≥3 distinct source dates each at exact same station and correct scheme, count of such sites and **original 117 outfalls** with ≥1 such station within2km;
3. among those, original EEA2022 benthic-site identity source overlap n≤68, retaining the station with 955m national/EEA coordinate anomaly as uncertain (not physically co-located);
4. do not call these locations hydraulically downstream: they were selected by 2km great-circle distance, and actual river source-reach direction must be verified separately.

A minimum for **new independent causal multi-compartment design** remains **≥20 true independent original outfall–receiver systems across ≥3 native H6 parent basins**, with BOTH actual 2022 water TN/TP chemistry and independent real 2022 macrozoobenthos sampling, plus ≥10 genuine untreated comparison reaches. If field source support fails, document `STOP_2022_SAME_WWTP_WATER_BENTHOS_SOURCE_SUPPORT`, don't fit or relax. This 225-location water-chemistry subanalysis on its own does not prove national WISA MZB source presence or site-to-outfall hydraulic connectivity.

Public H2O public quality query may require legal copyright acceptance in its JS UI, so no undocumented API, account login bypass or legally binding consent on user behalf. Using freely downloadable official EEA original public data instead retains provenance and avoids that uncertainty. National ECOweb benthos database explicitly restricts direct raw access by roles; unverified records must be requested via official environment information access, not invented.

## Output

Source-only original report JSON, compact 225 original PKINT measurement opportunity CSV, original full source checksum and stage counts, no full source zip rehosted. Reuse precisely original crosswalk artifact run #43865663130 and signed EEA2022 benthic station IDs main original source Gitblob 325f199e985d9862ffee5b0edb205c61a7be3833. Original main/chocho unchanged and this stays a draft until reviewed.
