# EEA Water Quality 2026 vs original biological EQR stations: measured TN/TP source gate

2026-11-15. Source-only availability audit, posthoc and observational. NOT a causal wastewater, biodiversity, food-chain or SEM mediator effect. Draft PR #91 stacked on PR #90.

## Immutable original source population
- Same European 454 H7×period and 232 unique H7. Source 2018/2022 EEA Biology field station original mapping 3,525 site-year rows, 2,485 native physical site IDs (source original GH run #43165326244).
- Genuine positive upstream annual *reported* WWTP N/P load at only 132 original station-years, spanning 56 H7-periods, 38 H7 and 31 H6. The 56 H7-periods contain 538 EQR site-year records, of which 132 actually had positive reported upstream source. Unmatched report means unknown wastewater exposure, not pristine control.
- Source match requires exact monitoringSiteIdentifier + monitoringSiteIdentifierScheme + exact year 2018 or 2022, RW river, W water matrix, total N EEA_31615-01-7 or total P CAS_7723-14-0, sample units, quality flags and at least 3 valid samples. Nitrate and orthophosphate are different determinands.

## Genuine measured-source findings and negative cases
- [run 43676663949](https://github.com/zuizui0223/gesui/actions/runs/43676663949): original n3525/2485 identity source verified; first-party [WISE_SOE].[latest].[Waterbase_T_WISE6_AggregatedData] source schema accessible; complex multi-predicate SQL returned HTTP 403.
- [run 43677094638](https://github.com/zuizui0223/gesui/actions/runs/43677094638): single-predicate SQL works and official TN/TP codes return authentic Water Quality records. Year/site combinations queried independently; do not treat API 403 as missing biology.
- [run 43677482455](https://github.com/zuizui0223/gesui/actions/runs/43677482455): 28 originally selected H7×periods /45 native station queries, strict same-station 2018/2022 TN AND TP yielded **0 H7 periods**. Source compatibility pilot, not exhaustive original population.
- [run 43677752344](https://github.com/zuizui0223/gesui/actions/runs/43677752344): a more relevant 32 original station-years, of which 24 genuinely had reported-positive upstream WWTP N/P, gave strict source TN 0, P 0 and both 0. 29 of 30 distinct station WQ queries completed; one stopped at the 1,000 record cap (unknown, not absence). No matched positive source chemistry was established.
- [run 43678020210](https://github.com/zuizui0223/gesui/actions/runs/43678020210): posthoc row-stage audit of eight source-positive stations with previously nonempty Water Quality records. Seven complete station datasets contained *historical* total N and P, but **none had a Water Quality record in target year 2018/2022**. Years present: Croatia sites 2003–2012 (two 2022 targets), Hungary 1987–2006/2007 (two 2018 targets), Italy 2007–2008 / 2012 / 2010 (three target years). Eighth site truncated, cannot interpret a nonmatch.
- All source receipts contain original SQL request / response SHA and completion flag. These are posthoc selected samples, NOT a measured 0/132 full source-positive population count.

## Why the two public Water Quality source tables matter
- [EEA official Waterbase - Water Quality ICM 2026](https://www.eea.europa.eu/en/datahub/datahubitem-view/fbf3717c-cd7b-4785-933a-d0cf510542e1) published July 7 2026 and lists 1900–2025 coverage, but this global coverage cannot establish any specific 2018 or 2022 EQR station has measurements.
- [EEA official metadata factsheet](https://sdi.eea.europa.eu/catalogue/datahub/api/records/82e194cc-2506-411c-b341-9eb50b6ad5e0/formatters/xsl-view?approved=true&language=eng&output=pdf): Part 1 DisaggregatedData (individual samples), Part 2 AggregatedData (annual site mean), AggregatedDataByWaterBody and spatial records; reported original values in one need not be in the other.
- [EEA WISE SoE data collection portal](https://water.europa.eu/freshwater/resources/wise-soe-data-collection) warns Discodata WISE-6 tables are not updated; the SQL latest table cannot be assumed to equal full 2026 publisher dataset without source freshness verification.
- Official 2026 direct download is approximately 4.3GB and not downloaded as part of this bounded API preflight.

## Biological inference rules
- Pre-response feasibility rule fixed: >=30 source-positive H7×period with same-site measured nutrients, >=20 original H7, >=10 H6, >=3 countries. Current confirmed support in finite audited pilots is insufficient; do not fit HNS/EQR/S50 concentration mediation.
- A real concentration of total N or P in the river does not directly identify wastewater origin or the fraction of municipal effluent delivered to a butterfly microhabitat.
- Preserve original PR90 static W vs HNS spatial association and its null annual reported N/P added-value, without translating biological source differences into direct resource or ecosystem-service flow.

## Disaggregated individual-sample check
- A further official source-only query is running in [Action #43678184779](https://github.com/zuizui0223/gesui/actions/runs/43678184779). This section should be updated only after actual observed sample outcome is verified; no assumed success or absence.

## Subsequent source-only outcome: 2026 WISE6 individual samples and direct publisher inventory

- [GitHub Actions run 43678184779](https://github.com/zuizui0223/gesui/actions/runs/43678184779), verified SUCCESS for running actual source feasibility workflow: the official [WISE_SOE].[latest].[Waterbase_T_WISE6_DisaggregatedData] table schema exists. Four original EEA physical monitoring station IDs from the 2018/2022 EQR cohort were individually queried; **all four returned no usable results array at the source API layer** (source response digest repeated, as archived), so 4-site source collection was INCOMPLETE. **Do not claim these four stations have no water-quality samples.** No same-year TN/TP value could be verified.
- [Direct 2026 publisher file inventory run 43678478810](https://github.com/zuizui0223/gesui/actions/runs/43678478810) SUCCESS for source-access diagnostics: the exact official source WebDAV path via simple PROPFIND responds HTTP301; HEAD and ranged GET return a Nextcloud HTML listing (28,429 bytes), not a source ZIP or SQLite table. Partition filenames and their source sizes are **not yet retrievable through that request form**; a 4.3GB direct-download corpus was deliberately not fetched.
- Fresh official [EEA datahub Water Quality ICM 2026](https://www.eea.europa.eu/en/datahub/datahubitem-view/fbf3717c-cd7b-4785-933a-d0cf510542e1) says the full dataset covers 1900–2025 and is split into aggregated/disaggregated parts. The Water.europa.eu WISE SoE portal notes the older Discodata WISE6 table collection may not be updated. Thus the old historical years seen on Discodata **must not be extrapolated as missing 2018/2022 measurements in the full 2026 published corpus**.
- All original raw/source SHA receipts remain accessible as GitHub Actions artifacts: aggregated source run 43678020210 digest sha256:28eea9aa160e2a4f80458640d4d1cd76f20da25dddb83cd74a002d39a1b8710a; disaggregated run 43678184779 artifact 14147800269 digest sha256:f4e20d7f0e38c48bdbacc3cca531b5870cf752ae26769a2498aec0d5f77a35e6.

### Updated STOP and scientific interpretation

Original site matching and load reporting are verified, but the required same-year water TNTP pairing remains **not established**. No chemical outcome model, mediation, adult butterfly attraction, or causal wastewater effect can be claimed. The next source source-quality step is to retrieve the publisher's 2026 Part 2 original annual-site AggregatedData directly (possibly subsetted with source-partition support), check actual publisher source SHA and native IDs, and only then apply the precommitted >=30 H7-period /20 H7 /10 H6 /3-country support threshold. Do not reintroduce old pre-2012 chemistry as if observed in 2018/2022.
