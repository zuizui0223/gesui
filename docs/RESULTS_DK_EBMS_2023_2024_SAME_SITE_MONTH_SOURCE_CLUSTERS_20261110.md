# Denmark eBMS 2023–2024: same-site same-calendar-month Lycaenidae source-distance sensitivity

**2026-11-10. Strictly POSTHOC.** Motivating 2022 Germany/Austria positive iNaturalist photography and 2023/24 Denmark count pattern were already observed. All new numbers below are standardized butterfly survey records, **NOT** a wastewater treatment experiment, confirmed effluent delivery, adult feeding, or population response.

## Verified sources and sample

- Denmark eBMS v1.4 Darwin Core Archive (GBIF caa38f94-ddeb-43e7-9a27-1ac46d7ddc5c), **CC-BY-NC 4.0**, original ZIP SHA256 `2cc64ec352754223db06a92eef5b01803968eac94bae34057c17f1abdb6da79b`.
- Original **15,606 Event rows** and **21,617 Occurrence rows**, replayed by event/date/site and original frozen taxonomy; **all 214 original site-year accepted-family totals and visit-day counts match** the previous independently audited panel.
- 107 public site-coordinate groups with at least two visits in both 2023 and 2024 and original positional uncertainty <=500m; 23 near (<=5km) vs 38 far (>=10km) from qualifying 2022-reported Danish EEA outfalls. No eligible <=2km exposure cohort. **Only 4 distinct registered outfalls underlie the 23 near site rows.**
- Retain only calendar months where the **same site** was surveyed in both years. Denominator is original distinct sampling dates, not number of original Event record fragments. Full same-site-month window and separate May–August sensitivity are both posthoc descriptive stress tests.

## Calendar-month matched original butterfly counts

| eBMS | Near <=5km: 2023 | Near <=5km: 2024 | Remote >=10km: 2023 | Remote >=10km: 2024 |
|---|---:|---:|---:|---:|
| Repeated site-coordinate rows | 23 | 23 | 38 | 38 |
| Source-matched original visit days | **95** | **134** | **372** | **346** |
| Lycaenidae individualCount | **73** | **20** | **386** | **437** |
| Lycaenidae per distinct visit day | **0.768** | **0.149** | **1.038** | **1.263** |
| Pieridae individualCount | 125 | 83 | 344 | 329 |
| Nymphalidae individualCount | 373 | 457 | 1,817 | 1,477 |

Near Lycaenidae 2024/2023 survey-day-rate ratio **0.1942**; remote **1.2172**. Near divided by remote relative annual rate ratio **0.1596**. For family-composition share (Lycaenidae relative to the 3 target families): 2023 near 73/198=36.9%, remote 386/730=52.9%; 2024 near 20/103=19.4%, remote 437/766=57.0%. **No Papilionidae were recorded** among these target families, making this effectively a comparison with Pieridae.

### Independent family contrasts

| Family | Near annual rate ratio | Remote annual rate ratio | Relative change |
|---|---:|---:|---:|
| **Lycaenidae** | 0.194 | 1.217 | **0.160** |
| Pieridae | 0.471 | 1.028 | 0.458 |
| Nymphalidae | 0.869 | 0.874 | 0.994 |
| Hesperiidae | 0.722 | 1.723 | 0.419 |

Fixed May–August matched-date sensitivity: Lycaenidae near ratio 0.194, remote 1.231, near/remote 0.158. The negative recorded relative change persists in this date-window sensitivity; it is not unique to Lycaenidae, because Hesperiidae also fall relatively, although Nymphalidae's relative change is near unity.

### Source-group pseudoreplication and instability

Near sites belong to only four registered outfall groups with 4, 7, 5 and 7 public site rows. Lycaenidae matched-month individual counts by source (2023 to 2024) are 0 to 1; 0 to 0; 44 to 17; 29 to 2. Only **two source clusters** had any Lycaenidae count in 2023, so the decline is not an independently repeated signal across four occupied populations. Omitting one outfall group at a time gives an observed near/remote relative annual rate-ratio range **0.048–0.237**, preserving the pooled direction but not implying significance or independence.

## Interpretation and limitations

**Supported:** descriptive temporal heterogeneity in standardized butterfly reporting, stronger Lycaenidae decline per site survey day in four source-proximate clusters than in sampled remote areas, persisting after matching months and in a May–August sensitivity. **Not supported:** actual wastewater-induced decline, direct site nutrient exposure, butterfly puddling/protein intake, larval host or ant mediation, or generalization to additional WWTPs. The EEA discharge source year 2022 is a *static* infrastructure proxy in both 2023 and 2024, not a measured change in treatment; near sites are spatially nonrandom, real river direction/habitat dose, exact survey dates/weather, land use, host distributions and within-month phenology are not controlled. Four near outfalls and only three shared 1-degree comparison regions cannot support a confirmatory test.

Original full-year PR #37 reports 2023 near/remote 74/387 and 2024 near/remote 34/494 Lycaenidae individualCount. The new matched-month counts are a distinct posthoc estimand; **both are retained without rewriting original results**. Next evidential threshold: geographically and hydrologically connected standardized monitoring at at least 10 genuinely independent source systems with matched receiving-water chemistry, host vegetation, flight season and non-source controls.

## Reproducibility

- [Original licensed monthly replay and complete aggregate result](https://github.com/zuizui0223/gesui/blob/analysis/dk-ebms-2023-24-season-aligned-families-20261110/results/ebms_denmark_2023_2024_WWTP_source_only/DENMARK_2023_2024_SAME_SITE_MONTH_LYCAENIDAE_POSTHOC_FAMILY_SOURCE_CLUSTER_QC.json)
- [Independent count/date/source-cluster verification](https://github.com/zuizui0223/gesui/blob/analysis/dk-ebms-2023-24-season-aligned-families-20261110/results/ebms_denmark_2023_2024_WWTP_source_only/DENMARK_2023_2024_SAME_MONTH_ORIGINAL_FAMILY_INDEPENDENT_REPLAY_QC.json)
- [GitHub Actions original v1.4 replay SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42762375670); [independent offline arithmetic audit SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42762613405).

**Keep separate Draft and no gesui main/chocho merges.** No original raw individual records, photo pixels or observer identifiers were redistributed in this new aggregate report.