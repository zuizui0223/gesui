# 2020 third-year measured N/P → aquatic EQR: pre-outcome source-support decision

**Data status: no numerical 2020 EQR opened.** This is an independently executed *source-only feasibility gate*, not an EQR prediction or an ecological result.

## Recovered original 2019–2020 public EEA sources

Original EEA UWWTD reporting-cycle 2019–2020 spatial files, version v08/r00, recovered from EEA's public Nextcloud/WebDAV under the official public dataset's catalogue path. Browser pages and the old /Tabular URL were insufficient; original official GPKG file bytes were downloaded with checksums.

| Original public source | Original records | Original SHA-256 |
| --- | ---: | --- |
| Treatment plants (7,159,808 bytes) | 21,448 | `13060c88a8d1d06f81b61b0a6886279db6051645a2b8d1f7d578ae8ee4e59300` |
| Discharge points (7,073,792 bytes) | 23,965 | `1b8b20e0ba9d106f659b0fd9716013b680a16a0e7ab2a779dac2bfb53a33d9a1` |

Original treatment fields include `uwwNDischargeMeasured` and `uwwPDischargeMeasured`, not modeled treated-water flow; initial complete historical public product has 4,710 plants with positive measured-reported N, 2,969 with positive measured-reported P, and 2,789 with both. **This public product mixes original official reporter codes with suffix 2020, 2022 and other codes**. A simple v08 dataset-label filter would have introduced other reporting cycles. All subsequent source cohorts explicitly require `repCode` ending in `_2020`, consistent full `repCode|uwwCode` plant/outfall key, a **single** associated point, direct freshwater (`FW`) and valid Europe study coordinates, both reported-measured N/P masses strictly positive, and no source inference from absent values.

These conditions produce **582** public original source-reported 2020 FW outlets. **521/582** connect to the pinned third-party HydroRIVERS EU 500m limit, using directed `NEXT_DOWN` paths ≤50km; same-reach outlet/receptor order is unknown and excluded. Each matched source retains 2020 reporter code, not necessarily a verified laboratory measurement date for 2020.

## Independently frozen biological site metadata (no scores)

Original EEA Biology 2026 archive was read **only for 2020 site IDs, assessment-method names, element labels, coordinates, country and reporting status**; no 2020 resultNormalisedEQRValue, resultEQRValue, status/score outcome columns were accessed.

- Original geolocated 2020 source-only site-element candidates: 8,624 (benthic 3,533; phytobenthos 4,024; fish 1,067).
- Directed river network candidate site elements after 500m snap: **7,837**; unique site keys: **4,282**.
- Geographic 3° blocks and deterministic original outcome-blind fivefold assignments are frozen in source-checked CSVs.

## Predeclared 2020 source support result — STOP for independent third-year EQR model

| Element | River-linked sites | Registered 2020 reported N/P mass upstream ≤50km | Source-exposed countries | Previously declared generalization threshold |
| --- | ---: | ---: | ---: | --- |
| Benthic macroinvertebrates | **3,175** | **17** (AT 12, HU 4, SI 1) | **3** | **FAIL:** compatible exposed ≥50 and countries ≥5 |
| Phytobenthos | 3,645 | 18 (AT 11, HU 6, SI 1) | 3 | FAIL |
| Fish | 1,017 | 4 (AT 4) | 1 | FAIL |

For benthic EQR 17/17 exposed have at least one 2018 country×method matching category, but just **1/17** matches a 2022 method category. No numerical EQR was viewed to select source/year/country or tune ecological model. The original benthic cohort also fails the predeclared basic model source support `≥20 exposed, ≥3 source countries, ≥5 total countries, ≥10 blocks` because only 17 are source-positive.

**Decision:** **Do not open 2020 EQR outcomes or compute a third-year heldout performance effect.** The year is a documented source-coverage failure under source-only predeclared thresholds, neither a failed nor successful 2020 ecological prediction. In particular, `2018/2022` original positive region-heldout +0.98–1.19% and +1.10–1.40% benthic MSE reductions remain unchanged and are not falsified by 2020 source incompleteness.

## Caveats and scientific next steps

The EEA GIS public product filters active and spatially correct sources, reports different years/countries unevenly, and omits ambiguous outlets and source-absent sites. `No registered upstream 2020 treatment plant measured in both nutrients` does **not** establish absence of WWTP pollution. Reported annual mass is not delivered nutrient load at a downstream biological station. The downloaded HydroRIVERS EU GeoParquet is a SHA-pinned third-party mirror, not an independently proven official HydroRIVERS topology.

Scientific goal remains **which wastewater-associated source features carry reproducible, geographically transportable information on river biological quality, and under which reporting/hydrological/ecological conditions**. A fourth independent source year or a broader independently measured source inventory must pass this same no-EQR source/method support gate **before** any new biological performance models or adaptation. Do not change the thresholds after learning the 2020 support failure to search for a positive result.

## Reproduction and provenance

- [Source-only official public 2019–2020 GPKG audit](https://github.com/zuizui0223/gesui/actions/runs/42445491034)
- [2020 biological monitoring-site and country/method freeze](https://github.com/zuizui0223/gesui/actions/runs/42445263038)
- [2020 hydrologically directed source-only river linkage and feasibility result](https://github.com/zuizui0223/gesui/actions/runs/42445631605)
- [Original 2018 temporal replication Draft PR #17](https://github.com/zuizui0223/gesui/pull/17)
- [New third-year source-only Draft PR #19](https://github.com/zuizui0223/gesui/pull/19)
- [Source-only predeclared eligibility Issue #18](https://github.com/zuizui0223/gesui/issues/18)

This stage remains **Draft**, with no automatic merge to main or the independent `chocho` manuscript.
