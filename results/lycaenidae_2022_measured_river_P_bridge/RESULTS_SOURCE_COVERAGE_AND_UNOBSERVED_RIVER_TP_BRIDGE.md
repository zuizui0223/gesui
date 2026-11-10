# 2022 real river-total-P × DE/AT Lycaenidae photo geography: measured-chemistry source and sampling-support audit

**Source-audited 2026-11-09. Status: no TP concentration or new wastewater-biota mechanism fitted.** This is an evidence-strengthening *data-availability* test, not a new proof of environmental phosphorus mediation. All original `gesui` and `chocho` ecological positives and negative-control results are unchanged.

## 1. Exact three-year original photo/source cells and original observed 2022 Peters TP sites

The previous DE/AT study froze originally EEA 2022 source-positive country×1° cells and independently ESA WorldCover2021 country land/urban/crop content, 2021 plant reporting effort, and original 2021/2022/2023 public photographed Lycaenidae/Pieridae/Papilionidae observation composition. The **strict balanced original sample is 42 country×1° cells, DE27/AT15**, with 22 original country×2° blocks.

The independently published Peters et al. (2026) original measured river nutrient study has an explicit 2022 TP site-table previously independently audited in `gesui` source-to-river PR #9. The original verifiable Git blob is **`a75a4cf140de64fd6038079a3fac6756a461f3f8`** at exact commit `6f037f65e9b948d6f8e324ad4eec9d8c03424b35`. **1,291** source-qualified 2022 measured TP river monitoring sites (at least four dated real 2022 TP samples in ≥3 distinct 2022 months). This new bridge audit read **ONLY** publicly reported site ID, TP monitoring 2022 sampling opportunity, station coordinates, source-registered upstream WWTP count, HydroRIVERS reach identifier, and original fixed photo geography; the numerical original 2022 TP values were not read, fitted or used for cell selection.

**Country assignment:** original Natural Earth 1:50m country polygon original SHA256 `3e458fc036ad0a66411f2c1e6cac49c5d7bfb81cb1123bc513b22511a2b7fdeb` (same as 2021 WorldCover land-cover audit), not a cross-border bounding-box approximation.

| Actual 2022 Peters monitored TP site support | Germany (DE) | Austria (AT) |
| --- | ---: | ---: |
| National TP sites after exact country polygon match | **68** | **3** |
| TP sites in the original strict 42 country×1° butterfly-photo cells | **35** | **3** |
| Original butterfly-photo cells with ≥2 2022 TP site measurements | **9** | **0** |

The **9** cells with ≥2 original independent water TP station measurements cover only **5** distinct original country×2° source blocks, **8** cells with at least one TP station with original hydrologically upstream registered qualifying WWTP source. The original prospective three-layer gate requires ≥20/42 photo cells with ≥2 distinct TP stations, ≥10 DE, ≥6 AT, ≥10 country×2° blocks and ≥10 tiles with original upstream source. **FAILED.** Do **not** fit a 42-cell photo-WWTP-*measured river TP* model on this geographically selective subset; no proper Austrian control.

[Source-only Peters TP geographic and official public country's polygon QC workflow SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42554476604) · [Full source-only station coverage JSON](2022_PETERS_ORIGINAL_TP_MONITORING_STATION_GEOGRAPHIC_SOURCE_SUPPORT_QC.json).

## 2. Newer official EEA 2026 Waterbase Water Quality ICM — original source availability verified

The EEA published *Waterbase - Water Quality ICM, 2026* on **July 7, 2026**, with 1900–2025 overall temporal coverage. The official data provider separates country-reported individual water sampling (`DisaggregatedData`) from country-reported yearly aggregates (`AggregatedData`) and provides source site coordinates in separate `SpatialObject_DerivedData` files. A particular country's absence from one reporting mode is **NOT** absence of water monitoring in all reporting modes.

### Original official 2026 ZIPs, SHA and record/column metadata

| Official original source ZIP | Compressed size | Uncompressed original CSV | Original SHA256 | Scope |
| --- | ---: | ---: | --- | --- |
| `Waterbase_T_WISE6_AggregatedData` | **199,878,007 bytes** | 2,702,939,173 bytes | `8492c92702524ce23ac488cba62631f277b1e9bd4e5ea8ece3bea37569307b7c` | Annual-site water quality, 32 CSV columns, `;` delimited |
| `Waterbase_S_WISE_SpatialObject_DerivedData` | **3,394,632 bytes** | 31,954,165 bytes | `afd64e1d5582492c09453f5a05aae0b1a9d9a38c8814778afdb5c07e43a1c974` | Public original country monitoring-site identifier schemes, lat/lon |
| `Waterbase_T_WISE6_DisaggregatedData` | **2,024,284,753 bytes** | **36,043,700,357 bytes** | FULL archive SHA **not computed; original bounded HTTP ranges validated** | Individual actual sample dates and observations; original single CSV; 29 original columns |

**Annual CSV source-only complete streaming audit:** **6,663,213** original annual rows checked without loading or using any numeric TP outcome. For **2022**, original country `AT` annual aggregate records **29,111**, all with source `parameterWaterBodyCategory='GW'` (groundwater), **no source `RW` river rows**, and original `DE` aggregate records in 2022 **0**. Austria's 2022 annual records here are primarily groundwater pesticide/metal measurements. **This does NOT imply Germany/Austria had no 2022 river TP measurements:** the **separate** EEA individual original `DisaggregatedData` is enormous and must be queried separately.

**Discodata public SQL API:** original documented station-metadata `SELECT TOP 1` attempts on `WISE_SOE/latest` and `v2r1` returned **HTTP 403** from this runtime, even though the original bulk ZIP public downloads returned HTTP 200. No national 2022 individual TP row counts were fabricated from a rejected query. Official EEA WISE freshwater page cautions that WISE-6 Discodata tables may lag the latest data.

**Bounded original individual ZIP-range audit (no original chemical records):** official HEAD reports compressed size 2,024,284,753 bytes; explicit HTTP `Range bytes=start-end` returns 206 Partial Content. Central-directory original ZIP metadata proves **a single original `waterbase_t_wise6_disaggregateddata.csv` of 36,043,700,357 uncompressed bytes**. The prefix 2 MiB HTTP range was decompressed *only for its exact semicolon-separated 29-column header*; required source fields confirmed:

`countryCode`, `monitoringSiteIdentifier`, `monitoringSiteIdentifierScheme`, `parameterWaterBodyCategory`, `observedPropertyDeterminandCode`, `observedPropertyDeterminandLabel`, `procedureAnalysedMatrix`, `resultUom`, `phenomenonTimeSamplingDate`, `phenomenonTimeReferenceYear`, `sampleIdentifier`, `resultObservedValue`, `resultQualityObservedValueBelowLOQ`, `procedureLOQValue` and other original metadata. **No numeric resultObservedValue, LOQ value, species record or private data read** in these checks. The full 36GB individual chemistry table was not decompressed/downloaded into the repository or analyzed.

Source access evidence:
- [Public SQL and ZIP availability source status SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42554162843).
- [Official annual/spatial ZIP original SHA/schema SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42554593278).
- [All 6.66M annual source row metadata / country TP support SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42555383185).
- [Official 2GB individual ZIP explicit byte-range central directory SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42555765234).
- [Official individual WISE6 source **29 exact fields** from first 2 MiB only SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42556062300).

## 3. Next defensible original data retrieval: 2022 DE/AT river total-P filtered WISE6 records, NOT whole 36GB

Original EEA WISE6 full observation data are **publicly available**, but the original individual dataset is 36GB uncompressed and the SQL endpoint returned HTTP 403 here. **This is source access/processing limitation, not ecological absence.** A narrow original national/year/matrix/parameter subsetting request from EEA via its official WISE SoE helpdesk or functional Discodata account should request ONLY:

- `countryCode` in `DE, AT`; `parameterWaterBodyCategory='RW'`; `observedPropertyDeterminandLabel='Total phosphorus'` with original observedPropertyDeterminandCode crosschecked before accepted, original sample dates **2022-01-01..2022-12-31**; original `procedureAnalysedMatrix` must indicate water and original `resultUom` explicitly retained to harmonize µg/L and mg/L safely.
- Original `monitoringSiteIdentifier` and `monitoringSiteIdentifierScheme`; `phenomenonTimeSamplingDate`, `sampleIdentifier`, `resultObservedValue`, `resultQualityObservedValueBelowLOQ`, `procedureLOQValue`, `resultObservationStatus`, optional sample depth for quality.
- Matching **public** original `Waterbase_S_WISE_SpatialObject_DerivedData` source original ID/scheme/lat/lon coordinates and source country. Do not reuse same station ID across schemes or deduplicate dated samples into false independent sites.
- Prior to a model, freeze source-only stations/countries and check ≥4 actual qualified 2022 TP sample dates in ≥3 distinct months per site, annual sampling support, and the already fixed 42-cell country-mask thresholds. If unavailable or below support, STOP. Only then read numerical river-TP concentrations and normalize original units/LOQ.

**Unsent data request and exact official-field SQL example** are in `docs/LYCAENIDAE_EEA_WISE6_2022_TP_FILTERED_EXTRACT_REQUEST_UNSENT.md`. No external email has been sent or authorization presumed.

## 4. Ecological claim boundary

Original `gesui` 2022 EEA source mass was publicly **reported measured annual N/P discharged tonnes** at WWTP outlet, not measured TP concentrations in the river or animal-substrate sodium/protein. Original 2021–2023 Lycaenidae photographic response has spatially persistent reporting composition; the positive 2022 P-source photo prediction becomes **worse** for predicting 2023 after 2022 butterfly history is accounted for. Our new 2022 real TP data source audit therefore supports **neither a new nutrient mechanism nor its absence**.

A three-way source→measured river TP→Lycaenidae reporting association is currently **NOT IDENTIFIABLE with properly country-matched and repeated sites**. Its first obstacle is linking independently measured TP sites to the same original geography; the best next step is retrieval of the already-public EEA individual 2022 DE/AT RW TP subset, not unblinded selection of favorable water-quality sites.

### References
- [Official EEA Waterbase Water Quality ICM, 2026](https://www.eea.europa.eu/en/datahub/datahubitem-view/fbf3717c-cd7b-4785-933a-d0cf510542e1).
- [Original EEA WISE SoE data collection and download guidance](https://water.europa.eu/freshwater/resources/wise-soe-data-collection).
- [Official WISE6 DisaggregatedData reporting field dictionary](https://dd.eionet.europa.eu/datasets/latest/WISE-SoE_WaterQualityICM/tables/DisaggregatedData).
- [Peters et al. 2026 original global river nutrient dataset paper](https://www.nature.com/articles/s41597-026-07625-1).
- [Previously frozen original 2022 source-aligned Peters river chemistry study, gesui PR9](https://github.com/zuizui0223/gesui/pull/9).
- [2021–2023 Lycaenidae photograph reporting historical negative controls, gesui PR27](https://github.com/zuizui0223/gesui/pull/27).

**All result claims are source limitations; preserve previously established aquatic and larval-resource positives unchanged.**
