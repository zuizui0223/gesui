# Independent exploratory cross-region Lycaenidae photographic opportunity × reported wastewater N/P (DE/AT)

2026-11-08. **Prospectively fixed before opening any 2022 butterfly photo observations in this independent 1° × 1° geographic cohort**. This is a deliberately new macroecological question after 2022 local river 1km / 5km paired photo-support failures and standard BMS sampling-source failures. It is not a sensitivity that rescues those studies, and it does NOT estimate WWTP effects on butterfly abundance, adult puddling, larval survival, ant interactions or direct organic-N/sodium availability.

## Ecology and estimand

For each broad 1° WGS84 cell that contains a verified **2022 EEA publicly reported, both positive measured N/P WWTP annual effluent mass** and a unique direct freshwater discharge point in country DE or AT: is the **share of public iNaturalist butterfly-family photographic observation records in Lycaenidae** versus Lycaenidae+Pieridae+Papilionidae related to source WWTP count or source reported P/N mass, after accounting for independent **2021 plant photo-record reporting effort**, geographic lat/lon and country?

**Observational reporting comparison only**. This does NOT show increased or decreased butterflies, establish actual biochemical nutrients or their transport to photos, confirm puddling or source-water exposure, or demonstrate causality. P and N source masses are collinear and must be tested separately. Source count is mainly a proxy of human habitation, historical monitoring and plant-source reporting coverage; its association with 2022 butterfly photographs may reflect effort, habitat and land use.

## A. Original source geography, frozen before any 2021 Plantae or 2022 butterfly outcome query

- Original [EEA 2022 DE+AT FW both measured N/P public source manifest](https://github.com/zuizui0223/gesui/blob/analysis/lycaenidae-de-at-public-event-discovery-20261108/results/lycaenidae_de_at_2022_source_only/EEA2022_DE_AT_PUBLIC_BOTH_NP_MEASURED_ONE_FW_OUTFALLS.csv), SHA-256 `926ab310904be048d5d36bbfebd5eabfe0b318279800142f2e2a86bee924f4ce`; original summary 630 DE and 526 AT.
- Previously frozen [source 1° WGS84 cells](https://github.com/zuizui0223/gesui/blob/analysis/lycaenidae-de-at-public-event-discovery-20261108/results/lycaenidae_de_at_2022_source_only/EEA2022_DE_AT_ORIGINAL_FW_OUTFALLS_IN_ONE_DEGREE_SOURCE_TILES.csv) contains **53 country×grid source-positive cells**: 32 DE and 21 AT. **Take all 53**, no replacement/selection from butterfly observations. Source intensity = `log1p(n_qualified_outfalls)`; additional separate sensitivity `log1p(source_P_tonnes_per_year)` and `log1p(source_N_tonnes_per_year)`, retaining original units and collinearity warning.
- Cell search uses exact 1° latitude/longitude **bounding rectangles** from source tiles, AND an official **iNaturalist Germany/Austria country place filter resolved independently by public `/places/autocomplete`**; no other-country records from transborder cells are allowed. Resolve unique official administrative country, else stop. A country×grid cell is not a WWTP catchment or river downstream exposure footprint. Original 1° rectangles vary in km² area and land fraction; control geographic latitude/longitude but those do not fully correct land-area bias.
- Do not treat zero source means zero sewage or zero human habitation. This estimand is explicitly **variation conditional on ≥1 qualified reported source**, not a true treated-vs-untreated comparison.
- Record source-only cell-level exact country, bbox, EEA original source count, original P/N tonnes, SHA, 2° geographical spatial block (`floor(lat/2), floor(lon/2)`) and cell country indicator; no iNaturalist observation data loaded at stage A.

## B. PRE-outcome independent 2021 plant-only photographic effort gate

Official iNaturalist v1 API: original `Plantae` kingdom (exact verified taxon ID, not a name-only family query), `year=2021`, `photos=true`, `geo=true`, `geoprivacy=open`, `taxon_geoprivacy=open`, `quality_grade=research`, `captive=false`, `acc_below=500`, country `place_id`, and the fixed 1° rectangle. Query **all 53** original source cells, no 2022 butterfly taxon/photo queries until ALL prior-year effort results are committed with original API request URLs/checksums.

- Minimum **40** 2021 public qualified plant photo observations per country×cell. Missing/failed API requests are NA and do NOT imply zero plants. Keep complete 53-row result and full source API receipts; no outcome-based cell replacement.
- To proceed require **≥25** eligible cells in total with **≥12 Germany and ≥8 Austria**; at least **3 distinct source-count bands** (1–4, 5–19, ≥20 qualifying outfalls) across the union, at least 2 bands per country, and ≥10 original distinct 2° spatial blocks.
- Any failure ends at stage B **before 2022 butterfly outcomes**.

## C. Frozen 2022 butterfly photographic reporting share (new exploratory response)

After B succeeds and both source and 2021 plant-effort eligibility are committed: query same fixed grid, country `place_id`, season `2022-04-01..2022-10-31` for Lycaenidae (`47923`), Pieridae (`48508`), Papilionidae (`47223`), including descendant species; **matching** photos/geo/geoprivacy/research-grade/noncaptive/accuracy filters. Original 2021 plant effort uses whole year, deliberately separate; 2022 family photo season is fixed **April–October**, not optimized by 2022 response. As a sensitivity only, add 2022 whole-year counts IF separately predeclared, but no discretionary season tuning. Do not count original taxon rank-only family observations as the whole family; confirmed taxon IDs include descendants.

- For each cell, counts per butterfly family and all-three denominator; photo fraction `Lycaenidae/(Lycaenidae+Pieridae+Papilionidae)`.
- Descriptive/adjusted ecological model support gate: **≥20 total 2022 three-family photos per cell, ≥25 response-supported cells total, ≥10 in each country and ≥10 distinct 2° spatial blocks**. If not met STOP before any regression; preserve raw source coverage.
- Primary descriptive rank association of source-count `log1p` and Lycaenidae photo-share within DE and AT separately and pooled, as well as country-equal effect. Fixed primary linear model: `logit((Lycaenidae_count+0.5)/(total+1)) = b0 + b1 log1p(n_EEA_both_NP_FW_outfalls) + b2 log1p(2021_Plant_photo_count) + b3 latitude + b4 longitude + country fixed effect`. **Unweighted at grid/country** to reduce citizen-science urban sampling hotspots. Rescaled standardized predictors for interpretable slope. This is exploratory association only; land area, land-use/urban fraction, climate, roads/observer and host plants may still confound.
- Mandatory negative controls/QA: compare coefficients after excluding highest source-count tile per country (source-influence), Spearman within each country, report all 3 separate source axes (count, P mass, N mass), and 2° block resample 999 draws with seed 20261108. No threshold changes after outcome. If source density and photo-share are associated, do not claim sewage nutrients are the mechanism or compare it to original chocho "resource expansion" as a proven mediation.
- Report the original 1km 24-pair and 5km 17-pair failed photo-support gates and source-only Norwegian, Swedish BMS gates alongside any exploratory macro association. Contradictions at different scales are scientifically plausible.

## D. Reproducibility and interpretation firewalls

- Do NOT read GBIF butterfly photos, user/observer personal data, photo pixels or individual 2022 butterfly IDs. Only API aggregated `total_results` with exact GET URL and body digest; do not interpret photos as abundance, specimens, abundance zero, or adult puddling.
- Complete source-only 2022 EEA CSV and 2021 non-butterfly effort check **must** be GitHub-committed before any 2022 family query. No 2022 EQR biology site outcomes are involved.
- A model fit can only support **geographical co-distribution of wastewater infrastructure with a Lycaenidae photographic reporting proportion**, conditional on the specified controls. There is **no data-driven attribution to effluent-fed adult behavior or larval fitness**.
- No `chocho` manuscript modifications and no merge of this new Draft without a new pre-outcome review and independent land-cover/urban validation.
