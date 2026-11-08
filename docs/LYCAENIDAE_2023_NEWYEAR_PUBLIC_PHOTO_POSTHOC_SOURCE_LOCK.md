# Independent 2023 photo-recording temporal source test (2026-11-08)

**Timing:** Developed after seeing 2022 WWTP-P × Lycaenidae positive 1° region photo-reporting results, 2021 ESA WorldCover land controls, and first historical **2021** butterfly photo diagnostic. All this test's rules are fixed **before opening any 2023 butterfly family counts**. This is a **posthoc new-year exploratory falsifier**, not a blinded confirmation and not a 2023 WWTP load study.

## Data source and identical cohort

- Original frozen 2022 EEA Waterbase public DE/AT UWWTD positive *reported measured* both N/P one-FW-outfall manifest (DE630/AT526) and 53 original country×1° source-positive areas.
- Preoriginal 2021 independent Plantae photo opportunity, 2021 ESA WorldCover categorical built/crop/country land area, and 2022 Lycaenidae/Pieridae/Papilionidae April–October photo-recording support already resulted in **43** 2022-supported original country×1° photo reporting cells after original land-cover QC (**DE28, AT15, 22 country×2° spatial groups**). These 43 are the only selection set here; no addition or replacement based on 2023 response.
- 2021 prior-year iNaturalist photographed three-family counts have been acquired for these exact 43 areas, no new 2021 API calls. Prior 2021 family photography ≥15 on 42/43 (DE27/AT15).
- New official iNaturalist public v1 **2023-04-01 through 2023-10-31**, fixed 43 cells, official country place IDs DE7207/AT8057, family-descendant taxon IDs Lycaenidae47923/Pieridae48508/Papilionidae47223, research-grade, georeferenced photo, public open/≤500m, noncaptive. Fetch only API `total_results`, the exact public query URL, and SHA256 response; **no photograph pixels, observer/person IDs, species records, or individual butterfly counts**.
- Query all 43 source/previously frozen land-QC 2022 cells regardless of 2023 photorecord counts. Minimum three-family **≥20** 2023 photo records per cell; require **≥30 cells** total **DE≥18 AT≥9** and **≥12** country×2° spatial blocks. Missing API query = NA, not absence. No retuning.
- A strict 2021–2022–2023 three-year descriptive balanced sample also requires ≥15 original 2021 photo records in each of these 2023-eligible cells. Report three-way support only; year trend is **photo reporting**, not population.

## Fixed 2023 geographic prediction after 2023 sample source gate

- Baseline (B) identical to 2022 posthoc WorldCover-adjusted heldout OLS: intercept, country, latitude, longitude, log(1+2021 Plantae photos), 2021 WorldCover built fraction, crop fraction, log(1+country terrestrial km²), all train-only standardized per left-out original country×2° spatial group.
- Primary predictive source covariate, separately: `log1p(2022 EEA reported P t/year)`, normalized *within country training folds only*; 2022 N and facility count separately as fixed sensitivity, never entered jointly due collinearity.
- Apply identical B and B+P to **2023** observed Lycaenidae log odds of public photographs versus Pieridae+Papilionidae. Compute group-heldout unweighted MSE, photo source relative improvement, country-specific heldout MSE on the 2023 photo-supported cells. Then, for same fixed 2023-eligible 1° cells, compute **2022** B vs B+P: compare two predictive associations without conflating different sample geometries. They are retrospective source-geography descriptions, not 2022+2023 source mass time series.
- *Historical community control model*: for 2023 response, augment both baseline and B+P with **2022 Lycaenidae photographic share logit** and **log(1+2022 three-family public photo record count)**, train-only scaling, then ask if 2022 P source contributes incremental heldout MSE predictive information beyond last year's photo-recording geography. As secondary model on three-year-balanced 2021-supported sample, include 2021 photographic share+effort only for comparison if target data supports. Do not posthoc select whichever has the better positive source gain.
- **Interpretation:** If historical butterfly photo share accounts for most P-source predictability, treat as evidence for persistent habitat/observer/reporting geography rather than acute wastewater nutrient action. If P-source remains predictive, it still does not prove sewage causation or organism benefits.
- **No source-registered 2023 annual P/N mass is claimed**, no environmental water substrate/organic N/protein or adult puddling evidence, no butterflies' individual counts measured.

## Frozen outcome-safety and traceability

Original `gesui` 1km/5km butterfly photographic source-matched sample support failures (1/24 and 4/17), Norway BMS and Flanders original location issues, `chocho` Lycaenidae 22/25 potential introduced-host resource geography, and `gesui` 2018/2022 river EQR predictive results remain separate. **Do not replace their source grain or predeclared status**. Source result files and all unsuccessful QA must be committed and reported.

This independent 2023 analysis is not a proper photo opportunity observation model, does not control direct land-use temporal change or observer identity, and cannot support biological wastewater effects.
