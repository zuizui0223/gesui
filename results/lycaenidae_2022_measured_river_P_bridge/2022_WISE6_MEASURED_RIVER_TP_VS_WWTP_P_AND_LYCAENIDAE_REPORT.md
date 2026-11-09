# Directly measured 2022 river total phosphorus does NOT support a landscape wastewater-P → Lycaenidae photo-recording pathway

**Completed November 9, 2026.** New, actual observed **river TP numeric measurements** linked to original publicly reported 2022 WWTP discharge-source geography and the previously frozen 2021–2023 Lycaenidae photographic-reporting source regions. **All results posthoc to previously known positive WWTP-photo association**; numerical TP extraction was independently source-gated before opening concentrations. The analysis is geography of *public photo-recording shares*, **not actual butterfly populations, adult mud-puddling or proven causal wastewater effects**.

## 1. New official source and genuine biological-source-aligned data support

**Official new source:** EEA Waterbase Water Quality ICM **2026** individual WISE-6 `DisaggregatedData` original public **2,024,284,753-byte** ZIP, full original ZIP SHA256 `2a2ef2a43a16ec068b7974a80ae5fca49e44697660a46ce3b91fb1813abe62d2`; original CSV `waterbase_t_wise6_disaggregateddata.csv`, **36,043,700,357 bytes** uncompressed. Read **102,569,057 individual source rows** in original order without expanding to disk. Source fields `countryCode`, original paired station ID/scheme, original year/sample date/sample identifier, `parameterWaterBodyCategory='RW'`, `observedPropertyDeterminandLabel='Total phosphorus'`, `observedPropertyDeterminandCode='CAS_7723-14-0'`, water matrix `procedureAnalysedMatrix='W'`, original `resultUom='mg{P}/L'`. Source station coordinates obtained from official `Waterbase_S_WISE_SpatialObject_DerivedData` original ZIP, joined on country+site ID+scheme and validated within original Natural Earth 1:50m national boundaries. No invented 2022 Germany/Austria annual aggregation source (2026 EEA original `AggregatedData` lacks same-year river TP for these countries).

| Real measured 2022 river TP source quantity | Germany | Austria |
| --- | ---: | ---: |
| Original sample rows satisfying 2022, river, water matrix and Total phosphorus | **3,690** | **1,281** |
| Original *date-source-eligible* sites (≥4 sample days in ≥3 months; original WISE6 station site key) | 242 | 100 |
| Original 42 three-year Lycaenidae photo regions: monitoring stations meeting original date metadata gate | **104** | **83** |
| Numeric TP site QC (below-LOQ censored without substitution, ≥4 quantified days and 3 months) | **103** | **81** |
| Country×1° source-positive photo regions with ≥2 numeric-quality WISE6 TP sites | **23** | **15** |

From the pre-fixed 42 original photo reporting cells DE27/AT15, **38 country×1° geographic cells (23 DE, 15 AT)** qualify in **21 independent country×2° source-blocks**. This **passes** the original prereading ≥20 (DE≥10, AT≥6) and ≥10 block source gate, unlike the earlier original Peters 2022 real TP source which supported only 9 DE country cells and 0 AT cells.

**TP numeric QA:** exactly **187** source-frozen qualifying sites were opened, **184** met numeric QC. There were **2,471 original 2022 TP source sample rows** at those sites, **53 below-LOQ records excluded**, no half-LOQ or zero substitution. Three preselected sites failed source quality (one 100% censored AT, one 54.5% censored AT, one 75% censored DE). Site-year measured TP = median of **original independent sample-date medians** in **mg P/L**, no weighting by number of repeated assays per day; country×1° cell TP = median of qualifying independent site-year medians. Among **38** eligible cells, measured TP cell medians range **0.0049–0.200 mg P/L**, overall unweighted median **0.0615 mg P/L**. These are actual river water concentrations, **not** reported WWTP annual P tonnes, and **not** nutrient levels measured where butterflies drink.

Source and code:
- [Original WISE6 full-ZIP source-only 2022 DE/AT river TP 187-station/daily-repeat support audit](https://github.com/zuizui0223/gesui/actions/runs/42572788975).
- [Original actual numeric TP source, 53 belowLOQ record exclusion, site-year medians and cell aggregation SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42574036643).
- [Original WISE6 numeric TP and LOQ site-year summaries, SHA-verified](2026_WISE6_2022_MEASURED_TP_SITE_YEAR_MEDIANS_AND_42_LYCAENIDAE_CELLS_QC.json).
- [Numerical outcome opening and regional testing rules committed BEFORE reading WISE6 numeric TP](../../docs/LYCAENIDAE_EEA2026_NUMERIC_TP_OPENING_AND_CROSS_LAYER_GEO_CV_LOCK.md).

## 2. Pre-numeric-locked 2022 river TP *environmental* predictions: WWTP P mass does not transfer

The first relationship is **actual in-river measured TP** at source-qualified country×1° photographic-landscape cells, with response `log1p(median 2022 observed river TP, mg P/L)`. Official EEA UWWTD original 2022 WWTP source variables are *geographic totals of positively measured-reported annual P/N effluent tonnes and number of uniquely verified freshwater outfalls in that broad cell*, NOT confirmed station-specific upstream discharge. Baseline is country, lat/lon, 2021 ESA WorldCover original built/crop/terrestrial area, all train-only scaled. Add each original WWTP predictor alone to the same baseline. Leave **one original country×2° geographic block out** at a time (21 blocks) with strictly train-only source scaling.

| Added original WWTP variable | Heldout baseline MSE | Heldout source-model MSE | Incremental MSE gain |
| --- | ---: | ---: | ---: |
| **2022 reported WWTP P tonnes/year** | 0.00099347 | 0.00109505 | **−10.22%** |
| 2022 reported WWTP N tonnes/year | 0.00099347 | 0.00110111 | **−10.83%** |
| Number of verified source outfalls | 0.00099347 | 0.00095778 | **+3.59%** |

An added predictor worsens out-of-region prediction for the P and N total reported mass models. Descriptive 999-sample country×2° prediction-error-block bootstrap on P gain **[−17.52%, −0.68%]**, and on outfall count **[−10.10%, +16.12%]**; the bootstrap intervals are about predictive-error heterogeneity, not randomized pollution effects.

Country-specific P-source prediction improvement: DE **−11.15%**, AT **−8.06%**. Country-specific 1° cell *unadjusted rank* correlations reported P-mass vs observed TP: DE **Spearman ρ −0.282** (23 cells), AT **ρ +0.604** (15 cells) — heterogeneous and uninformative for a universal physical dose model.

**Supported narrow conclusion:** A broad country×1° WWTP source reported P *mass total* does not provide transferable information about **actual in-river TP concentration** beyond the fixed geographic/WorldCover baseline in these source-eligible monitored areas. This is **not evidence that an actual hydrologically upstream WWTP load cannot change in-river TP**; local river connectivity/dilution/discharge must still be tested.

## 3. The actual river TP does NOT predict Lycaenidae photographic composition

Second relationship: the SAME 38 original 2022 measured-riverTP/source/photo 1° cells. Original 2022 response is public iNaturalist **photo reporting log odds** `log((Lycaenidae photos+.5)/(Pieridae+Papilionidae photos+.5))`, not species abundance, puddling or larval survival. Baseline country, lat/lon, WorldCover2021 actual built/crop/terrestrial area, 2021 independent Plantae photo reporting count, **2021 Lycaenidae photographic log odds AND 2021 all-three butterfly photo opportunity**. Every predictor standardized in training data only; leave 21 independent country×2° geographic source blocks out. No source-region selection by 2022 TP magnitude, photos or prediction gains.

| 2022 Lycaenidae photographic composition model | Baseline heldout MSE | Additional predictor heldout MSE | Incremental MSE gain |
| --- | ---: | ---: | ---: |
| Add **observed 2022 river TP** | 0.166399 | 0.183665 | **−10.38%** |
| Add original **2022 EEA reported WWTP P-mass geography** | 0.166399 | 0.153763 | **+7.59%** |
| With EEA P-source already in baseline: add actual 2022 observed river TP | 0.153763 | 0.175356 | **−14.04%** |

The original WWTP P-mass geographical predictor still has modest **photo reporting** predictive information on this 38-cell subset (+7.59%), but its descriptive 2° block error-bootstrap gain interval is wide **[−32.77%, +32.79%]**, and it has no supported measured river TP mechanism. Adding actual river TP makes photographic out-of-region prediction **worse**, both as a standalone candidate and after WWTP P mass is included. Descriptive error-bootstrap for the TP-only added gain **[−25.04%, −1.42%]** and for TP added beyond P-source **[−32.38%, −1.08%]**.

Unadjusted photo share log odds vs river TP differs by country: DE **Spearman ρ −0.223**, AT **ρ +0.418** (not a consistent positive mechanism). The **2023** Lycaenidae photographic report response with original 2021 history and 2022 TP likewise does not gain from measured river TP (**−5.71%** source incremental MSE), whereas the original 2022 WWTP P mass adds **+6.91%** photographic geographical prediction without 2022 butterfly-photo history; this is a *posthoc* temporal source-correlation comparison and **not 2023 measured nutrient dose**.

[Complete official observed river TP / WWTP source / Lycaenidae public photographic source region-heldout Action SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42574919214) · [Source-record predictions and country-wise audit JSON](2026_WISE6_ACTUAL_TP_2022_SOURCE_P_AND_LYCAENIDAE_PHOTO_REGION_BLOCK_CV.json).

## 4. What this shows, and what remains untested

**Empirically established advances:**
- A REAL independently reported **2022 environmental water chemistry** sampling set of 184 quantified TP monitoring site-years, sufficient to represent **38/42** original mixed Germany/Austria Lycaenidae photographic 1° source landscapes. The older Peters source could not do that.
- Official reported WWTP 2022 P-mass geography is not a substitute for actual in-river measured TP on these non-connected 1° cell regions.
- Actual measured river TP does not add positive out-of-region prediction of original source- and prior-year photo-controlled Lycaenidae photographic composition at this grain.
- Original previously reported WWTP-P × Lycaenidae **photographic recording geography** may largely reflect persistent landscape/source reporting, plant hosts/settlement and observation opportunity rather than identifiable biochemical eutrophication.

**Not tested and NOT refuted:**
- Physical effects of a specifically **upstream WWTP** on measured downstream TP after streamflow/dilution/connectivity matching; in this cross-layer model the WWTP source is just a total within the same broad 1° tile.
- Actual adult puddling, substrate protein/organic N or sodium, larval food quality, ant mutualism, survival, actual Lycaenidae density or population dynamics. TP is a river water concentration, not a field assay where butterflies feed.
- A true WWTP-induced P biochemistry mechanism: not identifiable from observational source-positive area/photo data and coarse noncausal regional correlations.
- Official individual WISE6 `resultObservationStatus` and expert QC metadata flags beyond LOQ/date/valid-unit filters have **not** been audited in this first analysis and must be checked before submission. EEA WISE Statistics applies additional automated and expert QC exclusions to some records; this source-derived model is **preliminary with respect to additional formal water-quality QC flags**. No flagged-sample issue is presumed present/absent without examining original flags.

**Future defensible independent physical test:** Instead of retrying macroecological photo models until a positive signal appears, spatially match each original WISE6 2022 TP sample monitoring station to the actual directed HydroRIVERS upstream EEA 2022 WWTP outfall network (including station-reach snap tolerance, uniqueness, hydrologic source-distance ≤50km and flow) and independently compare connected vs unconnected **real TP** station medians, with country/catchment and land-use controls. This would test a source→receptor link. Only after that would a properly effort-standardized independent butterfly behavioral response have a defensible causal mechanism candidate.

## 5. Existing results still separate

- Original `chocho` **22/25 (88%) Lycaenidae** increased *potential introduced-host larval resource geographic opportunity* — neither directly proven wastewater effect nor adult puddling.
- Original `gesui` 2018 and 2022 macroinvertebrate EQR incremental WWTP source predictive information — separate river biological response, not the present Lycaenidae outcome.
- Initial 2022 local WWTP source-aligned Lycaenidae photo effort fails at 1km (1/24 paired sources) and 5km (4/17 pairs); standard Norway BMS independent grid outlet coverage likewise too small for a source test. All failures are retained.
- Original 2022 DE/AT source-area Lycaenidae photo association was exploratory and posthoc region prediction +17.82% MSE improvement from reported WWTP P mass without independent river-chemistry/previous butterfly history; attenuation and reversal under historical photo composition and measured TP must be reported. Do not claim it is a wastewater-driven butterfly population increase.

### Original source references
- [EEA original Waterbase - Water Quality ICM 2026](https://www.eea.europa.eu/en/datahub/datahubitem-view/fbf3717c-cd7b-4785-933a-d0cf510542e1)
- [Eionet WISE6 individual sample/LOQ site key definitions](https://dd.eionet.europa.eu/datasets/latest/WISE-SoE_WaterQualityICM/tables/DisaggregatedData)
- [EEA WISE Statistics source expert QC description](https://sdi.eea.europa.eu/catalogue/srv/api/records/c6e353a2-c004-41eb-b0c1-5a2a1da5ef3b)
- [Original Peters 2022 WQ source geographic insufficiency, preceding PR #9](https://github.com/zuizui0223/gesui/pull/9)
- [Historical 2021–2023 Lycaenidae photographic falsification, PR #27](https://github.com/zuizui0223/gesui/pull/27)

**No numeric/mechanistic causal effects reported; all source original code and successful/failed gates are committed in this Draft branch, not `main`.**
