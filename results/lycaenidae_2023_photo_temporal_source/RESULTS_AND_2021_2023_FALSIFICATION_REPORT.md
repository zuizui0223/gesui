# 2021–2023 butterfly photographic-reporting composition and original EEA 2022 P-source geography: a three-year falsification

**2026-11-08; fully executed and original public API source hashes committed. ALL analysis here is explicitly POSTHOC relative to seeing positive 2022 wastewater-source geographical/photo reporting results.** No individual butterfly census, mud-puddling photo/behavior inspection, observer-level tracking, or direct nutrient delivery/water/substrate chemistry has been measured.

## 1. Original source-cohort design retained across three years

- Source original 2022 EEA Waterbase UWWTD Germany/Austria: **1,156** official WWTP plants with original both positive *measured-reported annual N/P discharge mass*, deduplicated to one direct FW outfall, 630 German, 526 Austrian; country×1° original source-positive cells selected without butterfly results, **53**.
- Independent 2021 Plantae public photo reporting and original ESA WorldCover2021 categorical built/cropland/terrestrial-area source quality screens yielded **43** 2022 butterfly-photo eligible country×1° source cells (DE28, AT15), **22** country×2° spatial validation blocks. The original analyses and positive heldout cross-region photo-reporting predictions are retained unchanged:
  - 2022 original 43-cell WorldCover built/crop/land-area-adjusted reported P-source incremental heldout MSE **+13.67%**; count-source effect **−1.52%**.
- Independent 2021 historical family-photo public API source QC: 42/43 cells have ≥15 original public family photographed records. (DE27/AT15). The 2021 source photo cohort was retrospectively assembled after positive 2022 results were observed; it is not a prospective negative-control clinical trial.
- New original **2023-04-01..2023-10-31** iNaturalist source (DE original `place_id=7207`, AT `place_id=8057`, Lycaenidae 47923, Pieridae 48508, Papilionidae 47223, georeferenced/photos/research grade/captive=false/position uncertainty≤500m) fetched at **exactly the same 43** preexisting 2022 original WorldCover-qualified/source-positive geographic cells, **129/129 family aggregate count APIs success**, no personal observer IDs or photo pixels. Before 2023 photo retrieval, the experiment's exact 43-cell source map, 2023 window, all filters, ≥20 family photographs per cell and geographic prediction specification were committed in `docs/LYCAENIDAE_2023_NEWYEAR_PUBLIC_PHOTO_POSTHOC_SOURCE_LOCK.md`.
- **42/43** 2023 photo cells meet the prespecified ≥20 family photographs cutoff (DE27, AT15) spanning **22 country×2° blocks**. All 42 also have historic 2021 photo coverage ≥15, giving one common 2021/2022/2023 42-cell public photo sample. Historical 2021 and 2023 source-qualified results are separate, and neither allows interpreting the EEA **2022** annual reported P-source as a known treatment introduced in 2022.

**Actual same 42 source-positive 1° cells — family photographed events (NOT actual butterfly individuals):**

| Observation year | Lycaenidae public photograph records | All Lycaenidae + Pieridae + Papilionidae photograph records |
| --- | ---: | ---: |
| 2021 | **2,237** | 4,356 |
| 2022 | **2,988** | 5,714 |
| 2023 | **2,926** | 5,780 |

Temporal differences in these opportunistic photo counts cannot be interpreted as population trends. The 2022–2023 Lycaenidae share log-odds correlation across the same country×1° source-positive cells is **Pearson r = +0.5877**, Spearman ρ = **+0.5856**. 2021–2022 corresponding Pearson r = +0.5143 on the same 42-cell source. These are cross-year spatial patterns in *photographic reporting*, not individual survival/fecundity.

Sources:
- [2023 original public three-family API source QC, 129 successful source requests](https://github.com/zuizui0223/gesui/actions/runs/42465810381).
- [Historical 2021 public photo source and paired prior-year heldout QC](https://github.com/zuizui0223/gesui/actions/runs/42465437203).
- [2023 new-year same-geography P-source heldout analysis SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42466025267).

## 2. Strict sample-consistent geographical 2° heldout comparison

Use **the same 42** 2022/2023 photo-supported 1° country cells and the exact same **22 country×2° geographic blocks**, each held out entirely in turn. All OLS transformations fitted on training-only data. Base predictors independent 2021 Plantae photographic reporting count, ESA WorldCover2021 **actual categorical** built-up/crop/country-terrestrial land area fractions, geographic lat/lon and country. Source addition **reported 2022 EEA P effluent mass**, within-country z(log1p(original reported mass)) fit from training only.

| Source/response year and baseline | 2° block-heldout baseline MSE | Baseline + original **2022** reported P-source MSE | Incremental MSE improvement |
| --- | ---: | ---: | ---: |
| Predict **2022** photographed Lycaenidae proportion (without previous-year butterfly photos) | 0.158674 | 0.142221 | **+10.37%** |
| Predict **2023** photographed Lycaenidae proportion (without previous-year butterfly photos) | 0.185150 | 0.175539 | **+5.19%** |
| Predict **2023** photographed Lycaenidae proportion **with original 2022 Lycaenidae photo log-odds AND total 2022 three-family photo reporting** already in baseline | 0.151373 | 0.162590 | **−7.41% (worse)** |

By country:
- 2023 no family-photo history adjustment: DE **+4.21%** incremental source P, AT **+7.30%**.
- 2023 after last-year's butterfly photographic composition/effort are in baseline: DE **−7.91%**, AT **−5.81%**. **Neither country improves** by adding the original 2022 P-source mass once 2022 photo geography is known.
- Original independent historical 2021 retrospective comparison: 2022 EEA P-source **worsens** region-out-of-sample prediction of the 2021 butterfly photo share by **−8.43%**, showing it is not a universal static predictor in every year.
- Original same-42 2022 photo response **with 2021 butterfly photo composition and target-group photo effort already in baseline**: P-source gain shrinks to **+3.70%** from +10.37%. This reinforces that previous-year butterfly-photo composition absorbs much of the initially positive P-source geographical information.

**Central scientific interpretation:** The 2022 source-reported P mass has reproducible, exploratory geographically transferable *association* with 2022 Lycaenidae **photographic reporting composition** (before last-year butterfly shares are included), but its added geographical predictive information weakens when the former year's photographic butterfly assemblage is modeled. In 2023, the same 2022 wastewater source adds **negative** prediction gain after 2022 butterfly photo geography is provided. This is more consistent with **persistent butterfly-community habitat and/or observation/reporter geography** than a uniquely 2022 biochemical wastewater-P mechanism. It cannot disprove an actual sewage ecological effect because source plants could have operated in 2021/2023, and unmeasured true nutrient delivery or butterfly dispersal is not observed.

## 3. Evidence hierarchy, robustness and limits

**Existing positive data-derived findings (preserved without overclaim):**
- `chocho` **22/25** Lycaenidae have expanded *potential larval-host resource geographic opportunity* under human-introduced plant host ranges. That is **not** measured butterfly increase, ant interaction, or source wastewater impact.
- `gesui` original 2018/2022 source-aligned benthic macroinvertebrate EQR additional heldout predictive information from reported upstream WWTP N/P remains a **separate aquatic organism / different grain** result with explicit detection/land use confounder limits.
- 2022 DE/AT iNaturalist photo composition source-P geographical association is a real **exploratory reporting pattern**, not evidence for wastewater causal biochemical impact.

**Competing and inconsistent falsifiers must all be reported:**
- Local 2022 1km/5km public butterfly photographic source-aligned effort support FAILED (1/24 and 4/17 matched WWTP source/no-registered-source river pairs).
- Standardized original 2022 Norway butterfly field grid source coverage near direct public both-positive measured N/P WWTP FW outlets FAILED (0≤2km / 2≤5km vs 63≥10km independent 70 Lucas grids); publicly coarsened Flanders 2022 locations cannot be matched at ≤500m.
- 2022 within-country reported P-source mass permutations yield stronger-than-random predictive reporting improvement, but **within 2° region** source shuffles yield similar gains in many draws and the 2022 landcover P-gain block-bootstrap 95% range crosses zero (approximately **−14.2%..+33.0%**).
- 2021/2022/2023 same-geography family photo records show substantial year-to-year compositional persistence; prior butterfly photographic composition/target-group recording opportunity materially explains original source-P gains.

**Unresolved confounding:** Original 2022 WWTP report and facility density correlated with human settlement, urban parks/roads, crop agriculture, incomplete direct butterfly survey effort, photographers and habitat/host plant distribution. Historical *public photograph* counts are not standardized survey visits or butterfly population census. 2022 P and N are strongly collinear and the mass year is not assured to be the chemical amount accessible to a given 2023 butterfly. The current 1° observations do not establish same-river downstream transport, or adult wetted-substrate/protein/Na feeding.

**Scientific claim boundary:** The macro-scale wastewater-reported P vs Lycaenidae photo reporting correlation is **at most landscape-level predictive observational geography**. The ecological mechanism "effluent supplies adult Lycaenidae nutrients" remains **untested**; do not frame 2023 negative historical-adjusted gain as experimental proof of no wastewater effect.

## Reproducibility

- Exact source-only original year and sampling-frame selection manifest SHA-pinned in both workflow results; source-only new-year panel original API receipts and taxon IDs included.
- Code files: `scripts/gesui_2023_lycaenidae_public_photo_source_original43.py`, `scripts/gesui_2023_lycaenidae_sourceP_geographic_year_check.py`; no new photographs/identifiable observers downloaded.
- Original WorldCover 2021 official categorical COG source and LandCover model: [Draft PR #25](https://github.com/zuizui0223/gesui/pull/25).
- Original historical 2021 diagnostic and exact 2022 response: [Draft PR #26](https://github.com/zuizui0223/gesui/pull/26).

**No main merge.** This report adds one new 2023 year beyond the original 2022 and historical 2021 and remains an explicit posthoc falsifier.
