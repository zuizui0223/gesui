# German–Austrian 2022 Lycaenidae photographic reporting: ESA WorldCover 2021 independent urban/crop/land-area falsification

**8 November 2026. STATUS: completed POSTHOC spatial predictive sensitivity, not causal ecological evidence.** The prior within-sample positive EEA source association and prior posthoc 44-cell MSE +17.82% P predictive gain were already viewed before designing this WorldCover covariate comparison. The independent environmental source itself was fetched and Git-committed **before this new WorldCover-adjusted model ran**, but this is not an independently preregistered butterfly outcome test.

## 1. Original ESA WorldCover categorical environmental source, separate from butterfly outcomes

- Original published [ESA WorldCover 2021 v200 10m categorical map](https://esa-worldcover.org/en/data-access), DOI [10.5281/zenodo.7254221](https://doi.org/10.5281/zenodo.7254221), read from original public AWS anonymous **COG** via HTTP-range windows, **not Terrascope WMS RGB visualization**. Original classes 40=cropland, 50=built-up, 10=tree, 30=grassland, 80=water, 0=NoData.
- Original [Natural Earth 1:50m sovereign/admin0 country polygons](https://github.com/nvkelso/natural-earth-vector/blob/master/geojson/ne_50m_admin_0_countries.geojson) country geometry; original GeoJSON SHA256 `3e458fc036ad0a66411f2c1e6cac49c5d7bfb81cb1123bc513b22511a2b7fdeb`. Only original country polygon land and valid nonwater original WorldCover classes contribute; cross-border areas are not counted as German/Austrian land.
- Full original EEA source-verified 2022 **53 country×1°** WWTP-source-positive cell sample (DE32/AT21) is preserved, with original 1,156 both measured-positive reported N/P uniquely one-FW-outfall source manifest SHA256 `926ab310904be048d5d36bbfebd5eabfe0b318279800142f2e2a86bee924f4ce`.
- **All 53 original source-region cells** successfully decoded from **11 unique original 3° COG raster tiles**, original CRS EPSG:4326. Original categorical nearest-neighbor class pixel window **512×512** per 1° cell; separate **256×256** original categorical resampled source sensitivity. This is a ~sub-kilometre systematic **sample of a 10m product**, **not** a 10m full-resolution pixel census, and sampling error is nonzero.
- Independent pre-model GIS-quality criteria: terrestrial country fraction ≥5% of original 1° cell, ≥512 country land classification sample points and ≤5% original NoData in country polygon; **48/53** pass (DE31/AT17). All 53 original source cells pass prespecified 512-versus-256 urban/crop classification stability check ≤10 percentage points. Five source regions fail country-land area, and are reported; missing source ≠ zero built-up.
- WorldCover source environmental dataset/QA was committed before reading original pre-existing 2022 butterfly photo-response results in this new model: [COG feasibility run](https://github.com/zuizui0223/gesui/actions/runs/42464131882), [all 53 independent original categorical landcover sources SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42464211489), [complete 53-cell source/QC JSON](ESA_WORLD_COVER_2021_DE_AT_ALL53_COUNTRY_MASKED_SOURCE_QC.json).

**Source attribution:** © ESA WorldCover project 2021 / Contains modified Copernicus Sentinel data (2021) processed by ESA WorldCover consortium. Natural Earth vector geographic boundaries © public domain.

## 2. Preserve original photo source population and GIS-supported subset

Original independent EEA2022×2021 Plantae reporting×2022 public iNaturalist Lycaenidae/Pieridae/Papilionidae study: **53** source-positive 1° cells, 52 pass independent prior-year 2021 plant photo reporting threshold, **44** pass 2022 butterfly three-family ≥20 photo-record gate (Germany29, Austria15) and 22 country×2° geographic regions. These are *photographic observations*, not butterfly census individuals, abundance or puddling.

The WorldCover source land area rule excludes only **one** of the original 44 2022 butterfly photo-eligible cells: Germany 47°–48°N, 13°–14°E, only about **0.94% of its source 1° cell** is German terrestrial land under Natural Earth50m and WorldCover2021 classification (much of the rest lies outside Germany). The original **44-cell posthoc signal** is retained in the older report; do not compare it directly to new 43-cell errors as if the sample were identical. **Same-sample comparison** uses **43 cells**, DE28/AT15, across the same 22 2° geographical block identifiers.

## 3. Geographically heldout predictions on SAME original 43 source-quality regions

Response: original observed log-odds of the proportion of Lycaenidae among three family **iNaturalist photograph reporting records**, photographed/research-grade/open GPS/≤500m, dated 2022 April–October. Country, geographic source selection, original prior-year plant photo effort data and photographic outcomes are unchanged in source IDs from original analysis. Train-only centered/scaled covariates, unweighted source cells, 2° geographic country-block-heldout OLS.

**Baseline M0 no WorldCover:** country indicator, 2021 Plantae photos, latitude and longitude.

**Baseline M0+WorldCover2021:** those same predictors **plus original 2021 country-clipped built-up land fraction, cropland fraction and log country terrestrial land area**. Source-added model = same baseline + each separate within-country-standardized log1p annual EEA original reported P mass, N mass or source count.

| Same-43-cell source prediction | Without independent land cover: heldout MSE improvement | With 2021 built/crop/country land: heldout improvement | New 512 source model MSE |
| --- | ---: | ---: | ---: |
| Source WWTP count | **+3.90%** | **−1.52%** | 0.158913 |
| Reported N t/year sum | +4.55% | **+8.03%** | 0.143962 |
| Reported P t/year sum | **+16.56%** | **+13.67%** | 0.135137 |

- **WorldCover-added baseline without WWTP** alone improves 2° heldout MSE from **0.195713→0.156538**, **+20.02%** relative to same-43-cell original non-landcover baseline; indicating real prior-year landcover/land-area structure is predictive of butterfly photographed reporting composition.
- Adding the original reported P source mass further improves heldout MSE **0.156538→0.135137** = **+13.67%** *conditional on country, 2021 Plantae photo effort, coordinates, 2021 built-up, 2021 cropland and country land-area* (posthoc).
- Regional breakdown for P source on same source-quality 43 cells: **DE28 +10.44%**, **AT15 +24.81%**.
- Original 256×256 WorldCover resampled-class sensitivity gives P source heldout MSE gain **+12.36%**, near 512 original **+13.67%**; N +7.07% (vs +8.03%); facility count −1.32% (vs −1.52%). Classification resolution change therefore does not reverse direction in this cohort.
- Important contrast: **WWTP count is no longer predictive after 2021 built/crop/land area control**. Treat any older facility-number association as plausibly explained by landscape/effort structure rather than source nutrient biology.

[Completed original WorldCover 2021 source-control OOS analysis GitHub Actions SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42464582745) | [full identical-cell 2° heldout P/N/count/256 results, original source QA, per-cell values](2022_LYCAENIDAE_ESAWC2021_DE_AT_URBAN_CROP_LANDAREA_POSTHOC_BLOCKCV.json).

## 4. Follow-up source-land area and source-P falsification, both favorable AND unfavorable

All these controls were newly fixed *before* this landcover-adjusted model executed, but as posthoc investigation after prior positive original EEA–photo associations.

- Original EEA P mass and independent 2021 WorldCover built share Spearman: **AT +0.568**, **DE +0.423**. P mass and 2021 cropland share: **AT +0.439**, **DE +0.210**. Independent land cover captures some, but not necessarily all, settlement / agricultural structure.
- Exclude strictly the country’s **one original largest P source-mass cell** without selecting on butterfly results (same resulting 41 source-quality cells): P-added heldout MSE improvement is still **+11.02%** over the landcover-adjusted no-WWTP baseline.
- **Country-preserving shuffle of reported P among regions** 199 fixed draws: observed P gain +13.67%; fraction null gains ≥ observed **0.005** (plus-one finite-sample adjustment); null median **−5.21%**. This tests geographic reporting-source alignment, *not randomized effluent exposure*.
- **Within 2° regional-block source-P shuffle** 199 fixed draws: fraction null gains ≥ observed **0.350**, null median **+11.28%**. **This is an important failure of within-2° spatial specificity**: much of the model’s predictive information is regional/geographic, not uniquely associated with particular neighboring 1° source pixels.
- Repeated 999 resamples of the existing 22-block prediction-error residuals: P incremental source gain **descriptive 95% interval −14.18% to +33.04%**; includes no prediction gain. This is not a new independent geographic experiment and not a conventional 95% causal or coefficient CI.
- Source original N and P annual masses were already strongly correlated across original photo-eligible cells (rank **ρ ≈ +0.888**), preventing clean source mechanism assignment. Reported annual WWTP P is not measured phosphate in a mud-puddling patch.

## 5. Defensible conclusion vs claims that remain unsupported

**New supported observational finding:** In 2022 across Germany/Austria’s *source-positive* EEA WWTP 1° geographic cells, reported original WWTP P-source geography provides a **posthoc incremental spatial prediction of Lycaenidae photographic reporting composition**, even after independent country-clipped **ESA WorldCover 2021 urban and cropland cover**, geographic country-land area, country, lat/lon and 2021 Plantae photographing effort. P has apparent larger predictive gain than facility count under this one fixed observational model; facility counts show negative additional gain once environment enters.

**Not yet supported:** A physiological **P nutrient effect on Lycaenidae**; that wastewater causes more adult butterflies, larval growth, ant occupancy or mud-puddling; that effluent reaches these land regions; that P vs N is biologically the mechanism; that photo fractions reflect actual abundance; or any geographic mechanism localized more finely than the 2° broad source-pattern signal. The restricted within-block shuffle, wide resample range, missing **butterfly-targeted observer effort** and true host-plant abundances/crop and built classification error remain major limitations.

**Recommended NEXT genuinely independent audit:** obtain **same geographic 2021 butterfly three-family photo counts** (source-fixed 53 cells and official iNat taxon descendants, original same date/accuracy/country filters), before fitting any new longitudinal composition model. If stable 2021 butterfly photo composition already predicts 2022 source association, this argues strongly for long-lived observer/host/urban geography instead of year-specific nutrient-response evidence. This is a new POSTHOC temporal-recording sensitivity, not 2018/2022 causal wastewater time series.

**Do not erase:** previous local 1km 1/24 and 5km 4/17 photo pair-support failures, Norway/Sweden/Flanders BMS near-outfall source failures, original positive aquatic 2018/2022 EQR predictive gains, or frozen `chocho` potential introduced-host resource expansion 22/25 Lycaenidae. None of those is the same ecological response/estimand.

## References and provenance

- [Original ESA WorldCover 2021 COG official access](https://esa-worldcover.org/en/data-access) and [ESA WorldCover dataset reference DOI](https://doi.org/10.5281/zenodo.7254221).
- [Original Natural Earth 50m country boundary source](https://github.com/nvkelso/natural-earth-vector/blob/master/geojson/ne_50m_admin_0_countries.geojson).
- [Prior 2022/2021 photograph–EEA exploratory source and 2° holdouts](https://github.com/zuizui0223/gesui/pull/24).
- [Full explicit posthoc landcover QA specification](../../docs/LYCAENIDAE_WORLDCOVER2021_INDEPENDENT_LANDCOVER_POSTHOC_FALSIFICATION_LOCK.md).
- Public source no restricted butterfly positions, identities, or pixels downloaded; only original total counts and source provenance.
