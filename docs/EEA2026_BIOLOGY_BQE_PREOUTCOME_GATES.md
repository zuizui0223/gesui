# Independent freshwater biological responses: EEA Waterbase Biology 2026 (preregistered source gate)

## Study question
Does **officially reported measured WWTP phosphorus/nitrogen discharge mass** have an observable spatial association with the ecological quality of **benthic macroinvertebrates, aquatic macrophytes/phytobenthos, and fish** across European rivers, after hydrology, country-specific classification and land use are controlled?

This is a separate biological response from the retrospective GBIF Plantae/EIV river pilot, whose original positive-source comparisons failed their strict trait/record effort gates. Do not choose EQR taxa, sites or reference years because they produce a favorable sign.

## Independently identified public monitoring source

- European Environment Agency [Waterbase – Biology 2026](https://www.eea.europa.eu/en/datahub/datahubitem-view/7bd881f0-59ca-4f0a-a43c-7e0cd0c66d6c), published **1 July 2026**, with reported data through **2024**. [Official source archive](https://sdi.eea.europa.eu/data/cc6cfbbb-58f1-4ad2-a1e1-f9f86734f0d6) (version identifier **eea_t_waterbase-biology-2026_p_1990-2024_v01_r00**). EQR data reported from 2012–2024 are normalized *biological quality element* indices, not numbers of unique species or individual abundance.
- The official EEA archive has spatial-object reference data; exact EQR site IDs, coordinates, missingness and classification-method fields are **not yet verified until the source inventory is run**. File checksum, size, dataset version, format, type and source output structure must be frozen before interpreting biological values. The current workflow never exports full raw records.
- Original WWTP source is [EEA UWWTD 2022 published measured N/P annual discharge](https://www.eea.europa.eu/en/datahub/datahubitem-view/6244937d-1c2c-47f5-bdf1-33ca01ff1715), filtered to unique direct freshwater (FW) outfalls, **2,969** previously mapped to source-hash-locked Europe HydroRIVERS v1-derived graph. [Source-only run](https://github.com/zuizui0223/gesui/actions/runs/42374105184).
- The river network is based on a third-party derived GeoParquet; matching all original EEA HydroRIVERS topology against the official vector remains a limitation, not proof of true effluent delivery.

## Before reading biological EQR results: eligibility gates

1. Keep **River** biological samples separated from lakes, transitional/coastal waters, canals/estuaries and data without verified waterbody category. Primary source loads are restricted to FW discharge points.
2. Use **2022** biological responses as the primary time match to official **2022** reported WWTP loads. A declared secondary window 2020–2024 may improve biological availability but is not a one-year causal before/after design.
3. Link BQE (benthic invertebrates, other aquatic flora/macrophytes/phytobenthos, fish) to river reach **using documented waterbody/site IDs and verified GIS**; nearest coordinates alone cannot prove source exposure. Flag missing coordinates, coarsened geolocation, duplicated IDs, monitored waterbody length, estuary transitions and unknown source quality. **A site cannot be matched to both source and control without a known river direction.**
4. Sum measured EEA annual N/P tonnes/year for unique direct-FW plant outlets **hydrologically upstream within 50km**. Same-reach facility allocation is excluded pending ordering within reach. Contrast physical downstream WWTP load as negative-control covariate. Source zero means *no registered eligible EEA measured-load WWTP* within that network, not an unpolluted river.
5. Independent environmental adjustment: river flow and upstream drainage area, geographic location, agriculture and population density when supported, waterbody type/size and spatial blocks. Country/biological-method class, original EQR calibration and reporting class must be controlled **without using the biological outcomes to tune source groups or select favorable countries**.

## Predeclared minimum source and quality support

For each BQE independently require **at least 100 unique river waterbodies/monitoring sites**, **at least 20 sites with measured upstream WWTP load**, **at least 5 independently reporting countries, at least 3 countries with a positive registered upstream source**, and **at least 10 spatial 3-degree blocks**, after dated 2022 alignment, correct EQR validity/QC and deduplication. Otherwise report not-testable as a data-coverage result. Do not change these bounds after seeing EQR signs. Biological sites may share waterbodies; the effective independent unit must be documented.

Only when gates pass, fit fixed baseline predicting normalized within-method EQR (if published as such) from country/method/bioregion and river Q/area/geography, with and without log1p(upstream WWTP reported P tonnes/year); run a separate N predictor. Evaluate within-country/bioregion heldout blocks with exact common denominators and consistent preprocessing. Compare to downstream-source negative controls, assess source-negative coverage and uncertainty. If EQR scales turn out noncomparable across countries or methods, restrict to standardized within-stratum contrasts or fail the cross-country test rather than treating raw scores as exchangeable.

## Ecological interpretation

A putative lower EQR in higher reported-load reaches could reflect sewage effects, but also agricultural P/N, channel modifications, other industrial discharges, wastewater treatment removal, stream temperature, flow regime and imperfect monitoring. An unchanged EQR would not prove no sewage effects at seasonal flow or individual species level.

**Strict claims boundary:** EQR is a nationally calibrated biological quality indicator; not population size, species richness, macroinvertebrate abundance or proof of local fish declines. The initial step is to verify genuine public measured EQR and its spatiotemporal topology. Only positive *out-of-region* evidence after negative controls, confounding, calibration and power checks could support a strong ecological association. Even this cannot establish causal WWTP discharge without experimental or quasi-experimental design.


## Exact numerical endpoint and model freeze (before observing EQR outcomes)

The original biological indicator has normalized EQR (higher generally means reference-like ecological condition, not number of species). Preserve reported scores: accept finite original EQR values >=0, do **not** arbitrarily clip >1 (flag out-of-usual-range), exclude nonnumeric/missing with counts. One unit is the original unique (country, monitoringSiteIdentifier, monitoringSiteIdentifierScheme, BQE code, reference year 2022). Repeated unit rows require separate source QC, never select the outcome-sign-favorable row.

Fixed independent response classes:
- **EEA_13-01-4** InvertebrateEQR_G (benthic invertebrates).
- **EEA_124-04-9** PhytobenthosEQR_E.
- **EEA_14-02-8** FishEQR_H, independently subject to the source support gate; no pooling with alternate fish EQR code.

The independent 2022 original EEA reported P and N loads, hydrologic upstream 50km source counts, downstream WWTP placebo, river mean discharge Q, upstream drainage area, observed country, original procedure classification system, geodetic lat/lon and deterministic region 3-degree fold labels are **source-only** and fixed before the EQR values are retrieved. Source-zero is zero registered eligible both-measured outlet in the source graph, not a sewage-free river.

Only BQE cohorts with >=100 unique monitoring sites, >=20 registered upstream source-positive sites, >=5 total reporting countries, >=3 positive-source countries and >=10 distinct 3-degree geographic blocks proceed to fitted EQR comparisons.

For every eligible BQE:
1. Response y = original station-year-2022 normalized EQR, never plant species richness or invertebrate abundance.
2. Fixed Ridge(alpha=10) five-fold regional holdout. Baseline M0 numeric lat/lon, log1p(river mean Q), log1p(drainage area); categorical countryCode and countryCode crossed with reported classification procedure. Missing method is a fixed UNKNOWN level. Fold-specific OneHotEncoder(handle_unknown=ignore) and StandardScaler are fit on **training only**. Original deterministic 3-degree block-to-fold map is stored from source-only sampling to prevent library-version instability.
3. P series M1 adds log1p(upstream reported WWTP measured annual P tonnes/year); M2 adds log1p(hydrologically downstream reported P load) negative control; M3 adds both. Independently repeat these fixed covariate comparisons using N. Do not choose whichever nutrient/sign gives better significance.
4. Report exactly the same withheld monitoring sites for M0 through M3, MSE and MAE, paired changes and fold/region support. Also report descriptive source-rank Spearman and block/country-stratified residual signs **without causal p-value claims**.
5. Record distinct waterbody IDs as effective grouping and report sensitivity to multiple monitoring sites in one waterbody. If country/method scales are noncomparable, restrict to within-class comparisons, not pooled unadjusted EQR interpretations.
6. Never declare downstream placebo improves/does not improve based on which outcome class is favorable; retain all classes and both N/P results.

A correlation or predictive increment alone cannot identify causal source N/P effects on aquatic community status: river engineering, farm/urban runoff, industrial discharge and monitor placement remain potential confounds.
