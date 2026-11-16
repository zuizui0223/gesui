# Frozen conditional ecological model after chemistry-first official 2026 Part1 source support

Date 2026-11-16. This protocol is specified BEFORE evaluating plant NUTRI or butterfly HNS response coefficients for the new official Water Quality site-first population. It does NOT replace the prior Part2-only STOP or the original PR #107 strict native EQR-station and facility-reporting STOP.

## Scientific question

Do official measured ambient river TN/TP at ANY true EEA Water Quality station within the original 232 European HydroBASINS Level7 polygons predict separately:
1. species-equal plant European EIV NUTRI (original 4,011 strict accepted vascular plant scoring) and
2. butterfly species-equal larval-host N-habitat affinity HNS (original 186 accepted butterfly numeric species scores),
after original country/period, land-cover, hydrography, GBIF effort, trait completeness and static modeled W geography?

This is a **cross-sectional catchment-scale association**, not a WWTP intervention, water→plant resource flux, leaf nutrient measurement, realized larval herbivory or estimated causal mediation.

## Pre-biology source and science gates

1. Part1 full original 96.6M sample source audit and fixed 2018/2022 country-matched WQ station-year geometry must FIRST pass the separate source-only >=30 H7period /20 H7/20 H6/3 countries threshold with **same WQ station-year both TN and TP**, original physical sample count>=3 each, correct LOQ and official source SHA. Missing or ambiguous source remains unknown, never zero.
2. Join the original immutable 454 H7period GBIF source panel SHA256 1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca and original 451/EIV numeric trait-eligible source v2 (full original source SHA gate). Join exactly on native H7 and the two frozen 3year periods. No EEA benthic EQR station matching is introduced.
3. **Stricter science gate after source support, BEFORE coefficients:** >=100 native H7period rows with original EIV NUTRI and annual WQ both nutrients, >=40 true original H6 parent regions, >=5 independent original country groups, >=10 numeric plant taxa and >=40% original plant species-key EIV coverage; no duplicate H7×period. If gate fails, emit status STOP_BIOLOGICAL_SAMPLE_SUPPORT without fitting effects; no reduction of thresholds after inspecting p-values.
4. WQ country must equal the original majority country for that H7period; WQ values were site-equal medians at H7period level, drawn from source-valid stations only. Report n_WQ stations as opportunity quality adjustment. A different EEA Biology assessment method is NOT included as a putative property of these WQ stations.

## Main response and model (only after gates)

Primary independent association:
- Plant_N_incidence mean of original EIV_NUTRI across original GBIF-detected scoreable vascular species at each H7×period (not biomass).
- Predictor: log1p(observed official ambient TN_mgN_L), per one source-cohort SD.
- Controls: country, period, static original modelled W, urban, crop, mean river flow, naturalout, H7 area, original nPlant/nAnimal/nButterfly GBIF recording opportunities, source numeric EIV scored-species count, plant trait coverage, log1p(number of official strictly TNTP-valid WQ stations). Same source cohort for all comparisons.
- H6-equal weighted least squares, H6-cluster 95 descriptive interval; independent genuine-H6-disjoint deterministic 5-fold Ridge(alpha10) with train-only numeric scaling, country and period dummies, baseline vs predictor addition evaluated by H6-equal MSE gain and resampling fixed out-of-fold squared errors at genuine H6 level.
- Primary classification: point estimate and heldout MSE interval. Positive in-sample slope without country-heldout transfer is NOT confirmatory.

Secondary fixed tests:
- Substitute log1p(TP_mgP_L) independently (do NOT interpret TN vs TP as independent causal mechanisms given nutrient co-correlation).
- Butterfly species-equal HNS ~ Plant_N_incidence + original W/region/GBIF effort/source trait coverage + TN, and separately TP; evaluate conditional coefficient and incremental H6 prediction. Plant NUTRI is a geographical covariate, NOT a measured mediator.
- Plant NUTRI→Butterfly HNS coefficient with observed TN/TP included vs original within-cohort baseline as a descriptive sensitivity, not an indirect causal effect.
- Country-out heldout as an exploratory external replication if minimum independent-country support (>=5 groups) passes.

## Output truth boundaries

- Official 2026 EEA Waterbase Water Quality ICM Part1 may be incomplete for some countries and years; selection is not random.
- Water station is somewhere within H7, not necessarily upstream of riparian GBIF plants or butterflies. TN/TP is ambient river concentration not wastewater origin or dose, GBIF is observational recording not ecological counts or actual larval feeding.
- EEA source-year=2018/2022 is only approximately aligned to GBIF 2017–19/2020–22.
- Retain all negative/null cases, no posthoc precision improvement by selecting one favorable country, nutrient or period, no merging qualitatively different indicators as one wastewater meta-effect.
- If source or science gate fails, stop and write the precise source opportunity; the older true aquatic EQR field results and parallel trait filtering story remain separate.
