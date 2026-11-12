# Repeated European watershed macroecology: benthic EQR, plant S50 and Lycaenidae temporal co-response

12 November 2026 | Source-locked original-data analysis, EXPLORATORY/POSTHOC | Draft unmerged

## Scientific question and original comparable scale

When actual benthic river quality changes between 2018 and 2022 at *exact identical monitoring stations and biological assessment methods*, do observational flora richness and relative butterfly-family reporting change in the same HydroBASINS Level7 river subbasin, in the corresponding 2017–2019 and 2020–2022 GBIF index periods?

This is a macro study of the SAME native Hydrological Basin×Time units, not an assumed causal serial sewage→benthos→plants→butterflies chain.

## Original reproducible source / strict support gate

- Entire original EU 2018/2022 benthic numeric EQR data from gesui PR #17/#14, and original authentic HydroBASINS level7 Europe polygon/source mapping from PR #62.
- Exact original 234 H7 subbasins with both original EQR year study sources, without selecting by response sign. 468 original basin-periods.
- Genuine original GBIF Plantae speciesKey facets calculate expected richness in exactly 50 indexed species-assigned records per original period, and original five accepted butterfly families (Lycaenidae, Nymphalidae, Pieridae, Papilionidae, Hesperiidae). **455/468** basin-periods have S50 and >=20 five-family records, and **223/234** basins have both periods.
- A different source-only audit froze **1,047 same original physical river station keys with unchanged assessment method and source country in both EQR years**, distributed across only 174 H7 basins and 7 countries. Within these, 170 basins had both period GBIF S50+five-family record support. One of these lacked the original source Animalia observer-opportunity recording covariate. It was excluded as genuine source NA (not zero), leaving **169 unique H7 basins and 1,023 exact paired river stations**.
- Source-exact SHA (before numerical effects): 234 source-reconciled candidate file SHA256 07fa9744cd4e9574587a7b25c4ed4fa7f3183e377cac61c28723a4f24a86c875, original unchanged-monitoring-station manifest SHA256 373a3a96ea790a7592cc4d10b0e4f9244fbdddab2c1af4901118792674874930. Original EEA2018 Git blob 0f30bfbc3dac13f6a182c8a910e9899c04778322; 2022 325f199e985d9862ffee5b0edb205c61a7be3833.
- **148/169 basins (87.57%) are in France**. Only 21 non-France H7 basins remain, which cannot support robust pan-European transfer.

Source-only 234 study support [Actions #43212724530](https://github.com/zuizui0223/gesui/actions/runs/43212724530), original same station/method study [#43204743959](https://github.com/zuizui0223/gesui/actions/runs/43204743959).

## Original results: direct within-H7 biodiversity-layer temporal cochange

Basin is the analytical unit. For each original H7 basin, average the exact matched physical-station EQR(2022)−EQR(2018) once, compare with S50(2020–22)−S50(2017–19) and the change in Lycaenidae log-odds out of five families. Spearman original basin ranks:

| Exact basin change association | Original Spearman rho |
|---|---:|
| Delta EQR vs Delta Plant speciesKey-S50 | **−0.05358**, descriptive p 0.489 |
| Delta EQR vs Delta Lycaenidae five-family log odds | **+0.03503**, descriptive p 0.651 |
| Delta Plant S50 vs Delta Lycaenidae log odds | **+0.01338**, descriptive p 0.863 |

Original source-controls baseline predicts the 169 H7 basin EQR changes with source country, log(area), static modeled HydroWASTE outlet-pressure density, original number of matched physical monitoring stations, Delta log indexed Animalia, Delta log indexed Plantae and Delta log five-family butterfly observer opportunities. All extra predictors were evaluated in the exact same rows using H7-ID disjoint 5-fold Ridge10, training-only standardization:

| Added GBIF observation-response layer | Fivefold out-of-basin EQR change MSE improvement relative baseline | Original 2k H7 basin bootstrap descriptive 95% |
|---|---:|---|
| Plant Delta S50 | **−0.278%** | [−3.121%,+1.849%] |
| Butterfly Delta Lycaenidae log odds | **−0.966%** | [−2.545%,+0.377%] |
| Both | **−1.204%** | [−4.419%,+1.290%] |

Partial HC3 descriptive regression effects: Plant Delta S50 beta (per original predictor SD) **−0.00856 EQR** 95 [−0.02512,+0.00801]; Lycaenidae log odds beta **−0.00100 EQR** 95 [−0.02295,+0.02095]. Neither excludes 0; country-equal total weight results also change magnitude/direction and include zero.

**Empirical conclusion:** no stable SHORT-TERM cochange was detected across these 169 historical, method-matched basin units. This is NOT evidence that cross-boundary insect/plant interactions are biologically absent. EEA source columns have possible study-year scaling differences even within unchanged methods, and opportunistic GBIF sample/recording changes are imperfect flora/butterfly indicators.

Main authentic outcomes and each original paired basin source CSV [Actions #43214039321 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43214039321). A completely separate second Python program reimplemented all three original Spearman rank correlations and all FOUR original 5fold CV MSE calculations exactly, with no new source download [Actions #43216686503 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43216686503). Independent audit is NUMERICAL, not independent field study.

## Confounds, evidence line, and next macro test

- HydroWASTE is worldwide registered/modelled outlet pressure to H7 *BASIN OUTLET* and does not measure an upstream dose actually delivered to any interior benthic monitoring station, nor actual 2018/2022 effluent N/P emissions. Static W cannot be a temporal difference exposure.
- The original GBIF botanical S50 is indexed species-key sampled richness, not field-established native flora, floral nectar, larval host or ecosystem service. Lycaenidae five-family reporting share is not butterfly abundance or actual pollination.
- Basin hashed CV is not contiguous region/country-independent ecological transfer.
- Country/method sampling is heavily biased toward France; a claim about all Europe or global systems cannot be made.
- Separate observational layer changes did not improve original EQR change predictions, but **spatial** cross-sectional source-gradient effects may still exist. These must be tested and reported as a distinct exploratory estimand, never rebranded as a supported serial causal chain.

Keep Draft. Never merge to gesui/main or chocho; retain all original positive/negative datasets and prior PR #62/63 results. No eBird.
