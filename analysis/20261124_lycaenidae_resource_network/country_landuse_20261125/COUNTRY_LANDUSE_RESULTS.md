# 2026-11-25: country-heldout terrestrial habitat-affinity correspondence after upstream and local landscape controls

Exploratory retrospective reanalysis of the original panel, not independent field replication. Rebuilt 207 paired H7 /81 H6 from original source archive hashes, the two accepted 94-species host-score rosters, and 4,011 plant EIV N scores. All 414 rows at each evidence cutoff replay previous archived plant and butterfly scores at 1e-12 tolerance. Restrict both periods to >=30 scored taxa outside all genera of documented hosts in both mapped catalogues: 205 H7 /81 H6; 1,876 eligible plant species keys. All endpoints use this fixed cohort. No new field samples or climate source. Existing upstream/local cropland and urban percentage only; climate data absent from the signed HydroATLAS subset. W is static modelled wastewater geography, NOT measured delivered dose. Plants and butterfly scores share a published EIV reference, even with documented host genera excluded. Not local actual feeding or abundance.

Primary H6 holdout fixed SHA256 modulo 5, independent sklearn ridge predictions max difference 3.51e-13; H6-equal weighted fitting, training-only standardization, alpha10. 5,000 fixed-loss parent bootstrap draws. For country holdout remove 13 H7 with mixed period country labels, remove all matching H6 from training and omit country effects (method retained). Source-only minimum 5 H7/3 H6 test and 30 H7/15 H6 training. Five eligible countries FR144 H7/34 H6, IT17/11, LT6/6, NL6/4, NO8/7 (181 H7, 62 distinct H6); other countries retained as potential training, all STOPs recorded. Equal-country loss ratio (not arithmetic mean of gains). Two evidence cutoffs separately, exact previous cohort. 252 result summary rows, 5,148 H6 loss rows. No p-value claims; descriptive fixed-loss resampling.

## Primary level 1–3, genus-excluded plant N added after W and upstream cropland/urban

| heldout country | old HNS gain % [H6 fixed-loss interval] | Clarke HNS |
|---|---|---|
| FR | 38.20 [29.25,45.81] | 36.91 [29.39,43.46] |
| IT | 41.65 [20.68,65.98] | 6.51 [-61.77,44.55] |
| LT | 53.38 [11.15,80.12] | 40.21 [-206.18,86.99] |
| NL | 18.74 [-21.24,37.32] | 27.24 [-11.66,41.03] |
| NO | 62.26 [50.25,69.84] | 64.74 [53.21,76.91] |
| Equal-country | 42.90 | 40.93 |

All five positive point estimates under both catalogues; only FR and NO intervals positive for both. Country resampling of 5 fixed country loss pairs: old [26.13,57.17], Clarke [21.56,56.19], omit-one-country score range old 36.41–52.37 and Clarke 29.08–45.48. Very small country count, overlapping training samples, not independent confirmatory intervals. Cutoff1–2 country equal old 42.64%, Clarke 41.22%. Local landuse substitution 50.52% old /48.08% Clarke. Original H6 gains reproduced old 22.74% [6.95,35.69] / Clarke 26.30% [7.82,40.97]; local landuse 24.21% [9.25,36.71] /28.36% [11.29,41.83]. Full plant N country equal old 44.55%/Clarke43.41%, but Clarke IT negative (−7.16%). No claim exclusion improves gains. Local landuse sensitivity not selected as new primary.

## Conditional landscape comparison, same 205 H7

| Added predictor | H6 old HNS | H6 Clarke HNS | country equal old | country equal Clarke |
|---|---:|---:|---:|---:|
| crop+urban after W | 6.58 | 5.11 | 19.34 | 21.97 |
| W after crop+urban | -1.53 | -2.81 | 2.88 | 5.66 |
| genus-excluded plant N after W and crop+urban | 22.74 | 26.30 | 42.90 | 40.93 |
| crop+urban after W and genus-excluded plant N | -2.58 | -5.59 | -3.65 | -4.02 |
| W after crop+urban and genus-excluded plant N | -3.37 | -3.44 | -0.80 | 1.94 |

The plant index offers predictive information beyond coarse land-use, not causal mediation. H6 W after landscape predicts plant N +15.79% [4.08,25.93], EQR +9.98% [1.12,17.30]; country equal +7.07% and -7.54% with fixed country resampling intervals crossing zero. Landuse after W country equal plant N +24.23% (uncertain), EQR -1.07% (uncertain). Country equal gain is not an estimate of a biological effect. Retain previous network nonreplication, FR measured TP-to-land null, and temporal/record-weighted nulls. No climate adjustment, no direct sewage effect, no cross-compartment causal link.

Replay: restore original archives and `gesui_fab_resource_consumer.py`, `gesui_rebuilt_network.csv` in source folder per SOURCE_MANIFEST.json. Copy committed filter_butterfly_scores.csv, filter_clarke_edges.csv and decompress filter_panel.csv.gz into prior output folder. `python gesui_filter_country_landuse.py SOURCE PRIOR OUTPUT`. All output derived locally; no specific remote replay yet.
