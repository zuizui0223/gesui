# Year-specific reported EEA upstream WWTP N/P × H7 HNS/EQR/PlantS50

A posthoc **predictive, not causal** continuation of `gesui` Draft PR #89. Rejoins a fixed reported-source-only upstream WWTP annual N/P mass index from previous original Actions #43222864099 to the fixed 454 H7×period biological outcomes in the existing signed HNS panel. Distinguish source-reporting coverage from physical wastewater-free controls.

- Locked source JSON SHA256: `012ad79feabc8c11e0d74f803730b5b7913b96a64a6af26fb0c97263d9aa7c26`.
- Locked HNS original 454-H7-period source CSV SHA256: `1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca`.
- Authenticated original count: 454 rows/232 H7/90 H6; positive upstream N/P reporting in only 56/454 rows (38 H7/31 H6); 398 no-positive-report rows are **not** verified no-wastewater controls.
- Compare nested modelled-W baseline with report-presence, then source N/P intensity separately. HNS source score is based on GBIF occurrence composition and larval host plants' *habitat* N affinity—not actual observed leaf N or local butterfly abundance. EQR not emergence mass. Modelled flow-scaled source mass not measured downstream TN/TP.
- `scripts/source_join.py` checks original data hash, full 454 row identity, duplicate keys, country-mode discrepancy and original flow continuity. `scripts/analyze_reported_np.py` uses source-locked prior `analysis/20261115_multilayer/scripts/multilayer.py` H6 blocked fitting; `scripts/independent_sklearn_qa.py` independently refits all 18 comparisons.
- `scripts/selection_sensitivity.py` explores France, non-France, fixed original two 3-year windows, and six-group speciesKey source coverage at least 80%. These are post-result descriptive diagnostics.
- To reproduce on GitHub Actions, download the original `gesui-H7-454-reported-2018-2022-NP-to-site-flow-source-only` artifact from run #43222864099. No original third-party full raw source data redistributed. `GitHub Actions` checks source SHA before computing.
- A reported 2018 or 2022 annual source index does not fully align with the three-year GBIF 2017–2019 or 2020–2022 observation windows; the source does not deliver butterfly-site treatment exposure.

Do not merge this Draft or promote across-units SEM causal arrows. Negative prediction and missingness are valid results. eBird/TinyFish not used.
