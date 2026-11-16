# Source-frozen conditional negative control: original upstream vs downstream reported WWTP N/P

- Study unit: fixed 454 European H7×period records, 232 H7, 90 authentic H6.
- Original outcomes: water EQR (site-mean ecological quality ratio), plant 50-record-rarefied GBIF taxon richness S50, butterfly host nutrient-affinity HNS (GBIF record-weighted; **not butterfly abundance**).
- Confounding source: reported upstream and downstream annual municipal N/P mass in the native same EEA field station network, modeled discharge-normalized. A 0 indicator means **no positive reported matched outfall**, NOT no sewage.
- Upstream vs downstream **nested** predictors in the same model, tested using 5-fold H6-group-disjoint Ridge α=10 with H6-equal MSE; exact original PR92 H6 folds and source covariates.
- Target estimates: 1) upstream incremental OOS gain after downstream covariate; 2) downstream incremental gain after upstream; 3) both vs fixed original source geography baseline. A source-reporting treatment effect is **not** identified.

## Original data dependencies

The source-frozen input files are already in GitHub artifacts and PR92 inherited code:

- Original 454 HNS/Plant/EQR original panel: `results/20261115_sync/hns_macro/H7_454_HNS_panel.csv`, SHA256 `1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca`.
- Original source-only upstream EEA actual annual reports: Actions #43222864099 and SHA256 `012ad79feabc8c11e0d74f803730b5b7913b96a64a6af26fb0c97263d9aa7c26`.
- Source-matched **downstream** directional negative-control original: Actions #43679382449 with source-only JSON and CSV; original EEA 2018/2022 source Git blobs pinned in PR92.
- Frozen source-layer numerical code: PR90 `analysis/20261115_reported_np_hns/scripts/source_join.py` and PR88 `analysis/20261115_multilayer/scripts/multilayer.py`, exact independent PR92 sklearn implementation `analysis/20261115_np_direction/scripts/independent_sklearn_qa.py`.

## Reproduce on GH runners

See `.github/workflows/gesui-20261115-np-conditional-placebo-H6.yml`: download the same signed original source artifacts and run 5 tests, 30 original Ridge 5fold predictions, 30 independent sklearn replicas, 45 archived source-MSE algebraic checks, and saved result comparison. The fixed H6 error bootstrap describes uncertainty conditional on frozen out-of-fold predictions; it does **not** refit the model in every bootstrap.

## Not claimed

No actual 2018/2022 receiver TN/TP measurements from exact original stations were added, no laboratory nutrient application, no independent untreated counterfactual, no EQR-based emergence biomass, no actual larval consumption/adult pudddling, and no fitted causal SEM. This is a posthoc stress test; don't present nominal bootstrap intervals as multiple-testing-confirmed intervention effects.
