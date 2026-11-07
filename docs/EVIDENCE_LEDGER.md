# 既存証拠台帳 — 元chocho PR #122の原本から引継ぎ

**数値は今回再推定していない。** 係数は標本・応答・説明変数が異なるため直接比較しない。

| ID | 結果 | 解析単位 | 読み方 |
|---|---|---|---|
| WWTP施設数25km | **β +0.1285**、国クラスタ95% CI [+0.0077,+0.2494] | 1,368港・52か国、候補外来／在来寄主のGBIF記録比 | 都市・植物索引量を入れると β+0.097、CI [−0.026,+0.220] |
| WWp（2015モデル下水発生量） | **β +0.05870**、CI [+0.01310,+0.10429] | 1,055港・42か国、WWTP・市街地・農地・記録量・国FE等調整 | 肯定的な**空間記録の関連**。R²増分+0.00171、Nの直接効果ではない |
| WWp−WWt（発生−計上処理量） | β −0.0673、CI [−0.1384,+0.0039] | 同都市調整標本 | 未処理下水の独立正効果は未支持 |
| Ellenberg N×記録港頻度 | 偏順位 **ρ +0.210** | 182植物種、国別普及・GBIF記録総数調整 | 正の栄養性指標と港記録頻度の対応。地元のN投入量でない |
| 全植物のN指標相対時系列 | **+0.1494**、島ブートストラップCI [+0.053,+0.252] | 欧州6島・同一27港／非港沿岸ペア（早期vs後期） | 植物記録構成差。観察120レコード等重みでは−0.094、80種等重みでは+0.136 |
| NIES国別未管理N×WWp | β −0.0324、CI [−0.0681,+0.0033] | 993港・38か国 | 正の衛生経路N増幅は未支持 |
| CIRAD人由来Nポテンシャル | β +0.02901、CI [−0.05163,+0.10966] | 1,041港・42か国 | 独立正効果は不明 |
| CIRAD農地Nポテンシャル | **β +0.03993**、CI [+0.01335,+0.06651] | 同1,041港、平滑化外来／在来記録対数比 | 農業Nの探索的な正関連。外来割合fractional-logit β+0.01421 [−0.03923,+0.06765]、高品質ラスターでは β−0.00402 [−0.06921,+0.06117] |
| SoilGrids試行 | 8沿岸画素の7つがゼロ／マスク | 4港／非港試験ペア | 当地の土壌N効果は推定不可 |

## 完全な元レポート
- [HydroWASTE・国固定効果](../archive/chocho_pr122/provenance/reviewer_defenses/results/global_hydrowaste_port_hosts_20261102/README.md)
- [Jonesモデル下水×WorldCover](../archive/chocho_pr122/provenance/reviewer_defenses/results/global_hydrowaste_port_hosts_20261102/JONES2015_WORLDCOVER2021_GLOBAL_WW_SOURCE_AUDIT_20261103.md)
- [NIESとCIRADの独立Nソース](../archive/chocho_pr122/provenance/reviewer_defenses/results/global_hydrowaste_port_hosts_20261102/NIES_CIRAD2022_SPATIAL_N_SOURCE_AUDIT_20261105.md)
- [農地Nの応答分母監査](../archive/chocho_pr122/provenance/reviewer_defenses/results/global_hydrowaste_port_hosts_20261102/CIRAD2022_CROPLAND_N_ESTIMAND_AUDIT_20261105.md)
- [植物全体のN背景](../archive/chocho_pr122/provenance/reviewer_defenses/results/port_nutrient_trait_pilot_20261102/ALL_FLORA_BACKGROUND_N_SOURCE_AUDIT_20261102.md)
- [記録と種数の希薄化](../archive/chocho_pr122/provenance/reviewer_defenses/results/port_nutrient_trait_pilot_20261102/RECORD_SPECIES_RAREFACTION_SENSITIVITY_20261102.md)

chochoの全球潜在寄主資源+54.9%、155種予測MSE改善4.12%、51島の寄主–チョウ類似ρ+0.300は独立した肯定結果。gesuiの「下水Nの因果効果」へ言い換えない。
