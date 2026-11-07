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


## gesui独立解析の新規結果台帳（2026-11-07、元chochoから移した結果とは別）

| 分類 | 出典 | 実測値 | 生態学上の判定 |
|---|---|---|---|
| 物理的な放流口↔河段の対応 | [河川ソースPR #3](https://github.com/zuizui0223/gesui/pull/3)／[成功ソース実行](https://github.com/zuizui0223/gesui/actions/runs/42364866919) | 欧州候補4,713/4,713 HYRIV_ID照合、4,497連続4河段、3,569非重複候補、16放流口。選定放流口の河道ずれ中央値25m・最大47m | **肯定的なGISデータ互換性の結果**。第三者GeoParquet未独立の公式原本検証。物質の到達や生物の移動を証明しない |
| 初回20地点の外来寄主非限定植物 | [植物初回PR #2](https://github.com/zuizui0223/gesui/pull/2) | 5か国ペア・20地点。処理水割合×放流口–施設植物種記録希薄化差ρ−0.150、n=9、p=0.700 | 初回の記録機会試験。地理的位置が異なり上下流・窒素曝露ではない |
| 新規河川植生の記録インベントリ | [指向性植物PR #4](https://github.com/zuizui0223/gesui/pull/4)／[原本取得](https://github.com/zuizui0223/gesui/actions/runs/42365305227) | 16放流口×3地点=48地点、1,402植物speciesKey、5,394地点×種の記録行。48/48地点に記録あり | **肯定的な新規データ構築・全植物記録の利用可能性**。すべて維管束植物とは限らず、記録＝植生被度でない |
| 放流口を挟んだ種数変化と処理水割合 | [PR #4](https://github.com/zuizui0223/gesui/pull/4)／[独立監査](https://github.com/zuizui0223/gesui/actions/runs/42365622805) | 厳格30記録に揃えるとρ+0.125874、n=12、p=0.69668。10記録感度ρ+0.173626、n=14、p=0.55277。高–低処理水割合差−0.05968（n高7・低5） | **明確な正の下水生態効果は得られず**。上下流の河段長・都市化・採集量・他排水の交絡あり |
| 実測窒素・リンとの生態的接続 | [Issue #5](https://github.com/zuizui0223/gesui/issues/5) | データ候補：Peters et al.2026 TN/TP・EEA Waterbase 2026 | **未取得・未照合・未検定**。河川の観測N濃度は施設放流N質量と区別する |

元chochoからの歴史的WWp・港・植物N指標の肯定的空間関連は上の新しいGIS／全Plantae結果と別の解析であり、再解析によって書き換えない。
