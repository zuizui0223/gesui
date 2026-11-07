# データソースと再現可能性

| 公開データ | 真の単位・生態学的意味 | 今回の扱い |
|---|---|---|
| [HydroWASTE v1](https://www.hydrosheds.org/products/hydrowaste) | 58,502 WWTP、位置・推定放流先・処理段階・排水量 | 施設・曝露候補であり局所N測定ではない |
| [HydroRIVERS v1](https://www.hydrosheds.org/products/hydrorivers) | 流域面積≥10km²又は平均流量≥0.1m³/sの河川網 | 小側溝までカバーすると仮定しない |
| [Jones 2015グリッド](https://doi.org/10.1594/PANGAEA.918731) | 下水発生/回収/処理/再利用、**百万m³/年**、5分格子 | WWpはモデル水量≠N負荷 |
| [NIES GNINPO-DSM](https://doi.org/10.17595/20260227.001) | 2005–2022国×年の衛生経路N、**kg N/年** | 国別値は港別局所Nでない |
| [CIRAD/Wassenaar](https://doi.org/10.18167/DVN1/AWQNSR) | 2022年約3km、人排泄N/農地Nの供給**ポテンシャル** | 水域へ実際に届いた量でない |
| ESA WorldCover 2021 | 10m土地被覆 | 2001–24植物記録と時期が一致しない |
| WPI 2019×GBIF | 港座標、2001–24植物**記録種数** | 現存量・定着・生物繁殖ではない |
| Ellenberg N | 主に中欧の植物の養分状態指標 | 窒素濃度や量ではない |

Jones原本ZIP SHA256：c413f7493fb5c2b5a07d3d974abf3b9af25a32fcadd8eb839cc3fe302838e5eb。元のソースロック・統計検証は archive/chocho_pr122/provenance/ 以下に保存。

## 移行が保証するもの
元PR #122の**変更された80件のテキストファイルのみ**を原本のままコピーしGit blob SHAを80/80検証。内訳 25旧workflow、21scripts、18tests、16audits。workflow YAMLは archive 内で無効。

**移行が保証しないもの：** 元のchocho/mainにある共通WPI/GBIF候補植物CSV等、GH Actions artifactに保存された港別の出力表、外部原データの自動ダウンロード、Python依存パッケージ、移植先でのテスト成功。旧スクリプトにはchochoの相対パスを使うものがあるため、コードとレポートを保存しただけでは完全に独立再実行できない。原本artifactは30日保持の場合があり、再現には固定元ソースからの再構築が必要。

## 次に必要な実測表
**キー：地点×観測期間×受理植物種ID。** 列：site_id、catchment_id、island_id、period、accepted_taxon_id、origin/native/alien/unknown、habitat_aquatic/riparian/terrestrial、occurrence_countとeffort、N_traitとsource、outfall_id、hydrologically_downstream、along_network_km、modeled_ww_volume、reported_wwtp_quality、**measured_effluent_N_kg_per_year_if_available**、cropland_N、built_fraction、population、climate、source_QC。

NoDataをゼロに置き換えない。海岸線で施設点と推定放流口の平均約6.5kmのずれを無視しない。GBIFのゼロ記録と真の不在は区別する。下水量に不明なN濃度を掛けて実測負荷として報告しない。
