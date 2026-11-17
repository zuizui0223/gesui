# 2026-11-17｜FR263 公開データによる水陸マクロ解析：source-only追加作業

状態：**新しい真正水陸間の生態効果量は未推定**。元263地点の実測TP・底生EQR、元73 H7と24 H6は保持した。

## 完了済み
- 元EEA2022 source Git blob SHA1 325f199e985d9862ffee5b0edb205c61a7be3833と元2018/22真正263河川地点IDのexact joinは263/263。元座標263組、原HYRIV_ID 262種類。保存済み source-only CSV：results/20261117_FR263_FIXED_ORIGINAL_EEA2022_STATION_COORDINATES_SOURCE_ONLY.csv
- 地点最近接距離の中央値10.80km。異なる地点で1km以内のペアは0、2km以内5ペア。流域独立性の証明ではない。保存済みJSON：results/20261117_FR263_GEO_DISTANCE_PRE_RASTER_QA.json
- 元2018/22水質観測日原本はGitHub Actions #43373061409、SHA256 2190ac905326f4f2f4f278eb0debb4a7bfda7c4c71347d4941b54cc19d3de264。元数値はActions #43373639391、カレンダー月一致の既存解析は#43378917089。季節調整でΔTP→ΔEQR予測が弱い既存結果を維持。
- Copernicus HR-VPP 10m年間TPRODは2018/22双方をカバーするが、公式の原画素ダウンロードには認証が必要。元263地点での品質画素数はまだ未取得。製品：https://land.copernicus.eu/en/products/vegetation/high-resolution-total-productivity

## 新たに追加した再現可能な原本検査（現段階で実行完了は未確認）
1. Official HydroRIVERS Europe v1原本（67,648,957バイト、SHA256 500da7d36ceee0aa4c82dd625f5c53d80ba2f7ac14ad022cd3a6601bf55068e5）を取得し、元263地点を元HYRIV_IDの河道ポリラインに限定して照合。地点からその線へ250m以内・±500m沿線の有効長900m以上か、正しいH7/H6母数を点検。コード analysis/20261117_FR263_riverland/fr263_true_hydrorivers_geometry_qa.py。実行定義 .github/workflows/gesui-20261117-FR263-official-HydroRIVERS-geometry-source.yml
2. 代替の公開Element84 Earth Search Sentinel-2 L2Aの原本候補を、H6ごとに元IDハッシュ最小の1地点×24（事前固定）、2018年と2022年の6–8月を対にして探索する。観測点を生物の符号で選ばず、同じ地点両年で各2か月以上のscene支持がある親H6が18/24以上か検査する。sceneは各月最大2、衛星タイル雲量40%以下。**この条件はシーンメタデータのため、雲なし実画素支持を保証しない**。コード analysis/20261117_FR263_riverland/fr263_s2_2018_2022_stac_24h6_source_pilot.py。実行定義 .github/workflows/gesui-20261117-FR263-24H6-S2-2018-2022-source-metadata.yml

## 事前固定した生態検定への条件
- 元の設計 docs/PROTOCOL_20261117_FR263_AT67_HRVPP_MACRO_RIVER_LAND_DISTANCE.md に準拠：水質地点前後の河道500m、横断距離帯20–100、100–250、250–500、500–1000m。
- 新しいsource gate候補は、同年で各帯100画素以上を持つ2018年・2022年ペアが120地点以上、異なる40 H7・15 H6以上。品質失敗・重複・欠測を逐次記録。失敗した条件を有意性に基づいて緩めない。
- 真正の年別衛星指標が揃って初めて、near−far植生勾配に対する真の河川TPと底生EQRの関連をそれぞれ推定、24のH6を境界に地域外予測と2018–22差分を検証。
- Sentinel-2 NDVIとCopernicus HR-VPP TPRODは別の応答指標。後付けの代替を原初の確認的分析に格上げしない。Austria67の前回Sentinel-2 10地点5/10 < 6/10 technical STOPを変更しない。
- TPの河岸植生との関連は、下水窒素の輸送経路そのものではない。France263での実測TPは下水特異的でなく、2018/22の親河川からのfluxも測っていない。蝶・現地寄主摂食・生態系サービスも未観測。

**研究判定：** 水域の真正地点情報はGO。元河道線・新しい公開Sentinel画素についての実行結果は現在未確認。新しい水陸係数はない。 
