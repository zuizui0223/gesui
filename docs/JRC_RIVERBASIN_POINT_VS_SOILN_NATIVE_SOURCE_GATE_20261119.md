# 2026-11-19｜JRC流域単位の陸への窒素入力×河川への点源N：新しいSOURCE-onlyルート

**状態：独立した事後探索検証の SOURCE-first。既知の水陸NUTS2共通2×2 STOP（25地域・6か国）を再開しない。** 旧269 NUTS2で点源/陸源の部分的な地理差、2016 EQRで点源固有の予測改善が不確実という結果を見てから開始するため、本解析は独立な確認研究ではない。

## 元データ

欧州JRCの公式 ArcGIS REST の以下は、**地理粒度が NUTS2 ではなく river basin** と表示されている。

- 入力土壌（EMEP2015 0.1° raster から地域内平均）：layer **5**、MEAN_mgsqmy、元フィールド Basin_3035_HydroID、Basin_3035_AreaSqKm
- 入力土壌（DAYCENT2010–19 1km raster から農地上平均）：鉱質N layer **10**、有機N layer **15**、生物固定N layer **20**。MEAN_kgha、Basin_3035_HydroID、Basin_3035_AreaSqKm
- 河川へモデル点源N（JRC GREEN 2014–18 家庭＋産業）：Losses_to_water_MIS layer **5**、PsN(t)、**BasinID**、AreaSqKm

配布元：[Inputs to soil](https://water.jrc.ec.europa.eu/arcgis/rest/services/inmap/MIS_Inputs_to_soil/MapServer)、[Losses to water](https://water.jrc.ec.europa.eu/arcgis/rest/services/inmap/Losses_to_water_MIS/MapServer)。

**最重要：JRC土壌投入流域の HydroID と JRC GREENの BasinID、元植物/蝶の HydroBASINS HYBAS_ID_L07 は同じIDではない。** 現段階で文字列や整数の一致による自動結合は禁止。別々の元Polygonで観測単位と包含関係が証明された場合だけ対応付ける。流域内の平均モデル沈着・施肥は河岸植物への到達量ではない。処理場排水Nの河川下流実曝露でもない。

## SOURCE-Gate 1 (成果と数値結果を見る前の固定)

1. 公式5層のJSON schema / Source REST URL / featureCount / ObjectID総数・重複/サンプル25件に対する真正native fieldの一致を確認する。取得不全・HTTP error・source fieldの違いがあればSTOP。全取得が難しい場合、ここで確認済み範囲だけを発表し、未取得の地域を0にしない。
2. 土壌4層を同じ HydroID ごとに正しく一意に対応できるかだけを見る。**候補として100個以上の異なるHydroIDが同時に存在**し、各層の元レイヤーに1,000 HydroID以上、3主要源の当初元フィールドに真正数値が存在することを続行の最低条件とする。これは源地理の取得機会の判定であり、結果の生物予測をまだ認めない。
3. 点源BasinIDと土壌HydroID間の**same river basin polygon**は次の別ゲートでgeometry overlayし、土壌側HydroIDと点源BasinIDの面積比・相互重複を算出する。geometryの未取得、片方が複数帰属、JRC流域と元GBIF H7が他流域の集合になっている場合は結合STOP。任意のpointN値を先に置いてから相関の良い近傍を選ばない。
4. いずれの測定年も2015/2010–19/2014–18で、植物2017–19と年が一致する因果介入ではない。ソース取得に成功しても原植物 EIV/HNSの数値を開かない。SOURCE-Gate2で実地理に結んだ独立H6/国/環境・trait採点率の支持数が必要。

## 先に固定する結論区分

- Native HydroIDが土壌レイヤー間で一意に対応し元source数と欠測率が十分なら、地域集計データの地理的分離の**source-only比較可能性**を示せる。
- 点源BasinIDやGBIF H7とのgeometry一致が不十分なら、陸モデル源と水点源を共通流域で比較できたとは言わない。
- 新しい元流域sourceの原APIが取得不能なら、**利用可能な既存269 NUTS2**と既存STOPのまま論文に戻し、元の閾値の緩和や相関値での方向選択はしない。

本研究が「水域は点源、陸域は面源」と主張できるのは、独立地域外の両生物応答に対する追加予測に成功してからのみ。今回の SOURCE-only 機会監査はその前の段階である。
