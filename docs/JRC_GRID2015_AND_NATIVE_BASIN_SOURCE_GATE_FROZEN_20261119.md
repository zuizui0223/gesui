# 2026-11-19｜JRC 0.1°大気N沈着・1km施肥N：原流域ポリゴンのSOURCE-only適格監査

## 位置づけ
前段の2018底生×2017–19植物の共通2×2は25 NUTS2・6国でSOURCE-STOP。これを維持する。本件はJRC元格子・流域別モデル集計という別の地理単位を使えるかという新たなSOURCE-only探索。以前の肯定・否定結果を知った後の事後研究で、独立した事前登録ではない。植物NUTRI・蝶HNS・水生EQRの数値は今回は一切開かない。

## 公的元データ
JRC [Inputs to soil MapServer](https://water.jrc.ec.europa.eu/arcgis/rest/services/inmap/MIS_Inputs_to_soil/MapServer)：EMEP2015大気N沈着の元ラスターモデルは0.1°×0.1° EPSG4326、mg N/m²/年、DAYCENT2010–19施肥/固定の元ラスターは農地1km EPSG3035、kg N/ha/年。元Raster表示はlayer 21–24。**サービスがレンダリング画像を提供することと、真の数値ラスターデータの直接ダウンロードを提供することは異なる**。その数値再利用可能性は未確認。

JRCが既に平均したriver basin Feature layer：大気N沈着ID5、鉱質肥料N ID10、有機肥料N ID15、N固定ID20。GREENの水へのモデル点源N（家庭＋産業・2014–18年平均）は [Losses to water MapServer](https://water.jrc.ec.europa.eu/arcgis/rest/services/inmap/Losses_to_water_MIS/MapServer) のID5。これらの元流域IDの見かけが同じ/異なるだけでは空間単位の同一性を決められない。

## 元SOURCE-onlyの凍結条件
1. 各Feature layerの公式RESTメタデータ、Polygon元ID、元OBJECTID、数値フィールド、全件ID inventory・countを別々に取得。途中切断、重複、無効な負荷をSTOPとし、原本スキーマ・ページ応答SHA・地域数を保存。
2. 4つの陸源river basinレイヤーを一意の元HydroIDで比較。大気・鉱質・有機Nの数値3源が同じネイティブ流域で揃うもの**100以上**。各指標に地域間変動が十分あることを別に記録し、欠測を0にしない。N固定は副次。どれか元レイヤー取得不可ならSOURCE-STOP。
3. 水への点源Nの元BasinIDと陸源の元HydroIDは、**IDの文字列・数値だけでは結合しない**。独立の公式polygon geometryで同一範囲を証明する前に源のpoint/diffuse対応を推定しない。
4. JRC native river basinとHydroBASINS H7/真正H6も、互いに別地理製品であり単純IDや重心最近傍で合わせない。

## 後続の生物地理GO/STOP（元研究では未実施）
元2017–19 EIV品質適格H7を元の451/454可否source-onlyで復元し、元GBIFに使ったH7真のWKTを使ってJRC流域ポリゴンと面積90%以上の一意な重複を要求。互いに連結するH6とJRC元流域の地理成分でfoldを固定。最低SOURCE支持は独立JRC河川流域/交差成分100以上、品質適格H7 100以上、真正H6 40以上、8か国以上、国内に十分な曝露変動を持つ国5以上、5foldの各fold独立成分4以上、さらに水点源と陸源が正式に同一地理と証明できること。この条件は生物の値を見る前に凍結。未達ならSTOP。結果を見て空間境界・閾値を変更しない。

## 禁止事項
JRCの水への拡散源Nを陸上植物への施肥と呼ばない。陸への平均N沈着/施肥を河道への実到達量と呼ばない。家庭＋産業のモデル点源をWWTP固有実測値と呼ばない。年代/単位の違う負荷源を足して同じ媒介効果にしない。EQRの栄養耐性種構成、植物/蝶の実際の寄主利用や摂食フラックスをこのSOURCE監査から推定しない。

main、草稿PR124、SOURCE PR125、追加PR126の履歴とSOURCE-STOPを変更しない。
