# 2026-11-19｜JRC native basin × 元H7 WKT：12流域のSOURCE-only地理QA

**解釈上の最重要点：ソースが取得できることは判明したが、元JRC GeoJSONの無効形状が多いため、独立GIS検証なしに「H7曝露取得成功」とは言えない。** 生物NUTRI・HNS・EQR数値は読んでいない。元行政区2×2はSOURCE-STOP継続。

## 元sourceと実行
- JRC公式の原入力：EMEP2015大気N沈着（ラスター0.1°）、DAYCENT2010–19鉱質・有機施肥/生物固定（元ラスター1km）、GREEN2014–18水へのモデル点源N。各公開ArcGIS native river-basinレイヤーの元件数はそれぞれ **74,016 / 12,338 / 12,338 / 12,338 / 71,779**。数字は独立した生物調査地域数ではない。
- [公式JRC元レイヤーmetadata取得の成功run #44445257924](https://github.com/zuizui0223/gesui/actions/runs/44445257924)。元全件を取得しようとした旧run #44444896204は最初の74,016件が元コードの技術安全上限50,000を超えて停止。元生物SOURCEゲートは変更していない。
- 2017–19年の元226品質合格H7から、元H7 IDをソートし均等順位で12流域を固定して選択。元GH Actions run #43593227162の2,808 facet中のGBIF request元WKT（component_SHA照合）を使い、ArcGIS各流域ポリゴンをbbox空間クエリで取得。源によって実測年・母数・単位が異なることを維持。
- [初回ポリゴン取得run #44445737039](https://github.com/zuizui0223/gesui/actions/runs/44445737039)は2番目のH7でJRC GeoJSONの無効ポリゴンを検出して停止。 [修復感度付き再実行run #44445917915](https://github.com/zuizui0223/gesui/actions/runs/44445917915)は60/60のArcGIS bbox queryを正常完了。[元JSON/CSV・SHA artifact #14492973226](https://github.com/zuizui0223/gesui/actions/runs/44445917915/artifacts/14492973226)。

## パイロット12流域で確認できたこと
- 元12 H7のsource_countryはFR 8、HR 1、NO 1、NL 1、HU 1：**5か国**のみのパイロット。国数の代表性を仮定してはいけない。
- 12 H7 × 5 model源で**60/60**の元ソース・ポリゴン属性取得が成功。
- GEOS make_validを施した幾何形状を使用した場合、**12/12のH7で大気N、鉱質/有機Nの3つの元モデル平均が流域面積の90%以上をカバー**。原大気Nと水点源Nのcoverは各12/12ほぼ100%、肥料3種も最小99.9879%。
- **重要な地理品質の問題**：12 H7のうち**10 H7で修復したJRC源ポリゴンが流域面積の5%超を占める**。この10 H7では元ジオメトリをそのまま用いる厳格なsource quality gateを通過したとは言えない。残る2 H7はこの感度フラグなし。JRC native GISの修復後面積が妥当か、別形式のEsriJSONや公式原ネイティブ面積で独立検査が必要。
- このpilotは元数値ラスターピクセルの取得ではなく、**JRCが元ラスターから集計したnative basin average**とそのジオメトリの取得。H7平均はネイティブ川流域平均を面積で加重した二次的近似で、元ラスターの直接面積平均に等しいとは限らない。

## 次のハード判定
元GeoJSONとArcGISネイティブEsriJSONの幾何学比較、修復法別の重なり面積、公式native AreaSqKmとの一致をSOURCE-onlyで確かめるまでは**正式なH7→JRC源のGOを出さない**。H7全226流域の本調査へ一括実行する前に、形状修復による曝露のバイアスの程度を評価する。生物EQR・NUTRI・HNSの値を利用して地図や除外地域を選ばない。

main・草稿PR124、先行PR125/126を変更しない。
