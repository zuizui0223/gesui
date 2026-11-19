# 2026-11-19｜JRC native流域の元EsriJSON×GeoJSON 3流域独立地理QA

**段階：SOURCE-only 技術検証。生物効果も全H7源支持も未推定。** [本番Actions #44446973908](https://github.com/zuizui0223/gesui/actions/runs/44446973908) SUCCESS。
[原本の符号付きリング/GeoJSON別アルゴリズム・元面積比較データ artifact #14492979186](https://github.com/zuizui0223/gesui/actions/runs/44446973908/artifacts/14492979186)。artifact ZIP SHA256 7811b873be488d7f2f8449d9dd8923389d4d5bd03e6b010702d2384c0c8fdf5c。

先の12-H7 SOURCE pilotでは、大気沈着と施肥3源は修復後12/12でH7の面積90%以上に対応していたが、10/12で元JRC GeoJSON不正形状が修復後面積の5%以上を占めていた。このため、そこから直接生物回帰に進めなかった。

今回は**初回パイロットで事前に観測された形状影響に基づいて**次の3 H7を技術診断として選択し（独立ランダム標本でも、前向き確認的検証でもない）、JRC元ArcGISレスポンスのEsriJSON Polygon rings、元GeoJSON→GEOS make_valid、元GeoJSON→buffer(0)の3手法でH7面積カバー率を別々に計算した：

- H7 2070012580：以前の修復不要群（HR）
- H7 2070021030：修復がH7面積約96.9%に関わる群（FR）
- H7 2070498990：修復がH7面積100%に関わる群（FR）

それぞれ2015大気N沈着、DAYCENT鉱質N施肥、GREEN水点源Nの3源を確かめ、**9/9の独立ソース値・元feature OID一致、3種の元形状解釈でH7面積カバー率1パーセントポイント以内に一致**。3流域でのcoverは各モデル99.9879–100%（修復後・独立Esri ringsとも同値）。各JRC native Polygonの公表AreaSqKmとmake_valid後のEPSG6933面積の差は**最大0.004%**。報告されたGeoJSONのinvalidフラグは、native地理単位が生物学的に不整合という意味ではなく、ポリゴンの複数リング等の表現上の問題である可能性がある。今回の3流域では元ring独立構築が形状補正後の結果と一致した。

**制限**：修復が必要だった残り8 H7の独立リング再確認は未完了。3流域は前の形状QA結果を知って選んだため代表性なし。元ラスター10km/1kmのピクセルを直接集計したのではなく、JRCが元格子から平均した各native river basinモデル値を元の流域面積で加重した二次的な近似。元JRC native basinと元HydroBASINS H7のGIS跨ぎはこのpilotで3 H7しか検証していない。水点源Nを陸域肥料と同じ物質経路の作用効果として比較したわけではない。

**科学的結論はまだ出せない。** 同じH7で点源と陸N投入の地理ソースが対応可能な例を示しただけで、元2018水×2017–19陸のNUTS2 2×2はSOURCE-STOPのまま。今後の新研究へ進むなら元226 H7について同様の**SOURCE-only**全数、国内変動、真正H6、独立空間foldを先に確認し、その後で独立した生物検定仕様を固定する。
