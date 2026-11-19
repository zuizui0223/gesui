# 2026-11-19｜元H7植物の流域とJRC流域窒素負荷源：真のPolygon被覆PASS後、解像度不足を特定

**結果の区別：初回のH7被覆・重複検証は12/12 PASS。しかしJRC source流域の解像度を別途厳しく点検すると2/12しかPASSせず、H7へ負荷を局所・独立な曝露として配ることはできない。** 生態的な植物NUTRIや蝶HNS、水域EQRを開いた新しい解析ではない。

## 真正source-only対象と先に固定した被覆ゲート

- 元2017–19の227 H7中、原451/454のEIV品質可否を用いた適格H7 **226**。新しいSOURCEは、そのIDをsha256で固定選択した先頭12 H7のみ。生物の指標値や近年の相関結果は選択に使っていない。
- original 2,808 GBIF facetに保存された真のWKTとcomponent SHAを復元。JRCの native point BasinID、river BasinID/soil HydroIDの全元属性9,739対応候補、元30/30一致Polygon、ArcGIS公式layer5の元PolygonからJRC候補とH7の実intersectionを計算。重心最近傍やNUTS2割り当てに戻らない。
- [元の固定pilot仕様](JRC_NATIVE_H7_PILOT_SOURCE_GEOGRAPHY_FROZEN_20261119.md)では、全H7ポリゴンに対しJRC native-sourceポリゴンunionで面積80%以上・重複ファクター1.05以下となるH7がpilot12中8以上、真正H6が4以上を条件とした。
- [初回の元pilot Run #44442417933](https://github.com/zuizui0223/gesui/actions/runs/44442417933)は**12/12 H7被覆PASS、重複は検出せず、真正H6が11**で、最初のSOURCE地理取得ゲート自体は通過。これは生態効果のGO判定ではない。

## 初回の100%被覆に隠れた空間粒度の問題

元pilotを見て、面積の大きいJRC native 流域がH7を丸ごと含むことで100%被覆になる問題があると分かった。そこで**結果を見た後の追加厳格診断**として、対応するnative polygonの公表元面積とH7真面積の比を連続的に出し、元native面積がH7の1×、2×、4×、10×以下であるもの**だけ**でH7面積を80%以上被覆できるか、同じ12 H7に再適用した。これは結果を見てから追加したため独立なpreregistrationではない。元の被覆PASSは取り消さず、空間解像度の別の否定材料として報告する。

[追加解像度 Run #44442798020 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/44442798020)／[元12 H7逐一source overlap/NativeIDと原sha artifact](https://github.com/zuizui0223/gesui/actions/runs/44442798020/artifacts/14492061011)。

| source-only証拠 | 観測 |
|---|---|
| 原12 H7の元ポリゴン被覆≥80% | **12/12** |
| native source面積がH7面積の4倍以下のものだけで被覆≥80% | **2/12** |
| 追加診断で決めた最低8 H7 | **不合格** |
| 優占JRC source BasinIDのdistinct数 | **11/12** |
| 主な被覆source Basinの面積/H7面積 | 最小 **0.791**、中央値（12流域の2中央の平均）**約18.49**、最大 **159.40倍** |
| 元植物NUTRI／蝶HNS／EQRの新しい数値分析 | **ゼロ** |

具体例として、H7 2070549930（面積1,352.4 km²）とH7 2070548410（面積2,848.8 km²）は、**同じ優占ネイティブBasinID 129487**に含まれている。その元JRCのregional pointN/soilNを二つのH7にコピーすれば、独立な源の複数反復として数える誤りになる。他のJRC area比が50倍以上や159倍というH7もあり、native候補の総数9,739に反して現在のH7位置では実際の露出変動の空間解像度が粗いケースが多い。

## 適切な解釈

1. JRC自体の9,739 SOURCE候補を比較した点源密度対陸N投入のlog-Pearson相関は大気沈着+0.089、鉱質−0.004、有機−0.017、生物固定+0.044という**記述的地理相関**であり、この数字は元のJRC native単位におけるsource geoを表す。国・生物の相関や介入効果はない。
2. 元H7の生物で「水は点源、陸は面源」と言うには、同じH7内の**独立した地理単位・真の曝露**に点源/陸源のsourceを割り当てる必要がある。大きなJRC nativeソースを複数H7に複製してもその独立性は得られない。
3. 本pilotはSOURCEの取得上限を可視化したものであり、元NUTS2比較（25地域/6国）や2016年水域単独の不確実な点源優位を後から肯定結果へ作り変えない。植物NUTRI/HNSやEQRを使った追加検定は未実施。
4. **次の新規独立研究での優先ソース**は、JRC/EMEPの実0.1°格子の大気N沈着およびDAYCENTの1km農業N投入を真H7へ面積で再集計し、水への処理場報告N/Pは元放流口の水文上流域・流量で結ぶこと。大気沈着の元sourceとは年や解像度・単位が異なり、現地植物葉Nの実曝露を直接測ってはいない。新しい生物模型はSOURCE-only適格標本を別途固定したときに限る。

このsource-only結果、元ポリゴンの一致、元API/元sourceSHA、初回成功、後追いの厳格解像度診断と不合格をすべて独立Draft PR #127へ保持。main・草稿PR #124・旧#125/#126はいずれも未変更。
