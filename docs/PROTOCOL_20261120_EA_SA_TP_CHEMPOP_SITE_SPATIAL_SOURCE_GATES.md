# 2026-11-20｜EA最終処理済み放流水TP採水地点→ChemPop底生調査地点のSOURCE-only地理照合

**事後探索の原source品質監査であり、処理場→河川→群集の因果解析ではない。** [2017/18Bath周辺SA地点監査](https://github.com/zuizui0223/gesui/actions/runs/44676173833) により、2017年46、2018年49のWater Company FINAL/TREATED EFFLUENT元地点でPhosphorus Total as P (0348) のmg/L観測が得られ、46地点が両年に存在。Bath50kmは英国内局所pilotであり全国推定ではない。

## 元資料と正しい曝露の意味

- EA新WQE v1から、同じ2017/18年Bath50kmの元水質コード0348、正式sample point type SA（water company final treated）、source sample material FINAL SEWAGE EFFLUENT、measurement unit MILLIGRAM PER LITREに限る。storm overflow SB、not-water-company UA、raw sewage RA、河川・湖沼を混ぜない。
- 補助的にsamplingPurpose WATER QUALITY UWWTD MONITORING DATAのみに限定した別集計を出す。source SA UWWTDと全SAを同じ物と主張せず、該当範囲のカバレッジを報告。
- ChemPop v2 EIDC元ZIP内の7地域分SITE matched location original CSVには原BIO_SITE_ID、BIO_EASTING、BIO_NORTHINGなどを含む。2010–18の元EA Benthic taxon・BMWP辞書とは別のサイト地理sourceである。水文/河川により放流が届いたという情報は含まない。
- 両original sourceの座標は英国のBNG Easting/Northing（m）として検査。真の座標変換・datum検証を伴わない地点の投影違い／欠測／(0,0)は source STOP。公開観測点の識別子・住所・座標の行単位転載は行わず、元source SHAと集約件数だけを公Gitへ保存。

## 研究対象母集団とソース品質

- 元ChemPop7地域1,519 BIO_SITE_ID（重複地域キーの同一位置矛盾検査）、もともとの同一データ系で2010–18 BMWP score-eligible 7,836原生物analysis/707地点をsource-onlyで再構築する。
- 元実際の観測採集日が2017または2018の適格サイト群を**予め別対象母集団**として数える。2010–16生物群集を2017/18下水の値で置き換えない。観測年2017/18と源の観測年が一致するものだけを本pilotの地理的機会に含める。
- SA点の座標が同じnative notationで異なる場合、差の件数と最大値を記録し、1地点の移動が100 mを超えるなら原地点結合をSTOP。複数地域/sourceで同じBIO_SITE_IDと異なる座標があるなら重複ゲートをSTOP。
- **距離は地理的 SCREEN に限る。** Bath50kmのsource STW地点と英国ChemPop original bio-monitoring siteをBNG直線の最近傍で比較。0.5/1/2/5/10 kmを計画固定し、各閾値で候補STW点数、候補ChemPopサイト数を、全1,519地点と同年2017/18 siteの双方について集計。両年観測の処理場別は別パネル。
- 最短距離だけで下流・上流を区別しない。候補siteが処理場の上流だったり支流側にあれば曝露は未定義。近い複数処理場のN/P濃度を地理的平均・合算して「源負荷」としない。

## 次の上流ルート判定に進む条件

源STWで2年にまたがる正式SA実測TPが検証されること（今回46 native pointsでsource opportunity成立）。真の同一年度ChemPop siteの地理的近接候補がある場合でも、**元流向を持つ河川ネットワークで source STW→receiving channel→downstream bio site の有向接続が確認されるまでモデルSTOP**。方法書の主研究規模条件（少なくとも30STW、50同年下流ChemPop site、100site-year、5独立hydrologic blocks）は位置だけで支持とみなさない。放流TP原液中濃度mg/Lを河川TP濃度や下水N/P質量にしない。

陰性・欠測も全て保存し、source positiveだけを選び直さない。主ChemPop観測結果（EDFに3河川化学種を加えるとEDF追加予測不支持）と混同しない。元EARaw/Wessex原ファイルのHTTP403と英国EEA2018両物質質量欠測は未解決である。
