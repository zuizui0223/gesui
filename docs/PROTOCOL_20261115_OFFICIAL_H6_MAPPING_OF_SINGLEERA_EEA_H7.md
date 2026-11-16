# 既存EEA単期H7を真正H6に割り当てるsource-only機会検査

**2026-11-15追補。真正EEA二期間H7が234、元研究232、追加可能2しかないことが #43723099758 で判明した後の別設計。前の二期間検定はSTOPのまま。**

## データ源

- #43165326244：EEA18/EEA22/EPA1819 original H7 station-geography source、6330 original rows、SHA256 `a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89`。対象は真のEEA18/22に限り、EPAをEU流域に接続しない。
- 真正原452? **元454** H7×期間、232 H7、90 native H6のパネル SHA256 `1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca`。source codeでは元period/H7/H6のみを使用する。原HNS,NUTRI数値やW数値は開かない。
- v2 original EIV品質採択結果451 H7×期間、元219 distinct H7 across two-period,86 native H6、うち56 H6 singleton。元Source-v2 #43718069015は175 named hybridを未得点として保持し、NUTRI scoresは1つも変更していない。
- HydroBASINS Europe official L7 `hybas_eu_lev07_v1c.zip` の元SHA256 `b9b810f4db271f0fa6b39fecf01784041b948944b3c3f296a393766388e86b44`、および公式Level6 source ZIPを取得。元の固定232 H7→H6結果と、最大面積重複≥95%・次点との差≥50%を1件残らず再一致検証する。公式L6 invalid geometryはmake_validするが、変更件数と元ZIP digestsを記録。PFAF文字列/H3/近傍距離で親IDを補うことは禁止。

## 新規H7の探索とソースゲート

- 2018のみEEA源H7=172、2022のみ=95（前の原本source gate）。本来の2期間ソースに入れ替えず、**単期EEA源が存在する別候補**と明記する。
- 本来の56 singleton H6に、元232にないH7が1つ以上重なるかを幾何でマッピング。重複H7、年違い、基準アウトカムでの選択禁止。
- 幾何のみの前段候補：異なるsingleton H6のうち**最低10 H6**で追加H7が1つ以上存在し、地理割当が重なり95%以上、既知EEA局の有効single-year presenceがある場合、**新規ソース取得を検討する機会**あり。10に届かなければ水域monitor registry由来の拡張ルートはSTOP。
- 候補を選ぶ場合は、1親H6に対し二期間EEA源を最優先、次いで元EEA物理局数の多い候補、H7数値IDの昇順で決定的に上位1件を選択し、候補となったすべてのIDとsource countsを保存する。結果を見て差し替えない。
- 一期EEA局のH7を取得しても、そこにGBIF植物・蝶両方のsource dataが2期間であるとは未検証。GBIF facet（2017–19/2020–22）×EIV 40% coverage×蝶数値host trait coverage、元W地理共変量、真のH6親、国と評価法を再取得できない限り、**研究のH6>=40合格とは扱わない**。単期EEAから水域の2期EQRを捏造しない。
- 元のH7期間2期を独立空間単位と数えない。別の新sourceコホートになるので、元の454観測結果と混ぜて既存W効果を流用しない。

## 範囲

今回のGIS source-only実行は**追加候補H7の地理的存在と監視局の有無**の確認まで。栄養嗜好NUTRI、蝶HNS、EQR、排水処理水フロー、回帰係数、生態系サービス効果、相互作用の解釈はいずれも新規推定しない。原本404? 454 cohortを保存する。結論不利でも出力を残す。main/chocho/PR #103は変更しない。
