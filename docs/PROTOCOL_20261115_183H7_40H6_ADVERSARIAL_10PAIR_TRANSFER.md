# 183 H7・40 H6：植物NUTRI↔蝶HNS地域内関係の反証検定（元結果後に固定）

**2026-11-15 / 事後 adversarial sensitivity。PR #104でH6 FE β=+0.1181、9,999 within-H6 residual-null p=0.0001 を既に観測した後の追加設計であり、独立した事前登録ではない。元の前段のSTOP、時間差分のゼロ、元の二期間水質データの制約は変えない。**

## 同一の原本と比較母集団

- 183 **異なる**本物のHydroBASINS H7、40真正H6に属する、元の2期間GBIF植物・蝶出現の均等平均1 H7 1行。
- 既存H7は173、source-onlyで後から追加した10 H7は10異なるsingletonだった親H6に各1つ。既存173 H7のうち、163は元から複数H7の30 H6、10は新設10 H6の旧H7。この分解を**コホート事実として検証**し、1 H7を2期に分けて地理的擬似反復にしない。
- source panel：`GENUINE_TERRESTRIAL_ONLY_183_H7_40_H6_PAIRED_H7_MEANS.csv` SHA256 `c0840b1c5c6e85342fea41d57dc40fa5ec295fabf6bcd6908a4072557db9a736`、原source-only QA schema `GESUI_ORIGINAL_NEW_183H7_40H6_TERRESTRIAL_SOURCE_COHORT_AND_TRAIT_REPLAY_V1`。前段の元454行/EIV/HNS復元誤差0を確認。
- 自然史の意味は、*植物NUTRI/EIV好栄養環境嗜好*と*蝶の既知幼虫寄主の固定HNS/生育地栄養嗜好*というGBIFで記録された形質構成の対応。現地栄養流量・葉N・幼虫摂食・蝶個体数ではない。
- 追加H7には2期のEEA水質観測が揃わない。ここでの解析は陸域出現記録2期間であり水陸の現地連鎖とは呼ばない。

## A. 同じ母集団の感度解析（絶対に母集団の変更を結果として使わない）

1. **全40 H6/183 H7：** H6 FEと13個の元固定source controls（W,crop,urban,naturalout,area,nPlant,nButterfly,nAnimal,plant_coverage,hns_coverage,hns_host_completeness,plant_scored_species,hns_scored_species）で主係数β +0.118076を再計算（fail fast）。
2. **元からの複数H7がある30 H6/163 H7：** 同じ13 controlsとH6 FEで植物NUTRI→蝶の**種均等HNS**係数を再計算、H6 cluster95と元H6内9999 Freedman–Lane残差置換p（種均等のみ）。これは**事後の30 H6解析**であって前段の40 H6ゲートを緩和した主要検定ではない。
3. 全183 H7における頑健性：事前列順を維持し、source controls全13、環境3つ `W,crop,urban` のみ、元H6 FEのみの**3本すべて**を並列比較。偶然都合のよい一つを主にしない。H6クラスタ95のみで新しい独立確認と主張しない。
4. 全183 H7の**元H6均等重み**と各H7均等重みの2種類を比較。新設10組にH6均等weightが過剰の可能性、source H6 size=2〜9による寄与を検査。
5. 元H6 FEにより除去された植物Nの変動に対して、残りのcontrolsが残差化でどれだけ削るかを加重VIF（H6-only残差分散÷H6+他controls残差分散）で報告。VIF>10またはデザイン条件数>100は警告、推定係数を原因と解釈しない。元条件数≈125は既知。

## B. 追加10真正H6を **未知地域ペア予測** として留保する検定

前段のsource-only選定結果を使い、10個の新しいH6ごとに「以前からのH7 1つ＋新たなH7 1つ」という**10ペア**が実際にあるか確認し、違えばstopする。

- **学習のみ：** 元から複数H7の30真正H6/163 H7。モデルは13 fixed source controls＋植物NのH6 FE回帰。train-only標準化。基準モデルは植物Nを抜くが、他13を保持。
- **テストのみ：** 新しい10 H6の各親で2つのH7の`new−old`差分。未知H6の固定切片はペア内差分で相殺し、学習時の係数を固定したまま`ΔY = Σγ ΔX + β_N ΔPlantN`を予測する。N係数や正負を見てペアを選び直さない。
- 10ペアのMSEを等重みにして`gain = 100(MSE_withoutN−MSE_withN)/MSE_withoutN`を主要予測値とする。10ペアに固定誤差2000回再標本化の記述95、負/ゼロ/正を全て残す。**10ペアしかないため不確実性は大きく、独立地域一般化を単純に主張しない。**
- 離散チェック：同じ10ペアの`ΔPlantN`、`ΔButterflyHNS`のPearson、Spearmanと、source W/effortのみの基準予測残差との相関（数値診断）。偶然相関を前段の公式同H6 permutation pと誤認しない。
- 対応する蝶の**記録重みHNS**は副次感度、主の種均等と一括平均しない。
- 先行の40 H6モデルはH6内**in-sample** SSE改善であり、今回の10 pair 外部予測と異なる。計算を混同しない。

## 成否判定・解釈

- A/Bとも正→同じH6内部の記録構成関連が古い30 H6に存在し、新しい10 H6の組内差の方向への転移も探索的に支持。
- Aのみ正、B不確実/負→新H6への転移確認なし、40 H6全体のin-sample効果だけを強調しない。
- 30 H6では0方向→新設10組の寄与や広域 source-variation 依存を疑い、全40 H6に基づく普遍性の主張を下げる。
- いずれも`W→Plant→Butterfly`因果媒介、実測局所WWTP投下窒素、真の植食者増加/バイオマス/サービスの証拠ではない。時系列では同H7のΔPlantN→ΔHNS同期がないという元結果を維持。

**操作：** main, chocho, PR #104に変更・マージなし。新Draftブランチで元source artifactを参照。元sourceと選定は一切書き換えない。
