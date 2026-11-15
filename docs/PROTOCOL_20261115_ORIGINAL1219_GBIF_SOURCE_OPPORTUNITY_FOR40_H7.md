# 267単期EEA H7のH6候補40件：真正11/12 GBIF原本だけで二期間の記録機会を判定

**2026-11-15。#43723429521 GIS-onlyで12真正H6に40新規H7が存在すると確認した後の二次source-only監査。既存PR #103の生物学的H6内回帰STOPを解除しない。**

## 実証ソース

1. 真正HydroBASINS L6/L7幾何監査 `#43723429521` の `EEA_SINGLEERA_NEW_H7_IN_ORIGINAL_SINGLETON_H6_SOURCE_REGISTRY.csv`、40候補／12 singleton H6。元研究232 H7親の100%一致、Level6 source SHA256 `6a4d9bcbadacbf6e32936970f84ec77c21809948085d79d1d4388964d607fa`、Level7 source SHA256 `b9b810f4db271f0fa6b39fecf01784041b948944b3c3f296a393766388e86b44`。
2. **独立に過去取得された真正** `#43199340731` の `hybas7-gbif1219-full-original-source-shard-0..3`、合計元1,219 H7×2 original periods×3 taxa groups（Plantae, Lycaenidae, Animalia）。元WKT source digest `47f597e0a798cd6155ee9fb69ce8916ed5629951d046d114feb9f667963310e4`。元期間2017–19、2020–22、GBIF index original source SHA per requestを保持。元のデータはsource-only `limit=0` countであり蝶の個体数や現地植物被度ではない。
3. 元EEA2018/22とHydroBASINS geometryの400? **501 H7**から出した候補しか選ばない。自然・生物的な応答符号やHNS得点との関連でサンプルを選択しない。

## 事前固定する二期間のGBIF機会条件

各候補H7の両方の期間について、**original 3 taxon-group API source requestがすべて成功**し、その原本インデックス件数が両期間とも `Plantae≥50 / Lycaenidae≥10 / Animalia≥50` を満たす。この数値は以前の1219 H7 source macro preflightで採用された同一の観測機会閾値であり、今の研究結果の方向とは無関係。

- 1 H6に複数候補があるときは、ソース源で上記機会がある候補のうち、**事前source順位**（EEA dual-era > single-era、元EEA物理監視局数多い順、H7数字ID昇順）で最初の1件を選ぶ。ただし「候補が1つも当該二期間のsource機会を満たさなければ失敗」とする。全候補も保存し、予測係数で差し替えない。
- 30既存H6＋候補H6が10以上で40に達する場合は「**次のspeciesKey facet/EIV/HNS源の収集へ進むだけの機会あり**」。40未満なら直前の固定研究課題はsource opportunity STOPとして保持。
- 成功したoriginal GBIF index件数はspeciesKey採点済みではないため、NUTRI対応率≥40%、各期の≥10採点植物種、全5科のHNS数値対応率、元の国・EEA法・W predictorと整合するかは**未確認**。ここで生態的係数をFitしない。
- 結果が不利なら閾値を下げたり、source外のH7と混ぜたり、2018のみのEEA観測を2022のEEA観測に代入しない。最低40 H6という前のsource内局所FE研究条件は維持する。

**科学的境界:** この監査は観測機会・元GBIF APIレスポンスに関する実在性を調べるだけ。葉N、植物資源量、成虫吸塩、幼虫摂食、処理場の実測N/P・蝶の個体数、因果SEMを推定しない。元US EPT/欧州EQRは依然別の測定単位。main・chocho・親Draftには変更なし。
