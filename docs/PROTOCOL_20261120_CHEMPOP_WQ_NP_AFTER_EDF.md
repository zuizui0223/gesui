# 2026-11-20｜ChemPop水質3指標での英国EDF→科均等BMWPの追加予測反証

**既存EDF→BMWP負結果を見た後の探索的・観測的反証仕様。** 追加水質のBMWP係数・予測符号を開く前に、原2010–18の水質・生物の採集キー／coverage source-onlyを完走した。[原N/P source-only #44670909657](https://github.com/zuizui0223/gesui/actions/runs/44670909657) では5物質候補とも7地域に対応し、生物採集site×dateペアで5,116〜5,391、約594〜607siteの水質測定summaryを取得できた。実際に同一生物原試料で完全な複数水質を揃える交差集合は、**この方法書に固定後に初めて検査**する。

## 主目的と三つの固定化学種

元英国2010–18のsite×sampleを固定した同じデータ上で「モデル化された地点下水混入率EDFが、河川の周辺時期の栄養塩を条件づけても底生科のBMWP構成と関連するか」を検定する。

- P系主指標：**Orthophosphate（リン酸塩、総リンTPではない）**
- 酸化N系主指標：**Nitrate（硝酸塩、総窒素TNではない）**
- 還元N系主指標：**Ammoniacal nitrogen（アンモニア性窒素）**
- 元候補にあったAmmoniaとNitriteは、同じ窒素プールへの依存・共線性の可能性が高く、主モデルへ結果の符号を選んで追加しない。将来の別源事前設定感度とする。元原票や未数値スコアを0補完しない。

ChemPop v2 EIDC公式326 ZIPから固定原ファイル名（7地域×3種）とoriginal site BIO_SITE_ID、BIO_SAMPLE_DATEを使う。対象原source：元2010–18で得点科3以上、付与率70%以上、地点EDF/固定共変量ありの7,836 original site×analysis samples。水質表のWQ_COUNT>0かつ WQ_MEAN_LOQ2 が非負・有限。各化学種は元BIO_SITE_ID×BIO_SAMPLE_DATEに重複ない1行であることを確認し、候補1種に複数値なら最大値・平均など都合よい値を後付けしない。

**注意：水質は元生物採集日に対応した期間集計であり、必ずしも採集当日実測値ではない。** LOQ1/2/3 の置換定義・原単位が統一されていない場合は因果や絶対量解釈を停止。主モデルは元数値LOQ2のlog1p変換を標本内またはtrain側だけで標準化、もとの mg/L など数値単位を原出典確認前に推測しない。

## 新しい源集合の固定ゲート

同じ元source original sampleの**全3物質が同じsite/dateに揃う**ことを要求する。最低500 original source samples・150 distinct sites・4 native region blocks（各10以上の独立sites）・200 distinct site×yearを事前条件とする。さらに全chem3の値が非負・有限かつ各標本内で分散を持ち、完全背景+EDF+3chem共線性のVIFが15以下でなければ単一多変量モデルはSOURCE/MODEL STOP。元生物原採集キーと栄養元まとめ統計のjoinにsite/dateの偽近接や翌年proxyを使わない。水質欠損とEDF/都市/農地/regionの差は元sourceだけで報告し、アウトカムを見て対象源を補完しない。

同じsite×yearに複数original sampleがある場合は真の元site×date単位を平均してからsite×yearへ集約。元BMWP出現科平均は同一対応した生物sampleのみで再計算。単一サイトの連続観測を「独立処理施設」や「実際に下水投入が年次変化」と呼ばない。

## 固定4モデルと検定の方向

主応答は前に使った2016公開英国BMWP scored familyの**出現科等重みCWM**のみ。主EDF曝露は静的EDF_MEANのlog1pをtrainのみ標準化する。背景共変量は以前のChemPop WLS/Ridgeと同一のurban、arable、河幅、標高、河床siltclay、観測年/月、得点率/科数、region、分析法。

モデル（同一 source matched sample/site×yearにすべて当てる）：

- B0：元背景
- B_E：背景 + EDF
- B_N：背景 + log1p(Orthophosphate,Nitrate,AmmoniacalN) 3種を同時
- B_NE：背景 + EDF + 3種

WLSのsite-cluster記述95をEDFについてB_E/B_NEで比較し、統一source同一標本での係数を報告。高いVIFやsource品質STOPを符号の都合で緩めない。

**一次予測**：7地域のうち一地域を完全留保、共有siteがあれば学習から除外、Ridge10/train-only scaling/categorical vocabulary、site均等損失、固定site予測誤差の地域等重みMSEと固定7地域2,000回再標本化95。2つの条件付き増分を区別する：
- EDF after NP：B_N → B_NE
- NP after EDF：B_E → B_NE

さらに B0→B_E と B0→B_Nを報告し、方向や予測性を片方だけ取捨選択しない。7地域のいずれも最低10 true sitesかつ他6地域で200site以上を事前必須、未達は対応foldをsource STOP。

## 科学的な意味の上限

EDFを追加しても説明力が残る場合は、**栄養塩summaryだけで置き換わらないモデル推定下水・河川環境の空間情報**であって、特定の毒物や下水起源のN/P供給の因果効果ではない。EDFを追加しても予測改善が消える場合も、下水由来の栄養媒介が証明されたわけではない。B_NEは水質指標が環境変動の結果・媒介・交絡やコライダーたり得るため、パス係数を因果的媒介率に変換しない。

生物の実個体数重み・目内シャッフルの別検証は [先行階層解析](RESULTS_20261120_CHEMPOP_BMWP_ORDER_TRUECOUNT.md)を採用し、この3種比較の結果を見て指数の種類や応答標本を変えない。mainへ自動マージしない。
