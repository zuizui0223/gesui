# 2026-11-19｜2016底生動物・水域単独点源N対陸域N投入：数値モデル固定（事後探索）

前段の2016年 SOURCE-onlyゲートが、73独立NUTS2、12か国、国内の点源N密度/大気N沈着の変動8か国、元5fold地域数15/16/14/13/15でPASSしたため、水域側の新たな比較を一つだけ実行できる。事前の2018×2017–19水陸共通2×2は25地域/6国でSOURCE-STOPしたまま。本解析でそのSTOPを取消さない。

**これは2016年の点源NとEQRに関する既存肯定結果を知った後の追加探索的調整であり、独立確認的な新データ検証ではない。** 2016 BQEの元の出現・位置・EQR妥当性・NUTS2 foldは既存原本のままで、標本やhyperparameterを新しい追加結果の符号で変更しない。

## 観測・尺度と出所

- Y=公式EEA2016年の底生無脊椎動物の元正規化EQRを、先行研究の元2016 SOURCE固定・独立検証済みNUTS2地域ごとに平均した**地域平均EQR**。元報告 [2016_EEA_RIVER_BENTHIC_EQR_JRC_GREEN_MODEL_POINTNP_NUTS2_HELDOUT_RESULTS.json] にある全73地域の原ソース数/同地域のEQR平均だけを使用する。元サイトのEQRを新たに取得したり1地点を独立NUTS2として疑似反復しない。地域平均EQRは有効な状態指数であり、水生昆虫の耐性形質・生物量ではない。
- P=元source GREENモデル2014–18平均の家庭＋産業由来点源N[t] / 元2016 GISCO NUTS2面積km²、log1p。2016年の個別WWTP実測量でも受水局のN濃度でもない。生態モデル結果元表の pointN density と独立 source-only台帳の原質量が同一か検証。
- D=EMEP2015大気N沈着とDAYCENT2010–19鉱質N施肥・有機N施肥のNUTS2平均を個別log1p。肥料は農地での平均投入、河川への拡散負荷ではない。元SOURCE-only 269地域のうち2016底生SOURCEで重なる73地域に限り欠測を0へ変換しない。生物固定Nは副次で、主ブロックには入れない。
- 単位＝**73独立NUTS2**、12か国。元SOURCEゲートで国≥8、NUTS2≥55、源内変動5国以上、5fold各≥8地域をすべてPASS済。これを変更しない。

## 凍結4モデルと予測評価

主背景 B0: GISCO NUTS2地域の重心緯度・経度、log1p面積、log1p(元BQEサイト数)、国固定効果。国別EQR評価法の細かい違いは原本の地域平均に残る。全モデルで同じ観測地域・同じ座標/国・同じ地域5foldを用い、1地域が学習と予測の双方に登場することを許さない。

- BP = B0 + P（点源N）。
- BD = B0 + D（大気沈着・鉱質/有機N）。
- BPD = B0 + P + D。

標準化と国ダミー辞書は**学習データだけ**で決定し、未知国は全ゼロdummyとする。Ridge α=10固定、切片非罰則、内側hyperparameter探索なし。元2016 NUTS2の事前固定5fold（元のsource manifest）をそのまま使用。各モデルのOOF予測と地域ごとの二乗誤差を保存。

主指標は地域等重みと**国等重みのheld-out MSE**（国ごとNUTS2等重み→国等重み）。副次としてサイト数等重みを表示するが、独立標本数は73。BPDを比較する固有追加改善は U_P=100(MSE_BD−MSE_BPD)/MSE_BD、U_D=100(MSE_BP−MSE_BPD)/MSE_BP、C_W=U_P−U_D。予測改善の符号は生態学的EQRの増減方向ではない。点源の標準化係数は探索的補助で表示し、因果N効果とは書かない。

国全体を復元再標本化する固定OOF損失ブートstrap3999回（seed=20261119）、各UとC_Wの95%記述区間を報告。固定予測誤差の再標本化でありモデル再学習の不確実性まで含まない。国単位leave-one-country-outは**予測結果の事後的な再集計**であって再学習ではなく、各国・各foldの損失を完全掲載する。

**水域で点源が相対的に優位とする強い判定**：国等重みの U_P≥1%かつ記述95%下限>0、かつC_W>0で記述95%下限>0、さらに点源Nを背景+陸N投入から予測するVIF<5。いずれかに不足したら点源が主因という結論は支持しない。陸N投入側が相対優位の場合も同じ計算と方向を全報告するが、当該源は河川への到達量ではないので「陸から水へ流れた因果効果」と言わない。全てが弱い場合は「追加情報が乏しい」とする。

## 生態解釈の上限

地域平均EQRと地域平均モデル源の関連は都市化・国別評価手法・農地/気候・水温・河川流量などの未調整交絡を含む。モデルは複合的な地域地理の説明力のテストであり、処理場の下水Nと蝶の寄主Nを結ぶ媒介・生態系サービス・負荷源の原因効果ではない。2016年モデル値と2018/2022年の実WWTP源報告も単純な独立再現年として統合しない。国/面積調整後のpoint-vs-land source相関の知見と整合的でも、水陸共通2×2のSTOPは解除しない。


## Source-vs-EQR eligibility implementation correction (after initial technical failure; before any new fit)

Initial code [GitHub Actions #44436319400](https://github.com/zuizui0223/gesui/actions/runs/44436319400) stopped before a model fit because it incorrectly required the **73 source-only candidate regions to equal 73 old numeric EQR-valid regions by identity**. Both totals equal 73 but they are not the same regions: the original historical numerical EQR summary has 6 Norwegian regions absent from JRC soil-N source, whereas the source-only JRC cohort has 5 Danish regions and Austria AT13 absent from historical EQR-valid regions. Therefore the true numerical outcome+source overlap is **67 regions in 11 countries**. This was discovered from preexisting **region ID and availability metadata**, not signs of the new four-model coefficients, and no four-model numeric effect was fitted in the failed run.

To honor the already fixed **minimum** conditions (>=55 region, >=8 countries, >=5 countries with 3 regions and both source variations, >=8 region per original fold), an explicit **pre-new-model source/validity compatibility audit** was added: analysis/20261119_point_land_N_source/eea2016_eqr_source_compatibility_gate.py. It will save the full set of eligible region IDs and reasons for exclusions *before the new EQR four-model fit*. Only if this audit confirms all predeclared minima will the new four-model script accept its exact region-ID list. No changes to period, metric, training fold, source models, coefficient thresholds or numerical model class are permitted. The original failed run remains available in the PR audit trail. This corrective stage is after previous numeric 2016 data were published and therefore remains exploratory, not prospectively preregistered.
