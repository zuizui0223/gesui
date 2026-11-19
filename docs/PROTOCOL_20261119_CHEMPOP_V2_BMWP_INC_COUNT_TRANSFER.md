# ChemPop v2：実下水混入proxyと底生動物のBMWP形質構成を検証する事後凍結仕様

2026-11-19。水生EQR/US EPT/欧州植物蝶の既結果を知った後の**事後的探索解析仕様**であり、独立の事前登録ではない。元解析は未実行。原本source qualityを新BMWPスロープを開く前に判定する。

## 元データの対象・独立単位

- EIDC ChemPop v2 DOI 10.5285/23f194b0-1e21-4f14-be9c-6eabdcc1feb7 の326ファイルのうち、7地域（cp_an,cp_mi,cp_ne,cp_nw,cp_so,cp_sw,cp_th）の**taxonabundance原票**、**sitevariables原票**、対応する水質・河川流量集計表のみ使用。旧版座標修正前のデータを混ぜない。
- 一次時期：**2010–2018年に原生物試料採集日があるもの**。過年度1965–2009は別の年代外挿・確認用（観測方法が違う可能性があり、主母集団に後付けで戻さない）。
- 元観測最小単位：地域×BIO_SITE_ID×BIO_SAMPLE_ID×BIO_ANALYSIS_ID×BIO_SAMPLE_DATE。複数taxon行に重複する同じ科は「出現」では一科一回、個体数重みでは正式個体数が確認できた場合のみ足す。1地点に複数採集日があればsite-cluster、1地点×年の複数採集を独立siteと数えない。
- 7地域とBIO_SITE_IDのID重複を検査し、元監視区域が重なり同地点が複数地域に再収録されるならサイト独立性ゲートをSTOPにする。元下水proxy EDF_MEANは現地で実測したN/P濃度・排出質量でない。site固定の値を年ごとの栄養変動として擬似反復しない。

## 原公式サイト変数と曝露の品質条件（先に固定）

- 主要勾配：\`EDF_MEAN\`（元LF2000-WQX由来のモデル推定 effluent share proxy の候補）。現地物質負荷ではない。分析前に**高値が下水混入率増大を意味するか、単位が0–1 fractionか0–100 percentageか、原モデルが評価した基準年**を原資料から確認する。未確認ならWと呼んで回帰せずSOURCE STOP。
- 背景変数の元候補：\`LC_URBAN_PER\`, \`LC_ARABLE_PER\`, \`LC_WOODLAND_PER\`, \`BIO_WIDTH\`, \`BIO_ALTITUDE\`, \`BIO_SILT_CLAY\`, \`BIO_SAND\`、河川流量参照、region、月。ソースの観測可能率/相関/単位をOUTCOME以前に開き、posthocに効いた変数だけを選ばない。
- 真の栄養塩を扱うときの各WQ original sample-dateは \`BIO_SITE_ID\` と \`BIO_SAMPLE_DATE\` を完全一致し、同一水質項目にLOQ1–3がある場合「複数手法の異なるLODルール」のまま比較。流量/水質期間集計と生物の単一採集日は完全同時観測でないため実測同日曝露と呼ばない。
- 個体数列 \`BIO_ABUNDANCE\` の生物学的単位を元metadataで確定。分類群の目/科/種グループ混合分類による重複計上を防ぐ。元値がlog classまたはカテゴリーのとき「raw individuals weighted CWM」と表現しない。

## 植物やUS EPTとの比較のための二つの応答

- 主応答：現地出現した**科**について、厳密文字列一致した第三者公開 UK biotic RData family \`BMWP\` スコアの**等重み平均**。原BMWPは通常10に近い方がより汚染に敏感な科。科得点辞書は2016公開旧版であることを必ず明記し、現在のEA 2026公式WHPT版スコアとは断定しない。
- 副応答：\`BIO_ABUNDANCE\` が真正数量の計数値として検証できた場合だけ、同じ得点可能科について出現計数を重みとした平均BMWP。出現重みの方法と完全に同じ適格標本だけで比較し、差（count−incidence）を報告。
- 絶対種多様性/現地絶滅、個体群動態、羽化昆虫・サービス量ではない。BMWP自体が水質分類目的で作られた先験的耐性スコアであり、現在の汚染下で独立に測定した生理学的耐性ではない。

## 先験的 source-only 支持条件（既存陽性を見た後の探索的設定）

- 原site/analysis/sample/date×分類群の参加元ID衝突0、各sample内での分類名重複が科へ正しく集約される。
- 適格sampleのうち、得点可能な出現科数 >=3、全出現科の得点付与割合 >=70%、site \`EDF_MEAN\` と urban/arable/width covariateがsource-valid、元採集年2010–18。
- 対象適格sample >=500、異なる真のsite >=200、少なくとも5地域でそれぞれ10以上siteかつEDFに5以上の異なる値。過年度データで不足を救済しない。
- 実測N/Pの同試料コホートのsource coverageが事前の200 sitesを下回る場合、EDF×水質×生物のSEM連鎖はSOURCE STOP。誤記値を0に補完しない。
- 先験的なスコア関係・空間スケールの欠測率と分布を署名済みSOURCE出力へ記録した後にのみ、次の生物回帰・帰無を開く。

## 上の全条件が通過した場合の検定順

1. **条件付き単層形質対応**：生物 BMWP incidence-CWM ~ +1 train-source SD EDF （EDF原スケールを事前に確定）+ land cover / river size / riverbed / month / region。重みはsite均等、未知siteを事前5fold spatial blockで留保。一次係数・siteクラスタ記述95、BH補正は事後的な多重比較を明記。
2. **出現組成か個体数頻度か**：全く同じSOURCE適格sampleでcount-CWMとincidence-CWMを比較。個体数重み付け差の区間を別に算出し、非有意だから個体数変化なしとは言わない。
3. **目内label帰無**：原科×目の正式分類対応を維持、科BMWPスコアだけ目内で999〜4,999回置換し、地理・source曝露・観測デザインを固定。EW／score欠測は変えない。これは完全な系統樹を用いたphylogenetic permutationではない。
4. **地域外予測**：少なくとも5原地域で別一地域を丸ごとholdoutし、該当地域に重複する元SITE_IDを学習から厳密に除外。B＝元地域/土地利用/河川物理/時期/effort、B+EDFのMSE差を地域を等しく数えて評価。別年代2010–18→過年度への外挿は別コホートで方法差が大きく主の代わりにしない。
5. **栄養の媒介性を反証**：正式な流域同日WQのTP/TNとEDFが同じサイト・年で得られた場合のみ、土木・農業・河川と採集方法を同時調整して +EDF after measured NP、逆 +NP after EDF のモデル比較。場所レベル共変動から下水起源P/Nの因果媒介を結論しない。

## 予測される科学的結論の分岐

- incidence/ラベル帰無/新地域外挿がすべて支持：英国の**モデル推定下水曝露地理と科レベル耐性形質の頑健な空間対応**。水陸共通の非因果的環境フィルタリング比較へ拡張。
- 係数だけ正／負、ラベル帰無・外挿不確実：過去のUS実測TP→EPT PTVと同様、**地理・監視方法・種プールの見かけの規則**として限定。
- このsourceゲートを一つでも通過せず：SOURCE STOP、空欄のまま、EQRを耐性形質へ解釈しない。

既存のEU底生EQR/US EPT/欧州植物HNSの肯定・否定・STOPを温存し、研究原稿やgesui/mainは変更しない。
