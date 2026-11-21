# 2026-11-21｜英国2018年の真正UWWTD放流地点とEA最終処理水TP観測・ChemPop底生調査を方向性で源照合

**SOURCE-only／事後探索。下水による底生BMWP変化やP/N流出量の因果効果ではない。**

## 真の原入力と再構築

- [HydroRIVERS/UWWTD/EA/ChemPop原SOURCE結合 Actions #44848946633](https://github.com/zuizui0223/gesui/actions/runs/44848946633) SUCCESS、original aggregate JSON artifact **14684837094**。
- UK5 EEA UWWTD v07r01 true 2018: 元1,799下水処理場／1,809原放流点／**1,498淡水FW放流点の真正座標**。施設ごとの放流点が唯一という条件、各元scientific member SHAを厳格確認。UK5の実測P・N元排出量は正値0件（欠測とゼロ排出を混同しない）。
- 2026年現在EA WQE Water Company FINAL/TREATED SAの元登録5,722採水地点。2018年英国FW原放流地点との距離250m以内で、UWWTD側から1 SA、SA側から1 UWWTD放流点、元施設も1放流点の**相互一意な位置候補930地点**。100mでは683、500mでは1,026（こちらは感度のみ）。この候補は元行政permit native ID一致を実証したものではない。
- 930候補を2018年EA WQE本当のTP0348/FINAL SEWAGE EFFLUENT/SA/mg/L数値sourceで確認し**360元SA採水地点に2018年TP実測あり**。そのうち歴史的EEA報告**放流地点そのもの**がHydroRIVERS500m以内の河道へスナップしたのは292。元BIO2018は366 source-qualified site（354が500m以内に河道スナップ）。
- HydroRIVERS v1 EU official SHA 500da7d36ceee0aa4c82dd625f5c53d80ba2f7ac14ad022cd3a6601bf55068e5、2018年の元UWWTD放流地点からBIOまで原NEXT_DOWN reach距離50km以内の上流候補・逆流向placebo・同reach順序不明をすべて出す。原放流点座標とSA採水位置を入れ替えない。

## 元の2018年だけの結果（河川網の方向を使用）

| 元SOURCE監査（2018年） | 件数 |
|---|---:|
| 上流側UWWTD放流地点×下流側BIO地点のreachペア | 182 |
| **上流側の歴史的放流点候補を持つ異なるBIO調査地点** | **93** |
| 上流側の原UWWTD放流地点で同年TP測定SAと位置対応したsource | 112 |
| 逆方向：BIOの下流側にUWWTD放流点がある対照BIO地点 | 124 |
| 同一reachでUWWTD放流点とBIOの前後順序不明の地点 | 52 |
| 元BIO調査日と真のEA SA TP採水日があり、**BIO調査前90日以内にTP測定あり**のBIO地点 | **92** |
| BIO元site×date192組のうちTP前90日 | **188** |
| BIO元site×date192組のうちTP後90日 | 190 |
| TP測定が同暦日で時刻順が不明のBIO日付 | 16 |

原2018 BIO source地点数は366。同じ年だけなので上流候補93地点＝93地点×年であり、前固定の50 distinct sites（通過）かつ100 site×year（**未通過**）を同時に満たさない。TP前90日でも92地点であり基準は通らない。特にBIO調査後90日にもほぼ同数のTP測定があり、監視頻度だけでは処理水Pの生物への影響を示せない。

## 主張の上限と次段のSTOP

以前の全国のEA処理水**採水地点**をHydroRIVERS上流と見なした2018年のBIO地点は119。元2018年UWWTD歴史的な**放流地点**に位置候補が相互一意で、2018年TP実測がある条件へ絞り直すと93地点。分母・識別条件が違うため119→93の差を生物学的な負効果や脱落率と読まない。

今回の真正sourceで通ったのは歴史的放流地点座標・2018EA原TP採水地点の位置候補・原2018TP実測・河川網上の幾何的方向・真のBIO_SAMPLE_DATEだけ。**permit IDで同一施設と認証したわけではなく、2018英国UWWTDのN/P元実測質量は欠測、EA TP mg/Lをkg/yearに変換する流量も欠測**。ほかの事業所や農業由来栄養、河道工学、monitoring site選択の交絡も未制御。下水源TPの実際の受水域到達を示す化学検証、原因別物質収支、BMI/FD/生物量への因果経路はSOURCE STOP。

英国PRTRの都市下水処理場年度N/P報告にも同一施設・同年の両物質値がそろうケースは0（2017年0、2018年0）。施設側の緯度経度とUWWTD歴史的アウトレット、EA SA現行採水点の3者間native facility IDを推測でつないではいけない。

関連方法：[original 2018 source strict 250m reciprocal ID geo protocol](PROTOCOL_20261121_UWWTD2018_OUTLET_TO_EA_TP_TO_BIO_DIRECTION.md)、[UK UWWTD元放流地点データ監査](RESULTS_20261121_UK_UWWTD2018_FW_OUTLETS_TO_EA_SA.md)、[英国 PRTR source mass](RESULTS_20261121_UKPRTR_TRUE_2017_18_WWTW_PN_SOURCE.md)。公開GitHubには個別元化学値／座標／施設ID／生物応答を転載していない。PR #137 Draftでsource QAを固定し、main/PR135/PR136を変更しない。
