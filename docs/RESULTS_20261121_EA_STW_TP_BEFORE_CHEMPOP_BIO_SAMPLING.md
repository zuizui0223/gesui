# 2026-11-21｜処理済み下水TP採水が底生動物の調査より先だったか：同一原採集日・上流河道のSOURCE監査

**SOURCE-only再解析。環境影響の回帰や因果推論ではない。**
[元採水日照合の実行済みActions #44838235844](https://github.com/zuizui0223/gesui/actions/runs/44838235844) 成功、成果JSON artifact ID **14678789421**。原EA Water Quality Explorer 2026 API、EIDC ChemPop v2 official 326 original files、UK 2016 BMWP採点辞書のSOURCE品質ゲート、HydroRIVERS v1 EU original source SHA256 **500da7d36ceee0aa4c82dd625f5c53d80ba2f7ac14ad022cd3a6601bf55068e5**を使用。生物BMWP値・実個体数・栄養による生物効果係数は計算していない。

## SOURCEと年月日の原本ゲート

- Water Company FINAL/TREATED sewage effluent (native SA) の実測Total Phosphorus as P、code 0348、**単位mg/L**、原sample material FINAL SEWAGE EFFLUENT、numeric非負値。年2017/2018。最終処理水の**採水場所**であり、施設の正式な物理放流口とは限らない。年間総リン kg/yrの測定ではない。
- 元WQEは2026現行登録5,722 SAを原英国BNG座標で確認し、HydroRIVERS v1 NEXT_DOWN 50km上流候補・逆方向placebo・同reach順序不明を**最初**に選択。その後元TP実際の採水があった2017/18のSAを確認。元同年TP監視SA点 **2017年281/2018年293**、そのうち275地点が両年。年次の真正SOURCE測定 **8,753件**。2026現在のSA登録は2017/18年の施設稼働・処理場本体数と同値でない。
- 対応するEIDC ChemPop v2 original BIO 2010–18品質合格 **7,836 source samples／707 BIO地点**。2017年389 BIO地点、2018年366、両年和集合463。今回のHydroRIVERS upstream reach条件に真の元SOURCE測定年が対応するBIOは**2017年118地点／2018年119地点、和集合148異なる地点・延べ237地点×年**。元source-lockedの数字を再現できなければSTOP。
- 元**BIO_SAMPLE_DATE の日付のみ**を使用（調査年から月を推測しない）。調査対象の上流SAの**phenomenonTime**も正式APIの元観測日付から読む。同年だけで適当に値を付与しない。上流SOURCEとして選ばれた元BIO 237地点×年に、BIO元採集日 **510 source site×date**。重複した同一BIO日付や同じ年の複数original analysesを独立な処理場や生物地点には数えない。
- BIO当日とSA当日が等しい場合は採水時刻の前後が不明であるため**「先行したTP」として扱わず**、same-dayとして別指標にした。観測日の前30/90/365日、**後90日**を日付の元データから別々に判定。2017年の365日前は**2016年のTP監視履歴を取得していない左側打切り**を含み、2018年の365日前と対称な観測窓ではない。

## 調査前後のSOURCE監視記録の有無（同一元source 237 BIO地点×年）

| 元TP測定の時点（任意の上流候補SAで一件以上） | 2017年BIO site×year（118中） | 2018年BIO site×year（119中） | 異なるBIO siteの和集合 | 合計BIO site×year | **全BIO採集日に条件が満たされたsite×year** |
|---|---:|---:|---:|---:|---:|
| **調査前30日以内**（同暦日は除外） | **116** | **117** | **145** | **233** | 192 |
| **調査前90日以内（一次）** | **118** | **118** | **147** | **236** | **233** |
| 調査前365日以内（2017年は左側打切り） | 118 | 118 | 147 | 236 | 235 |
| **調査後90日以内（監視頻度の対照）** | 117 | 119 | 147 | 236 | 234 |
| 同じ暦日（時刻順不明のため先行に入れない） | 20 | 28 | 45 | 48 | 7 |

元510 distinct BIO source site×dateに対応するTP先行測定のあるBIO日付は、前30日 **460/510**、前90日 **505/510**、前365日 **507/510**。調査後90日にも **506/510**、同じ暦日には57/510でTP監視記録が存在。前後90日ほぼ同じという監視頻度の性質を保持し、**時間的先行＝負荷の変化・機構の因果確認**と混同しない。2017年の2016年SOURCE歴は欠測なので前365日窓の欠測・打切りをゼロ源と扱わない。

### source-onlyの品質判定

- 元事後手続きで固定した前90日の **50異なるBIO site／100 site×year** は、観測された **147 site／236 site×year** で数値的に通過。さらに保守的な「そのsite×yearの**すべてのBIO採集日**に前90日TP測定が存在」でも **233/237 site×year**。元SOURCE流向/日付の標本機会を検証するには十分。
- しかし主SOURCE候補237のうち236で前90日TPの採水があり、後90日でも236あるため、**「先行TP測定あり/なし」というバイナリはほぼ飽和**。これを単独のWTP曝露のあり/なしとしてBMWP回帰に入れれば、源別の真正な差・施設効果を識別できず、監視頻度・源選択の内生性を過小評価する。**原河川BMWP値の新しい回帰は実行しない**。
- 逆流向placebo BIO地点は元source条件付きで2017年170、2018年176、同reach未確定61/71。これらを**下水非曝露地点**に置換しない（ほかの上流SOURCEがあり得る）。HydroRIVERSは小流路・側溝を欠落させ、EA official SAの採水位置は実際のpermit-outfallとは限らない。

## 統合研究への意味：水質とSOURCEはまだ別の測定層

元ChemPop主結果[PR #135](https://github.com/zuizui0223/gesui/pull/135)は、591の同じ実河川地点・2,624 site×yearに対してモデルEDFと出現科均等BMWPの空間相関が、原河川orthophosphate+nitrate+ammoniacalNの3水質窓集計を調整すると追加地域外予測できなくなる（EDF after water chemistry MSE−0.519%、0/7地域で改善、LOQ1/2/3・591siteのみでも未支持）。**この結果は「下水由来P/Nが河川生物へ因果媒介した」の証明ではない。**

今回のTP原SOURCEが生物調査以前に「測定されていた」という時間の機会は実際に確保できたが、未解決なのは**物理放流口→同一受水河道の特定、処理水の年度実測流量、下水源別リン／窒素質量、現地河川源別N/P受水、処理条件の外生的な変更**。濃度mg/Lだけから年間排出量kg/yrを捏造しない。現時点では、**source-firstを強くしただけでBMIの新しい生物学的因果結果は増えていない**。

## 再現・保存

- 元実行: [GitHub Actions #44838235844](https://github.com/zuizui0223/gesui/actions/runs/44838235844)（旧hard-coded source gate判定フィールドには方法上の誤りがあり、観測件数の結果は本表の通り）。後続実行では**計算から導く実際の90日source合否**に置換した。修正後実行を採用版とする。
- 元データ・元コードのSOURCE-first [日時条件プロトコル](PROTOCOL_20261120_UK_STW_TP_BEFORE_BIO_SAMPLE_DATE_SOURCE_ONLY.md)、[ソース計算](../analysis/20261120_wessex_wwtp/ea_SA_TP_upstream_before_ChemPop_bio_dates.py)。元原データ化学・生物・地理行や原station IDは公Gitへ再配布せず、集約数値と原source SHAのみ公開。

研究原稿・gesui/main・ほかのPRは変更しない。結果の符号や有意性を見てTP源の時間閾値や空間閾値を緩めることをしない。
