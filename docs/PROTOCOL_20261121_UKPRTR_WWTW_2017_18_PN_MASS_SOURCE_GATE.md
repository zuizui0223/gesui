# 2026-11-21｜英国PRTR：2017・2018年下水処理場の施設別P/N年間排出量を本物の原報告から取得可能か

**この段階は結果を見た後の探索的source-first調査であり、処理場の河川生物への因果影響モデルではない。** 先行の水域結果は、2017/18年のEA final-treated SA採水点の元TP mg/L測定と、元ChemPop適格BIOのHydroRIVERS上流方向・採集日前90日という条件を満たすことが多い（147 unique BIO sites、236 site×year）。しかし本当の物理outfall点・実際の年間P質量・河川中の源別化学到達は未確認。

## 本来の新しいsource：UK Pollutant Release and Transfer Register

- 英国政府の原[UK PRTR 2017/18 CSV](https://www.gov.uk/government/publications/download-industrial-emissions-and-waste-transfer-uk-prtr-data)（2026年公開アーカイブ）を年別に独立確認。2017は `uk_prtr_dataset_2017.csv`、2018は `uk_prtr_dataset_2018.csv`、両方原URL・SHA・元列名・対象年度・元sourceラベルでsource-lock。
- PRTRは施設／permit単位の **年間の排出量**（単位は必ず元sourceから読み、通常kg/yかtonnes/yのいずれかを勝手に仮定しない）。その施設が水道会社WWTWか、民間産業の工場、産業廃水処理施設等を、原の活動code、説明、source文書から識別する。
- [EA公式 Pollution Inventory reporting](https://www.gov.uk/government/collections/pollution-inventory-reporting) は下水処理場が原則 **100,000 population equivalents以上**の場合に報告対象となる。大規模施設の一部のみであり、元EA SAの全5,722採水点とは別母集団。PRTR不記載はP排出ゼロを意味しない。
- 元[DEFRA PRTR documentation](https://www.gov.uk/guidance/uk-pollutant-release-and-transfer-register-prtr-data-sets)では報告量の由来として M=Measured、C=Calculated、E=Estimated、U=Unspecified、正式測定方式の下位codeを定義。**Mを独立field検証できない限り「すべて実測負荷」と書かない**。Mだけの相対的coverageを事前に数え、C/Eを感度・独立の列で管理。
- PRTRは許可区域の空気/水/土壌への各放出や、施設外の下水への移送がある。**receiving water releaseのみ候補**。`Offsite transfer in wastewater`・`Release to sewers`・`Land/Air`等は源の河川への直接放出ではなく除外。更に施設年度の総排出と、物理outfallの個別水路への排出は別。
- P/Nの名称は元正式コードを識別し、`Total phosphorus`／`Phosphorus compounds`／`Total nitrogen`など名称が違えば推測で統合しない。リン酸塩濃度（orthophosphate）、総リン、アンモニア性N、全窒素を同質の指標としない。年別sourceハッシュ・行重複・missingness・LOD/thresholdの別を必ずQA。

## 段階別source gate（生物BMWPを開く前に固定）

1. 両年の元CSVのGET真正200、元byte量・SHA、元source schema/年度、category・活動code、pollutant names、unit、release mediumを検査。公式元アーカイブにyear/plant/quantity/method/mediumのどれかが無い場合は該当部分SOURCE STOP。異なるsource年の数字を都合よく補完しない。
2. 英国のWWTW施設に限定し、水への**年間**PとNの値が同一年度・同一真正施設IDで存在する数を数える。`measured (M)`／`calculated (C)`／`estimated (E)`は原値に従い分け、測定方法不明なら別statusとする。水域放出と管渠送出の境界を保持。高い数値だけ選ばない。
3. 真のWWTW facility originが見つかっても、EA SA source point採水IDとの全国crosswalkには、施設名・許可ID・同じ年の活動/地理の独立照合が必要。河川に沿う自然流下のNEXT_DOWNに**施設の主建屋位置をそのまま排水の物理outfall位置として置かない**。許可区域あたり複数outfalls/工業排水の可能性を記録する。
4. さらにBIO 2017/18と同じ上流 reachへ流れることが真の原outfallで示され、局所の測定河川P/N chemistryと採集時点の対応が成立して初めて、源→河川の条件付き予測へ進む。下水への曝露の有無だけの回帰（90日source monitoringは147/148 sitesでほぼ飽和）はSTOP。
5. 原施設源の年排出量が取得できても、2017/18だけで河川BIO地点×年に実験的に外生変動があるとはみなさない。施設選択、農業・都市・河床等交絡を調整し、真正逆方向placebo・地理外留保、実測水質を追加した反証を事前に固定し直す。

## 既存英国源と独立したsourceカテゴリー

- UK PRTR（施設年間報告mass）はEA Water Quality Explorer SA（元最終処理水TP **mg/L**）ともChemPop（元受水域orthophosphate/nitrate/ammoniacalN期間集計）とも異なる。
- 元EEA UWWTD2018『N/P両方を実測正値＋一意淡水放流口』1,889施設の厳格source cohortには英国報告元が存在しなかった。UK PRTRに正しく別年・別制度の英国annual massがあれば相補的sourceとなるが、EEAの欠損を別制度の施設コードで同一sourceと錯覚しない。
- Wessex WaterとDEFRA Article15の直ファイルはActions環境でHTTP403だった。正式なUK PRTR CSVの取得が成功するまではPRTR由来の係数・元件数は未確定。

元許可台帳のConditional Licenceに関わる原個票や各施設原データの公Git転載は行わない。公式PRTR自体の利用条件も元ファイル/元情報から検証し、公Gitには**集約source-coverage、code、SHAと判定・方法書だけ**を保存する。gesui/mainおよびPR #135・#136の既存結果を上書きしない。
