# 2026-11-20｜英国処理場の元最終処理水（TP濃度）を、河川底生ChemPopへ結びつけるためのSOURCE first監査

**本資料は放流元SOURCEが存在し取得できるかの検証であり、STWがBMWPを変えた新しい生物学的効果量や因果連鎖の証明ではない。**

## 既知の英国底生群集側の因果境界

- 原ChemPop v2 2010–18年の707site、3,884site×yearではモデル化された静的EDF（処理水の推定河川寄与）と出現科均等BMWPの条件付き空間勾配 **β−0.0911 [−0.1336,−0.0486]**、7地域留保EDF追加MSE改善+4.52% [0.59,8.35]。
- 同じ591site・2,624site×yearでOrthophosphate、Nitrate、Ammoniacal Nitrogenという原**河川水質期間集計3種**を追加すると EDF係数β+.0071 [−.0454,+.0596]、**EDF after 3chemの地域外MSE変化−0.519% [−.761,−.253]**、0/7地域改善。水質3種はEDFの後に+10.559% [5.062,14.290]、6/7地域改善。この結論はLOQ1/LOQ3と591site平均でも保持。**3化学種は現地にある濃度であり、下水から排出されたN/Pの発生源別実測量ではない**。元結果はgesui Draft PR #135に保持。

## 実際の下水処理場SOURCE探索で新しく分かったこと

1. **Wessex Water公式2010–19年の処理後排水実測原CSVと処理場地理XLSX**：公開カタログのファイル識別子はあるが、GitHub Actionsから原4 URLが全てHTTP403。数値・原列を読めず SOURCE ACCESS STOP。Wessex公式の元日別処理場流量も掲載が更新により停止されており、濃度からkg/yrへの変換は不可。[正式source receipt](https://github.com/zuizui0223/gesui/actions/runs/44674766035)
2. **EEA UWWTD2018年の実測N/P両方正＋淡水一意出口1,889施設**：17の適格報告コードに英国なし。施設N/Pの両方が測定正値というEU適格母集団内で英国sourceが0であり、英国STWが存在しない・排出ゼロであることを意味しない。[元EEA/source-frozen監査](https://github.com/zuizui0223/gesui/blob/main/results/eea2018_independent_source_frozen_20261108/2018_MEASURED_SOURCE_BIOLOGY_COHORT_QC.json)
3. **EA新WQE公式API (v1)**：旧WQAは2025年12月廃止。GET official determinand codelistから0348 Phosphorus Total as P、6008 Orthophosphate : Load as P、6003 Ammoniacal Nitrogen : Load as Nを真正確認。[API実行成功 #44675194578](https://github.com/zuizui0223/gesui/actions/runs/44675194578)
4. **観測から得られるものを識別**：追加の2017–18年 limited-source分析ではコード0348の「最終処理された水道会社の下水、Water Company、SA」から、**TP濃度 mg/Lの真正field/lab記録**を取得。Bath周辺50kmの2017年・2018年はページ末尾まで受信したが、英国全国の全源量を調べたものではない。公式6008/6003 loadのコードは存在するが同範囲の2017–18の限定検索ではHTTP204＝記録不在。負荷の無い施設としてゼロ値を代入しない。[小規模loadの原ソース #44675442028](https://github.com/zuizui0223/gesui/actions/runs/44675442028)
5. **地理パイロットの元の粒度**：[Bath50km・2017/18 Water Company最終放流TP SA 原件数取得 #44676173833](https://github.com/zuizui0223/gesui/actions/runs/44676173833)。2017はすべての0348の1,581観測中 **水道会社の最終処理水TP 932記録／46 native samplingPoint.notation／783 native地点×採水日**、うちUWWTDモニタリング目的540記録。2018は全1,726元観測中**最終処理水TP 1,066記録／49 native sampling points／874 native地点×採水日**、うちUWWTDモニタリング目的550記録。**46点の採水地点IDは両年で繰り返す**。46地点は「異なるSTW施設数」と検証できていないため、正式採水点数として記述する。source page cap 3,000/年に達することなく年ごとのソース取得が全て完了、元のWater Companyとnon-water company、storm、河川・湖沼を区分。公式原source地点名・番号・原濃度行は公Gitに再配布しない。
6. **ライセンス／登録の制約**：英国EA排水許可台帳にはsource座標・許可日時等があるが、公開条件の一部にはEnvironment Agency Conditional Licenceによる再利用制約がある。permit limitは「測定処理場負荷」ではなく施設排出許可の制限値。無許可の公Gitへの転載・代替の非公式ミラーでの回避はしない。EA Water Quality Explorer自体は公式Open Government Licence。
7. **英政府2026年3月のリン回復調査**はリン濃度・生態状態の回復を評価しているが、実際の処理場介入→河川→生物形質を直接因果的に実証した研究ではない、という出版物自身の限定を保持。

## 次のSOURCE連結（地理的近さと流下を混同しない）

Bath50kmの2017/18 EA SA TP元採水点の正式英国BNG(easting/northing)と、ChemPop v2原7地域の BIO_SITE_ID/BIO_EASTING/BIO_NORTHINGを同じ元データから取得し、同年（2017または2018）のsource-eligibleなChemPop地点に限って0.5/1/2/5/10kmの直線距離SCREENを行う。これは本質的には**SOURCE供給と受水候補の機会**の数え上げであり、河川上流から下流への到達があると証明しない。

同一河川の正式ネットワーク流向・上流下流順序、STWの正式permit/outfall identity、同年P放流水濃度、下流生物サイトが結合できるまで **STW→BMWPの新たな因果回帰はSTOP**。下水の元量（kg/yr）は実測流量がない限り計算できない。到達源でない地点をゼロ曝露として扱わない。

関連プロトコル：PROTOCOL_20261120_TRUE_WWTP_TO_BENTHIC_2010_18_SOURCE_LINK.md、PROTOCOL_20261120_EA_SA_TP_CHEMPOP_SITE_SPATIAL_SOURCE_GATES.md。全部Source-only成果物、元EFX/TPの正負結果・アクセス制限/未成功を削除せず、gesui/mainと別の生物系原稿は触らない。
