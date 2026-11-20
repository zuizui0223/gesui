# 2026-11-20｜最終処理水の原実測TPとChemPop生物地点の真正同年地理照合：Bath50kmでは必要数に不足

**SOURCE ONLY。厳格な原source IDsと座標による地理的候補スクリーニングであり、処理場の栄養塩が下流へ届いたこと・生物学的反応・下水由来P実測排出質量の推定ではない。**

## 元EA処理後排水とChemPop生物試料の真正照合

- [2017/18英国EA Water Company SA原TP全ページ取得 #44676173833](https://github.com/zuizui0223/gesui/actions/runs/44676173833)：Bath周辺50km、元TP (Phosphorus Total as P, code 0348) mg/L、採水地点type SA（Water Company final treated sewage effluent）、元sample material FINAL SEWAGE EFFLUENT。2017は46 native SA測定地点・932原records・783 source site×date、2018は49 native SA測定地点・1,066原records・874 site×date。**46 native sampling-point notationsが両年に反復**。2017のUWWTD monitored SA40地点、2018は42地点。採水地点数をそのまま処理施設数へ解釈しない。Storm SB、非WaterCompany UA、河川・湖沼等は別に集計し、source値0とみなさない。
- [新EA source×ChemPop真英BNG座標監査 Actions #44676597375 成功](https://github.com/zuizui0223/gesui/actions/runs/44676597375)。EA側は原samplingPoint.notation、samplingPoint.easting/northingで元同一点内の位置変動最大**0m**（重複source geometryドリフト100m停止を通過）。ChemPop EIDC元release v2 ZIPの7地域matched-location原CSVにあるBIO_SITE_ID/BIO_EASTING/BIO_NORTHINGを読んで**1,519真正生物調査site全てに元BNG位置**あり。2000年代の原論文とは別データが再収録され得るため、EA採水地点とChemPop原siteを別レイヤーにして独立性を勝手に仮定しない。
- 元BMWP source-only適格は2010–18の**7,836 original analyses／707 BIO site**で検証済み。今回必要な「同じ採水年だけ」の eligible ChemPop は2017年 **389 BIO site**、2018年**366 BIO site**、いずれかの年に真正source適格なもの463 BIO site。同じsiteの2017と2018を独立した別施設数に数えない。

## 元地点のBNG直線距離による機会SCREEN（河川流向ではない）

| EA正式SA TP採水地点までの元地理直線距離 | 全1,519 original ChemPop site（2017/18 source pool） | **2017年のsource適格 ChemPop BIO site** | **2018年のsource適格 ChemPop BIO site** |
|---|---:|---:|---:|
| ≤500 m | 2 | **0** | **0** |
| ≤1 km | 8 | **1** | **1** |
| ≤2 km | 17 | **2** | **2** |
| ≤5 km | 30 | **3** | **3** |
| **≤10 km** | **40** | **5** | **4** |

最初の全site 40は2017/18年の生物試料が存在することを要求していないため、**過大なsource機会**の参考のみ。主sourceスクリーンは2017年5・2018年4の生物site。元水道会社UWWTD目的だけに制限しても同年10km以内の候補は2017年5 site、2018年4 siteのままだった。Bath50km範囲以外にある英国処理場の全曝露を検証していないため、「全国で処理場近くに生物がいない」とは言えない。

**重要な停止判定：** 元事後SOURCE事前下限は「30異なる処理場原採水地点・≥2測定年・5独立hydrologic blocks・真に有向到達するChemPopの異なる**50site/100site-year**」であった。Bathパイロットだけでは(1) Water Company SA pointの年代供給は十分(46/49)、(2)同年source適格生物siteの地理的近接スクリーンは5/4に過ぎず、(3)同じ河川の**流向・下流の到達**を全く実証していない、(4)排出源のannual massは流量不在で未計算。よって**SOURCE GATE STOP／英国ChemPopでこのBath50kmだけからSTW→BMWPの新回帰を実施してはいけない**。

## 真の次のデータ取得に向けた仕様

- EA公式新Water Quality Explorerの全国「Water Company final-treated (SA)」原採水点APIで、SA地点に本当に絞れているかを同じAPIのSB storm-overflow陰性フィルタと比較する**source contract test**を追加。APIが未知paramを無視する場合、全国source code一覧とみなさずSTOP。
- もしSAのnationwide inventoryが実際に機能すれば、全国7地域のoriginal ChemPop同年BIO地点の元座標と同じBNGで処理場samplingPoint source候補を抽出し、元permit/outfall位置・流向河川ネットワークの独立確認に進む。曝露陰性と未登録・未計測を混同しない。
- 年間下水由来TP/TN質量（kg/yr）が必要なら、総リン濃度mg/Lとは別に原実測水量または公式報告排出massを追加取得する必要がある。旧EEA2018の実測NP両方正一意淡水出口1,889には英国の適格報告コードがない。Wessex原CSVの直接ファイルはHTTP403で元値不取得、元WRC日流量sourceも更新停止。いずれも隠して補完しない。
- 許可台帳は限定ライセンスであり公式許可制限値をそのまま観測源濃度・年質量に読み替えない。EA Water Quality Explorer原観測はOGLだが、元個票や住所を無制限に公Gitへ転載せずsource aggregate/QA/SHAのみ保存。

**生物学的主結果は不変：** 英国ChemPopではモデルEDF→BMWPの元係数は負、しかし河川orthophosphate・nitrate・ammoniacalNを同一site×dateで条件づけたらEDF after 3chemの別地域予測は−0.52% (0/7)、これは別既存PR #135の結果であり、EA STWのsourceデータを加えた新しい確認実験ではない。次は**sourceを実際に近隣/同年の河川有向経路で到達させられるか**に集中し、追加の肯定が出るまで実装条件を緩めない。

ソース方法書：PROTOCOL_20261120_TRUE_WWTP_TO_BENTHIC_2010_18_SOURCE_LINK.md、PROTOCOL_20261120_EA_SA_TP_CHEMPOP_SITE_SPATIAL_SOURCE_GATES.md。作成物はDraftで、gesui/main、米国EPT、欧州陸域原稿は変更しない。
