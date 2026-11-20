# 2026-11-20｜英国下水処理場の原放流情報→受水域→底生動物の同一地点・同一時期の接続

**SOURCE FIRST／結果を見た後の探索的データ取得プロトコル。** 真のSTW処理施設のN/P測定と英国ChemPopのBMWP変化が因果的に結びついたと主張するものではない。未取得データに都合のよいゼロ値を入れない。

## 統合する目的

英国ChemPop v2の同じ591地点で、EDF（静的な処理水混入予測）と底生科等重みBMWPの負の関連は、orthophosphate、nitrate、ammoniacal nitrogenの**河川水質期間集計**を追加すると消える。7地域外予測ではEDF after 3chemはMSE−0.519%、水質3chem after EDFは+10.559%。LOQ1/2/3や591site平均でもEDFの追加予測は不支持。ここからは「実際の処理場で出たP/Nが受水域に届くか」を源から独立に照合する必要がある。

## 現時点で使えるものと不可のもの

1. **Wessex Water公式原STW最終放流水質2010–2017、2018–2019** のCSVとWRC地理XLSXはデータページ上で公開が確認されたが、GitHub Actionsからの原ファイルGETはいずれも**HTTP 403**。取得・原値・現地STW名を検証したとは主張できない。別の非公式ミラーで制限を回避しない。
2. **英国STW流量**について、Wessex公式の元日別放流流量ファイルが更新のため掲載停止。源内濃度×観測流量の厳密な時間積分ができず、**元の処理場別年間N/P kg/yrの計算はSTOP**。
3. **EEAの2018年UWWTD真正「N/P両方測定正値＋一意淡水出口」1,889施設**の品質対象に英国の報告コードは0。これは英国に下水排出がない意味ではなく、適格な両物質の報告母集団での欠測。英国ChemPop地点へ2018年EU報告質量を無条件に転用しない。
4. **英国EA新Water Quality Explorer v1**：2025年12月に旧測定APIが終了。現行POSTの /data/observation、/data/sampling-point とGET /codelist/determinand が使用可能。水質原観測には native samplingPoint.notation、phenomenonTime、sampleMaterialType、samplingPurpose、determinant notation、result、unit がある。2017/18年Bath周辺で、元分析項目0348 Phosphorus Total as P に対して**最終処理済み下水のWater Company (SA)**、samplingPurpose UWWTD monitoring、単位mg/Lが原記録で確認された。014...など未確認コードを絶対に推測しない。公式の6008 orthophosphate P load、6003 ammoniacal N loadの項目辞書はあるが、その年・範囲の元load記録は0件の限定検索だった。
5. **EA正式採水地点の区分**：SA=water company FINAL/TREATED effluent、SB=water company storm tank/overflow、UA=not water company final/treated effluent、RA=crude sewage. **主STW曝露候補はSAに限定**。公定サンプル用途もWATER QUALITY UWWTD MONITORING DATA等の適切なものが明示されたもの。雨天時越流、未処理下水、農場／工場／非水道会社放流、河川そのものは別系列。
6. **EA排水許可台帳**の元個票はEnvironment Agency Conditional Licenceで再配布や利用目的に制限がある場合がある。公Gitへの個票転載やアクセス制限回避は禁止。許可上限（permit limit）は施設の実測平均値・kg/yrではない。

## SOURCE-ONLYステージ1：水道会社最終放流の年代・項目・単位

- 確定コード0348 Phosphorus Total as Pを主候補とし、EA WQEの2017/18各年について、Bathを中心とする50kmの元公開観測を元date/coordinates/determinandで**上限付きページ送り**する。各元source rowにSA、UWWTD目的、FINAL SEWAGE EFFLUENT、mg/Lかどうかを独立判定。地点数／年、同一点が両年で再現する数、月カバレッジ、検出下限・欠測数を記録。
- 水道会社以外のUA、嵐/越流SB、未処理RA、河川観測、湖沼観測を候補から除外し、その数も別に記録。最初の250/20件だけでゼロを宣言しない。ページ上限到達時は SOURCE CENSUS INCOMPLETE とし全国の無記録とは言わない。
- 補助項目6008と6003はコードの公開だけ確認済み。原実測値とunit/source sampleがそろわないなら kg/day・kg/yearと呼ばない。濃度mg/Lを年度平均P massと乗算しない。
- 将来の施設年度原水質曝露の最低データ条件：**30異なる正式SA採水地点、2年の実測source、5以上の独立河川流域ブロック**。今回の狭いBath pilotはこの条件を自動的に満たす設計ではない。

## SOURCE-ONLYステージ2：真の上流放流水→ChemPop受水地点の地理的経路

- ChemPop v2の元1,519siteの BIO_EASTING、BIO_NORTHING と、WQE SA samplingPoint easting/northingをEPSG:27700で扱う。site ID、試料日、元source official type/versionで対応の誤差・重複を検査。
- 同じ元河川網（HydroRIVERS等）の流向に両方を≤500mでスナップし、正式 NEXT_DOWN / 河川下流方向の**有向到達経路**があるときのみ「上流STW→下流の底生site」とする。同一河道segment内での出口順序が不明なら「不明」で除外。単純な半径・隣接流域・上流語だけは不可。
- 論文主要source gate：真正SA outletが**30 sites、≥2測定年、少なくとも5独立流域ブロック**、それに接続するChemPop生物地点が**50 distinct BIO_SITE_IDかつ≥100 site×year**。stream-connected陽性が不足ならモデルSOURCE STOP。放流水の年代は2010–2018の生物観測年と同じ年または固定された前年だけで、2018年施設値を2000年群集へ転用しない。
- 下流のSTWを上流と誤判定しないplacebo方向、河川の上流/下流未登録源、農地/都市と河道工学、流量や河床を事前固定の対照に含める。施設側実測ゼロと非報告・未取得を混同しない。

## 因果モデルの許可条件

- 一次BMWPは元ChemPop出現科等重みの平均、個体数加重は元sourceで実際に数えられた「Actual abundance」方式のみ。元source qualityと全陰性・不確実な転送結果を消さない。
- 施設年度の真の測定P・N**質量**が得られない場合は、STW原最終放流の**濃度**のみの別解析候補として保留。河川内濃度とは流量・希釈過程で異なる。
- 空間経路・時点両方成立した場合にのみ、地理・農地・都市・河床・採集法→下水施設source濃度／質量→同地点河川water-chemistry→BMWPを対照付きで比較。真の因果媒介を認めるには施設の変更時期や外生性の検討が別途必要で、単なる回帰の係数消失を「栄養媒介の証明」にしない。

必要な一次文献：EA新WQE 2025/26 APIの公式announce、EA WQE source licence OGL、Wessex原データ・流量withdrawn notices、EEA UWWTD2018正式元品質監査、ChemPop v2。元原票を公Gitに再配布せず、集約JSON・固定SHA・アクセス可否・SOURCE STOPを残す。gesui/mainや他研究原稿は変更しない。
