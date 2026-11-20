# 2026-11-20｜EA全国最終処理水SA登録5,722点と英国ChemPop生物調査：地理的候補の真の上限

**これはsource-onlyの地理的候補スクリーニング。登録2026年の最終処理水採水点を2017/18年の実測・稼働地点、処理場の施設本体、実際の排水口、上流流下到達、WWTP年度N/P質量とみなしてはいけない。**

## 1. 以前の局所パイロット（不十分）

2017/18 Bath50kmのEA公式Water Quality Explorer APIから、元Water Company FINAL/TREATED EFFLUENT標本（type SA）、元 TP＝0348 Phosphorus, Total as P、元unit MILLIGRAM PER LITREを抽出：2017は46 SA採水点/932元TP濃度観測行、2018は49 SA採水点/1,066元観測行、46元SAコードが両年に存在。原同じ年のChemPop scored family source eligible siteは2017年389、2018年366。Bath50km由来SAから10km以内の同年BIO sourceは2017年5、2018年4だけであり、距離のみでは下流曝露を識別できない。元元source run #44676173833/#44676597375。

## 2. 全国の正式 Water Company SA 台帳を完全取得したSOURCE陽性

[全国正式5,722 Water Company final-treated (SA) 登録点の原本ページ全取得＋BNG原座標チェック Actions #44677688277](https://github.com/zuizui0223/gesui/actions/runs/44677688277) 成功。API /data/sampling-point のPOSTで真のsamplingPointType=SAフィルタを使い（storm overflow SBを使った陰性フィルタ契約 run #44677052429）、250件×23ページで元資料最後まで取得。current official totalItems＝**5,722**。全5,722点のoriginal asWKTを英国BNGメートル座標として正常解釈。欠測・unparsable original WKT 0、異なるnative samplingPoint.notationの重複0。元SA地点のstatus code O＝5,060、C＝662。ただしO/Cの公式意味と2017/18年の存在/稼働は検証前なので、「5060の稼働中処理場」と呼ばない。

元ChemPop v2原ZIP original biology SITE matching-locations sourceより**1,519 BIO_SITE_IDの座標がすべて存在**。元生物の2010–18年SOURCE品質合格7,836試料／707地点、同じ対象で2017年**389 site**、2018年**366 site**、いずれかの年**463 site**（同地点の年を独立site数へ足し算しない）。

**全国現在登録SA点とのBNG直線距離候補**（同年のTP実測STWという意味はまだない）：

| 現在登録SAへの最大直線距離 | 元2017年生物SOURCE 389 siteのうち | 元2018年生物SOURCE 366 siteのうち | 2017/18いずれかに観測463 siteのうち |
|---|---:|---:|---:|
| 500m | 54 | 59 | 73 |
| 1km | 136 | 134 | 171 |
| 2km | 237 | 235 | 293 |
| 5km | 380 | 357 | 450 |
| 10km | 389 | 366 | 463 |

近接しているのは同じ川とは限らず、河川源の投棄位置までの正しい距離ではない。絶対ゼロSTW sourceを仮定しない。 Bath局所パイロットから全国へ広げると**地理的には2km以内だけで2年和集合293 BIO site**が存在し、原水路SOURCE連結を試す意味があるが、これは曝露陽性の実数ではない。

## 3. EAの真の2017/18年総リン測定地点を原APIで確認する方法

[EA WQE 2026 native pointNotationの一点絞り契約 #44678002358](https://github.com/zuizui0223/gesui/actions/runs/44678002358) では、2018の2つの正式SA TP採水点に対してそれぞれ独立24件が返り、元sourceの別SA点は混入なし。架空のnative source notationは正式HTTP400拒否（既存sourceかどうかについてAPIが区別する）。この事実は「旧台帳のどんな採水地点でも当該年の測定が存在する」という意味ではない。

さらに[EA原APIの複数native点のqueryフィルタ契約 #44678384278](https://github.com/zuizui0223/gesui/actions/runs/44678384278)で、comma-separated pointNotation＝2独立STWの各24件が結合された**48件を別ID混入なく正確に回収**。一方、pointNotationパラメータを二度繰り返すと1点24件しか返らず**source網羅性を損なう**。正式取得はカンマ連結で最大100 native sampling point codes/回、server skip/limitでページ送り、sourceかchemistryが未報告なら0でなく missingとする。

[同一生物観測年のTP原実測sourceのみの全国近接スクリーニング](../analysis/20261120_wessex_wwtp/ea_WQE_national_SA_actual2017_18_TP_chemistry_bio_screen.py) は、現在登録SA中のChemPopの2017/18 SOURCE適格463siteから直線2km以内にある**367元SA地点**（原BNG座標と元samplingPoint.notation）を、毎回最大40地点ずつ同一元年で直接問い合わせる設計。2017/18の元phenomenonTime・FINAL SEWAGE EFFLUENT・source type SA・0348・単位mg/L・非負numeric TPで絞り、元データの真の過年度SOURCE対象施設数/同年近接BIO site数を数える。符号やBMWP生物係数を見て母集団を選ぶことはしない。**この検証が終了するまでは2017/18の真の近接元source nを未確定として扱う。**

## 4. 真の下水→群集の因果ルートのSOURCE STOPはまだ残る

今後も必要なのは、(i) original SA観測点と本来のSTW outlet/permitのsource ID・同年度の測定値の同定、(ii) Hydrological networkでそのoutletから下流のChemPop BIO siteまで**有向の水路到達**が証明されること、(iii) 真の同年度処理場N/P質量に流量と濃度測定の時間積分、または正式源質量報告があること、(iv) 非下流や別支流のnegative controlである。

元既存EEA2018実測両NP正値・淡水一意出口1,889には英国のsource report codeが0。 Wessex Water原2010–19のCC BY CSV/XLSXはこのActions環境からHTTP403で元dataを閲覧できず、元WRC日流量も更新で公開停止。EA水質APIの源は処理水中のTP濃度mg/Lに限り、2017/18年リン酸塩P負荷code6008、ammoniacalN load code6003はBath限定小規模searchで原数値が無かった。ゼロ負荷を推測してはいけない。

**主張の上限：** 全国の地理的機会は十分ある可能性が示せたが、実際の2017/18に測定されたSTWの流下P/N→科BMWPの関係は未実証。現時点の最も強い生物結果はPR #135内の英国河川3栄養関連水質summaryを調整するとEDFの独立的なBMWP予測増分がなくなるというsource-verified観測結果である（原同一591sites、EDF after chemistry MSE−0.519%、LOQ/single-site頑健）。

公開GitへEA original native samplingPoint ID、座標、TP原観測行、生物元試料表やBMIの正式辞書をそのまま転載しない。元出典・ソースIDのSHAと集計結果を保持する。Main・他原稿・他生物リポは変更しない。
