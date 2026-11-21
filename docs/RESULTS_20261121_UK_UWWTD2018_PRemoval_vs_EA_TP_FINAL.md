# 2026-11-21｜英国2018：リン除去処理場の報告フラグと実測最終処理水TP濃度の関連は、施設規模と地域の調整後に残らなかった

**別原資料をつないだ事後探索・観測的な処理場source比較。UWWTD処理場のnative施設IDとEA採水地点のpermit/native ID一致ではなく、相互一意の位置対応を使っているため、元処理場のリン除去の因果効果と呼ばない。**

## 元データ・支持条件

- [EEA UWWTD2018英国原3科学表](https://sdi.eea.europa.eu/datashare/s/mGmXg3zcEqKQLTG/download)：内部reporter UK5・元2018年、1,799施設／1,809放流地点／1,498淡水FW放流地点の原座標。全原放流点とplantsのoriginal SHAを過去成功runと照合。正式[Eionet UWWTPs辞書](https://dd.eionet.europa.eu/datasets/latest/UWWTDArt15/tables/UWWTPs)ではリン除去フラグPRemoval、窒素除去NRemovalはbooleanの処理特性。元全施設のPRemoval=1は418/1799、NRemoval=1は40/1799。
- EA公式Water Quality Explorer新APIの処理済み水道会社最終処理水SA採水地点。UWWTD元FW放流点とSAが250m以内で双方相互一意、かつ元UWWTD施設に原放流地点1つという地理候補は**930**。このうち**360**SAで真の2018年総リンTP as P、mg/L、FINAL SEWAGE EFFLUENTの非負数値観測が存在。原native source code・TP原個票は公開Gitに転載していない。
- [原2018処理フラグ↔TP比較 Actions #44850033349](https://github.com/zuizui0223/gesui/actions/runs/44850033349)、元集約JSON artifact **14684493300**。360地点はPRemoval=0が113、PRemoval=1が247。TPは同一SA採水地点の独立日ごとに実測値の中央値を取り、それを2018年各地点1つの年間中央値へ要約。各SAが別施設の独立単位であることは行政IDで未確認。
- 元支持条件（PRemoval各30地点、少なくとも3元NUTS地域群で各10site、元cap正値が各群70%以上、両sourceフラグ真偽が存在）はすべて通過。360の全地点で元容量population equivalentの正値があり、元9つの3文字NUTS地域に対象が分布。2018 TP測定が2か月以上あるSAは354/360。元観測に基づく方法と係数、失敗した出力実装runは監査として残す。

## 主データ：PRemoval=1グループの方が単純中央値は高かった

| 元source 2018 SA 1地点＝2018 TP中央値 | PRemoval=0 | PRemoval=1 |
|---|---:|---:|
| 元source対象地点 | **113** | **247** |
| 最終処理水TPの中央値（mg/L） | **0.652** | **0.965** |
| source第25パーセンタイル（mg/L） | 0.493 | 0.589 |
| source第75パーセンタイル（mg/L） | 0.9705 | 1.2275 |

全sourceの中央値比PRemoval1/PRemoval0=**1.4801**、差は**+0.313 mg/L**。源観測地点を群内等重みで2,000回再標本化した中央値差の記述区間 **[+0.205, +0.4295] mg/L**、比の記述区間 **[1.299, 1.760]**。これは除去処理がリンを増加させる証拠ではない。処理導入は地域・処理施設の規模・規制基準・流入負荷などに依存し、監視される地点の選び方も異なる。

## 事前固定の元source背景調整

同じ360 source地点に対し log(1+原TP mg/L) を目的変数、元PRemoval flag、NRemoval flag、元処理場容量population equivalentのlog1p、元uwwNUTS3地域の固定効果を説明に含む観測OLS、HC3記述95区間。

- **調整後PRemoval係数 β=+0.000597**、HC3記述95 **[−0.09519,+0.09639]**、名目 p=**0.9903**。ゼロを明確にまたぐため、単純中央値の差をリン除去の独立した同年化学効果としては扱えない。
- 本研究は現地処理前水中のP質量や処理効率、許可のP基準、施設ごとの源流入N/P負荷やsource実流量を計測していない。調整した地域・施設規模による偏りを部分的に減らすものの、処理の非ランダム導入や地理的な施設/採水点ミスマッチ、監視努力の偏りを完全に除去しない。
- 元2018 UK UWWTPSのN/P discharge *measured/calculated/estimated* 計6報告欄は**全1,799施設で空欄**、ゼロ放流とみなさない。したがってTP mg/L×流量で年kg質量に換算する根拠はない。UKPRTRは制度・施設IDが別で、P/N同一施設・同年報告0。源負荷から河川生物BMWPへの因果媒介はSTOP。

## 強い選択過程の可能性

元2018 UWWTP 1,799施設のPRemoval1は418（約23%）だが、地理候補かつ実TP観測がある360地点では247（約69%）。これは**TP監視・対象施設の選択がリン除去処理と独立でない可能性**を示す。ただし両集団には施設規模や淡水放流の有無、2018観測の有無、geographical one-to-one matchといった異なる選択条件が入っており、この比だけで欠測機構の原因や全国の処理場の代表性を断定しない。追加source-selected分母（930）におけるPRemoval内訳を別途QAし、単純比較を因果的な処理効果へ昇格させない。

## 水域の統合研究への結論

この元source比較の新しさは、処理方式PRemoval→処理水中TPという**発生源側の独立した化学特性**を元2018 UWWTD×EA WQEで接続できたことであり、単純比較では予測と逆、地域・施設容量・NRemovalを調整すると明瞭な差がないという**反証的な結果を得た**こと。これは既存ChemPopの河川orthophosphate/nitrate/ammoniacalN→科等重みBMWPの予測情報がEDF proxyを上回るという別母集団の事実とは、まだ一本の因果経路で結べない。処理水の実際のoutfall到達、同一source streamの測定、栄養輸送、河川BMI変化を同一河川網・時期で直接推定する必要がある。

事前方法書：[2018 UWWTD PRemoval vs EA TP源品質・調整](PROTOCOL_20261121_UK_UWWTD2018_PREM0_REMOVAL_TO_EA_TP_SOURCE.md)。別PR #135には英国底生のモデルEDF/river chemistry結果、PR #136にはEA SA source network/time quality監査。本稿はPR #137 Draftの新独立source成果であり、mainや既存原結果に変更を入れない。
