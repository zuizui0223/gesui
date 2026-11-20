# 2026-11-20｜全国SA全5,722点を流向で先に絞る、処理水TP観測から底生調査へのSOURCE-only反証

**事後探索・完全source-onlyの必要な設計修正。既知の水域EDF→BMWP、3河川化学種の条件付け結果と、前回の47?施設/79地点×年の地理・流向源QAを知った後に固定。**

## なぜ前回の「79地点×年」では母集団の上限ではないのか

以前の流向source gate [Actions #44679903035](https://github.com/zuizui0223/gesui/actions/runs/44679903035) では、EA現行Water Company FINAL/TREATED EFFLUENT type SA の全5,722地点を最初にBIO原siteとの**直線2km以内**に絞り、そこからのみ2017/18実TP 0348 mg/Lを取得した。その67（2017）／71（2018）source pointからHydroRIVERS v1の流下経路 ≤50km を評価：原BIOの上流となる候補36（2017）、43（2018）、和集合51 distinct BIO／79 source BIO site×year。これは「その事前2kmで拾ったSA地点だけ」という条件付き結果である。

一方で、河川に沿った50kmの真正上流STWはBIOから直線2km以内にあるとは限らない。従って**最初に2kmでSTWを切るのは、定義された上流曝露源を除外する選択バイアス**になり得る。ここでは効果や記録BMIに関係なく、現在の全国5,722SAを元英BNGでHydroRIVERSにスナップし、流向50km以内の到達候補を構築してから、**同年2017/18に実際にTP測定のある採水点のみ**を利用する。これが科学的なsource sensitivityであり、前のQCの合格基準を結果に合わせて緩めるものではない。

## 元母集団・停止条件を変更しない

- Source: 2026年WQE全国 /data/sampling-point 5,722 unique Water Company FINAL/TREATED EFFLUENT type SA、原native notation・元asWKT英BNG原座標、23 original pages、source type strict SA vs storm SB original negative control。**現行登録は2017/18に測定を意味しない。**
- Biology source: ChemPop v2 EIDC original site×specimen source key、2010–18 score eligible 7,836 original assays/707 BIO sites、うち2017年389、2018年366、union463 original BIO sites。BMWP値やEDF効果、原生物個体数・種の形質応答はこの段階で**一切開かない**。
- River: 正式HydroSHEDS HydroRIVERS v1 EU original ZIP SHA256 **500da7d36ceee0aa4c82dd625f5c53d80ba2f7ac14ad022cd3a6601bf55068e5**、NEXT_DOWN有向reach。BNGへproj、最大500m snap、同距離1m内で複数reachが競合すると除外。原サンプル点が実際の**許可されたoutfall**と一致しているとは仮定しない。
- Source eligibility: 同年の元WQE /data/observation POST、source year2017/18、determinand0348 Phosphorus Total as P、Water Company SA、material FINAL SEWAGE EFFLUENT、unit MILLIGRAM PER LITRE、nonnegative numerical result。コンマ区切りpointNotationが正しく複数native codeを選ぶ契約確認済み、繰返しパラメータは禁止。元個票・濃度や地点住所をGitへ公開しない。
- ネットワークsource候補は、同年のBIO siteと、SA pointの位置するreachが **(A)≤50km upstream in NEXT_DOWN**, **(B)≤50km downstream negative-control**, **(C)same HydroRIVERS reach order UNKNOWN** のいずれかにあるもの（source-level preselection）。別支流/50km外、un-snapped、wadeable small reaches missingは**測定ゼロ**とは扱わない。Aのみを取得してBを省く選択は禁止。
- まず**元全登録SAのネットワーク位置関係の件数**を固定し、次に元同年TPを取得した採水点に限定して出現SOURCEを再分類する。候補の数が多くWQE APIが制限する場合、元source-selected batchを分割し、全件取得できなければSOURCE CENSUS INCOMPLETEとして報告。年や施設を結果の都合で選ばない。
- 検定の研究規模基準は旧版から固定：50 distinct true BIO site、100 site×year、最低5 independent hydrological groups。今回はdistinct/hydrologicalの情報からSOURCE-supportのみを判断し、満たしても**真のoutfall IDと放流量kg/yearの未識別**のため生態モデルの因果認可は別途STOP。

## 公式EA仕様で実装後に見つかった検証上の注意点

- 公式[EA Water Quality Archive移行のお知らせ（2025-11-27）](https://dsp.environment.data.gov.uk/support/announcements/275811447/1156939777)によると、**/data/observation の最大ページ長は250**、一回の検索のsampling point上限100、ページ送りは skip を使う。旧 source-only 実装で limit=2500 を要求し len(rows)<2500 で完了扱いした部分には、**APIが250に制限した結果を全件と誤認する潜在的な欠落**があった。正しい limit=250 で全sourceバッチを再実行して初めて採用する。未修正版のActions runはSTOP/診断履歴として残し、過小なsource observation countは用いない。
- WQE SAの**採水地点**と元permitの物理的な**放流口**が別な可能性は行政要件にもある。[GOV.UKの排出監視ガイダンス](https://www.gov.uk/government/publications/waste-water-treatment-works-treatment-monitoring-and-compliance-limits/waste-water-treatment-works-treatment-monitoring-and-compliance-limits)では複数outletとcombined sample pointがあり得る。[事業許可の管理システム](https://www.gov.uk/guidance/develop-a-management-system-environmental-permits?lang=en)も排出monitoring pointsとsurface-water outletを別項目として扱う。流下の候補を出しただけで実際の物理outfall→河川の受水と認めない。
- 生物観測年と処理水のTP測定年を別々に確認する。2017のsourceを2018のBIO地点に後付けで付与しない。
## 繰り返してはいけない表現

HydroRIVERS reachが上流にあるだけでは、WQE SA「採水点」が原処理場の物理放流口と一致すると保証されない。HydroRIVERS v1は10km²/0.1cmsより小さい流路を落としており、支流/感潮/施設敷地内の水路/下流方向の不確実性は残る。BWMPの見かけの地理傾向をSTW由来の因果N/P選別に格上げしない。処理水TP mg/Lは質量kg/yrの代わりにならない。下流比較のplaceboと同reach曖昧、別branch=陰性と決めつけない。

**前の2kmパイロット79site×yearを削除しない。対象とするsource universeが変わった別の感度設定として併記する。** 過去の陽性探索後の新source選択のため確認的推論ではない。mainにはマージしない。
