# 2026-11-21｜2018年英国UWWTDのリン除去処理とEA最終処理水中の実測総リンの条件付き比較（原source-first）

**事後探索的・処理場単位の観測比較。論文の処理場→河川生態への因果連鎖の最終証明ではない。** Eionet公式[UWWTD UWWTPs data dictionary](https://dd.eionet.europa.eu/datasets/latest/UWWTDArt15/tables/UWWTPs)では原 uwwPRemoval と uwwNRemoval は boolean の処理方式指標（濃度・排出量の実測値ではない）。UWWTPS の原WWTP報告UK5-2018施設1,799のうちPRemoval正の施設418、NRemoval正40。実測N/P annual dischargeはMeasured/Calculated/Estimated 6欄いずれもUK5元原表で全空欄。したがって総リンkg/年や除去率を計算しない。

## 2018年原母集団とソース選択（TPの値を見る前に固定）

1. 原2018 EEA UWWTD UK5（ReportPeriod repReportedPeriod=2018）UWWTPS+DischargePointsの原3科学表の固定SHAを確認。元plant native uwwCodeが一意、原元処理場のdcp outletが**1つだけ**、元outletの受水域FWかつ原英国WGS84座標妥当であるものに限定。PRemovalが明確に0または1であることをsourceで検査。元WWTP施設活動状態・capacity、uwwLoadEnteringUWWTPは元列で保存し、欠測率を先に監査。
2. EEA2018物理FW放流地点からEA WQE Water Company final-treated SA（2026現在registry）をEPSG27700で主**250m以内の相互一意位置**（UWWTD放流点→ただ1SA、SA→ただ1UWWTD FW放流点）だけ採択。前のsource auditではこの位置候補は930あり。100/500mは感度であり結果の符号によって選ばない。これは原permit/nativeID一致の確認ではなく**地理的source-pair候補**。
3. 本当に2018年のEA SAに決定項目0348 Phosphorus Total as P, FINAL SEWAGE EFFLUENT, Water Company SA, original unit MILLIGRAM PER LITRE, nonnegative finite TP measurementがある場合だけ生物化学試料として採択。元同一SA/同暦日で複数独立分析行があれば元source重複・LOQ/数量を監査し、反復数で擬似標本数を膨らませず**SA原採水点1地点の2018年中のTP median**を主応答とする。元TP値を公Gitへ再配布しない。
4. 元source cohortを応答TPの符号・大きさを見る前に固定。処理方式PRemoval1/0が**各30独立SA採水地点以上**、源の元region（uwwNUTSを機械的に地域群として使う場合にはそのコードの実在とカバレッジ）**3以上**、原WWTP capacityの正値が**各群で70%以上**ある場合のみ調整分析へ進む。元の実際の処理場IDがWQE SAコードと一致しない点をsource identification limitationとして明記。

## 固定比較と交絡制御

- Primary: source等重みの2018 year-median effluent total-P mg/L の **PRemoval1 vs PRemoval0** の中央値、中央値比、比の95% site-bootstrap記述区間（2,000 draws）を報告。比較は大規模処理場のsource indexの観測的違いでありP除去による介入効果ではない。
- 同じsource地点において y=log1p(year-median TP mg/L)、PRemoval + log1p(uwwCapacity) + source大域BNG East/North + original uwwNUTS地域群 + native NRemoval を調整する。欠測を0/平均補完せず、train/testingにまたぐ施設重複なし。2018年には元source地点1行なのでsite cluster=siteが標本独立と同値でないことを明記。
- 2017/18の別年再現・別流域留保の転送は元UWWTD sourceで2017の独立報告が確認できないためこの研究の主張へ加えない。地域間外挿が成立しなければ「どの地域でも効果」や「普遍的」を主張しない。原元sourceのbetaや単純中央値が望む方向と逆になっても同じreportに残す。
- 感度として250m主に対し100mと500mのsource-only相互一意状況を並べる。未元確認の施設同一性を強引に広げて肯定を探さない。真の個別管理/permit IDとUWWTD uwwCodeの対応は別sourceゲートであり、相互一意地理source同士の近さでは満たせない。

## 関係づけるものと絶対に主張できないもの

この測定は**処理場にP除去処理が報告された場所では、処理後放流水中の総リン濃度が低いか**という同じ2018年の処理方式と化学指標の対応を検証できる。これはChemPopの英国河川3栄養関連水質→底生BMWPの予測結果とは**別の監視母集団**であり、元同一施設の真のoutfallとBIOがdownstreamで物質輸送することを示さない。処理施設の選択は処理規模・地理・許可の内生性があり、施設P除去の因果効果ではない。TP濃度は年間放出kg/yearではない。元のP/N排出量の欠測は解決せず、源別化学到達や水陸資源循環への因果連鎖も未推定。

原施設・地理・TP元source rowsや生物BMWP値は公開GitHubへ転載しない。source SHAと集約監査のみPR #137 Draftとして残す。gesui/main/PR135/PR136の原結果を変更しない。
