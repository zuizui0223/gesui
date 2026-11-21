# 2026-11-21 英国2018年の正式UWWTD放流地点→EA最終処理水TP採水点→ChemPop底生調査の方向性：source-only確証条件

**事後探索的な源照合であり、現在の本命は源座標をすり替えず正しい過去の河川放流点を使うこと。** 2018年元UWWTD英国報告には1,799下水処理場・1,809原放流地点・座標があるFW（淡水）放流地点1,498。これとは別のWQE水道会社最終処理水SA現行台帳5,722点に対し、元放流地点との距離250m以内で1,032/1498放流地点にSA候補があり、そのうち938地点は採水点がひとつだけ。さらに元UWWTD側も放流地点が一つの施設に限定した932地点が地理的に一意候補となる。ただし**元WQE SA native採水IDとUWWTDの排水口IDが直接一致する証明ではない**。

## 既存sourceを固定し、結果は後から標本調整しない

- **2018出発点**：EEA 2017–18 UWWTD v07r01のoriginal ReportPeriod (`UK5`, reference 2018)、UWWTPS.csv、DischargePoints.csv、既存独立検証の三元SHAを必須。2018年英国の淡水放流点のみ（dcpWaterBodyType FW）、元UWWTD処理場IDとのjoin有効、元UWWTD施設に放流点が一つだけ（当年UK5全植物1787）。元EEAの2018年UK下水処理場に実測N/P共に正値は0＝**mass量の検証不可能**。SOURCE定義には原報告discharge-point lat/lonを使用。
- **水質確認先**：EA WQE API Water Company final-treated type SA、正式登録のsource BNGと元UWWTD放流点と地理比較。主唯一結合は**250m以内で (A)UWWTD放流点が1つのWQE SAとだけ近接 (B)そのSAに250m以内のUWWTD候補放流点が1つだけ**という相互一意条件。施設側放流点も1件でなければ除外し、名前が似ているだけの自動結合をしない。距離100mと500mは後付け感度であり主250mの合否を置換しない。現在のSA登録が2018年に稼働・実測という意味はないので、**EA本当の2018年TP0348 mg/L 原採集記録があるSAだけ**をSOURCEとして採用。欠測は0排出扱いしない。
- **河川ネットワーク**：公式HydroRIVERS EU v1 ZIP SHA256 500da7d36ceee0aa4c82dd625f5c53d80ba2f7ac14ad022cd3a6601bf55068e5、UK元WGS84 discharge pointをEPSG27700へ正しく変換し、唯一の最寄りreach500m内にスナップ。ChemPop2018年の原SOURCE適格366siteも同じBNGで独立スナップ。`NEXT_DOWN` reach距離50km以内の出発OUTLET→下流BIO reachをupstream candidateとして数え、反対方向BIO→OUTLETをdownstream placebo、同reachで真の順序未知はunknownにする。しきい値100/250/500mや10/25/50kmは別の地理的感度であり、元主500m/50kmは維持。
- **生物採集日の元の条件**：同一ChemPop2018年BIO_SAMPLE_DATE以前90日、同2018年の対応SAから本物のTP0348元測定があったかを集計。生物と化学採水が同暦日だけならintraday順序不明として別。前90日と後90日両方の監視記録の比較を必ず保持。正式原BIO採集日が追跡できないならSOURCE STOP。

## SOURCEの解釈上限と完全STOP

- 主SOURCE品質最小は以前のデータ探索と同じ**50 BIO site、100 BIO site×year、5独立hydrological blocks**を基本候補とする。ただし2018の1年だけでは同じBIO site×year=site数なので100地点が必要。新しく該当するhistorical UWWTD physical outlet subsetだけで100サイトに届かなければ元の2017+2018 SOURCE合計の研究規模に届いていても2018源確証ではない。元source監視地点が一施設を代表するとは断定しない。
- 2018 EEA報告のdcp pointは歴史的な**放流地点**であり、EAの現在のSAは施設の**処理水採水地点**。その中間を紐付ける許可IDが一致するかまだ未検証なので、今回の250m相互一意で絞った関係は**地理的・時間的に整合的な候補**でしかない。HydroRIVERSは10km²未満等の小川と工場内の配管/下水路を省略する。無スナップを0曝露にしない。
- 水域へ届いた源P/Nの実測kg/yrや下水からの発生源別化学測定が無い。元UWWTD UKの2018年N/P実測排出欄は空なので、EAのTP mg/LとPRTR施設年質量を足す・補完する・別IDを同一施設として統合することはできない。源時空条件が通っても**下水によるBMWPの因果効果モデルは許可しない**。

全ての元出力は原施設ID/元放流地点ID/元採水ID/元座標/元化学値を公Gitへ再配布せず、元SHAと集約SOURCE品質のみ。既存PR #135/#136、別PR #137（この文書）の結果・source STOPを保持し、mainへの自動マージは禁止。
