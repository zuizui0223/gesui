# 2026-11-19｜フランスROSEAU/SILLAGE汚泥施用経路C：元サイトとCSV取得監査（結果）

**現時点の結論：公開輸出フォームと2017–22年の年選択・公式POST処理先は実在する。しかし2018年の元CSVの行数・真正施用自治体数は未確認。植物・蝶の生態効果も未計算。**

## 確認済みの公的事実
公式 https://www.assainissement.developpement-durable.gouv.fr/pages/data/basededonneeboue.php はROSEAU/SILLAGEの「boues urbaines épandues」のCSV/GeoJSON輸出を案内。元汚泥の生産乾燥重量・処理方式・品質・施用先自治体（communes d'épandage）を公開対象に挙げ、2,000人口当量未満の処理場は信頼性が低いと注意する。

[GitHub Actions元GETフォーム監査 #44455724925](https://github.com/zuizui0223/gesui/actions/runs/44455724925)で政府の3公式ホスト全てHTTP200。元年のselectには2015–2024があり、**元対象2017,2018,2019,2020,2021,2022の6年間すべてを選択できた**。bassin・tranche d'obligation・formatは必須（CSV/GeoJSON選択あり）。これは年次の選択肢であり、各年に実際の汚泥施用が記録されていることを保証しない。

政府自身の公開JS [exportbasededonneeBoue_min.js](https://www.assainissement.developpement-durable.gouv.fr/public/js/exportbasededonneeBoue_min.js) は元SHA256 f4cf1e2fe84e0ccf608ece6c5d363426960e2b730e21900d48265eba4a78f0e5。元JSのgetDataExport関数に、政府正式のpublic POST /controllers/GetDataExportBoueHandler.php、フォームキー bassin/region/annee/format/tranche/departement、CSVファイルとして保存するロジックがある。原公開フォームのALL指定を用いる以外のパラメータ推測・権限回避は行わない。

## 2018年実CSV取得の最初の試行
[2018年元CSVの取得・行構成監査 #44456283313](https://github.com/zuizui0223/gesui/actions/runs/44456283313)はプログラムの単体試験を通過したが、政府のwwwホストへのHTTPS接続がタイムアウトしたため **SOURCE_ACCESS_STOP**。このSTOPは欠測値・施用ゼロではない。元CSV未取得のため、施用先INSEEコード、実施用量や2018適格データの存在は**未確認**。元URLとJSの検証は成功済みであり、アクセス問題と源記録がないことを混同しない。

別の公式ホスト（wwwなし、www2）についても前のGETフォーム監査でHTTP200を確認済み。JSを同一元SHA・同じpublic controllerで照合し、**2018年のまま**再試行する最終SOURCE-onlyコードとワークフローをPR129に追加。これでも実CSVを取得できなければACCESS_STOPを保存して、生物のH7を公開元データで照合できたとは言わない。

## 汚泥が本当に植物のいる場所に届くために必要な項目
- 発生元処理場コード（STEU）、汚泥の発生量と処理・販売先
- 施用先の自治体コードINSEE（処理場の位置ではない）
- 実際の施用年・施用量・肥料化後製品の種別と窒素率
- 計画と実績の区別、処理場2000人口当量以上などの収録品質、畑の具体的な圃場境界が一般公開されるかの判定

自治体への施用実績が得られたとしても、H7に住む蝶の寄主植物への曝露ではなく**面積重複からの曝露機会**にとどまる。台帳の未記載自治体を「無汚泥」と解釈しない。元の水陸2×2のSOURCE-STOPは継続、DAYCENTの一般有機施肥Nモデルを下水汚泥施用と呼び換えない。

**main、草稿PR124、既存点面源のPR125/126/128は未変更。**
