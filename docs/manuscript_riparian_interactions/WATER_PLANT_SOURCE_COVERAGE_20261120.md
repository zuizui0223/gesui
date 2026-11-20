# 水・植物接続の原資料監査（2026-11-20）

新しい統計解析・効果量計算は実施していない。公開原著・補足資料の取得と観測単位の確認を行った。取得失敗はデータ不存在と区別する。

## 結論

**Allen et al. 2025 は同じ処理湿地の水と河岸植物を追う有力な接続資料。ただし公開の行データは未取得で、論文はデータを要請に応じ提供すると明記する。Hall et al. 2015 は広域の河岸植物同位体資料だが、水質との同時対応表ではない。** どちらにも植食者の観測は確認できず、全連鎖の実証資料として登録してはならない。

## 1. Allen et al. 2025：オーストラリア南部

[原著 DOI: 10.1016/j.scitotenv.2025.178816](https://doi.org/10.1016/j.scitotenv.2025.178816)、Science of the Total Environment 969:178816、巻の発行日2025-03-15。著者は Allen, Farrell, Wong, Huang, Mosley。CC BY published version を[Adelaide公開PDF](https://digital.library.adelaide.edu.au/server/api/core/bitstreams/b6852c43-b530-4c5c-b481-80e603af7cdc/content)から取得・本文抽出して確認した（方法pp.2–4、制約・Data availability p.8付近）。

| 項目 | 原著で確認した内容／限界 |
|---|---|
| 排出源・接続 | Mount Barker WWTP → Laratinga constructed wetland → 満水時 Mount Barker Creek。湿地内の植物接続と自然河川の植物接続を同一視しない |
| 時間 | 2022年8月、冬季。本文は2回の採取イベント前の雨量を記すが、正確な日付・各試料へのイベント割当は未確認 |
| 水・堆積物 | B0, B1, B2, B3a, B3b, Creek上流US・下流DSの7地点。各地点4反復を方法に記載。図説明ではB3aが2点の場合もあり、予定反復数を全変数の有効nに流用しない |
| 植物・土壌 | 各basinの4場所。岸から3m以内をriparian、3m超をterrestrial。新葉を採取。B3aは植物・土壌を採らない。US/DSの植物も揃うとは確認できない |
| 接続単位 | basin/場所・媒体。4場所は施設4反復ではない。独立した個体4株か、混合試料か、同株の複数葉かは本文から確定しない。植物種名と種を対応させた比較も未確認 |
| 観測媒体 | 水のNH4-N/NO3-Nと同位体、懸濁物、堆積物、河岸/陸域の植物・土壌同位体。植物δ15Nを葉窒素量・栄養価・食草供給量の代用にしない |
| 独立性 | 1施設・1湿地の連続basin。上流から下流の処理・混合・分別と水理条件が変化。複数の独立下水影響研究に分割しない |

### ファイルとスキーマの到達点

- [機関レコード](https://digital.library.adelaide.edu.au/items/6348fa0a-b954-4d7f-95a1-a94eee108cb0)の通常HTMLは今回の直接取得で404となったが、[bundle API](https://digital.library.adelaide.edu.au/server/api/core/bundles/45836238-0daa-4668-b2ee-797d40fc9b2e/bitstreams)経由で原著PDF実体を取得できた。
- ファイル名：`9916951434401831_12299344970001831_Isotopically resolved fate and processes of nitrogen in a constructed.pdf`。このORIGINAL bundleで確認できたのは原著PDFであり、行データではない。
- 本文が参照する補足の内容：Table S1先行研究、S2分析項目、S3分析法、S5各媒体の集計（植物はS5d/e）、Fig.S1位置など。**これらは本文中の参照による内容把握で、補足実ファイルのスキーマ検査完了を意味しない。**
- [出版社ページ](https://www.sciencedirect.com/science/article/pii/S0048969725004516)は今回403。[Elsevier XML API](https://api.elsevier.com/content/article/PII:S0048969725004516?httpAccept=text/xml)は200だが2065文字の書誌情報のみで、補足リンクを含まなかった。架空のmmcファイル名は登録しない。
- 原著のData availabilityは “Data will be made available on request.”。機関レコードの「Supplementary data」表示を公開rawデータと誤認しない。今回、観測行・列名・欠測コード・植物種・個体ID・採取イベントIDは取得していない。

**利用判定：** 水→河岸植物の同位体接続の機構的ケースとして採用候補。湿地と自然河川を分層する。原データ取得前は、個体レベル媒介解析や全段階共分散の推定には投入不可。原データ依頼時は地点・日付・媒体・個体/混合単位・植物種・濃度/同位体・反復対応表を必要項目とする（連絡は未送信）。

## 2. Hall et al. 2015：米国Utah

[原著 DOI: 10.1890/ES15-00319.1](https://doi.org/10.1890/ES15-00319.1)、Ecosphere 6(10):200、2015-10-28公開。[著者公開PDF](https://www.ehleringer.net/uploads/3/1/8/3/31835701/467.pdf)を取得して方法・Table 1・補足へのリンクを確認した。

| 項目 | 確認内容／限界 |
|---|---|
| 階層 | Salt Lake・Cache・Heberの3流域、20河川等、255地点、819葉試料。819独立地点ではない |
| 採取単位 | 約30mの河岸トランセクト。水路2m以内または水中に発根した植物。種ごと5–15葉、複数個体があれば最大5個体を混合するため、行は個体ではなく地点×種の混合試料に相当 |
| 時間 | Salt Lakeの多くは2013年9月末–10月初め、Cache/Heberは2014年7月末。Salt Lake追加・再採取も2014年7月。著者は時期を統合したが、再利用では同地点再訪の識別が必要 |
| 植物 | 11種、葉δ15N・δ13C・C/N等。N固定種は対象外。Phragmites australisとTypha latifoliaは著者の後続空間比較から除外。全819試料を最終空間解析nと解釈しない |
| 水 | 同時水サンプルの濃度・同位体表は本文方法にも公開補足にも確認しなかった。植物に現れたN源・過程の指標であり、河川水濃度→葉への対応傾斜ではない |
| 下水源 | 補足はJordan Riverへの4 POTWs、Mill Creek末端の処理場、Provo農業区間のPOTWを記載。Emigration等の浄化槽、農業・動物・交通・下水漏出も存在。施設名/放流量/処理法の対照表ではない |

### 補足実体の検査

[補足DOI](https://doi.org/10.1890/ES15-00319.1.sm)は[Figshare collection 3308607](https://figshare.com/collections/3308607)へ転送する。[collection API](https://api.figshare.com/v2/collections/3308607/articles)で1記事を確認し、[article API 3564279](https://api.figshare.com/v2/articles/3564279)から以下2ファイルを取得して内容を確認した。補足版 DOI は [10.6084/m9.figshare.3564279.v1](https://doi.org/10.6084/m9.figshare.3564279.v1)、Figshare掲載2016-08-10。採取年と保管庫掲載年を分ける。

| ファイル | bytes | MD5 | 内容 |
|---|---:|---|---|
| [appendix-A.htm](https://ndownloader.figshare.com/files/5638659) | 3381 | b288769318701e598bdfb25e06dde41d | 補足の案内HTML |
| [appendix-A.pdf](https://ndownloader.figshare.com/files/5638665) | 390004 | 92e4d9ee6195c1e26a5fd69403cd2903 | 7ページ、Table A1と追加図 |

Table A1列は watershed、stream/river、watershed area(km²)、traffic density、montane/suburban reach comments、urban/agricultural reach comments。A1–A3図は植物δ13C・標高分布等。**この公開補足に819試料の生物行データCSV/Excelはない。** Figshareのdatasetという分類だけでrawと判定しない。植物種一覧は原著Table 1にあり、同種対応確認には利用できるが、葉化学の地点行を復元しない。

**利用判定：** 混合起源の人為N→河岸植物の広域比較・異質性を示す候補。下水特異効果、水→植物の同一地点媒介、植物→植食者としては未適格。追加公開rawの存在自体を否定せず、今回確認した補足にはない、と記録する。

## 統合時のゲート

1. Allenのwetland basinとHallの河川地点を同じ独立単位としてプールしない。施設/研究・水文単位・植物種・採取時期を保つ。
2. 水・植物同位体の並行変化、葉窒素/C:N、植物量、植食量は別アウトカム。δ15N上昇を食資源増加・植食者増加と読み替えない。
3. Allen2022/2023は同じLaratinga系の先行研究として原著が引用する。別論文であることだけで独立研究数を増やさない。
4. この2研究から植食者・蝶・成虫タンパク利用へ実線を延長しない。資料探索の到達点は接続可能性と欠落の確定である。
