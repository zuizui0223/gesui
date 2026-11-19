# 2026-11-19｜ドイツの下水汚泥「陸域への資源分配」の年次変動：2016–2022年の独立公表統計

## 結論（SOURCE-only、非因果）

全国の農業向け汚泥乾物は2017年 **306,230 t → 2022年 233,065 t（−23.892%）**と減少。一方、ザールラント州由来は **2,140 t → 3,502 t（+63.645%）**と増加した。全国はザールラント州を含むため、独立比較として州と**州を除いた残りのドイツ**を分けると、農業処理比率の2017→22年変化は、州 **+8.238ポイント**、残り **−4.135ポイント**、差 **+12.374ポイント**。これは**行政ソースの資源処分経路が異質で時間変化する**という記述的証拠に限られ、施用農地の所在・生物曝露・河川放流量ではない。

## 公式元表と粒度

- **全国**：ドイツ連邦統計局Destatis [GENESIS-Online 32214-0001](https://genesis.destatis.de/datenbank/online/statistic/32214/table/32214-0001)「Entsorgter Klärschlamm ... Deutschland, Jahre, Entsorgungsarten」。2016–22年の農地利用、景観利用、その他物質利用、単独・混合・不詳焼却、その他の乾物質量を転記。**最新版公開表の年次修正値**を用い、2017年の旧年鑑値（311,905 t）と混在させない。
- **ザールラント州**：州統計局発行 [Statistisches Jahrbuch Saarland 2023](https://www.statistischebibliothek.de/mir/servlets/MCRFileNodeServlet/SLHeft_derivate_00012170/SL%20Jahrbuch%202023.pdf) 印刷p.143、表「Klärschlammentsorgung aus der biologischen Abwasserbehandlung 2016 bis 2022」、農業、景観、熱処理、その他の列を年別転記。ウェブ上のテキスト索引を取得したが**PDFページ画像の再取得は403で不可**、バイナリSHA未確保。旧国家統計と異なる版の公表なので厳密な系列としての統一・最新改訂は残課題。
- **単位は乾物 t/年**。州の計数は基本的に**処理場所在地**に帰属。土壌への施用先農地の実所在・施用量・N/P濃度・放流口と対応する施設IDを示すものではない。

## 取得した現時点の年度別 SOURCE 系列

| 年 | ドイツ全国・農業利用 t乾物 | ザールラント州・農業利用 t乾物 | 全国の農業処理割合 | ザールラント州割合 |
|---|---:|---:|---:|---:|
| 2016 | 447,922 | 5,089 | 24.9% | 26.8% |
| 2017 | 306,230 | 2,140 | 17.9% | 11.5% |
| 2018 | 277,601 | 3,883 | 15.9% | 21.6% |
| 2019 | 277,278 | 4,538 | 16.0% | 25.0% |
| 2020 | 263,161 | 3,776 | 15.1% | 20.2% |
| 2021 | 226,747 | 2,658 | 13.2% | 14.6% |
| 2022 | 233,065 | 3,502 | 13.9% | 19.7% |

ドイツ全国の**熱処理総割合**（単独＋混焼＋不詳）は2017年約69.7%から2022年約80.0%に増えた。これは汚泥の別の処分経路への振替であり、**河川に排出されるN/Pが増えたことを意味しない**。

## 事前の地理・独立性・生態的証拠の規則

1. 集計上の**ドイツ全国はザールラントを包含する**。全国対州の単純な独立2標本検定をしない。比較グループが必要なら「ザールラント vs 全国−ザールラント」で記述。
2. 2016–22年の7年は独立した7政策実験ではない。特定の政策効果、置換因果、相関の有意検定、または複数州への普遍性はここでは主張しない。2017年の制度変更の影響は別の反事実と他州を必要とする。
3. *統計上の農地利用*と個別**施用済み圃場の位置・日時・量**は違う。今探索中の[Saarland公式WFS Klärschlamm_Verfahren](https://geoportal.saarland.de/arcgis/services/Internet/Abfallbeseitigung/MapServer/WFSServer?request=GetCapabilities&VERSION=2.0.0&SERVICE=WFS)には真正なレイヤー名が存在するが、施用年・実績か計画・処理場IDの属性内容は未確認。
4. 同じ施設・年の**水側の実測排水TN/TP年間負荷**と**陸側の汚泥由来N/P施用量**がそろうまでは、二経路のN/P配分率、植物NUTRI、蝶HNS、水域EQRの因果的媒介推定はSTOP。2018年Feurs単一施設の実施用43台帳行と混ぜない。

## 再現性

- `analysis/20261119_de_sludge_saarland/DE_GERMANY_SAARLAND_2016TO22_SLUDGE_ORIGIN_ALLOCATION_SOURCE.csv`：7年間・二地域の公的表転記。
- `analysis/20261119_de_sludge_saarland/de_sludge_origin_2016to22_SOURCE_trend.py`：元処理量の全経路合計チェック、年次欠測禁止、その他ドイツとの差分による重複排除。
- `tests/test_de_sludge_origin_2016to22_source_trend.py`：8つのfail-closed unit tests。
- `.github/workflows/gesui-de-official-sludge-2016to22-origin-trend.yml`：GitHub Actionsで独立再計算と結果JSONの保存。

**判定：2016–2022年の発生処理場側の汚泥処理経路の年次変動 SOURCE=GO。圃場曝露の地理・水側同年N/P・植物や昆虫への移行=STOP。** main・原稿PR124・Feurs SOURCE PR129は変更しない。
