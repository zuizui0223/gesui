# 2026-11-21｜英国PRTR真正原CSVを2017/18原報告年で統合：大規模下水処理場のP/N年間報告量が同じ施設で揃わない

**原source-only監査完了**：公式CSV取得・source schema run [#44847278163](https://github.com/zuizui0223/gesui/actions/runs/44847278163) SUCCESS、実際のReportingYearで2017/18の都市下水処理場のwater-release報告を数えた [#44847500143](https://github.com/zuizui0223/gesui/actions/runs/44847500143) SUCCESS、管轄当局/施設地理を追加確認した [#44847721292](https://github.com/zuizui0223/gesui/actions/runs/44847721292) SUCCESS。採用原集約JSON artifact ID **14684330631**。生物BMWPの係数・施設間のP/N濃度とその影響・個別河川下水源の因果影響は**計算していない**。

## 1. 新しい2017/18英国の年間質量sourceが取得できた

- [DEFRA UK Pollutant Release and Transfer Register 2017、2018、2019年元CSV](https://www.gov.uk/government/publications/download-industrial-emissions-and-waste-transfer-uk-prtr-data) は正式 GOV.UK assets から全てHTTP200で取得成功。各元CSVのsha256・元schema・当該年の真のReportingYear・施設ID・活動code・経路・汚染物質・方法code・TotalQuantityを監査。元のfacility-level値、施設ID、住所や原座標の行は公Gitへ再配布していない。
- 元PRTRはUK全域の施設・permit単位の**年間排出**。元E-PRTR定義でwaterへの`TotalQuantity`は **kg/year（事故放出量を含む）**。`Quantity`は別の施設外移送／廃棄物量でありその代用としない。`MethodBasisCode`が **M=measured、C=calculated、E=estimated** なので、全量を「実測P/N」と呼ばない。総水量や源別河川への到達フラックスとは異なる。[EEA E-PRTR metadata](https://www.eea.europa.eu/data-and-maps/data/member-states-reporting-art-7-under-the-european-pollutant-release-and-transfer-register-e-prtr-regulation-23/database-structure-and-use-information/eprtr_database_metadata_v11.pdf-1/file)、[GOV.UK PRTR方法](https://www.gov.uk/guidance/uk-pollutant-release-and-transfer-register-prtr-data-sets)。
- 元の都市下水処理場は`AnnexIActivityCode=5.(f)`かつ`NACEMainEconomicActivityCode=37.00`（Sewerage）の**両方**を主対象とする。5.(g)は独立した工業排水処理場なので別原母集団。元sourceの総窒素を`TOTAL - NITROGEN`、総リンを`TOTAL - PHOSPHORUS`、排出媒体`MediumCode=WATER`で厳密に区別。アンモニア`NH3`、総有機炭素TOC、移送先下水道、AIR等を総リン/総窒素に含めない。都市下水処理場主活動の報告閾値は100,000 population equivalentsで、ここに載らない小処理場を水域N/Pゼロとしない。

## 2. 重要な原年度の齟齬：ファイル名で年を割り当てない

| 正式CSVファイル名 | 収録行数 | 内部ReportingYearの多数派 | 内部2017年行 | 内部2018年行 |
|---|---:|---|---:|---:|
| **uk_prtr_dataset_2017.csv** | **6,258** | 2016年（5,012行） | 912 | 0 |
| **uk_prtr_dataset_2018.csv** | **5,965** | 2017年（5,478行） | 5,478 | 274 |
| **uk_prtr_dataset_2019.csv** | 原sourceで監査済 | 2018年（5,059行） | 81 | 5,059 |

公式のファイル名の年をそのまま測定年に使うと2017年と2018年の生物試料・施設源量を誤接続する。対象年の定義は**CSV内部のReportingYearだけ**とした。2017/18のwater-release P/N候補の原キーは`(ReportingYear,NationalID,PollutantCode,MediumCode,Designation)`で監査。元3ファイル間でこの同一キーに**6件の分類・年間値などの矛盾**があり、該当原キーは集計から停止・除外してsource conflictとして残した。反対に、他ファイルに値が無い施設へ勝手に補完しない。PRTRの一部施設が漏れていても0放流ではない。

## 3. 大規模都市下水処理場：水への年間P/N量の原sourceカバレッジ

母集団：`AnnexI 5.(f)`と`NACE 37.00`が一致、`MediumCode WATER`、元`TotalQuantity>0`、対象年は真正`ReportingYear`。

| 元施設レベルの水域放出報告 | ReportingYear 2017 | ReportingYear 2018 |
|---|---:|---:|
| **TOTAL - NITROGEN の正の年間報告を持つ異なる施設** | **124** | **118** |
| うち窒素の報告がM（測定）方式の原行 | 42 | 42 |
| 窒素の報告がC（計算）方式の原行 | 82 | 65 |
| 窒素の報告がE（推計）方式の原行 | 0 | 12 |
| **TOTAL - PHOSPHORUS の正の年間報告を持つ異なる施設** | **0** | **8** |
| うち総リンM方式の原行 | 0 | 8 |
| **同一施設・同一報告年に正のN/P両方がそろう施設** | **0** | **0** |

なお2018年の総窒素は119元原行で118 unique NationalID、元MethodBasisの42+65+12は原行数119であって118施設を分類したことではない。年の同一施設の二重行を質量として加算しない。2017年の原窒素124施設は全124行。既存2018 EEA UWWTDの両物質実測正値＋一意淡水出口の集合に英国報告元が無かったことと独立な国別/制度別資料である。

### 施設活動コードの緩い代替の確認

`AnnexI=5.(f)`だがNACE37.00ではない下水処理施設候補の水域総窒素は2017年11件・2018年12件、総リンは2017年1件（C）・2018年0件。同じ施設で両方正値がそろう候補は依然0。独立工業廃水処理 `5.(g)` は2017/18に総窒素M報告各1件、総リン0で主母集団へ混ぜない。

### イングランドとの空間的互換性を厳格化

PRTRの対象は英国全域。元`CountryID`はこれらのsource行では空欄であり、**英国全域の124/118をイングランドChemPopの施設数と呼ばない**。主厳格都市下水施設の`CompetentAuthorityPartyName`は、2017窒素**EA 103／NRW 7／SEPA14**、2018窒素**EA98／NRW7／SEPA14**（2018は元119行でEA98、NRW7、SEPA14）。総リン2018年8施設はすべて**EA管轄**、元NUTS等でも英国内の場所を確認。もともと報告したPRTRの位置座標は**施設の位置**であり物理放流口と同値ではない。

## 4. 厳密な論文・モデル解釈の上限と次段のsource GO/NO-GO

- **取得できた肯定的な事実**：英国大規模処理場には施設年度の水域への年間総窒素massのsourceが多く、2017/18の同年報告と元方法（M/C/E）が区別可能。英国正式EA監視原TP mg/Lが2017/18に存在することも別SOURCE監査で確認。水源に沿ってBIOを探索できる事後的マクロ研究への**データ取得機会**は成立する。
- **停止すべき検定**：同一施設の総窒素・総リンmassを用いたP/N比・通年源別曝露・施設→受水化学→底生BMWPの因果媒介は、正のP/N両年両物質報告施設が**0**なのでSOURCE STOP。PRTRなしの施設を排出0として埋めない。元PRTRの施設座標は正式河川放流口ではなく、SAの実測採水点もpermit-outfallと一致する保証がない。真の下流輸送を直接実証したと主張できない。
- **検討価値のある独立ルート**：Environment Agencyの[Pollution Inventory](https://www.data.gov.uk/dataset/cfd94301-a2f2-48a2-9915-e477ca6d8b7e/pollution-inventory)は別の年次大規模sourceがあり得るが、公開原データは**EA Conditional Licence等の利用制約**が示される。利用権/原sourceを取得できるまで原個票を公Gitへ配布せず、異なるレポート制度を同じ年度outfall massとして継ぎ合わせない。[英国2026 Water B1c公式統計](https://www.gov.uk/government/publications/environmental-indicator-framework-theme-b-water/environmental-indicator-framework-theme-b-water)は全処理場の河川への推計P/アンモニア負荷を示すが**5年ごとの全国集計**であり、2017/18個別施設の値や日付に戻せない。
- **最小の今後の原source検証**：正式EA SA samplingPoint.nativeCodeをPRTR NationalID施設の物理outfall/permitに独立リンクさせる。元PRTR facility geomとの単純最近傍は監査候補の数にすぎない。PRTR年次NはM/C/Eが混在し、源のNなら河川中Nitrate/AmmoniacalNとは異なる化学種・期間である。原sourceを条件付けた観測モデルは主既存のChemPop EDFと河川化学の負・正の両結果を保存した後、別地域外留保/逆向きplacebo/農業交絡で事前定義し直す必要がある。

今回の原資料はsource-only PR #137（PR #136の上にDraft）として提出・保持し、gesui/main・他の解析PRには変更を加えない。
