# 2026-11-19｜ザールラント州 汚泥施用レイヤー WFS SOURCE-only 実アクセス結果

[公式GDI-DE/GovDataデータセット](https://data.gov.de/suche/daten/klarschlammverfahren-sl-klarschlammrichtlinie-86-278-eg)では下水汚泥の農地施用面が描画されると記述。実際の属性・施用年・処理場・実施用区分はこの宣伝文句だけでは分からない。

## 実ネットワークSOURCE監査

**[GitHub Actions #44522249816 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/44522249816)**。2026-11-19 12:20 UTCの公式ドメインへのGETのみで、元metadataレスポンスを再配布せずSHAで固定した。

| 対象公式URL | HTTP | Content-Type | bytes | SHA256 |
|---|---:|---|---:|---|
| WFS 2.0.0 GetCapabilities | 200 | text/xml;charset=UTF-8 | 15,006 | f6b7ecbabbb2b038c148c45c4ae7ada46b666c2dcf52eb26c90e4d1dc0cf5ea5 |
| OGC Features collection landing | 200 | application/xml | 27,293 | f8974a80aaa16fec44bb949f77173e24b66cbc7910a52140af917baa4010468c |

- **WFS FeatureType は2個**、うち1個は正式な `Abfallbeseitigung:Klaerschlamm_Verfahren`（Title `Klaerschlamm_Verfahren`）。元WFS XMLのNameから完全一致で採用し、名前を推測していない。
- WFS GetCapabilities の段階では**施用した年、実施用／許可済み・予定、処理場ID、乾物/N/P、実圃場属性の存在を確認できない**。他のレイヤー名を流用しない。
- OGC collection landing は `application/xml` だったため、初回監査の単純JSON／HTML二分岐ではレイヤー属性を読めない。WebサーバがXMLを返したのは**エラーではない**が、これを実施用の証拠とはしない。
- 元ポータルの「データ公開開始2019年」「最終メタデータ更新2026年」は**実際の施用年ではない**。
- 実GIS `GetFeature`・個別農地・生物データは**一切取得していない**。匿名化地理であっても、年次と実績の可否を明らかにする前に曝露分類しない。

## 次の段階

公式のWFS `OperationsMetadata` が **DescribeFeatureType** を広告している場合だけ、確認済みName `Abfallbeseitigung:Klaerschlamm_Verfahren` を使い、標準操作で**XSD属性名のみ**取得する。これも元の行政記録が「実施用」だという証明にはならない。フィールド年・処理場ID・施用イベント、N/P質量の定義が欠ければ SOURCE-STOP、代わりに2016–22年のDestatis全国/ザールラント州**発生源側処分経路パネル**を記述的に使用する。

- 今の実行コード：`analysis/20261119_de_sludge_saarland/sl_saarland_official_ogc_SOURCE_probe.py`
- Schema-only次段：`analysis/20261119_de_sludge_saarland/sl_saarland_official_wfs_describe_feature_SOURCE_gate.py`
- 独立原結果：GitHub Actions #44522249816 artifact `gesui-DE-Saarland-sludge-spreading-official-OGC-source-gate` (artifact 14528656644)。

**研究判定：WFS layer exists / realized field-year, actual fields, N/P, STEU identity NOT identified / H6 independent biological replicates 0 / ecological model STOP.**
