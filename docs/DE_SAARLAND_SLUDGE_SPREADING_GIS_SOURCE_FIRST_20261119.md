# 2026-11-19｜フランスFeursの1処理場事例からのマクロ拡張：ドイツ・ザールラント州の公開汚泥施用GIS SOURCE gate

## 一次資料から確認済み

ドイツ公式公開データGovDataに[Klärschlammverfahren SL (Klärschlammrichtlinie 86/278/EG)](https://data.gov.de/suche/daten/klarschlammverfahren-sl-klarschlammrichtlinie-86-278-eg)が掲載され、「汚泥が施用された農地の面」を図化した GIS と説明。公開開始2019年、対象地域ザールラント州。公式[MetaVer SLの該当OGC API Featuresレコード](https://md.umwelt.info/dataset/metaver-sl/222f0723-3436-aaca-5f7d-778779f1d407)には地域bbox約経度6.355–7.379、緯度49.111–49.623、ライセンスCC BY 4.0、次の公式データ提供先がある。

- https://geoportal.saarland.de/spatial-objects/321/collections/Abfallbeseitigung:Klaerschlamm_Verfahren
- https://geoportal.saarland.de/arcgis/services/Internet/Abfallbeseitigung/MapServer/WFSServer?request=GetCapabilities&VERSION=2.0.0&SERVICE=WFS

これは「実際の農地施用を表す地理レイヤーかもしれない」という**SOURCE候補**であり、2017–22年の実施用台帳だと確定したわけではない。公式ポータルの公開開始年と個々の施用年は別。GISの特徴量/ポリゴン数・時点・汚泥量・実際の受入処理場数をまだ報告しない。

## SOURCE-only 検定順序（生物の結果を見ない）

1. URLは上記公式ホスト固定。**メタデータ2点だけGET**、status・bytes・SHA256・WFS FeatureType名とタイトルを受領。認証回避・推測APIの総当たりなし。サイトエラーはACCESS-STOPであり「施用ゼロ」と解釈しない。
2. 実際のレイヤー名が公式WFS XMLで確認できた場合だけ、DescribeFeatureType等で原属性schemaを監査。**年次列、実施用対計画／許可、STEUに結びつくID、施用N/Pの量と単位、ジオメトリ、ライセンス**をSOURCEだけで判定。
3. 年2017–2019および2020–2022の両期間に実際の施用を支持できるか、元STEUまたは独立な施設/キャンペーンが識別できるかを検証。「施用の承認／利用可能区域」だけなら生物曝露と区別してSOURCE-STOP。
4. SOURCEの実圃場・排出河川口が両方識別された場合のみ同年流域H6/H7と合わせる。元フランス研究と同じH7名だけをドイツへ適用するのでなく、真正な対応流域・期間・国別の評価法と地理独立性を監査。報告されていない区画を未施用対照に使わない。
5. **水域EQR / 陸域NUTRI・蝶HNSをこのSOURCE監査では読まない**。独立H6≥10、実施用支持H7≥20、実施用自治体≥30など既存の閾値を勝手に減らさない。生物モデルと因果的経路はSTOP。

## 研究の意義と水域との接点

Feurs 2018の実施用は**43台帳行/8自治体/6か月**に記録されたが、**独立な元処理場×年は1つだけ**。スケールの異なる各結果を関連と因果の同義として後付け統合するのではなく、同じ処理場を共通起点にして**処理水→河川（A）と汚泥→農地（C）**を別々の source edge としてつなぎ、年・月・真の施用点・N/Pの実測質量がそろうかを検証する。

ドイツ公開レイヤーが真正な実施用ポリゴンと年次を持てば、単一施設に偏らない空間・時間反復の可能性が生じるが、植物や蝶に下水栄養が移行したことにはならない。

## 実装・監査範囲

- `analysis/20261119_de_sludge_saarland/sl_saarland_official_ogc_SOURCE_probe.py`: 公式metadataのGET（非GETFeature）
- `tests/test_de_saarland_ogc_sludge_source_probe.py`: 4ケースのoffline fail-closed tests
- `.github/workflows/gesui-de-saarland-sludge-ogc-metadata-only.yml`: 公式メタデータのaudit JSONだけ保存。元地理レイヤーを無確認で複製しない。

現段階の判定：**公開GIS metadata candidate VERIFIED / layer schema・施用年・処理場接続・生物曝露 NOT VERIFIED / biotic model STOP**。元のフランスSOURCE解析・草稿PR124・mainは変更せず、別PRで調査する。
