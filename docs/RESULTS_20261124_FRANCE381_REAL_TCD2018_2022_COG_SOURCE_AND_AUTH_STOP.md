# 2018–2022年の河畔樹冠変化：France381 genuine EEA × official Copernicus 10m TCD source-only 実測可能性判定

**2026-11-24 JST｜結論：SOURCEの原本・年度・物理的site→タイル・品質レイヤーは全381地点で一致。しかし原COG画素への匿名HTTPアクセスは2018/2022両年403。実樹冠変化・生物との相関はまだ計算していない。** これは源データへの取得権限のSTOPであり、河畔植生変化が存在しなかったことや底生生物への影響がないことを示す生物学的nullではない。

## なぜ必要か・元の実際の反証結果

[前段の新しい時系列結果PR #149](https://github.com/zuizui0223/gesui/pull/149)は、**フランスの本当に同じ381水生生物評価方法のriver monitoring site、14つの2018原3度地理ブロック、固定の5地理fold**で、2018年Riparian Zones土地被覆が2018/2022両年の横断的EQRの空間差を予測する一方、同じsiteの2018年EQRを調整すると2022年EQR予測改善は**−1.540% site／−2.640% 14地域均等 [−6.440,+0.345]**だった。これらの値は探索的だが真正原EEA+Copernicus site source、NumPy独立再計算で再現済み。

**残る本当の問い**：2018〜2022年の河畔土地/樹冠そのものに変化があり、その差が同じsiteの`ΔEQR=EQR2022−EQR2018`に対応するか。2018年固定地図だけから「河畔植生変化が生物の動態を変えた」と言うのは不可能。静的LCLU2018は空間的環境文脈しか表さない。

## 年次地理源の選択と厳密な非同等性

- 公式[Copernicus Riparian Zones](https://land.copernicus.eu/en/products/riparian-zones)は現段階で公開済みの**2012/2018土地被覆・2012–2018変化**が中心。[2026-08-31制作契約](https://land.copernicus.eu/en/news/new-contract-signed-for-future-production-and-consolidation-of-the-riparian-and-coastal-zones-products)が2021、2024、2027 status追加を対象とする。**「RZ2018とRZ2022で同じ56土地クラスの変化を取れる」と捏造しない。**
- 代わりに、全年2018/2022に原本がある**樹冠被覆率TCD**を独立の比較する場合、公式[CLMS Tree Cover Density 2018–present, 10m](https://land.copernicus.eu/en/products/high-resolution-layer-forests-and-tree-cover/tree-cover-density-2018-present-raster-10-m-europe-yearly)をprimary元sourceに固定した。source DOI **10.2909/e677441e-fb94-431c-b4f9-304f10e4dfd8**。TCD0–100は空中から見た**樹冠面積被覆の割合**であり、実落葉量、植物バイオマス、生物種多様性、2018/2022の草地/湿地増減、葉のP/N、成虫蝶の寄主、河川羽化昆虫、下水の影響とは異なる測定量。
- 別の公式[HR-VPP Total Productivity](https://land.copernicus.eu/en/products/vegetation/high-resolution-total-productivity) 2017–年次10mはPPIの季節積算であり、原植物現存量ではない。元primary TCDのsource STOPを「別指標に替えれば有意になる」と変更することはしない。

## 公式原データ・日付・投影タイルを原本で全件監査

[CDSE STAC公式](https://stac.dataspace.copernicus.eu/v1/)は最初の`/collections`レスポンスが10件の**ページネーション**だった。TCDは「STACに存在しない」わけではなく、正しいネイティブ**`clms_vlcc_tree-cover-density_europe_10m_yearly_v1`コレクション**とsource itemが実際にHTTP200で取得できる。最初の一覧だけを見て不存在と誤認しない。

| signed source実行 | 真正データ品質結果 |
|---|---|
| [Official 3 service transport #45204778028 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45204778028) | OData TCD **6,314** items、元COG CSV HEAD200／Range206、STAC collection list一枚目の10件（ページネーション） |
| [Official full TCD catalogue original 3,059,396 bytes #45205351185 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45205351185) | **2018～2024年の各年902 COG tiles**、Europe EPSG:3035の原同一Grid ID **2018 897 /2022 897**、**同IDの両年ペア897**。海外territories他投影は年5枚で別母集団。元CSV SHA256 `5a6f0dc9f67ec7e3671f5984e9cf0425988e23dbb59ce1f6a83a18dc8661d6c9` |
| [STAC exact original 2018+2022 COG asset header #45206093035 SUCCESS SOURCE QA but DATA ACCESS STOP](https://github.com/zuizui0223/gesui/actions/runs/45206093035) | フランスのネイティブ確定タイル `E36N28 EPSG3035`で両年Source真正itemと`data: image/tiff; ...cloud-optimized`資産が登録。**2018＆2022とも素の原TIFFへのByte-range HTTPアクセスが403**。source STAC `data`は`auth:refs=["s3"]`、元zip productは`auth:refs=["oidc"]`を明示。元の画素値は**一切取得していない**。 |
| [2018/2022真正France EEA source exact381 site→TCD/TCDCL年次タイル #45206831800 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45206831800) | source RZ2018×EEA2018 EQR対象ID×EEA2022 original same site/method/reach：再訪 source 418→安定409→France **381 site**／source5fold 122/58/48/33/120を生物EQR数値未使用で復元。公式年度TCDとTCDCL品質両方とも381/381地点の±500m範囲で**同EPSG3035原タイルに原商品登録あり**。必要な100kmタイルは**47**、うち374地点は1 tileで500m圏全域を覆い7地点では境界の両側2 tileが必要。 |

TCD full CSV source metadataはフランスの真2018/2022原タイル2018 EPSG:3035 897/2022 897を保持。品質補助`Tree Cover Density Confidence Layer`の元CSVは2018 897 /2022 897、native SHA256 **`347053a15fe8e3cf5df7fa261daf7ae0cf9dcab4cf1fe5f81237a28f735b4222`**、3,191,375バイト。**quality layerの2018/2022画素も未取得**。目視国勢地図や色づけWMSサムネイルをそのまま生TCD percentに変換しない（定量ソースの厳格性を維持する）。

元2018 EEA river source SHA256 `2a1a72d43656aca588a832dfae13860054c01a651e76bf25f984c5a374709b02`、original 2022 EEA Git blob `325f199e985d9862ffee5b0edb205c61a7be3833`を確認。Copernicus RZ2018 original=2,639 source river candidate中978真正地図適格、同じ2018数値EQR原siteは976、2022再訪418、同一評価手法・同一HYRIV・座標一致409、FR381。対象を**TCDの値やEQR変化の符号で選んでいない**。site IDs/元coords、生物生データ、元GeoTIFFはpublicリポに追加していない。

### SOURCE GO / STOP の結論

| 検証した主張 | 科学的判定 |
|---|---|
| 2018/2022両年10m TCDが出版社のデータとして存在し、同じ欧州100kmタイル単位で対応 | **SOURCE_META_GO** |
| 真のFrance381 EEA再調査地点で2018/2022両年TCDと公式Confidence layerのネイティブカタログ領域被覆 | **SOURCE_GEOGRAPHY_GO 381/381** |
| 両年の真正TCD画素値とConfidential/Quality画素を使用したsource測定 | **SOURCE_PIXEL_STOP (HTTP403両年)** |
| 381siteの`ΔTCD2018→2022`近岸・外側の樹冠変化の定量値 | **NOT_MEASURED** |
| `ΔTCD近岸−外側`追加による`ΔEQR`の空間分割予測の検証 | **NOT_FITTED/SOURCE_STOP** |
| 河畔の本当の2018→2022植生変化から底生生物動態へ因果的に伝播した | **NOT_IDENTIFIED** |
| 下水由来N/Pが河畔の樹冠・河川EQRや生物間相互作用を変えた | **NOT_IDENTIFIED** |

**次段のreal sourceピクセル取得条件**：Copernicus Data Spaceの本当のS3/OIDC認証手段を用意し、2018+2022 exact source COGと**同年のTCD Confidence Layer**から、本当の0–250m/250–500m 10m画素を読み込み、NDV=255などnative no-data/品質・2018–2022時差に伴う季節サンプリング差を監査する。画像の画素値は元ソースから、元フランス381サイトの正確なriver source座標に対して出す。**original pixel source GOまでは2022 EQRやΔEQRに対する新たなモデル係数を表示しない。** 途中まで進んだ公式meta/geometry source GOは、完全な生物学的分析実行GOとは違う。

### 照合コード・認証を要するscope

1. `scripts/gesui_france381_CLMS_TCD2018_2022_catalogue_source_probe_20261124.py`
2. `scripts/gesui_CLMS_TCD2018_2022_full_provider_COG_inventory_20261124.py`
3. `scripts/gesui_france_CLMS_TCD_original_2year_STAC_COG_header_20261124.py`
4. `scripts/gesui_FR381_TCD2018_2022_quality_tile_matched_source_20261124.py`
5. `scripts/audit_FR381_TCD2018_2022_signed_source_GO_STOP_20261124.py` は4つの独立集約署名artifactを再読込し、source COUNT、年グリッド、質マスク、匿名403、source-only no-EQRことを検証。**最後のSOURCE GO/STOP独立監査が完走していれば別途リンクを追記。**

このSource GO/STOPは生態学的な否定結果でなく、必要な本物のsource raster pixelsを取得できない**認証/可用性STOP**である。非公開のS3 credentialsをGithubにコミットしたり、画像を推測して誤った使い方をしない。元のPR149（2018河畔景観の断面予測は正、EQR2018調整後のEQR2022への追加予測は負/不確実）を変更しない。
