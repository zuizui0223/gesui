# France381 河畔の2018→2022年樹冠変化：公式GeoServer公開WMS/WCSと認証済S3の数値アクセス検証

**2026-11-24 JST | 実験種別：SOURCE ONLY | EEA底生動物数値・元Copernicus TCD画素・2018–2022樹冠変化は未取得／未解析。** [親PR #150](https://github.com/zuizui0223/gesui/pull/150)で見つかった匿名CDSE元TIFFの2018/2022両年HTTP403に対して、公開OGCおよび公式S3認証の正規アクセス手順だけを追加監査した。

## 研究の最初の固定状態は不変

- 実際のEEA 2018→2022同じ物理河川モニタリング地点、BQE method/HydroRIVERS reach一致・元座標一致のフランス **381 site**、14の元2018 3度地域ブロック、元5地理foldを事前固定。元Riparian Zones LCLU2018の0–250m/250–500m環境; 2018/22の真正EQRを用いた以前の時系列テストでは、既知の2018 EQRを入れて2022 EQRを予測するとき河畔土地被覆の追加改善は **−1.540% site、−2.640% region [−6.440,+0.345]**。土地被覆2018は土地被覆変化・実樹冠の2022年値ではない。
- 正規のCopernicus HRL VLCC 10m年次Tree Cover Density TCD 2018と2022の本当の元公的商品識別子がある。source DOI `10.2909/e677441e-fb94-431c-b4f9-304f10e4dfd8`、native `clms_vlcc_tree-cover-density_europe_10m_yearly_v1`。元の公開COGカタログには2018/2022年の各**902タイル**、EPSG3035 native grid**897 matching tile IDs**。公式Confidence layerも両年sourceあり。フランスの同じ381 site周囲500mには、両年の正規TCDとTCDCL商品のcatalogue footprintが完全に存在：**381/381、必要native 47 tiles**。**これはラスタの画素値が読めたことではない。**

## 新たな別経路A：公開GeoVille GeoServerは両年レイヤーを掲載、数値APIは今回の実行環境で全てタイムアウト

公開[source layer preview](https://geoserver.vlcc.geoville.com/geoserver/web/wicket/bookmarkable/org.geoserver.web.demo.MapPreviewPage?0=&filter=false)に、実際に`HRL_TCF:TCD_S2018`と`HRL_TCF:TCD_S2022`という二つのレイヤーが登録されている。ただしこのWMS製品名`CLMS_HRL_TreeCoverDensity_S2018_R10m`と、現在の公式S3/STACネイティブCOG`CLMS_HRLVLCC_TCD_S2018_R10m_..._V01_R00`が同じsource processing versionかは確認できていない。

[真正GeoVille公開OGC源監査 GitHub Actions #45239758429 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45239758429)、監査結果 artifact **14849943610**。

| 正規公開サービスの要請 | 原APIの結果 |
|---|---|
| GeoVille `WMS 1.3.0 GetCapabilities` | **NETWORK_ERROR / timed out** |
| GeoVille `WCS 2.0.1 GetCapabilities` | **NETWORK_ERROR / timed out** |
| GeoVille `WCS 1.0.0 GetCapabilities` | **NETWORK_ERROR / timed out** |
| `WMS GetFeatureInfo` 2018の元レイヤー | **NETWORK_ERROR / timed out** |
| `WMS GetFeatureInfo` 2022の元レイヤー | **NETWORK_ERROR / timed out** |

そのため **WCSから本当のTCD GeoTIFF画素値が取れることを証明できず**、GetMap PNGの配色を数値％に変換する処理は全く行っていない。WMS名の登録があることと、その元数値がこのruntimeから読めることは異なる。これは「すべての利用者からGeoServerが永続的にダウンしている」という主張ではない。

## 新たな別経路B：公式Copernicus Data Space S3の認証・安全な試験入口

[Copernicus公式S3 docs](https://documentation.dataspace.copernicus.eu/APIs/S3.html)は、EOData元ソースへのアクセスを署名済S3 URLと利用者発行のaccess/secret keyによって認可する。元の2018/2022 TIFFへの匿名HTTP `Range`で**両年403**が返ることはPR150の真正原本で既に実測済みである。これは回避するものではなく、正規の資格情報を使用する。

今回**保護されたGitHub Actions Secretsを使える場合だけ実行するsource-only読取スクリプト**をリポへ追加：
`scripts/gesui_FR381_CDSE_authorized_S3_TCD_TCDCL_2018_2022_preflight_20261124.py`
および
`.github/workflows/gesui-FR381-TCD-CDSE-S3-auth-source-preflight-20261124.yml`。

特徴は以下の通り：
- GitHub Actions`CDSE_S3_ACCESS_KEY`、`CDSE_S3_SECRET_KEY`の二つの**repository Secrets**を利用。ログやsource JSONに実資格情報を表示しない。
- 元の2公的カタログCSVのSHA256を再検証して、**実在する年・元物理グリッドE36N28と製品名から**CDSE S3内の`eodata`真正ファイルを特定。2018/2022 TCD 10mおよび同年`TCDCL`Confidence layerの**4つのTIFFについてhead 32 bytesだけのS3 `Range`でnative TIFFヘッダー**を検証。サイトの樹冠画素データや生物EQRは読まない。
- 認証情報が無い場合は安全に`STOP_NO_CDSE_S3_SECRET_CONFIGURED`で停止し、欠損を0や「植生変化なし」に置き換えない。
- User account credentials are not passed via chat or committed; Github secret store only.
- 本workflowは**まだ実ラスタ成功と誤認してはならない**。

[実際のauth-optional source proof Actions #45240408046 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45240408046)、原本source aggregate artifact **14850882118**：
**`source_status=STOP_NO_CDSE_S3_SECRET_CONFIGURED`**。これは「Workflowの失敗」ではなく、「sourceの鍵が存在しないため実数値の測定を行わず、正規に停止した」という安全な実行結果である。

## Source claim gates

| 進捗・表明できる内容 | 判定 |
|---|---|
| 2018/2022同じフランス381水生調査地点と2018公的RZに対応 | **SOURCE_SITE_GO** |
| 2年分のCLMS TCDとTCDCL元10m同一投影タイル・品質商品の存在 | **SOURCE_META_GO 381/381** |
| 元TIFFの匿名CDSE APIで直接取得する | **SOURCE_STOP_HTTP403** |
| 別GeoVille公開WMS/WCSから実raw numerical band値取得を確認 | **SOURCE_STOP_NETWORK_TIMEOUT** |
| 正規CDSE S3でsource keyを使った実TIFF 4ファイルの先頭を取得 | **SOURCE_STOP_MISSING_REPO_SECRETS** |
| 真の2018/22品質画素判定、近岸と外側の樹冠被覆変化 | **NOT MEASURED / no ecological regression** |
| 本当の河畔植生変化が底生無脊椎動物EQRの2018→2022変化に影響 | **NOT IDENTIFIED** |
| WWTP由来P/Nの河畔植生吸収や水域–陸域相互作用 | **NOT MEASURED** |

## 再開条件を限定

CDSEアカウントを利用可能な者が[公式S3 access credential発行手順](https://documentation.dataspace.copernicus.eu/APIs/S3.html)に従いaccess key/secretを取得し、GitHub repository `zuizui0223/gesui` の **Settings → Secrets and variables → Actions → Repository secrets** に**それぞれの名前で**保護して登録。本人がworkflowをRunするだけで、本当に同一EPSG3035 source tileの2018/2022 TCDとTCDCLが各32bytes取得できるかを安全に確認できる（no authentication bypass）。その次のSource stageでは初めて両年の本当の10m樹冠割合/品質数値を381siteの元近岸・遠方円環に測定する。そこがGOになるまで新たな生態学的`ΔTCD→ΔEQR`係数やp値は**表示しない**。

## 原結果に対する意味

本件は残存のsourceアクセスボトルネックを正確に特定し、認証を安全に整備するための研究準備が進んだ一方で、**前段のPR149の空間的肯定結果/2018EQR調整後の2022年追加予測欠如の科学的証拠を変えない**。TCD source metadataと元公的RZのclass-codeは植物 biomassではない。WWTP負荷が元のFrance381では両年report positive14/409のみだったことも変更できない。source year mismatch or lack of causal delivery path remains a limit.

No numerical TCD raster or third-party original protected products were copied into Git; this branch is Draft based on PR150, with `main` and WWTP source STOP unchanged.
