# フランス376同一河川地点：真正2018→2022年次土地被覆変化は底生動物EQR変化の予測に追加的な情報を与えるか

**2026-11-24｜実際の公開2018/2022年次10m衛星土地被覆クラス画像・同一河川での2回の底生生物状態を使用。SOURCE/地域分割/品質を生物結果の前に固定した探索的・非因果分析。**

## 結論と位置づけ

**結果：2018–2022年の近岸0–250mに対する外側250–500mの樹木分類割合の相対的変化を、同じEEA水生生物の`ΔEQR2022−2018`へ追加しても**、元2018年EQR/地理/河川規模/元2018年のRiparian Zones土地被覆からの追加予測改善は**小さく、不確実**だった。地理的に未使用の5foldで、original376 sites地点均等**+0.474% MSE improvement**、**14地域等重み+0.805%**、元14地域の固定OOF損失からの記述的95%区間**[−1.228,+2.901]%**。同じ固定sourceの外側250–500m樹木クラス変化のみを追加しても地点+0.510%、地域**−0.044%** [−2.401,+2.510]で明確でない。

この結果を「河畔樹木・森林には水生群集を支える機能がない」「河畔木本植生の衰退が下水で起きなかった」「汚水起源の生態効果なし」とは**解釈しない**。本データは年次**衛星AI推定の樹木vs非樹木のclass label**であり、実樹冠密度％・植物の生体量や落葉供給・水→陸の生物量フラックスではない。2018→22のクラス変化とEQR変化は同じ期間であり、時間的先行性・河道の連続的曝露・逆因果を識別しない。

## Genuine sourceと旧TCDアクセスSOURCE STOPは別

1. 前のCopernicus [Tree Cover Density 10m yearly 2018/2022, PR #150](https://github.com/zuizui0223/gesui/pull/150)は両年の元CLMS native GeoTIFFがCDSE S3/OIDC認証を要し、匿名HTTP403。**TCD原画像の数値源アクセスSTOPは取り消さない。** 本研究が使ったのは別製品（植物現存量ではない9分類land cover）。
2. 新たに開いた独立公開ソース：**Impact Observatory / Esri / Microsoft `io-lulc-annual-v02` 10m Annual Land Use Land Cover 9 classes V2**、全球2017–2023、**CC BY 4.0**。[Microsoft Planetary Computer STAC original source](https://planetarycomputer.microsoft.com/api/stac/v1/collections/io-lulc-annual-v02)、[原著作側 Maps for Good](https://docs.impactobservatory.com/lulc-maps/maps-for-good.html)、[AWS open registry](https://registry.opendata.aws/io-lulc/)。本当のnative **2018 `31U-2018`と2022 `31U-2022`** itemにpublic SASを付け、2018/2022原10m COGの真TIFF headerを**HTTP206**取得。[Genuine source-first #45244690575 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45244690575)、original source metadata artifact ID **14853750935**。元class code：0 NoData、1 Water、**2 Trees**、4 Flooded Vegetation、5 Crops、7 Built、8 Bare、9 Snow、10 Clouds、11 Rangeland。独立年画像のclassifier errorやoutcome-dependent data-qualityの選択を避けるため、元生物2018/2022 EQR数値を開く前に全検定指定。
3. **実際の同じ物理河川地点**：原EEA2018底生BQE EQR2,634サイト original SHA256 `2a1a72d43656aca588a832dfae13860054c01a651e76bf25f984c5a374709b02`、original2022 EEA source2,475サイト git blob `325f199e985d9862ffee5b0edb205c61a7be3833`、原2018 Copernicus RiparianZones LCLU exact site/geo eligible978、元2018数値EQRと合致976、元2022再訪418、評価method/HYRIV/元座標一致409、うち**フランス381サイト**。これらは[前段PR #149](https://github.com/zuizui0223/gesui/pull/149)の生物数値を見ずに固定済み、元5 geographical source folds **122/58/48/33/120**と14 original3度ブロックを維持。
4. Original STAC MGRS年次同じUTM源 item検索を実行。[source-only France381 site mapping #45245728321 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45245728321)、artifact **14852638642**。ここで**データソース自身の注意が必要**：国際日付変更線の元MGRS `01U`、`60U`についてSTAC `bbox=[-180,47.96,180,56.04]`という**全世界幅の経度表示**が返るため、WGS bboxだけのマッチではFrance siteをPacific source tileに誤割当した。生物EQRを見る前にsite longitudeから**厳密なUTM Zone EPSG32630/32631/32632**を要求して再計算した。True mainland France381 = native30T40,30U8,31T256,31U13,32T42,32U22、**381/381**二年V2 source supercell対応。これを含む誤tile事例・修正前後のソースQAsは全てActions runとして残す。

## 事前固定した品質条件と本当の生データ画素

本来の直接測定TCDではなく、実年度別のSentinel2-derived**人工知能分類クラス**を真10m UTMで使用した。元TIFFの同じnative UTM image grid 2018 vs2022で**CRS・geotransform・width/height・解像度10m一致**を求め、河川観測点から0–250mと250–500mの非重複円環を元投影のピクセル中心距離で定義。**両年の共通画素**だけを比較し、年を跨いで NoData(0)/Clouds(10)・Water(1)・Snow(9)を除外、真の「樹木(2)」面積割合を計算した。これにより水位で濡れた面積や欠測を「樹木の消失」へ誤変換しない。

事前閾値：元近岸>=200・遠方>=400の真正paired land pixels; 残存cloud/no-data除外<=5%; 最終全地区>=300実地点、元5foldの**それぞれ>=30**、元3deg region>=10。そのsource gateは応答EQRを開く前に凍結。[pre-outcome protocol](https://github.com/zuizui0223/gesui/blob/analysis/fr381-io-lulc-2018-22-original-source-20261124/docs/PROTOCOL_20261124_FR381_IO_LULC_V2_2018_2022_YEAR_CHANGE.md)。

- original fixed SHA site one per foldのsource pixel pilot **5サイト**：raw sourceの両年画像はすべて真正の10m class codeとして読み取れ、同年度 native grid一致。[Pilot #45246664786 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45246664786)、artifact **14854437254**。5中**2サイト**は外側250–500mのNoData/Cloud lost fraction **20.045%・7.854%**のため、元の5%条件FAIL。先にあった失敗でsite/半径/欠測条件を緩めず、source full censusへ進めた。source品質FAILは生態のnullではない。
- [真正全381地点のoriginal native10m class pixel quality #45247194137 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45247194137)、原public tile年次source一対一。固定source gatePASS **376/381 = 98.69%**、5不適格（元cloud/nodata quality fail, not biological effect), original 5fold eligible **121,58,46,33,118**、**14** original 3deg region残存、全source最低閾値合格。原ソースのquality+site category source artifact **14854692199**、site IDはGitHub artifactではSHA256ダイジェストで識別し、元EEA site dataの実座標・元raw raster本体はGitHubにコピーしていない。
- 品質適格376地点の2018→2022 source `Trees`面積率の変化記述：近岸0–250m **−1.517パーセントポイント**、外側250–500m **−1.616ポイント**、近岸−外側の相対変化平均 **+0.099ポイント**。各地点の近−外変化差の10–90%点 **[−4.381,+4.506]ポイント**、範囲[−24.117,+32.160]ポイント。これは**衛星で推定された分類クラス率の差**であり、現地植生減少の原因や地上のバイオマス変化を示すものではない。

## Genuine fixed 5 geographical heldout actual biological test（**探索的・因果ではない**）

真の同一France siteで測った`EQR2022−EQR2018`平均**−0.02117**（この376-site source eligible subsetだけの記述値）。元生物数値を開く前に前段original source GOが成立していることをプログラムで検証、SOURCE failならモデルを一切fitしない。

固定のgeo5fold、original site equal and 14block equal, 2018元地理foldを変更しない、train-only numeric standardization、Ridge `alpha10`。全予測対象は同じ376source site、入力baseline B = **実測EQR2018 + lat/lon + log river Q + log watershed area + 2018 original RZ outer-ring natural-like・urban・cropland class composition**。

| 2018→2022年の`ΔEQR`を予測（実独立geo heldout） | 元baseline MSE | 追加後MSE | source site均等MSE改善 | 14region均等改善 | 14region固定OOF loss bootstrap95（**記述区間**） |
|---|---:|---:|---:|---:|---|
| **Primary：∆Trees近岸−外側**をBに追加 | 0.01790636 | **0.01782146** | **+0.474%** | **+0.805%** | **[−1.228,+2.901]%** |
| Secondary source: **∆Trees外側のみ**をBに追加 | 0.01790636 | 0.01781503 | +0.510% | **−0.044%** | [−2.401,+2.510]% |

14地域のうち予測MSE改善したのはprimary**9/14**、secondary**10/14**。どちらも元地域ブート記述区間は0を含む。これは「増えた樹木がEQRを改善した」「減った樹木がEQRを悪化させた」という**回帰係数の符号や因果**を示すものではなく、**樹木クラス変化という時変説明変数の追加予測力の点推定が小さい**ことのみ示す。新しい`Δtree` sourceはclass-label誤差を持ち、真の植生ストレス応答を測ったわけではない。

[Genuine source-locked EEA2018/2022 delta-EQR original numeric Actions #45247849723 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45247849723)、derived aggregate artifact ID **14853634027**。このフランス376 EQRのbaseline MSE=0.01790636は元PR149フランス381のbaseline 0.01752151とは対象source n/response異なるので直接比較しない。

## Independent numerical audit

別スクリプトで元`sklearn.Pipeline`/`Ridge`を使わずに、`NumPy`の正規方程式 `β=(X'X+10I)^{-1}X'(y-ybar)` + 元fold train-only mean/popSDで再実装した。元source 376サイトと元EQR2018/2022のSHAを再監査し、元14regionのbootstrap2000（固定seed20261124）を元の結果を見ずに再計算する。

**[Independent real-source/5fold 376-site NumPy audit Actions #45248258976 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45248258976)**、signed independent QA artifact **14853409644**：primary +0.47415% /region +0.80548% [−1.22837,+2.90103]、secondary +0.51008% /region−0.04444% [−2.40143,+2.51011]を含め、**全5fold MSE・14ブロック誤差・区間が一致**。ここでの独立性は計算経路が異なることであり、別の時期の独立現地調査ではない。

## 先行研究ギャップへの節度ある解釈

- 元の[EEA RiparianZones same-site 2018/22 PR149](https://github.com/zuizui0223/gesui/pull/149)では、元固定2018年の土地被覆は2018/22同一siteの横断生態状態に地理的予測情報を持つが、2018年状態を知ると**2022年生態状態の追加予測に支持なし**。今回ようやく**2018→22年実画像クラス変化**という時変候補を同じsiteに合わせることができたが、**予測改善は小さく14地域の不確実性が大きい**。
- 相関の機構：Tree class changes may affect shading, litter, erosion or habitat; but no source measures true tree physiology, biomass, leaf N, aquatic subsidy or sewage exposure dose. Other land development, hydrology, climate and EQR method/detection changes can jointly influence satellite land class changes and EQR. Class changes may reflect mapping error as well as real changes. **Temporal synchrony cannot prove causal lead**.
- 元の2018国全体の2018 RZ mapに対する国別ソースreporting selection、2018/22同一河川でのWWTP P報告が両年陽性**14/409**だけであること、国内France381 eligibleだけの解析の限界を維持。**WWTP-specific SOURCE STOP [PR #139](https://github.com/zuizui0223/gesui/pull/139)**を変更しない。
- Original sources Attribution: Impact Observatory/Esri/Microsoft annual LULC V2 (CC BY 4.0), Copernicus RZ2018 56 class map, EEA Waterbase Biology 2018 and2022 official EQR, original HydroRIVERS matched site registry. No raw copyrighted third-party GeoTIFF, EEA individual site coordinates or private credentials committed to GitHub.

### いま採用できる狭い科学的主張
> The same 376 French river-monitoring sites had genuine satellite-estimated annual tree-class changes over 2018–2022, yet incorporating the difference between near-river and outer-riparian tree-class changes produced only a small, geographically uncertain improvement in predicting changes in benthic macroinvertebrate ecological status beyond baseline site condition and landscape context.

「陸→水への生物機能伝播の欠如」や「生態系サービス非存在」ではなく、**静的landscape associationを実時系列の作用機構へ変換するには、真の変動・分類精度・源物質・時系列の実測が別途必要**という反証的な境界。これは単独でNature級/トップ級の一般法則ではない。探索的な同国site cohort上の差分統計と予測に限定。

**すべてoriginal GitHub Draft branch review-only。メインブランチ/WWTP causal STOP/原他結果を変更しない。**
