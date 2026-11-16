# 凍結仕様 v1｜図4「河川からの距離に伴う水域―陸域の形質対応の変化」

**作成：2026-11-16。元451 H7×期間における正・負・STOPを見た後に立案した事後反証計画であり、独立した事前登録ではない。** 結果を見て距離帯や閾値を変更しない。既存 main PR #106 + #112 の結論は変更しない。元PR #110/111のnullと源取得制約も保存する。

## 0. 単位と表現を取り違えない

- 同一欧州451 H7×period（232 native H7、90 native H6）は**生物条件EQR、植物EIV NUTRI、蝶の潜在寄主HNS、モデル排水関連地理W**をまとめた同一流域×期間の分析母集団。
- EQRは**現場で測定された底生動物の生態学的品質比率**であり、河川の実測TN/TP「水質濃度」ではない。EIV NUTRIは植物の生育地栄養嗜好、蝶HNSは潜在寄主の栄養嗜好であって、植物葉N、幼虫の実摂食、実際の排水由来Nではない。
- 欧州2018/2022**処理場による年間報告排出質量**→底生EQRの小さな地域外予測増分は、元451 H7×periodの完全同一標本ではなく、別のEEA Biology監視局コホート（2018 n2634、2022 n2475、重複サイト1167）である。**同じ451流域で4層の放流水化学まで揃っているとは書かない**。
- 既知の同一451空間リンク：PlantN→Butterfly incidence HNS β+0.1956 (95% [.1446,.2465])・H6留保予測+36.56% ([15.34,49.51])；W→PlantN β+.2152 ([.1171,.3134])・+18.51%；Aquatic EQR←PlantN β−.03446 ([−.05969,−.00923])、H6留保MSE改善+2.838% ([−4.698,+9.804])、**地域外追加予測不確実**；Butterfly HNS←Aquatic EQR after PlantN β+.000597 ([−.0257,+.0269])、予測−1.57%。ΔPlantNとΔButterfly HNSは同一H7時間同期null、同method継続のΔAquatic EQR←ΔPlantNもβ+.0121 ([−.00354,+.02772])でnull。
- **β・記録重み・種均等・MSE改善の数値を直接縦軸に並べて「減衰係数」とはしない。** 本図の直接比較には**全く同じH7×期間コホート、同じ真H6、同じ共変量、同じ対照モデルからの偏相関rとH6留保MSE改善**を再計算して初めて採用する。

## 1. 研究質問・2つの競合仮説

**問い：** EQRと植物NUTRIの流域内条件付き関連は、河川近傍で強く、遠隔で弱いか？（動機＝水陸境界で観測指標の空間粒度が食い違う可能性）

- **H`river`（河川近傍に局在）：** 同じH7×periodの植物NUTRIについて、0–250 m（主近傍）でEQRとの負の対応が1–5 km（主遠方）より**明瞭に負**、同じH6留保でも近傍Nは追加予測に寄与。これは河畔環境選別と整合するが、水質からの栄養移動・下水起源を特定しない。
- **H`landscape`（流域全体に共有する環境・観測努力）：** EQR–植物Nの偏相関や予測は距離帯でほとんど変わらず、または調査努力・道路近接・国×method・土地利用で消える。既存の同一流域空間的な植物–蝶の関係は残っても、水際からの移送・距離減衰は未支持。
- **第三の結果（距離帯の観測が不足）**：距離と効果は識別不能、図4に距離別傾向を描かない。肯定・否定のどちらとも扱わない。

## 2. 真正原資料・追加取得・位置の検証

- 元のH7×period 454/232/90パネルSHA256 `1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca`。元植物数値EIV accepted vascular speciesKey 4,011、元NUTRI-scorable H7×period 451、元EIV SHA256 `ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba`。元GBIF 2,808 facet JSON SHA256 `251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4` は**H7×period×PlantaeのspeciesKeyごとの件数集計**で、**各観測の緯度・経度・GBIF occurrenceKey/coordinateUncertaintyを保存していない**。したがって**元facetの集計値だけで河川距離を復元してはならない**。
- 距離モデルを動かす前に、同じGBIF occurrence source universeを原本source時と一致する**original polygons WKT・2017–19/2020–22・taxonKey Plantae6・hasCoordinate=true・hasGeospatialIssue=false・occurrenceStatus=PRESENT**で再取得する。現行API結果の観測時点の違いがあるため、既存H7×period countとGBIF source endpoint `count` のずれを記録・SHA監査し、同じ科学的標本だと仮定しない。**個別occurrenceKey・speciesKey・座標・座標不確実性・記録日・datasetKey・basisOfRecord・observerが必要**。選択は科学的アウトカムを開く前のGBIF source-only source gateで固定。
- GPS点を持たない原本標本・GBIF不確実度欠測・基準を超える座標精度の観測は**NA**（河川から遠い0ではない）。主：coordinateUncertaintyInMetersが明示され≤100m、感度≤250m；原本GBIF地理問題なし。国・年・種を含むidentical eventKeyや完全重複位置/日付/種/observerを元レコードにより重複検査。座標がよくない元sourceを安易に増補しない。
- 河道：公式 **HydroRIVERS Europe** の真の線データ・バージョン・ZIP SHAをまず記録し、HydroBASINS v1c L7の真の232ポリゴン（原ZIP SHA `b9b810f4db271f0fa6b39fecf01784041b948944b3c3f296a393766388e86b44`）に属する河道のみを用いる。**川からの距離**は、元流域内で最も近いHydroRIVERSの河道線へのWGS84座標からの測地線距離（メートル単位）を、適切な局所等距離投影または幾何学的に精度検証した方法で計算。河道に入らない小支流はHydroRIVERSの縮尺・最小河道の閾値で欠けうる。道路・農地・都市縁の距離との相関を記録。
- **厳格な橋渡し副検定**：同じH7内に複数の川があるため、EQR局が所在する**同じHydroRIVERS河段／流下ネットワーク**まで固定できたH7のみで検定。最短「いずれかの河道」距離は**河岸環境proxy**であり、EQR局の河道から栄養や物質が移動した証拠ではない。
- 後年の追加GBIF記録数、無座標件数、国別・H6別採取努力、GBIF dataset/basisOfRecordの偏り、距離帯間のハビタット面積・道路の密度を**モデル前に監査**し、都合のよい地点だけ採らない。

## 3. 固定距離帯・測定アウトカム

- 主3帯（互いに非重複）：**近傍：0–250 m、移行：250–1000 m、遠方：1000–5000 m**（メートル、0は河道線上の幾何距離）。5kmより遠い記録は本比較に使わないが**各H7×periodの除外割合を別表に記録**。buffer面積とGBIF人間記録密度は異なるので、面積が同じと扱わない。
- 各帯・H7×periodの観測結果は、同じ元GBIF厳密受理植物種の**EIV NUTRI種均等 community-weighted mean（CWM）**。補助：観測GBIF record数重みのCWM。植物個体数・被度ではない。種が見つからない帯をNUTRI=0にしない。
- **主の固定努力対比**：各帯・H7×periodの数値EIV採点植物観測から**50 record**を重複なしに抽出、seed 20261116、事前に指定した100 Monte Carlo反復。それぞれの希薄化のscoreを計算し、GBIF記録数を揃える（修正後の追加試験）。近距離と遠距離は**同じH7×periodのpaired unit**。MC分位は**現地独立研究のCIではない**。
- 植物数値スコア付き種と総植物種両方のcoverageに注意：全数site×bandで元GBIF個別speciesKeyのユニーク数に対する有効NUTRI accepted numeric-scored speciesKey率≥**40%**、得点種≥**10種**、得点付きGBIF記録≥**50件**。中間帯も同じ判定だが、近傍と遠方が主要 paired unit；移行帯欠測で遠近比較を不必要に落とさない。
- EEA生物EQRの**現地site→original H7×period mean**は元パネルで固定し、近傍植物点からWQ/EQR現地siteまでの実距離も別に監査。EQR自体は2026 WQ測定TN/TPではない。

## 4. 原本source-only標本支持ゲート（植物NUTRI/EQR効果値を開く前）

**主の近傍 vs 遠方比較が成立する最低条件：**
- 同一H7×periodで近傍・遠方それぞれ得点付き50 records・得点種10種・採点種coverage40%以上を満たす、**異なるH7×period ≥100**；
- **異なる真正H6 ≥40、元country_holdout区分 ≥5、独立H7 ≥80**。元期間は2017–2019/2020–2022を混同しない。両期間があるH7は最大2H7×periodで、予測分割は同じH6ごとに行う；
- 完全原本451 H7×periodからなぜ除かれたか（座標不明、50記録未満、スコア欠測、川データなし、2帯片方欠損など）を年/国/H6別に**生物値を開く前に**公開；
- 1件でも原source SHAやWKT/H7 IDの整合性が壊れ、座標重複や元個体観測に異常がある場合はSTOP。元451全件を合格したとみなさない。

**未達時：** `STOP_RIVER_DISTANCE_PAIR_SUPPORT`。図4は地理的なサンプル制限図のみ（効果を描かない）。生態学的境界減衰に関する主張は**既存の同一H7との空間対応差**にとどめ、距離減衰を「実証」と書かない。

## 5. モデルと予測は同一コホート／同じ推定量で比較

主の水陸指標＝`EQR_field-site_mean ← PlantN_bin`。実証前に全モデル共通のソースは国×assessment method・期間、モデル静的W、crop、urban、流量、地形/面積、元EEA EQR source effort、植物GBIF(dataset数/記録努力)、各距離帯のNUTRI score coverage、coordinate uncertainty、土地利用/道路への距離と元source固定13 controlsの識別可能な非共線項目を指定。

- 同じpaired合格H7×period・同じH6、同じ調整Xと重みで、N_near, N_middle, N_farの**各1SDに対するH6クラスタWLS β・記述95%**、環境/努力残差に対する**加重偏相関r**を全帯報告。
- **主要の直接比較統計量：`D = r_far − r_near`**（EQRとNUTRIの負の関連が河川近傍で強いならD>0）。各距離帯の相関の差を**同じpaired H7×period/true H6についてH6クラスタ2,000回再標本化して共分散も保持**。両側95%区間がゼロを除外しD>0、かつr_near<0（同一評価法）を、**距離による空間的結びつきの減衰支持**の必要条件とする。近傍だけp<.05で遠方は非有意、という判定は禁止。
- 主の予測診断：各bandを共通ベースモデルに1つ追加したときの**元H6 disjoint 5fold Ridge alpha=10、H6等重みheldout MSE改善**と**fixed-error H6再標本化記述95%**を一枚に載せる。対比の予測増分差を同一foldの各H6損失の**paired差**としてbootstrapする（元の異なるH7系の数字をそのまま比較しない）。
- 100回のGBIF記録50件希薄化すべてで、beta、D、追加予測・そのH6区間を評価。主判定の頑健性は`D>0かつ95%下限>0`となる反復の**割合**と、MC分位を報告する。2,000回のH6 bootstrapは固定re-fitまたは固定予測損失のどちらか明記し、MCのみで独立確認pとしない。
- **河川特異性の偽対照**：同じ閾値・同じGBIFデータで河道に近いことと**道路に近いこと**の予測の違いを事後に選ばず比較し、県/国・urban/crop・アクセス・記録source effortに条件付ける。道路、川、農地の共配置と「水陸生態学的移送」を分離できなければ**条件付き河畔組成の対応**と表現する。
- データを開いてから効果方向と上流下流の定義を変更しない。下水由来WWTP N/Pの質量は局所EQRとGBIF観測点へ実際到達したか未確認なので**因果媒介SEMを実施しない**。

## 6. 蝶について

蝶の個別GBIF座標は元facet集計には保存されていない。全5科・両期間で同じsource query/GBIF speciesKey個別記録から真正位置が取れれば、同じ距離帯の蝶種均等HNSの**補助的解析**を行う。主とは別の事前source gate：各band ≥20個別score付き蝶レコード、≥5受理HNS種、同じpair H7×period ≥100/true H6≥40/5country。未達ならSTOP。蝶の成虫を河岸で見たことは幼虫の現地寄主利用・下水からの吸塩・腐肉利用を意味しない。

## 7. 図4と論文の言い分け

- **図4が支持**：タイトル例「Riparian proximity predicts stronger spatial correspondence between aquatic condition and plant habitat-affinity composition」。これは**近い植物観測点の環境適性と水生EQR**が対応したことであり、処理水Nの移送速度・距離減衰率を測ったわけではない。
- **図4が不支持**：タイトル例「Watershed-scale plant–butterfly trait covariation persists without a demonstrable riparian gradient」。水陸の対応が薄いことを「生態系隔離の実証」とは書かず、**空間分解能と観測過程に依存する条件付き関連**として議論。
- **データ不足**：「HydroRIVERS distance test not identifiable under coordinate quality and spatial source support」。残りの陸域・水生結果は論文でそれぞれの推定母集団をラベル付き比較、直接の境界減衰の断言を避ける。

**実務：** 元PR #112が2026-11-16にmainへマージ済み。新規解析はmain起点のDraftブランチのみで、main・chocho・既存PR #110,#111を変更しない。最初に真正source-only feasibility（GBIF個別座標再取得・HydroRIVERS本物の線・H7地理/努力）を実施。元GBIF facet aggregateだけでモデルを走らせない。停止を不利な結果として隠さない。


## 第1波のsource-only実測（固定仕様後、結果を見る前の確認）

2026-11-16付け **Actions #43825788157 SUCCESS**：original source facet `251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4`から真H6・元country_holdoutが各8群異なる8 H7×両期間（計16 source query）を、GBIF生物アウトカムを見ずにSHA選定。現行GBIFへ原ソースのWKT/年/taxon6/座標あり/geospatial問題なし/在を再発行し、**16/16 GET成功**。各query先頭最大300件、実際に計**4,654個別植物GBIF records**の元sourceフィールドschemaを確認、WGS84緯度経度4,654、coordinateUncertainty≤100m**2,991（64.27%）**。**先頭ページの質監査**であり、元sourceの全植物speciesKey facet、11,511種、4百万全記録における構成割合ではない。元索引作成時点と今のGBIF生息地record indexが完全に同じと断言できない。EIV NUTRI・蝶HNS・EQR・Wの数値応答は開いていない。

2026-11-16付け **Actions #43826262099 SUCCESS**：HydroSHEDS公式 [HydroRIVERS v1 Europe / Middle East Shapefile](https://data.hydrosheds.org/file/HydroRIVERS/HydroRIVERS_v10_eu_shp.zip) のソースZIP **67,648,957B**, SHA256 **`500da7d36ceee0aa4c82dd625f5c53d80ba2f7ac14ad022cd3a6601bf55068e5`**。線データ**938,544 reach**、元 native `HYRIV_ID`, `NEXT_DOWN`, `HYBAS_L12`, `ORD_STRA`, `DIS_AV_CMS` が含まれる。WGS84 line geometryで河川距離と流下方向のsource実装が可能となった。原本zipも全ライン座標も本リポへ転載せず、短い署名付きsource QAだけを保存。**まだGBIF個体pointを線へ投影した距離帯やEIV平均、EQRとの偏相関を計算したわけではない。**

**HydroRIVERS v1の観測対象制限**：公式仕様上、河道は集水面積概ね10km²以上、または推定自然流量0.1m³/s以上のネットワークが基本。小さな支流・排水溝・都市暗渠まで列挙した流路ではないため、実際の水際距離より**遠く推定する方向の誤差**が起こりうる。主解析とともに水系密度／小支流欠損の地域差を記述しなければならない。

**次のsource-onlyゲート**は、元8 H7×2期で「個別GBIF source point + 真のHydroRIVERS reach」から距離帯の記録数・受理EIV得点種数**だけ**を算出し、近傍/遠方の≥50得点記録・≥10得点種・≥40%対応が各地域で可能かを確認。その後、元全451母集団のsource completenessを調べる。現在は`RIVER_DISTANCE_MODEL_GATE=NOT_FITTED`。
