# 下水マクロ生態学：全球の正のGBIF記録相関と、標準化された植物群集の実測栄養塩勾配を分けて検証する

**2026-11-10 / 証拠統合版 / POSTHOC; actual sewage→plants ecological causation NOT established.**

## 研究全体の問いと今回の検証

「下水発生量が多い地域では生物が多いのか？」を測定尺度で分ける。世界の港とモデル下水発生量の正の相関は実在するが、GBIFの観測記録や未希薄化種数が高いことは植生区画内の真の植物種多様性増加を意味しない。今回の目標は **既存の記録相関と独立した標準化植生調査から同じ正の傾向が見られるか** を検証することであり、個別のシジミチョウ・放流口の局所解析は主役にしない。

## 1. 既知の肯定的マクロ相関：モデル下水発生量と生物記録

World Port Index 1,055港・42か国、PANGAEA Jones et al.モデルWWp（2015年の**下水発生量**、港25km内）、HydroWASTE処理場数、2021年ESA WorldCover市街地/農地、国固定効果および地理統計を再現。港5km内の2011–24年 GBIF全Plantaeデータを独立に取得。以下は全 **1,055港・42か国** の未希薄化GBIF記録に対する回帰係数（国クラスタ95%CI）。

| 全植物のGBIF応答 | WWp回帰係数β | 国クラスタ95%CI |
|---|---:|---:|
| **植物の記録件数** | **+0.31556** | **[+0.18161,+0.44951]** |
| **植物の記録種ID数** | **+0.22884** | **[+0.10938,+0.34831]** |

2指標の正の係数は42か国を1国ずつ除外しても維持。植物の種IDと十分な記録を持つ**752港・38か国**に限定したモデルでも未希薄化記録種数β=+0.17567 [0.06685,0.28449]、全植物記録数β=+0.29411 [0.12871,0.45952]。ただし**50 GBIF記録にそろえた希薄化種数**は β=+0.01256 [−0.00088,+0.02599] と小さく、不確実性区間が0をまたぐ。

同じ752港で**Animalia（動物界）**記録件数β+0.37115 [0.06924,0.67306]、**Fungi（菌類）**β+0.38012 [0.13648,0.62377]も正。非植物の記録件数を共変量に入れると植物記録種数の係数は+0.17567→+0.04117 [−0.0140,+0.0963]へ縮小する。**下水発生量に対応して複数分類群のGBIF記録が集中する**ため、記録努力・地域観測集積を強い代替説明として扱う。動物・菌類は『下水の影響を絶対受けない陰性対照』ではなく、共変量調整だけで因果を除去できない。

**重要な定義：** WWpは2015年のモデル下水**生成量**であり河川/土壌へ実際に出たN・P・Na・有機物の量ではない。GBIF記録は種数や植物量の完全な標準化現地調査ではない。2015年のWWpと2011–24年の記録には年代のずれと都市・港バイアスがある。

## 2. 全球95,104植生区画 sPlotOpen — 原本実測地理の適格性を監査し、STOP

実際の植生区画の候補として、sPlotOpen （95,104植生調査区、42,677維管束植物分類群、1888–2015年）を調査。公開GitHubの2024年ミラー版**24,975,505 byte**ヘッダーのみをGit LFSから取得し、原本sha256=`b20e21c70bf1c7fc2d79db403797e018bacbc38f74bb851dcf28825bc131977a`を検証。公式2021年原本とのバイト完全一致は未検証。対象の湿地・人為植生は元のsPlotOpen選抜設計から除外されるため、そもそも下水影響の強い水辺/都市植生を代表しにくい。

測定値を見る前に固定した、2000–2015年、全維管束植物を記録した区画、1–1000m²調査面積、元の位置精度≤1000mという条件で、**8,311植生区画・14か国**のみ適格。うち5,257区画はオーストラリア、1,615区画は米国。元の全球1,055 WPI港から25km内には**12調査区**、50km内24区画、100km内150区画しかない。固定した≥15か国という最低品質条件を満たさず、**Jones2015のWWpを区画地点に割当てて植物種数を検証するステージには進まない**。このゲート失敗は生態学的下水効果の否定ではなく、現地調査区の地理・生息場所が下水発生量比較に不向きというデータ設計の限界。

- [sPlotOpen original source-only plot QA JSON](https://github.com/zuizui0223/gesui/blob/analysis/global-splotopen-2015-WWp-plots-feasibility-20261110/results/global_splotopen_WWp2015_20261110/ORIGINAL_SPLOTOPEN95K_2000_15_PLOT_SOURCE_ONLY_GEOGRAPHY_QC.json) — [GitHub Actions original LFS header SHA validation SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42832317681).

## 3. 独立した標準化された植物群集：米国 EPA 2021年NWCA湿地

U.S. EPA National Wetland Condition Assessment 2021 の公開公式4 CSVから、**確率標本設計の湿地初回訪問**、調査重み有効、緯度経度有効などを条件に920地点を抽出。植生指標`H_ALL`（植物全体のShannon指数）、`H_NAT`, `H_ALIEN`, `FQAI_ALL`等と実測水質の全窒素`NTL_RESULT`（mg N/L）、全リン`PTL_RESULT`（µg P/L、解析時mg P/Lに換算）、湿地タイプ、9生態地域を**元の`SITE_ID + VISIT_NO`で568地点結合**。このSOURCE ONLY段階ではHやTN/TPの数値を開かず地点数/州/生態地域を確認。

固定した品質フラグ・正の実測TN/TP・H_ALL非欠損条件で**499地点・48州・9生態地域**が実際に解析可能。湿地の内訳は内陸369、潮汐130。元の2021測定TN・TPは**下水由来と同定されていない現地水質**である。

### 3a. 実測栄養塩と植物群集（すべて独立した現地測定）

| 指標：元の湿地全植物Shannon多様性 H_ALL | 結果 | 95%CI（州クラスタ） |
|---|---:|---:|
| 単純TN vs H_ALL順位相関 | ρ **−0.30062** | — |
| 単純TP vs H_ALL順位相関 | ρ **−0.29651** | — |
| **TNのみ、湿地クラス・生態地域・緯度を調整**（TN10倍当たり） | **β=−0.20251** | **[−0.33955,−0.06546]** |
| TN・TP同時調整の**TN**係数 | −0.19051 | [−0.39357,+0.01254] |
| TN・TP同時調整の**TP**係数 | −0.01330 | [−0.18341,+0.15681] |
| 元のNWCA調査設計重みを使ったTN・TP同時調整のTN | −0.25644 | [−0.63334,+0.12045] |

主解析のTN・TP同時投入と面積重み調整では信頼区間が0をまたぐ。したがって米国湿地全体で栄養塩が実際に植物多様性を減らすという確証的・因果的結論ではない。だが**標準化された植物群集からWWpによる一律な正の生態学的応答は確認されていない**。

### 3b. 後から確認した湿地水文タイプの違い（POSTHOC）

| 2021年の固定済み湿地区分 | 元の実測TN・TP・H_ALL品質適格N | TNとTPを同時調整したTN10倍当たりH_ALL係数 | 州クラスタ95%CI |
|---|---:|---:|---|
| **内陸湿地** | **369** | **−0.26923** | **[−0.48386,−0.05460]** |
| **潮汐湿地** | **130** | +0.02055 | [−0.27054,+0.31164] |

内陸湿地ではTNの負方向が同時調整後も残るが、潮汐湿地では信頼区間が広い。**この水文タイプ別比較は499地点の全体結果を見た後の事後的感度分析**であり、内陸と潮汐間で係数が有意に異なる正式な交互作用検定は実行していない。別途全499地点の電気伝導度`COND_RESULT`も調整したところ、TN係数−0.15634 [−0.35596,+0.04327]に弱まる。イオン濃度・塩分環境の共変動は代替説明として重要であるが、電気伝導度は下水固有の識別マーカーではない。

### 3c. 原本と独立再計算

- [EPA2021 original vegetation/water/hydrology SITE_ID source-only gate](https://github.com/zuizui0223/gesui/blob/analysis/us-nwca2021-wetland-plant-chemistry-source-20261110/results/us_NWCA2021_plant_TNTP_20261110/EPA_NWCA2021_ORIGINAL_VEG_WATERCHEM_SITEONLY_GATE.json) — [original source original bytes SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42833429529).
- [499 actual wetland plots measured TN/TP × plant H_ALL/FQAI state-cluster models](https://github.com/zuizui0223/gesui/blob/analysis/us-nwca2021-wetland-plant-chemistry-source-20261110/results/us_NWCA2021_plant_TNTP_20261110/EPA_NWCA2021_REAL_PLANT_SHANNON_TNTP_WETLANDS_STATE_ECOREGION_QA.json) — [official data model SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42834025457).
- [Independent NumPy estimator original 568 matched source sites, 499 quality visits and 48-state exact covariances](https://github.com/zuizui0223/gesui/blob/analysis/us-nwca2021-wetland-plant-chemistry-source-20261110/results/us_NWCA2021_plant_TNTP_20261110/EPA_NWCA2021_TNTP_PLANT_SHANNON_499SITE_INDEPENDENT_SOURCE_COV_QA.json) — [independent source arithmetic SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42834918331).
- [Inland vs tidal and original measured specific-conductivity posthoc sensitivity](https://github.com/zuizui0223/gesui/blob/analysis/us-nwca2021-wetland-plant-chemistry-source-20261110/results/us_NWCA2021_plant_TNTP_20261110/EPA_NWCA2021_INLAND_TIDAL_CONDUCTIVITY_TNTP_SHANNON_POSTHOC_QC.json) — [original source sensitivity SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42835309877).

## 4. 一本の科学的結論：正の記録相関と栄養塩の生態学的効果を混同しない

**確かに肯定的結果がある：** 42か国のモデル下水発生量と全Plantae/Animalia/Fungiの記録件数が正に相関し、全1,055港の全Plantae未希薄化記録種数とも正に相関する。一方、観測記録件数を50にそろえた希薄化植物種数はWWpと強い正ではなく、米国の実測湿地群集ではTN・TPとShannon多様性の単純順位相関はむしろ負であった。これらは**別の場所・別の曝露変数・別の生物応答であり、1つの横断的な『下水が生物多様性を増減させる』係数に統合できない**。

**現在優先すべきマクロ生態学研究は『下水発生量に伴う観測記録集中と真の群集応答がどこで分岐するか』** である。そのためには、同じ標準化植物調査に**実際のWWTP排水受水域への水文接続、源別の実測N/P濃度、流量、調査努力、植生区画**を結び付ける必要がある。今回のNWCA水質N/Pは下水源を特定せず、sPlotOpenは圧倒的な地理・生息場所不一致がある。

2018/2022年のEQR × 公式実測WWTP N/P報告質量が与える予測改善（2018 P+0.982%/N+1.190%、2022 P+1.103%/N+1.397%）は**別の水生無脊椎動物の予測情報に関する肯定結果**として保持し、EQR正の生物学的効果とは称さない。既存のシジミチョウローカルPRは主線に戻さない。**全研究の主張は観測と生態学的測定尺度を区別して大きくまとめる。**

**原本PDF・CSVs・植生区画個別ファイルは公開コミットしない。下流の下水因果と現地の水質は未接続。Draft PRのまま、`gesui/main`や`chocho`原稿にはマージしない。**