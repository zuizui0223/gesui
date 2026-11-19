# 2026-11-19 欧州底生動物 EQR→分類群別形質解析：原本ソース検証（Source-only）

**結論：EEAの2018/2022 EQRから分類群の出現・個体数・耐性CWMを逆算することはできない。別の真正原データは取得可能であり、英国EA Biosysは正式な同一site→sample→taxonスキーマの実取得・ZIP検査に成功した。ただし下水特異的な形質選別モデルはこの段階では未実施。**

## 新たに確定した元資料とデータ境界

1. [EEA Waterbase Biology 2026](https://www.eea.europa.eu/en/datahub/datahubitem-view/7bd881f0-59ca-4f0a-a43c-7e0cd0c66d6c) は正規化EQRなどの生物状態指数を提供するが、地点×分類群の出現計数表を提供しない。既存のEEA2018/22上流WWTP報告N/P→底生EQR地域外予測約+1%と、分類群耐性変化は別の推定対象。
2. [EA England Biosys original bulk downloads](https://environment.data.gov.uk/ecology/explorer/downloads/) は実ファイルを検証済み。公式 `INV_OPEN_DATA.zip` **62,373,479 bytes**、`OPEN_DATA_TAXON_INFO.zip` **825,708 bytes**、HTTP200、サーバーのbyte-range 206と完全本体ダウンロード・ZIP整合検証に成功。[公式実ファイル原本schema検査 Actions #44543020599](https://github.com/zuizui0223/gesui/actions/runs/44543020599)。
3. Biosys ZIPの正式CSV3表：**INV_OPEN_DATA_SITE.csv**（SITE_ID、CATCHMENT、WFD_WATERBODY_ID、地形・河川環境）、**INV_OPEN_DATA_METRICS.csv**（SITE_ID、SAMPLE_ID、ANALYSIS_ID、SAMPLE_DATE、SAMPLE_TYPE、BMWP/ASPT/WHPT sample index）、**INV_OPEN_DATA_TAXA.csv**（ANALYSIS_ID、TAXON_LIST_ITEM_KEY、TOTAL_ABUNDANCE、ABUNDANCE_CATEGORY）。結合用独立辞書 **OPEN_DATA_TAXON_INFO.csv**（TAXON_LIST_ITEM_KEY、TAXON_RANK、PARENT_TLIK、PREFERRED_TAXON_NAME）。実ZIPの内部非圧縮量：SITE 約7.74 MB、METRICS 約137 MB、TAXA 約806 MB、辞書 約5.56 MB。行数、join合格数、実個体数の可否はまだ数値監査していない。ファイルがあるだけで生物結果とは言わない。
4. [ChemPop v2 2025 EIDC 原本 DOI](https://doi.org/10.5285/23f194b0-1e21-4f14-be9c-6eabdcc1feb7)：**英国1,519河川サイト、1965–2018年**、taxon abundance、41水質指標、モデル化下水寄与・流域環境を1つの公開データ製品内に収録。公式metadata v2ハッシュ検査成功。公式配布indexから本物の **/data/23f194b0-1e21-4f14-be9c-6eabdcc1feb7.zip** へのリンク確認。初回のread timeout後、bounded downloadで本物のZIP形式は確認されたが、最初の意図的ファイル数上限300に達しschema未検査としてSOURCE STOP、無理に好ましい生物モデルを開いていない。より一般的な全memberとCSV列の限定走査に修正して再実行。旧版 DOI 10.5285/faa526cc-52f5-4468-97bb-295660bcea34 はv2に置換され、旧誤座標版を使わない。
5. [Pit et al. (2025) Pan-European Macroinvertebrate Zenodo](https://zenodo.org/records/15268354) はSampleID・分類群・IndCount・国・座標・方法などの原2010–2015年200万件超のレコードをもつ別データ。Zenodo公式APIのrecord/filesメタデータ取得に成功。原公開CSV約804MBのdownload/サンプル粒度は未開封。**2018/22 EEΑ WWTP源と同じ日時の直接曝露実験ではない**。

## BMWP/WHPTの耐性形質の正式な扱い

BiosysのMETRICSにある **BMWP_ASPT/WHPT_ASPTは試料単位の合成生物品質指数**で、各分類群の耐性スコアと同値ではない。別TAXON_INFOは階級・親名の辞書であって、個別BMWPスコアを含まない。別の公開 [Rob Briers, CRAN biotic](https://github.com/robbriers/biotic) Rパッケージには公開されたBMWP/WHPT family score `indextable` があるため、固定原コミット(08cea3f...)からの source-only辞書取得と公式EA nomenclatureの厳密照合を追加。旧指数定義がEA2018年修正後のすべての指数と一致すると前提にせず、行・階級・スコア欠測/複合科を源で検証する。

## この時点で言えることと禁止すること

- **データ結合の実施可能性（schema）は肯定**：同一の実EA Biosys原ZIPからsite×analysis×taxonと公式分類辞書のキーが存在する。成功は「同一試料の観測による分類群構成と説明変数」を作るための**source availability**に関するもの。
- **生物学的効果量は未計算**：この新規EA/Bachiller-Jareno/Pitソースに対して、下水投入→耐性種、TP→耐性形質、豊富さ、密度や多様性の新しい回帰係数をまだ推定していない。
- **WWTP定量曝露は未識別**：EA Biosys ZIPだけでは年別処理場報告N/Pの上流到達量も実測植物・蝶との因果媒介も作れない。モデル化下水比率がChemPopに含まれるとしても、地点×サンプル×年の元キーを検証するまではsource→biota因果対応と言わない。
- **年代・分類群・独立性**：EA BiosysとChemPopは英国の同じ監視網を部分共有する可能性があり、元SampleIDを照合するまで「独立2研究」と数えない。PanEU欧州データにも同じ国家監視が再収録され得る。解析の地域間メタ研究のnを単純に足さない。

## 2026-11-19 実行完了：ChemPop v2全ZIPと英国の分類群採点辞書

**重要な更新：旧「300ファイル超でSOURCE STOP」は、原ファイルが326件であったための監査側の意図的ファイル数制限です。その後全326ファイルのZip・CSVスキーマ検証に成功して解除されました。** 元データを改変したのではなくファイル一覧の監査上限だけを変更しました。[真正元ZIP検証成功 #44544328596](https://github.com/zuizui0223/gesui/actions/runs/44544328596)。

- **7地域の底生動物分類群原票**：各ファイルは正式 `BIO_SITE_ID,BIO_SAMPLE_ID,BIO_SAMPLE_VERSION,BIO_SAMPLE_DATE,BIO_ANALYSIS_ID,BIO_ANALYSIS_METHOD,BIO_ORD_NAME,BIO_FAM_NAME,BIO_SPGRP_NAME,BIO_ABUNDANCE`。地域ごとの同じ試料の「科」「目」別構成を検証する直接の原入力。
- **7地域のsitevariables**：`BIO_SITE_ID` と `EDF_MEAN/EDF_STD/EDF_Q90/EDF_Q95`、`LC_URBAN_PER`、`LC_ARABLE_PER`、流路幅・底質・地形など。EDFは推定下水混入・希釈の**モデルproxy**、単位や高値方向・モデル期間は元資料で検証するまで未固定。実測N/Pの質量ではない。
- **287種の地域×水質項目集計表**：生物site+sample dateに対応するWQ各項目の最小・最大・中央値・平均・LOQ階級・観測数。生物採集日に取得した元試料と完全同時の単発水質実測ではなく集計窓の適合性が必要。
- **第三者公開された英国BMWP/WHPT形質辞書**：ソース原コミット `robbriers/biotic@08cea3f06eec93703bb10b53d1015d8fe8c626c4` のR `indextable` に得点つきBMWP分類群92件。別のEA Taxon Infoにおける異なる1,793科名のうち88の正式科名が完全一致。**88/1793は元の観測試料に占めるスコア付与率ではない**。試料で何件対応するかはChemPopのBIO_FAM_NAMEに当てて初めて判定する。公開された旧スコア辞書はEAの2026公式最新スコア表の代替とは限らない。
- 新しい原ボリュームZIPの中にあるBIO_SITE_ID・BIO_ANALYSIS_IDと原site variablesの整合、試料の正の`BIO_ABUNDANCE`、BMWPカバー率・EDF分散・時期一致を、同一元ZIP上で**source-only**で検証中。生物の係数は未計算。

## 2026-11-19 英国EA Biosys全原資料の結合監査（完成）

[同一EA公式元ZIPでの全原票結合監査 #44544801464](https://github.com/zuizui0223/gesui/actions/runs/44544801464) 成功。**30,752原河川site、281,255異なるANALYSIS_ID、7,590,163分類群原票行、3,069 TAXON_LIST_ITEM_KEY、614 distinct CATCHMENT**。SITE→METRICS→TAXAの原ANALYSIS_IDの照合率100%、SITE_ID/SAMPLE_ID conflict 0。年代1965–2026年。正の数値`TOTAL_ABUNDANCE`は**6,786,884/7,590,163行=89.4%**、残る803,279は未記載/非正値でありゼロを「個体数なし」と無条件に補完しない。原の試料品質・河川曝露・正式形質辞書との結合モデルは未実施。2018/22 EEA EQRのサンプルと同一の実試料ではない可能性を維持。

## 次に具体的に検証する項目

1. ChemPop v2元ZIP：構成ファイルのCSV列・原source 1,519サイト・元水質とモデル下水のどの列が共存するか、試料単位でのモデル更新有無。SOURCE/TRAIT両方が同じsource単位で測られるかを見てからモデル開始。
2. EA Biosys：正式raw TAXA行で `ANALYSIS_ID`→`METRICS`、`SITE_ID`→`SITE`、`TAXON_LIST_ITEM_KEY`→`TAXON_INFO` の原join重複・カバレッジ・真の出現と計数の割合を先に固定する。検定と追加サンプルの基準は [source-only方法書](PROTOCOL_20261119_EU_BENTHIC_TAXON_SOURCE_ONLY.md) に記載。
3. 形質辞書は公開family-level BMWP/WHPT等の**分類群ごとの**スコア源を独立確認し、原合成BMWP_ASPT値を個々の原分類群スコアとして付与しない。
4. 上記を満たしてはじめて incidence / count CWM、科/目内ラベル帰無、別地域・時期の予測を実行。比較する水・陸の元結果と共通のN/P環境処理効果量のメタ解析にはしない。

**新しいsource QAと報告は専用Draftブランチで保存し、元EQRの正・負両結果および gesui/main・原稿・他PRは変更しない。**
