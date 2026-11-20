# 2026-11-19 欧州底生動物の分類群原票への移行：source-onlyと反証可能な解析規約

**現状：原本取得可能性の段階。形質・水質モデルを新たに推定したと主張しない。** 生物群集の分類群等重み（incidence-CWM）をEEAの総合EQRから逆算してはいけない。

## なぜ既存のWWTP→EQRから直接進めないか

EEA 2026 Waterbase Biology の河川底生無脊椎動物アウトカムは国・手法で較正したEQR値で、原分類群の出現・個体数ではない。既存の2018年2,634、2022年2,475河川地点、報告された上流N/P追加によるEQR予測改善1%前後は、分類群・属内・感受性の変化の証明ではない。

**確認した別データ源（独立の標本であり、元EEA EQR地点を自動で補完しない）：**

| ソース | サンプル×分類群の機会 | 曝露/特徴 | 識別すべきこと |
|---|---|---|---|
| Pit et al. (2025), Pan-European Freshwater Macroinvertebrate Dataset, Zenodo [10.5281/zenodo.15268354](https://zenodo.org/records/15268354) | 2010–2015年、200万件超、Taxon/TaxonID/IndCount/SampleID、方法・日時・座標 | 多国・広域・分類群を統合できる | 2018/22 WWTP源と時間一致しない；直接N/P源を付与できるとは限らない |
| Bachiller-Jareno et al. (2025), England ChemPop v2 [10.5285/23f194b0-1e21-4f14-be9c-6eabdcc1feb7](https://doi.org/10.5285/23f194b0-1e21-4f14-be9c-6eabdcc1feb7) | 1965–2018年、1,519サイト、分類群数/原計数・41水質項目、河川環境 | 同じデータ製品にモデル化 sewage effluent exposure、土地利用、地形・流量あり | site×sample×yearの元連結、静的/時変曝露、測定TP/TNの完全性、モデル下水地理との分離 |
| EA England Biosys official [downloads](https://environment.data.gov.uk/ecology/explorer/downloads/) | Freshwater river macroinvertebrate surveys + official Biosys Taxon Info | 別公式原票の同定・独立調査系の候補 | 生物分類群とBMWP/LIFE等の正式スコア、ChemPop原試料との重複と独立性 |

**特にChemPop v2**：旧版 DOI 10.5285/faa526cc-52f5-4468-97bb-295660bcea34 は撤回／v2が置換。v2の座標誤りの修正を保持し、旧版や他試料との重複を独立研究として数えない。

## Phase 0: Source-only（実装済み）

- EIDC v2 metadata の公開UUID/版・公開日・ライセンスを検証。ChemPopデータ本体とサポートファイルはHEAD + 高々4KiBのRangeのみで問い合わせ、巨大生物原票を無断でraw Git履歴へ配布しない。
- Zenodo公式Record id・files大小とchecksumを確認。大容量CSVをこのsource-probeでは取得しない。
- EA Biosys 公式2ZIPの提供/アクセス状態、content-length、応答タイプを監査し、アクセス不可は**SOURCE STOP**とする。
- 曖昧な EEA-BQE 指標値から分類群の存在/不在/耐性を推測しない。
- 全リモート応答のヘッダと有限バイトのSHAをsource receipt JSONへ保存する。バイオアウトカムや結果符号は開かない。

## Phase 1: ChemPop v2 実データ取得可能なら最初に行うゲート

対象の生態単位は **原調査の site×date×sample×taxon**。少なくとも以下のすべてを、耐性・環境のモデルを回す前にsource-onlyで数える。

1. 原データの正式ファイル一覧、ライセンス、バージョン、download SHA、site_id/sample_id/date/taxon/rank/abundance units と曝露の対応元情報。元列名を実際のschemaで検証し、推測列名・日付の補完は禁止。
2. 1サイトに複数の調査日/試料があるとき、複数計数行の重複、subsample、同じ family/種の混在同定、季節、計数精度を調べる。raw abundanceと ordinal abundanceを合算しない。
3. **同一試料・同一調査期**の水質（実測TP/TN、環境）とモデル化された effluent exposure を分離し、真の下水由来N/P到達量と呼ばない。サイト固定のモデル下水比率が多年度同じ場合、年次の dose-response としない。
4. 目標：少なくとも200品質合格のsite×periodと80独立実地点、4以上の環境/地理ブロック、少なくとも2段階の元曝露値の実変動、元EQR/処理場報告と別の原出現記録があること。原site/sample/habitat/source-score での欠測が曝露の高低で偏る場合は先に表にする。基準は**後付け源選択防止のため、アウトカムを見る前に固定する探索的設定**。
5. 正式な汚濁耐性形質はoriginal BMWP/WHPT等の taxon/family scoring lookup を**公開原表の同定/版検証後**に付与。既存EQRや観測地点のBMWP総合値を分類群名札とみなさない。科と種が混在したsampleは重複集計を禁止し、科以上の同定行は別感度へ。
6. Score coverage：主に出現分類群の得点可の割合が70%以上、少なくとも3 scoring taxa per sample。原始の報告された個体数が正・実計数である試料のみ count-CWM を計算。countがordinalなら等重みincidenceだけ、count-CWMはSTOP。

## Phase 2: 検定（Phase 1合格時のみ）

- ① incidence-CWM score（各現地出現分類群等重み）、count-CWM（元正式個体数がある場合）の各sample環境・処理場proxy関連を比較。差の係数も95%区間付きで報告する。標本数の見かけ上の急増を防ぎ、siteとcatchmentの両方でcluster。
- ② 形質名札の帰無：公式 taxon/family scoring のみを同一科/目／系統に近い階層内で入替。原出現と環境・missingness固定。正式系統樹ではないことを明記。
- ③ 別catchment/別regional authority holdout、さらに別年の源連続性が確保できれば別年holdout。trainのみ標準化、真のsite同一性がtrain/testを跨がないよう厳密除外。
- ④ 年比較は断面の異なる観測サイト vs 真の同地点の曝露変化・分類群変化を厳格に分ける。
- 水質の実測TP/TNとモデル下水比率が共存する場合は、**背景調整＋実測TP/TN → ＋下水proxy**、逆順の順序も固定して差分の外挿予測を出す。実際の処理場放流量を参照できない場合、下水特異的な選別・栄養を媒介とする因果効果はSTOP。
- 因果を主張する際には農業負荷、汚濁物質、流量、地理、土地利用、河床特性・季節の交絡を構造化し、source-onlyで真の上流/下流を構築できない負の対照は「存在する」と仮定しない。

## 論文との境界

Pan-EU→England ChemPop→EEA EQRは標本・年代・曝露が異なるため、サイト数/年をそのまま独立研究数としてメタ解析しない。**共通の見かけの耐性方向**よりも、独立データの地域外予測・ラベル帰無・曝露の確かさを比較する。データが不足すれば具体的な SOURCE STOPを結果として報告し、次に外部サンプルを探す。

eBird、TinyFishは使用しない。元原稿と gesui/main は変更しない。
