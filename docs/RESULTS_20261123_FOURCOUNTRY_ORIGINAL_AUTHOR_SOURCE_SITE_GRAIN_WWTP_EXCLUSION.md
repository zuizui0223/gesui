# 4か国外部細粒土砂–EPT解析の原資料による調査設計監査：独立地点性と下水曝露のSOURCE STOP

**日付 2026-11-23 JST｜科学的判断：独立地点単位の推論STOP／下水特異的な一般化STOP。** 本資料は、既存[PR #146](https://github.com/zuizui0223/gesui/pull/146)および[PR #147](https://github.com/zuizui0223/gesui/pull/147)の実数値を再推定して肯定を探すものではなく、元の出版資料のsampling designとsource inclusion/exclusionを真の原本で検証した追加監査である。

## 原資料の真正性

- 元研究：McKenzie et al., *Global Change Biology* 30:e17084 (2024), [doi:10.1111/gcb.17084](https://doi.org/10.1111/gcb.17084)。この研究は既に4か国で細粒土砂に対する群集応答やEPT指標の地域差を解析した。
- 元データ：Loughborough repository [doi:10.17028/rd.lboro.24615570.v1](https://doi.org/10.17028/rd.lboro.24615570.v1)、**CC BY-NC 4.0**。
- データ原本 `MathersEtAl-2023-PairedVisualFineSediment-Data.xlsx`，**4,640,542 bytes**，MD5 `1e4d6e37b73c663406e5567dd79d9452`，SHA256 `a02573bef0e11ea19c9465a78ebecc130fb9c39d76d35a79863c2efd5debcd48`。
- 今回新たに取得した**出版者の調査方法原本** `MathersEtAl-2023-PairedVisualFineSediment-DataSourceSummaries.docx`，Figshare original file id **43299384**，**279,236 bytes**，publisher MD5 **`08d81e4d97af34053b453a1c1e4f8c7d`**，computed SHA256 **`20dd6fc28888e6e28efeb5a9ce65e8e32fc01a3cf55d9fb6251868c556fff1bf`**。
- [原DOCXの構造・方法記述・付属表のSOURCE ONLY監査 Actions #45084627476 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45084627476)，signed aggregate artifact **14777734573**。[原DOCX×6,491行のXLSXの識別子構造を照合 #45084789187 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45084789187)，artifact **14778488184**。[原文の標本除外条件・反復季節・各NZ元研究の標本数を機械的に再照合 #45084925062 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45084925062)，artifact **14779131206**。
- 原著作物DOCX/XLSXそのものはGitHubにコミットしていない。再取得元URLとSHA、解析コード、限定的な構造情報、検証判定だけを残す。

## 実際の原調査設計：国ごとに別の母集団である

| 国 | 真正Excelのsource行 | 同じ元出版資料に明記された源データの性格 | **下水・標本独立性**に関する境界 |
|---|---:|---|---|
| オーストラリア | 3,099 | NSW州の**複数の長期河川モニタリング事業**（Monitoring River Health 1994–2000; Sustainable Rivers Audit 2002–2012; Monitoring Evaluation and Reporting 2006–継続）を統合。春・秋の採集を選択、AusRivASのreference条件判定を満たす地点のみ採用 | 原元sourceは**reference条件へのフィルタ済み**。国全体の全河川・全ストレスレベルのランダム標本ではない。元XLSXの`ANALYSIS_ID`は数値連番で、真のphysical site id不明 |
| ブラジル | 232 | Cerrado地域の複数研究を統合。各sample reachで11 transectに分けた採集は**1つのcomposite sample**にまとめて評価。相対的に汚染されていない場所を選択した旨が原文にある | 11 transectを11独立地点として数えていないのは妥当。しかしreaches across different studies / source periodsの独立性は原XLSXだけでは確認不可 |
| ニュージーランド | 342 | Hunter/Death, Lange, Magbanua, Matthaei, Wagenhoffなど**5つの異なる研究**を統合。元論文の採集法・subsampling個体数が異なる。Magbanuaのbefore/after研究では介入**前のみを採用** | 元出版者の方法表で5系列のn合計は**352**、配布Excelは342。10件の差は原方法表のMagbanua **36** vs ExcelのMagbanua **26**に一致。ただし10件の個別理由は不明 |
| 英国 | 2,818 | Environment Agencyの行政モニタリングと別の農村部学術研究を統合。春・秋に採集。学術sourceの説明では**同じ調査地点を春と秋に採集**すると明記 | **行政データ選別では、下水処理場放流口の下流にある地点を除外**。*英国全標本が下水を受けないとまでは主張しない*が、少なくともこの明示的source部分を下水曝露解析の標本として扱えない |

いずれも**単独の同質な観測設計ではない**。原研究はデータ取得の段階で他の共存するストレス因子の影響を小さくするために地点を選別している。source covariate rangeと生物群集の条件付き母集団は、独立な「全球河川の代表標本」と同一ではない。元のデータ6,491件は**paired biological sample `ANALYSIS_ID`**であり、独立流域やsite countsではない。

## ニュージーランドのANALYSIS_IDのsource-code構造（推測でsite数へ変換しない）

真正原XLSXのNZ 342個のIDは「先頭の研究名＋区分＋…」形式。以下は**単に先頭2区分が同じか**で機械的に分類したコードパターンであり、完全なphysical site IDではない。

| ID先頭の原sourceラベル | 真正source観測行 | 異なる先頭2区分コード | 2回以上現れるコードに属する行 | 最大同コード反復 |
|---|---:|---:|---:|---:|
| `Death` | 206 | 80 | 201 | 3 |
| `Lange` | 43 | 29 | 23 | 5 |
| `Magbanua` | 26 | 14 | 24 | 2 |
| `Matthaei` | 24 | 12 | 24 | 2 |
| `Wagenhoff2008` | 43 | 43 | 0 | 1 |
| **合計** | **342** | **178** | **272** | — |

**先に作った初版の集計「136接頭グループ／315行」は正しくない。** 2分割だけの`Wagenhoff2008.X`をすべて`Wagenhoff2008`一群へ落とした誤処理が原因。正しく原表の先頭2区分を使うと**178コード、272行**。ただし、`Lange.Becks.Creek…`等で最初の2区分が同じでも**物理地点が同じではない**可能性がある。逆に`Death.GRES1.42417`, `Death.GRES1.42759`, `Death.GRES1.43131`のように別時点を示唆する構造もあるが、publishersはこの数値suffixの意味を文書で明示していない。**178は観測地点数の推定値ではない**。

NZを基点に元source期の構成から独立な空間holdoutを作ることは、同じ物理siteを確実にまとめる実鍵がないため**GOできない**。研究名/地点名と日時の厳密なcrosswalkを元データ著者から確認したときだけ再検証可。

## 下水研究に含めてよいことと、含められないこと

1. **不可**：この4か国比較を「下水処理場由来の富栄養化がEPTを減らす」または「下水で水域→陸域の相互作用が変わる」証拠にしない。元データには原本のPOTW放流量、N/P曝露量、河岸植物・蝶や羽化昆虫のflowが**存在しない**。英国行政sampleはむしろWWTP下流を除いている。元の下水曝露直接推定**SOURCE STOP** [PR #139](https://github.com/zuizui0223/gesui/pull/139)を堅持。
2. **条件付き可能**：PR #146/147で取得した公開データ真値・EPT正しい分類照合・4か国外予測のpoint-estimate・正則化感度・within-/between-country数学分解は**監査済み記述統計**として引用可能。ただし本件は元論文が既に答えた細粒土砂-EPT群集関係の再解析であり、完全な新しい生態学的法則ではない。
3. **推論STOP**：サイト単位bootstrap・独立site分割・4か国全標本のp値・6,491物理河川への外挿。国単位のholdoutは同じ国のtrain/testを避ける設計なので同じ国のphysical siteのtrain/test漏出は起きないが、**各国のheldout MSEとその不確実性が独立地点数を表すわけではない**。
4. **標本選択を因果機構にしない**：元データの地域・source差（UK行政のWWTP除外、AU reference-site filtering, NZ before-only）を適切な層別条件として表示する。ただし元データだけから除外された地点の反実仮想は推定できない。環境ストレスによる群集差とsource-selection processの区別が必要。

## 独立データでEPAの一般化限界を確認できたか

**NO GO**。EPA NRSAの`XEMBED`は河床粒径の組成によって定義上100%となる特殊な埋没度、4か国の元曝露はchannel bed <2mmの目視細粒被覆率で、測定量そのものが違う。またEPAの「2008/13/18適格source標本にHUC8が含まれていたか」という履歴は4か国にはなく、検証対象の交互作用を別sourceで定義できない。したがって元4か国の総MSE改善だけでPR #145（EPA slope × historical source group）を独立再現したとは書かない。

## 技術QAと進行判断

- 原DOCXのsource-only監査 [#45084627476](https://github.com/zuizui0223/gesui/actions/runs/45084627476) → 原DOCX×XLSXのsource-code照合 [#45084789187](https://github.com/zuizui0223/gesui/actions/runs/45084789187) → 原文のUK sewage排除・季節反復・NZ研究別採集の**原本機械検証 [#45084925062](https://github.com/zuizui0223/gesui/actions/runs/45084925062)**、すべてSUCCESS。
- 区分とsource-code数に関する初版の誤集計は今回の**系列別厳密パーサー監査**で訂正。独立コード・原MD5、数値source行数、4か国のdataset metadataを再現。
- **優先順位：4か国の同一データへの更なる有意差探索はSTOP。** 今後必要なのは、(a)原著者側のphysical site・sampling occasion crosswalk付き真正メタデータ、あるいは (b)本来の下水研究に戻り本当のPOTW effluent/nutrient loadと水域・河畔・陸生群集を同じ流域・期間で結合する原データの取得。未取得の環境計測を仮定で代替しない。

Source reports are limited to metadata, numerical summaries and few very short source contexts; full original proprietary DOCX/XLSX are not committed.
