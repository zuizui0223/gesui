# EPA 2013–14同一訪問データ：並列経路と直列経路の統計的比較（2026-11-15）

## 新しい問い

**河川実測TN/TPと底生EPT割合の対応は、同時に測定した生物膜Chla・AFDMで説明されるか。それとも河床・河畔構造の予測情報が支配的か。**

本解析は**既存NRSAの事後探索**であり、因果SEM・下水処理による間接効果・事前登録実験ではない。既に結果を見たPR #79の同一元データを異なる推定量で再構成する。差があれば双方保持し、結論を更新する。係数の符号だけで解析対象を絞らない。

## 実際に同一UIDで観測できるもの

米国EPA NRSA2013–14、元PR #79で厳密化済み**1,094訪問・693 HUC8・47州**。2013/14 source-first semantic gate original Actions #43560769667、既存スクリプト `scripts/gesui_NRSA1314_independent_cycle_numeric.py` をそのまま再実行し、下記の一時CSVを再構成：
`$RUNNER_TEMP/gesui1314_independent_numeric/AUDIT_ONLY_NRSA1314_AREAL_AFDM_SOURCE_PANEL.csv`。
`UID`は本来のEPA同一訪問IDで重複禁止。

| 段階 | 実測の指標 | 本解析の意味 |
|---|---|---|
| 資源環境 | 実測logTN・logTP | **河川水中**、下水由来N/Pとして分離されていない |
| 生物膜 | log10面積当たりChla（µg/cm²）、AFDM（mg/cm²） | 立存量/色素、NPPではない |
| 生息地フィルター | XEMBED、樹冠・草本・河幅等6項目 | 基質・河畔環境、下水以外の攪乱も含む |
| 消費者群集 | 底生EPTの**カウント割合** | 絶対個体数・羽化量・体重・送粉ではない |

## 統計的なpiecewise path decomposition

同じ全1094行・同じカテゴリーと共変量・同じHUC8均等重みで次の線形射影を適合させる。

1. `z(Chla) ~ z(TN) + z(TP) + geography/season/state/eco9 + six PHab`
2. `z(AFDM) ~ z(TN) + z(TP) + identical controls`
3. `EPT_fraction ~ z(TN) + z(TP) + z(Chla) + z(AFDM) + identical controls`
4. `EPT_fraction ~ z(TN) + z(TP) + identical controls`（総条件付き係数）

各栄養塩の総条件付き係数は**同一測定の線形代数上の恒等式**：
`c_total = c_direct + a_Chla*b_Chla + a_AFDM*b_AFDM`。

これは「媒介」を発見したという意味ではない。時間順序・曝露源・共通原因の制御が未証明だから、係数積は**観測された条件付き関連の代数分解**としか呼ばない。

- HUC8均等WLS、HUC8クラスタ頑健95%、500回HUC8単位再標本化の係数積95%区間。
- 物理環境と生物膜が共通する誤差・環境履歴を持ち、栄養塩・生物膜・EPTが同時測定である点は常に明示。
- 両期NRSAを独立サイトとして擬似的にプールしない。別の2018–19期にはAFDMの面積単位が非互換。

## 競合する予測モデル

同じHUC8完全留保5分割、train-only数値標準化・州/生態地域ダミーを使い、アウトカムEPT割合の未知HUC8予測を比較する。

- M0：geography/season/river conductivity + 実測TN/TP。
- M1：M0 + 6物理環境。
- M2：M1 + 実測Chla/AFDM。

HUC8均等MSE、同じ保持予測残差をHUC8で2000再標本化した記述95%改善を報告する。

主張の判断：M1の改善が安定し、M2の追加改善が不確実なら**資源媒介の強い普遍的予測は支持されない**。ただし統計的媒介に非有意でも真の因果媒介が不存在とは言えない。

## 機構と一般化のための比較プラン

1. 米国の同一地点piecewise解析：資源パスの係数積と物理環境予測を分ける。
2. 元PR #83の2008・2013・2018期でEPT付着形質と耐性値：物理環境による組成選別の外挿を別検証する。
3. 欧州PR #96：同一H7流域×期間EQR・S50・蝶HNSの共有環境応答・H6残差共分散（パス矢印の因果推定はしない）。
4. 共通概念は**資源量と選別フィルターが別の生物階層にどう映るか**。別国・別生物層の生係数を1つのSEMに混ぜない。

## 技術と状態

新コード：`analysis/20261115_piecewise_path/observational_paths.py`。新Actionで真正6 EPAテーブルをSHA検証して元の1,094行を再構築し、統計的経路・HUC8外挿を出力する。旧解析結果に見える正負を全て保持。**Actionが成功するまで新しい係数・有意性を主張しない。**
