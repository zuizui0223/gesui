# 2026-11-17 | 原本2018年 河川底生EQR×Copernicus河畔土地被覆：地域内の予測と国際外挿の差

**Status: source-authenticated, exploratory / posthoc observational reanalysis. No causal water→land transfer, actual plant productivity or biological interaction was measured. Earlier Austria Sentinel 2 QA STOP and FR263 HRVPP status unchanged.**

## 元の共通地点で実際に完了したこと
- 2018年オリジナルEEA淡水底生観測地点 **2,639**（元有効数値EQRは2,634）。元source全件を、真2018年Copernicus Riparian Zones土地被覆ポリゴン（CODE_1_18）と同じsite_keyで対応。母集団を有意性で選んでいない。
- 観測地点を中心とする **0–250m内円**と**250–500m外円環**の両方で70%以上が本物の分類済みポリゴンになった **986地点**。非水域土地面積が各帯10%未満だった8地点と原数値EQR対応なし2地点を除き、真EQR+陸域面積分数のモデルは**976地点/13国/49元3度空間ブロック**。元上流P報告正は169地点（source-only 986中は175）。これは正規の下水なし対照ではない。
- 陸域変数：natural=元CODE_1_18の3(森林),4(草地),5(低木),7(湿地)/分類済み非水域面積。nearとfarの差平均 **+0.02970**（モデル976）。2018年の1時点、土地被覆の**面積構成**であって植生生産性・現地バイオマス・栄養吸収量ではない。
- 元物理的河道ポリラインとの側方距離で切った帯ではない。流域面でなく**水質・底生動物観測地点の近傍同心円**である。正確な水際の選別とは異なる。

## 水域生物への地域内予測：元5地理的空間foldを完全保持
地点（site）重みと、元49地域の等重みを区別する。Ridge10、学習foldのみ標準化、国・生物評価法、緯度経度・流量・集水域面積を基準モデルとして再現。
| 比較 | 地点重み MSE改善 | 49ブロック等重み MSE改善 [固定誤差記述区間] |
|---|---:|---:|
| 基準→外円環の自然・都市・農地分数 | **+3.4228%** | **+9.4341% [1.876,17.005]** |
| 土地被覆→+内外自然分数差 | **−0.2485%** | **−0.3492% [−1.257,0.493]** |
| 土地被覆→+上流報告P | **+1.3632%** | **+1.2410% [−0.703,3.158]** |
| 土地被覆→+下流報告P placebo | +0.1143% | +0.1167% [−1.685,1.853] |
| 土地被覆+上流P+自然差→+交互作用 | +0.2504% | +0.0964% [−0.381,0.715] |
| 土地被覆→+上流報告N | +1.4237% | +2.2723% [−0.295,4.486] |

重要：異なる国・異なるRZ地理カバレッジの条件で、距離依存の**追加**情報は原5foldでは検出されなかった。上流報告P自体の地域等重み区間もゼロをまたぐ。「下水で自然被覆が増えた」「水から陸に窒素・リンが移った」「EQRが生物多様性増加した」という結論は不可。

## 強い観測ソース選択バイアス（独立再集計）
| 元2018年2,639地点 | RZマップ両帯合格 | 不合格 | 合格割合 |
|---|---:|---:|---:|
| 上流P源の報告が正 | 175 | 29 | 85.784% |
| それ以外（報告なし≠排水なし） | 811 | 1,624 | 33.306% |

合格率比 **2.576**、単純オッズ比 **12.084**。国別層化した共通オッズ比（9つの計算可能な国）**8.761**。これはRZカバレッジ・報告制度・地域選択の非無作為性であって、排水による植生変化ではない。非水域面積≥10%かつ両帯70%の元source-only n978には、内円の自然被覆比平均0.57152、外円0.54239、内−外+0.02913（49地域等重み +0.02278、地域固定平均再標本化記述95% [+0.00328,+0.04464]）。こちらもsource-onlyの近傍の構成差。

## 国を完全に留保する予測（実行済、厳しい事後ストレス検査）
元同じ976地点、13国について国を丸ごと訓練から外す。小さい国（SE2、SI7、HR8）も事前に含め、結果に合わせて除外しない。国別固定誤差の等重み集計と記述95%区間を示す。
| 予測変数の追加 | 地点重み MSE改善 | 13国等重み MSE改善 [固定誤差95] | 改善した国数 |
|---|---:|---:|---:|
| 基準→遠方土地被覆 | +2.399% | **+1.860% [−3.27,+10.12]** | 6/13 |
| 土地被覆→+上流報告P | +3.239% | +2.123% [−5.98,+7.18] | 8/13 |
| 土地被覆→+内−外自然差 | +0.055% | **−0.169% [−1.44,+1.14]** | 8/13 |
| 土地被覆→+上流報告N | +3.518% | +2.101% [−5.96,+6.77] | 8/13 |

従って**最初の地域内foldの+9.43%を欧州一般的な強い生態系間の結合として示すのは不正確**。国間外挿では符号・予測の利得が揃わず、河畔の内−外自然被覆差は独立したEQR情報として弱い。過去の仏FR73実測TP→底生EQRと、欧州全域の陸域植物→蝶寄主の肯定結果は**異なる観測単位**であり、両者を媒介SEMにつなげない。

## 出典・実行・監査
- Official 2018 RZ full2639 source signed [Actions #44004357868](https://github.com/zuizui0223/gesui/actions/runs/44004357868), artifact14293433193, ZIP SHA256 00b1b72b3812d9c541740d3542da82cf4717265ed014df3a6852508865ce3f6b, JSON SHA256 ff60c807be99468ae95e0c8cd6fabbb60a70a5e55aa3bf51f12dc5f67cdb9f9b.
- Actual EQR × RZ same-station original spatial 5fold [Actions #44005855019](https://github.com/zuizui0223/gesui/actions/runs/44005855019) artifact14293429762, ZIP SHA256 a898195e9136a5d801b2ad1e29aee06c15a613ea8638d15aa9c5ea0d2a12cb32, JSON SHA256 5bf7f983e99ffaf7abdb1fb2eef6058e5560184ed4a5a6b91921429a63030fc0. Site weighted MSE confirmed separately by original 5 fold n=220,178,121,242,215 weighted average.
- Country-purge model [Actions #44006787936](https://github.com/zuizui0223/gesui/actions/runs/44006787936) artifact14295425313, ZIP SHA256 38da0db1d7c170c13a409f5630f84ef0e7d1b2c6b3aef2a6de0885357f18642e, JSON SHA256 00d29f8d86cf0a58cd27d456bee87e8ddb490b60f80371d02c2599f1b9119257.
- Original program: analysis/20261117_RZ2018/eea2018_RZ2018_benthic_EQR_geoblock_model_if_fullsource.py and country QA companion script, run on Draft PR #123. No main merge or confirmed aquatic-to-terrestrial nutrient flux.

## Decision and next non-sign-chasing follow-up
This **confirms an authentic paired cross-domain macro source panel and a within-region landcover–aquatic association** but shows no robust river-localized extra signal and weak out-of-country prediction. Do not switch distance windows or classes seeking significance. The next distinct dataset, not a posthoc redefinition of this one, should be genuine 2018+2022 10m *vegetation productivity* at the same water sites / actual river reaches, with predetermined satellite QA and matching, or measured hydrologically connected emergent aquatic invertebrate biomass; otherwise keep present inference as spatially specific environmental co-variation rather than propagating sewage.
