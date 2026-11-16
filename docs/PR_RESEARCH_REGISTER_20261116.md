# 全PR研究登録簿（2026-11-16 snapshot）

[**全104 PRのCSV登録簿**](PR_RESEARCH_REGISTER_20261116.csv) — PR #2〜#110のうち実在する104件を、2026-11-16のGitHub REST APIから取得して記録したものです。PRのURL、実際のopen/merged状態、作成・更新日時、元ブランチ名、**40桁のHEADコミットSHA**、解析単位、暫定判定、結果文書とそれが今回の統合前のmainに存在したかを1行1PRで保持しています。欠番のPR（issueなど）は勝手に補完しません。

## 分類の読み方

| 判定 | 件数 | 意味 |
|---|---:|---|
| 肯定 | 10 | 検証した特定の限定的な関連・予測増分が正または所定の予測方向に再現 |
| 否定 | 17 | 特定の統計的支持/予測再現が不成立、または期待された傾向が頑健でない |
| STOP | 13 | 原資料、地理・時期・実測曝露・標本独立性の事前ゲート未達で効果モデルを実行しない |
| 混合 | 43 | 肯定の一部と否定/地域・感度差・予測不能が併存 |
| 空欄 | 21 | 現在の登録時点で結論を原本に基づき確定せず、推測しない |

これらは**研究内容の暫定分類**であり、GitHubのopen/closed/merged状態とは別です。最終的な論文用の結論やFDR補正はこの欄だけで確定しません。「source-onlyでデータが見つかった」ことを排水の生態学的正の効果とは分類しません。

**結果文書欄は出所を明示できた39件で記載**し、それ以外の65件は空欄としました。結果文書を「存在しない」と認定した意味ではありません。今後の監査時に元ブランチの原文とSHAを照合したうえで補完します。

重要な点：GitのHEAD **SHAだけでは、ブランチ/PR参照の全削除後にもオブジェクトが永久に残る保証はありません**。本登録簿は追跡の索引であり、長期保存を確実にするにはコミットへの保持参照（タグなど）または証拠資産のバイト一致したmainへのコピーが必要です。そのため、今回は下記の主要な水生側原資産**218ファイル**を実際にmainへ選択収録する計画とし、各オリジナルGit blob SHAを `provenance/GESUI_20261116_AQUATIC_PR_SOURCE_BLOB_MANIFEST.json`で検証します。残るPRのブランチや証拠を消したわけではありません。

## 今回の統合選択：水生側の2系統

| 系統 | 元PR | 内容・保存位置 |
|---|---|---|
| **欧州・実測排出N/Pと底生EQR** | [#14](https://github.com/zuizui0223/gesui/pull/14), [#16](https://github.com/zuizui0223/gesui/pull/16), [#17](https://github.com/zuizui0223/gesui/pull/17) | [2022 EQR/N/P実測報告](../results/eea2022_biology_bqe_20261108/RESULTS.md)、[2022欠測・交絡反証](../results/eea2022_bqe_falsification_20261108/RESULTS.md)、[独立報告年2018再現・負の対照・跨年転送](../results/eea2018_independent_temporal_bqe_20261108/RESULTS.md) |
| 欧州・2020の未成立と2016の異曝露 | [#19](https://github.com/zuizui0223/gesui/pull/19), [#56](https://github.com/zuizui0223/gesui/pull/56), [#57](https://github.com/zuizui0223/gesui/pull/57) | [2020原本取得のSTOP](../results/eea2020_source_only_20261108/RESULTS_AND_PREOUTCOME_STOP.md)、[2016 JRCモデル由来point N/P（実測WWTP年間負荷ではない）](RESULTS_EEA2016_REAL_BENTHIC_EQR_JRC_GREEN_NUTS2_MACRO_20261111.md)、[地域・時間反証](RESULTS_BENTHIC_REGIONAL_TEMPORAL_POSTHOC_20261111.md) |
| 欧州・河川**実測濃度TP**と底生EQR | [#70](https://github.com/zuizui0223/gesui/pull/70), [#71](https://github.com/zuizui0223/gesui/pull/71), [#72](https://github.com/zuizui0223/gesui/pull/72) | [フランス273観測局の空間関連](RESULTS_20261113_FRANCE_I2M2_MEASURED_2022_RIVER_TP_EQR_273SITES_74H7_THREE_LAYERS.md)、[同一263地点の空間対時間の非同期](RESULTS_20261113_FR_REAL_2018_2022_TP_I2M2_263_SAME_SITES_SPATIAL_TEMPORAL_DECOUPLING_3BIOTA.md)、[同月比較の反証](RESULTS_20261113_FR263_SAME_SITE_TP_2018_2022_COMMON_MONTHS_I2M2_THREE_LAYERS.md)。**河川TP濃度は下水処理場由来の質量ではない** |
| **米国EPA河川・3調査期の基礎** | [#73–#83](https://github.com/zuizui0223/gesui/pull/83) | [真正NRSA3期のHUC8非重複地域外予測](RESULTS_20261114_NRSA_THREECYCLE_EPT_AMBIENTNP_PHAB_HUC8_EXTERNAL_TRANSPORT.md)、[3期のEPT形質同定・厳格支持](RESULTS_20261114_NRSA_3ERA_RARE100_AND_EPT_FUNCTIONAL_TRAITS.md)、同一訪問Chla/AFDM、河床・光環境、欠測・群集等重み/個体重みの検証。一部EPA POTW報告ゲート失敗や生物膜媒介nullを含む |
| 米国EPTの経路null・形質の上積み | [#97](https://github.com/zuizui0223/gesui/pull/97), [#99](https://github.com/zuizui0223/gesui/pull/99), [#100](https://github.com/zuizui0223/gesui/pull/100) | [同一訪問N/P→Chla/AFDM→EPTの統計的パスnull](RESULTS_20261115_NRSA1314_PIECEWISE_PATH_VS_PHAB.md)、出現分類群の形質×環境、州/水文地域外の感度と結果。元[#98のFourth-corner EPT](RESULTS_20261115_EPT_FOURTHCORNER_TRAIT_INCIDENCE.md)はPR #106ですでにmainへ収録済み |

## 方法論・結論の厳格な境界

2018・2022年のEEA実測**報告された処理場の年間N/P排出質量**が底生EQRの未知地域予測へ小さな追加情報を持つことは、「同じ年に河川の底生生物へ下水のN/Pが届き、そのため増えた/減った」ことではありません。異なるEEA年の生物調査局が重複し、国別生物EQRの評価法が異なり、排水以外の農業・都市化・河川形態の交絡が残ります。

米国NRSAの**河川内の実測TN/TP**とPHab河床物理環境を使った群集形質モデルは、欧州のWWTP出水施設から報告されたN/P質量とは**曝露と地理単位が違う**ため、ひとつの「下水効果量」として統合メタ解析しません。

このページとCSVは、[統合前に保存された水陸系の肯定・否定台帳](RESULTS_20261116_VERIFIED_MAIN_INTEGRATION.md)、[従来の証拠台帳](EVIDENCE_LEDGER.md)、原PRブランチのSOURCE_RECEIPTに並ぶ索引です。文書だけで独立数値再計算が成功したと主張せず、元の原本取得・数値監査の成功Actionsを各RESULTS文書に残します。
