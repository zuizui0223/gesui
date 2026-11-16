# 全PRの検証結果登録簿（2026-11-16）

本登録簿はGithubの全PRの元head SHAと実際の変更ファイル一覧を読み、単位・判定・元結果文書・mainに同一Git blobが存在するかを分けて記録したもの。

判定が空欄なのは未精査であり、肯定・否定・停止のいずれとも推測していない。main READMEに言及があっても結果ファイルの統合とは数えない。

- PR全数: **105**。
- 結論精査済み: **23**、未判定: **82**。
- PR結果文書の異なるパス: **102**、mainに元blobが一致しない結果文書パス: **88**。

全PRの行別テーブル: [原本PR登録簿CSV](../provenance/PR_ALL_RESEARCH_REGISTRY_20261116.csv)。詳細な元PR変更文書とGit blob SHA: [doc manifest JSON](../provenance/PR_SOURCE_FILE_DOC_MANIFEST_20261116.json)。

## 主な精査済みPR

| PR | 単位 | 判定 | head SHA（短縮） | 元結果文書のmain同一数 |
|---|---|---|---|---|
| [#14](https://github.com/zuizui0223/gesui/pull/14) | EEA生物観測局×2022 | 肯定＋制限 | 625af95fc12e | 0/1 |
| [#16](https://github.com/zuizui0223/gesui/pull/16) | EEA生物観測局×2022 | 肯定・否定混合／欠測 | 651281a175d9 | 0/1 |
| [#17](https://github.com/zuizui0223/gesui/pull/17) | EEA生物観測局×2018/2022 | 肯定＋地点依存 | 69e140643e8b | 0/1 |
| [#81](https://github.com/zuizui0223/gesui/pull/81) | US NRSA UID×調査期/HUC8 | 肯定＋制限 | 3e206a79b03d | 0/1 |
| [#82](https://github.com/zuizui0223/gesui/pull/82) | US NRSA UID×調査期/HUC8 | 肯定＋感度 | 9452a92128f1 | 0/1 |
| [#83](https://github.com/zuizui0223/gesui/pull/83) | US NRSA UID×調査期/HUC8 | 肯定＋感度 | 5c85b83f0526 | 0/1 |
| [#90](https://github.com/zuizui0223/gesui/pull/90) | EU H7×期間 | 追加予測不確実 | b502cb7f6af9 | 1/1 |
| [#92](https://github.com/zuizui0223/gesui/pull/92) | EU H7×期間 | 方向対照で不確実 | e3ee4f124f0c | 0/1 |
| [#93](https://github.com/zuizui0223/gesui/pull/93) | EU H7×期間 | 方向対照で不確実 | b557e8da768c | 0/1 |
| [#95](https://github.com/zuizui0223/gesui/pull/95) | EU H7×期間 | 肯定＋得点感度 | 853ae25e4e0a | 1/6 |
| [#96](https://github.com/zuizui0223/gesui/pull/96) | EU H7×期間 | 層間残差は不確実 | d21d4599ef74 | 0/1 |
| [#97](https://github.com/zuizui0223/gesui/pull/97) | US NRSA同一訪問/HUC8 | 間接経路不確実 | f226dd5cd643 | 0/1 |
| [#98](https://github.com/zuizui0223/gesui/pull/98) | US NRSA UID×3調査期/HUC8 | 肯定＋感度 | bd51607d23e7 | 1/1 |
| [#101](https://github.com/zuizui0223/gesui/pull/101) | EU H7×期間 | 肯定＋科別注意 | f4aa1fc94764 | 1/1 |
| [#102](https://github.com/zuizui0223/gesui/pull/102) | EU H7×期間 | 肯定＋時間null | f6c3c6b09d9a | 1/1 |
| [#103](https://github.com/zuizui0223/gesui/pull/103) | EU H7×期間 | 肯定＋共通6形質STOP | a90dfaab1ecf | 2/2 |
| [#104](https://github.com/zuizui0223/gesui/pull/104) | EU H7内異なるH7/真H6 | 肯定＋水質未確認 | ae740666246a | 0/0 |
| [#105](https://github.com/zuizui0223/gesui/pull/105) | EU H7内異なるH7/真H6 | 属間正・属内不確実 | 25b46200013d | 1/1 |
| [#106](https://github.com/zuizui0223/gesui/pull/106) | 研究原本単位別 | main統合済み | 4ec510cf087e | 14/14 |
| [#107](https://github.com/zuizui0223/gesui/pull/107) | EEA生物観測局×年/H7 | 水質媒介source STOP | 22497c297f4b | 0/4 |
| [#108](https://github.com/zuizui0223/gesui/pull/108) | EU H7×期間 | 方向特異性は不確実 | 64c1f730e3b0 | 0/1 |
| [#109](https://github.com/zuizui0223/gesui/pull/109) | EU H7×期間 | 方向維持＋外挿感度 | 2a360515298e | 0/1 |
| [#110](https://github.com/zuizui0223/gesui/pull/110) | EU H7×期間 | 空間正＋時間外挿不確実 | a4015e5cffbb | 0/2 |

## 同じ平面で比較する際の単位

EUは真正HydroBASINS H7×元の2期間の454行（有効植物NUTRI451行）。水生EQR・植物NUTRI・蝶HNS・WWTP報告N/Pを同じ流域と期間へ揃える。EEA実測負荷と底生EQRは元観測局からH7に再集計するが独立観測数は増えない。
米国NRSAの3期はUID×survey-cycle、地域外検証はHUC8単位。欧州のモデルWWTP地理Wと米国の河床埋没度XEMBEDは別の曝露であり同一の生態学的効果量としてプールしない。

GitHubにPRのhead SHAを記録しても、元ブランチ削除後の恒久的blob保存までは保証しない。重要な結果ファイル・コードは元Git blobを選択的に統合する。第三者のEEA/GBIF/EPA原アーカイブは本リポに再配布しない。
