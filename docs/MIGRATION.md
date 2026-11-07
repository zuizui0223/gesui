# 2026-11-07 引継ぎ記録

**元：** [zuizui0223/chocho PR #122](https://github.com/zuizui0223/chocho/pull/122) のコミット f5ba867cce723767b31fa536312fa5da6333a661、ブランチ analysis/nutrient-port-trait-pilot-20261102。

**先：** [zuizui0223/gesui](https://github.com/zuizui0223/gesui) の archive/chocho_pr122/ 以下、元の相対パスのまま保存。**80/80 Git blob SHA一致**（[全件マニフェスト](../provenance/CHOCHO_PR122_FILE_MANIFEST.json)）。元ファイルを削除、上書き、再解析、論文本文へ反映していない。

- **chochoの責任範囲：** 全球のチョウ寄主資源、外来寄主、島・港、風、AIS航路、海上Circuitscapeと有効抵抗、チョウ類の分布拡大。既存GEB成果は保存。
- **gesuiの責任範囲：** 処理場・排水量、放流先と河川接続、栄養塩、人為負荷と植物相の配置変化、淡水・河岸・河口・沿岸の多分類群解析。島や港を必須としない。

旧workflowを新しい .github/workflows として有効化しない（過去のGBIF・ラスター外部APIを不要に繰り返さない）。旧コードは前提データをchochoに依存するので、**保存＝完全に再現実行可能**ではない。今後の水文的な主解析は独立して実装・検証する。
