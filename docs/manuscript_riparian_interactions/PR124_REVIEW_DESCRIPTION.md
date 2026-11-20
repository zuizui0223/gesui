## 問いと現在の結論

下水関連地理と蝶の寄主資源特性の対応を、人口集積・景観条件から区別できるか。水域と陸域の共変動と、水域内の化学・処理場指標の予測情報を、標本と推定対象を分けて統合する。

最後の固定人口感度解析では、同じ437流域×期間（225小流域・88親流域）のW係数が人口調整前+0.08154［0.04337,0.11970］から調整後+0.03099［−0.01576,0.07775］へ弱まった。人口調整後のW追加予測改善は0.98%［−3.56,7.12］で、固定した2条件を満たさなかった。人口集積から独立した下水作用は未確定であり、下水の無影響や人口の因果作用を示したとはしない。

元451行の人口未調整の肯定結果（β+0.07750、予測改善7.98%）も保持する。437行では人口調整前から予測区間がゼロを含むため、標本差と人口調整の影響を分ける。旧モデルにも耕地・都市・面積・流量は入っていた。

## レビュー用成果物

- [要旨と最終3図・日本語キャプション・PNG/PDF](https://github.com/zuizui0223/gesui/blob/docs/wastewater-riparian-interactions-outline/docs/manuscript_riparian_interactions/REVIEW_PACKAGE_ABSTRACT_THREE_FIGURES_20261120.md)
- [最後の人口感度：仕様・取得元・品質判定・全数値](https://github.com/zuizui0223/gesui/blob/docs/wastewater-riparian-interactions-outline/docs/manuscript_riparian_interactions/POPULATION_FINAL_SENSITIVITY_20261120.md)
- [更新本文](https://github.com/zuizui0223/gesui/blob/docs/wastewater-riparian-interactions-outline/docs/manuscript_riparian_interactions/MANUSCRIPT_FULL_20261119.md)
- [各経路のスケールと支持状態](https://github.com/zuizui0223/gesui/blob/docs/wastewater-riparian-interactions-outline/docs/manuscript_riparian_interactions/PATH_EVIDENCE_SCALES_20261120.md)

図1は蝶の人口感度、図2はフランスの水陸の空間対応と時間差、図3は英国・オーストリアの水域内の参照証拠。水域結果を蝶の曝露を保証する陽性対照とはせず、河川経路の不支持を「川を通らない」と解釈しない。従来のパス図と広域資料監査は補足として保持する。

## 検証と限界

- 仕様をコミット223da1aに固定後、GHS-POP R2023A・2020年・1 km格子を上流512構成流域へ境界画素の面積按分で集約。元232 H7のうち7 H7を固定した上流面積整合条件で除外。欠測人口をゼロにしない。元格子はGitへ再配布せず、URL・SHA256・派生表を保存。
- 別人口製品だが、HydroWASTEの補完と共通の国勢調査情報を持つ可能性がある。静的な上流人口で、局所活動や年間変化を完全に調整したものではない。
- 固定応答・共変量・親流域留保法を維持。仕様と入力のハッシュ、ネットワーク・面積、画素按分・欠測の合成検証、係数の別計算、MSE算術、前後の同一W尺度を確認。結果後の応答・科・モデル探索をしない。
- 水陸共変動の150相関を別実装で照合済み。既存139行の13経路再現、源監査と反証も保持。全129 PRを独立研究として扱わず、新しい全鎖メタ解析や因果媒介推定は行わない。
- 3図の出所と出力ハッシュを保存し、表示確認済み。オーストリア図は引用集約に区間がないため点推定のみ。異なる標本のMSE改善を同じ効果量として統合しない。
- Virginia汚泥資料はこの論文の生態解析から除外し、別の野外研究候補として保管。成虫栄養利用・汚泥は未測定の考察仮説。

生態解析の追加は終了し、指導教員確認用の要旨・3図を用意した。投稿規格への最終調整前。指導教員への送信はしていない。Draftを維持し、mainにはマージしない。
