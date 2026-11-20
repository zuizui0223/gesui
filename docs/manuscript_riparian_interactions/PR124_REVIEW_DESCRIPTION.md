## 問いと現在の結論

下水関連圧力に沿って、水域と陸域の生物指標は対応するか。その対応は共通の景観条件や植物資源との関係で説明できるか。水陸の並列応答を比較する新規解析と、物質連鎖を問う既存の経路解析を分けて統合する。

フランス73小流域・24親流域の新規解析では、水生昆虫の低い生態状態と植物の栄養立地構成が空間的に対応した（観察条件調整後r=.514、主要4比較の近似同時95%区間[.147,.880]）。景観調整後は弱まり、下水指標Wの独自の説明寄与と2期間の同調変化は不確実だった。空間的対応を同調変化や下水由来物質の移行とは解釈しない。

水域内の予測情報と植物–潜在寄主構成の対応は支持される。一方、フランス共通139 H7×期間ではW→TP→植物→シジミHNSの係数積が0.00176［−0.0302,0.0377］で不確実。河岸の植物吸収・幼虫摂食・成虫タンパク質利用を実測した全鎖の証明には至っていない。

## 変更内容

- フランス73流域の期間平均・差分で、水域と陸域の残差相関と共分散を並列に推定。蝶全体を主解析、両期間適格67流域のシジミを補足とし、3,000回の親流域再標本化、全結果CSV、入力ハッシュ、比較図、本文2.8・3.6を追加。

- 本文、4パネルのパス図、経路別スケール表を更新。下水→河川環境→植物資源→植食者の既存経路図を保持し、底生EQR/BMWPは水域側の別応答として置く。新規の並列比較とは推定対象を区別する。
- 129 PRを17のレビュー系列に整理。初期の局所植物、EEA底生、独墺写真・河川追試、デンマーク蝶、米国湿地について主要結果と依存関係を追記。全効果量の抽出・全標本ID照合の完了とはしない。
- 肯定結果と追試の限界を併記。独墺写真の44セルでの予測改善17.82%を保持し、実測TPによる接続の不支持と区別。広域植物の属内関連を保持し、局所18対の形質適格性不足・欠測偏りも記載。
- 蝶5科→シジミの図は比較の構成を示す。シジミへの仮説上の関心はそれ以前にあり、対象選定の歴史を5科比較から始まったものとして描かない。
- フランス以外の公開一次資料も調査。北米・南豪州の水–植物、アフリカ・中東・欧州の処理水灌漑、一般栄養添加研究を、曝露経路と測定段階で分ける。Arizonaには全段階を扱う候補があるが、対応原表は未取得。一般N/Pメタ解析間の原研究重複も確認。
- 外部資料の監査を本文の方法2.7・結果3.5・考察4.5へ反映。同一標本の経路、比較可能な部分経路の効果、定性的な機構比較を分け、植物量・植物相・葉の栄養価・茎の取り込み・同位体を混同しない。ヤナギ原報告の取得・測定部位確認も収録。
- 物理的生息環境、水供給、汚泥農地施用を代替経路・条件として保持し、異なるデータの肯定結果から全鎖係数を合成しない。

## 検証した範囲

- 新規共変動解析：固定入力・一意キー・同一水域地点の差分・静的Wを確認。4組すべて3,000/3,000再標本でフルランク。合成データと共分散分解恒等式、別実装による全相関点推定、図表示を検証。気候は未調整、Wは時間変化を持たず、景観の説明は標本内分解。

- フランス139行：既存13経路と3,000回の同時再標本化を再現し、係数・区間の一致を確認。共分散・反復値を保存。
- ドイツ6放流口：固定PRの月対応TP表・植物表から既存対比を再構成。元全採水記録の再取得ではない。
- EEA2018/2022：固定地点別CSVのハッシュとキーを確認し、共通1,167地点、同じ評価法1,148地点、報告源状態の変化11地点を再確認。モデルの新規推定ではない。
- 18対の植物形質条件、デンマーク公開集計、米国湿地の源条件・標本所属を、それぞれ保存済み表から確認。何を原データから再現したか、何を原報告から引用したかを資料ごとに明示。
- 外部データは取得済み原本・掲載集計・原表未取得を分ける。15設計の測定範囲図を作成・表示確認したが、15独立研究や網羅的な世界レビューとはしない。

新しい全鎖メタ解析・因果媒介推定は実行していない。既存解析の再実行と監査を「解析なし」と一括せず、独立した新規生態研究とも数えない。

## レビュー対象

- [新規の水陸共変動解析・図・再現コード](https://github.com/zuizui0223/gesui/blob/docs/wastewater-riparian-interactions-outline/docs/manuscript_riparian_interactions/WATER_LAND_COVARIATION_20261120.md)

- [論文本文](https://github.com/zuizui0223/gesui/blob/docs/wastewater-riparian-interactions-outline/docs/manuscript_riparian_interactions/MANUSCRIPT_FULL_20261119.md)
- [図1 A–Dと関連資料](https://github.com/zuizui0223/gesui/blob/docs/wastewater-riparian-interactions-outline/docs/manuscript_riparian_interactions/README.md)
- [経路別のスケール・支持状態](https://github.com/zuizui0223/gesui/blob/docs/wastewater-riparian-interactions-outline/docs/manuscript_riparian_interactions/PATH_EVIDENCE_SCALES_20261120.md)
- [既存結果の依存関係](https://github.com/zuizui0223/gesui/blob/docs/wastewater-riparian-interactions-outline/docs/manuscript_riparian_interactions/PR_DEPENDENCY_INTEGRATION_20261120.md)
- [広域統合の設計](https://github.com/zuizui0223/gesui/blob/docs/wastewater-riparian-interactions-outline/docs/manuscript_riparian_interactions/CHAIN_INTEGRATION_PROTOCOL_20261120.md)
- [達成状況と未完事項](https://github.com/zuizui0223/gesui/blob/docs/wastewater-riparian-interactions-outline/docs/manuscript_riparian_interactions/GOAL_EVIDENCE_AUDIT_20261120.md)

原表の対応・比較可能な分散の取得、未確認元結果の監査、最終的な本文・図との整合確認は継続中。Draftを維持し、mainへはマージしない。
