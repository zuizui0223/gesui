## Posthoc secondary amendment after six-trait STOP — 2026-11-15

**この追補は厳格な6形質同時比較が41 H7×期間・16 H6でSTOPした後の新しい事後設計であり、前段のSTOPを取り消すものではない。**

元の6形質同時モデルを推定しない。独立した二形質の共通植物種検定を**5対全てに実施**する。

- 元の11,511 speciesKey、元454 H7×期間、元のEIV source SHAと、既存4,011 accepted plant N scoresを固定する。NUTRIの比較対照は原本のLIGHT/TEMPE/MOIST/REACT/SALINの**5つ全部**。
- **対ごとに厳格な共通形質種集合**：元植物speciesKeyのうち、NUTRIと当該対照の両方に数値がある種。site×periodごとに両CWMを同じ植物種で計算する。異なる指標間で別々の植物種集合を使った単純比較は禁止。
- 各対で元植物記録数≥50、共通数値種≥10、元全Plantae speciesKeyに対して共通形質種率≥40%、さらに群全体≥100 H7×期間、≥40真正H6、≥5 country_holdout groups の解析前ゲートを満たした場合のみ分析する。
- 分析コホートは元451 EIV-NUTRI supported H7×期間の部分集合。対ごとに標本が異なりうるので、対をまたいだβやMSE改善の大小を「Nの特異性のランキング」として使わない。各対で**同一サイト・同一得点植物・同一モデル式**を比較する。
- 蝶の種均等HNSが主応答。固定元W/geography/effort + period + country×method + pair-specific coverage + pair-specific log(scored_species)をベースとする。βのH6-cluster CI、H6完全非重複5fold Ridge alpha10予測、2000回固定H6予測損失再標本化を用いる。
- **各対の主反証：** `NUTRI` を当該対照CWMのあるモデルへ追加（N|other）したβ、H6-out MSE改善。反対向き `other|N` も同じコホートで報告。5つのN|other主検定にBH5、他の5対にもBH5を明記。Nと対照の群集CWM間相関、collinearity VIF_N>10は警告。
- 5対の全解析が良好でも、全5対照を**同時に**調整したNの独立性は実証されない。全6共通種の厳格STOPは維持する。
- 同一H7の時間差分のゼロ、PR #96のEQR↔HNS不確実、下水の局所Nの未測定、記録組成≠生物個体数、EIV habitat preferenceと幼虫真の摂食は異なる、を維持する。main/chocho/親PRは変更しない。
