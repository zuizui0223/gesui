# 183 H7 / 40 H6：植物NUTRIの属間・属内成分分解（事後）

**2026-11-15。全種ラベル入替p=0.0012、同属内のみ入替p=0.6864という結果を見た後の事後分解。独立検証・事前登録ではない。**

## 原本と数値の恒等式

- 元の183 distinct H7 / 40 genuine H6、2期間GBIF種出現のCWMを合算した既存パネルSHA256 `c0840b1c5c6e85342fea41d57dc40fa5ec295fabf6bcd6908a4072557db9a736`。元のH6固定効果と13 source controls、植物4,011 scored accepted vascular taxa/正式な`gbif_canonical_name`/NUTRI scoreは既存原本からのみ使用。
- EIVスコア付き種sの属gについて、厳密に`N_s = mean(N|genus=g) + (N_s−mean(N|genus=g))`。属平均は**今回の真正4,011種の数値scoreの単純平均**。元研究にない属名や数値は補完しない。属当たり1 scored speciesなら属内差は0。
- 各H7の原本2期間における得点種均等の重み`a_{i,s}`を使用。`PlantN_total_i = Σ_s a_is N_s`を再構成し、`PlantN_between_genus_i = Σ_s a_is mean(N|g_s)`、`PlantN_within_genus_i = Σ_s a_is (N_s−mean(N|g_s))`と定義する。
- 全183行で`total == between + within`、元source plant_n_incidenceを絶対誤差≤1e−10で厳密再現しない場合STOP。元H6/H7/observational HNS/13 source controlsは変更しない。

## 主な検査（どの数値を見ても全部報告）

1. H6固定効果・13 controlsで調整した蝶の種均等HNSと、植物NUTRI *total*、*between genus*、*within genus* それぞれの加重偏相関と基準モデルからのin-sample partial weighted SSE改善。全40 H6 /183 H7、H6均等重み。
2. 同じ設計にbetweenとwithinを両方入れた**多変量**H6クラスタWLSの個別標準化β、記述95%区間、モデルrank、condition number、変数間相関と部分共線性を報告。非同定ならSTOP、指標を一つ省略して有利なモデルにしない。
3. betweenだけモデルとbetween+withinモデルの**追加in-sample weighted SSE増分**を数値として報告。「追加説明できる率」の統計であり、**地理的heldout予測・個体数・因果の説明率ではない**。
4. 元の全NUTRI WLS β+0.118076との加算分解は、固定の元total NUTRI residualを用いて回帰係数分子の共分散をbetweenとwithinに分ける**数学的寄与**としてのみ行う。独立した「媒介割合」「選択圧割合」「進化効果」にはしない。
5. 副次として蝶記録重みHNSで同様の偏相関を計算して全結果保持。主の蝶種均等HNSを上書きしない。
6. 後から属を選択・除外しない。例示する場合、全属単独・複数得点種の数、属内N分散の四分位点はsource-onlyで記述し、予測可能な高効果属を見て選ばない。

## 解釈限界
属間成分が主・属内成分の追加情報が乏しい場合、先のラベル入替と整合し「大きな属レベルの生息地嗜好の地理的構造」が主な記録組成関連である可能性が高まる。ただし系統共通祖先・属そのものへの選別と属に保守的なN形質を切り分けられない。属内スコア差の検出力は原本種数・変動量で制限される。

元のモデル化排水関連地理Wは同一H6で明瞭ではなく、実測局所WWTP放流水N、葉組織N、実際の幼虫寄主使用・植食・昆虫個体数、上位捕食者、サービス、生物エネルギーフラックスや因果媒介は測っていない。元の時間変化の非同期と10 H6外部予測の不確実性は変更しない。PR #105上で検証し、main/chocho/親PRを編集しない。
