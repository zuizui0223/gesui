# 独立蝶モニタリングでの属内HNS×植物NUTRI検証：ソースゲート監査（2026-11-22）

## 結論：独立の栄養環境反復は未成立（SOURCE STOP）

PR #138 の元のGBIF 183 H7×2期間における「植物NUTRIが高い小流域ほど同属のHNS高値種が相対的に記録されやすい」という探索的推定を、**異なる観察努力制度で直接検証**するため、元の実データを2系列調べた。両方で種別蝶データは見つかったが、**同地点の植物種別Ellenberg栄養値（EIV N）／実測土壌栄養値が確認できないこと**と、もともとの>=10蝶属の情報要件を満たさないことから、独立の植物NUTRI–属内HNSの新規係数は推定しなかった。継ぎ足し・代替曝露・閾値緩和による肯定値の探索は行わない。

| 条件 | デンマーク eBMS 2025 | ドイツ Sturm et al. (2026) |
|---|---|---|
| 第三者の観測単位 | eBMS 2025の調査地点・個体記録の既存集計 | 原始 13,718 観測行、187地点、原始 780地点×月の記録単位 |
| 既存の種名→数値HNS | 63受容蝶種のうち49種 | 原始110分類名のうち51受容蝶種 |
| HNSが異なる同属種2種以上 | 7属 | 7属 |
| 調査区間で実際に種比較が可能 | near/remote両群に各種記録：5属（同一地点の共存を意味しない） | 同地点で2同属種以上記録：6属 |
| 記載済み寄主のHNS採点率≥80%でも残る属 | 1属 | 1属（Aglais） |
| 同地点の植物NUTRI | **未確認** | **未収録**（植生の被度・植物種数は収録） |
| >=10情報属／植物NUTRI同地点の両条件 | **不成立** | **不成立** |
| NUTRIに対する独立推定 | **STOP** | **STOP** |

**ご注意：元の≥10蝶属・≥20 H6という基準はもともとH6を持つマクロ解析を実行するための必要条件である。別デザインでの妥当な検出力を直接定義するものではない。** 同じ基準を保守的な比較不能判定の参考にした。真の独立調査で新しい設計を事前に決めること自体は可能だが、肯定値を得るための事後的な再設定はしない。

## 1. デンマーク eBMS 2025

元データ：gesui [Draft #43](https://github.com/zuizui0223/gesui/pull/43) の既存の[Actions #42786494070](https://github.com/zuizui0223/gesui/actions/runs/42786494070)、artifact `13749857350` 内の受容63種別・群別実個体記録集計。元のeBMS v1.4 archive SHA256 `2cc64ec352754223db06a92eef5b01803968eac94bae34057c17f1abdb6da79b`、利用表示 CC BY-NC 4.0。元303蝶HNS表の受容・スコアあり・別名競合なしを**同じGBIF speciesKeyで照合**し、49/63種で使用可能。

同属HNS差のある7属は `Aglais, Aricia, Coenonympha, Lycaena, Pieris, Thymelicus, Vanessa`。ただし元の既存 near<=5km・remote>=10km 両群に比較種各々の記録があるものは5属 `Aglais, Lycaena, Pieris, Thymelicus, Vanessa`。HNS採点率80%以上を適用すると `Aglais` 以外は比較が成立しない。

**本質的なストップ理由：2022年登録下水処理施設からの距離は、蝶の2025年の実際の栄養曝露・河川下流への物質移送・寄主栄養・植生EIV Nではない。** さらに元の[Draft #42](https://github.com/zuizui0223/gesui/pull/42) では水路沿いの正しい受水域方向の実験単位が2つのoutfall/reachにしかならず元の物理曝露ゲートは不成立。距離群をNUTRIとして読み替えてはならない。near/remoteも同地点内の同属種交代ではない。元の登録施設の近接正負関連と、栄養環境に沿ったHNS選別は異なる仮説。

## 2. ドイツ草地 Sturm et al. (2026) 独立原始データ

公開CC BY 4.0 原データ：[LineSturm/2026_grassworks_butterflies](https://github.com/LineSturm/2026_grassworks_butterflies)、論文 DOI [10.1111/1365-2664.70513](https://doi.org/10.1111/1365-2664.70513)、データ DOI [10.5281/zenodo.20763921](https://doi.org/10.5281/zenodo.20763921)。関連するPANGAEA個体観測 DOI [10.1594/PANGAEA.982396](https://doi.org/10.1594/PANGAEA.982396)。

GitHubから**原ファイルを直接照合した**。Gitのblob SHA（SHA256とは異なる）は `data/raw/data_butterfly_raw.csv` = `b63c9611c527da85b992a25ce43dca216ae357ef`、`data/raw/data_site_environment.csv` = `0f7daaea51f506b85594de18144a4b11bc3fb0a8`、`data/raw/data_vegetation_variables.csv` = `94221b422d48ce7489eed053d54665dbec69a882`。

原始 13,718 観測行・187地点・110異なる `butterfly.name` ラベル（未特定名等を含む）。`cf == TRUE` は286行あり、**この不確実同定行と数値HNSのない種を除いて**11,128記録行・51厳密受容種が残った。187地点は蝶・環境・植生ファイルで地点IDが完全一致。780「地点×月」は掲載されている出現記録の組合せ数で、全観察実施日数・観察者の時間・非検出を含む訪問数ではない。

同属種2種以上でHNSが異なる7属の詳細：

| 蝶属 | 採点済み同属種数 | 1種以上記録された地点 | 2種以上が同じ地点で記録された地点 |
|---|---:|---:|---:|
| Pieris | 3 | 182 | 163 |
| Aglais | 2 | 135 | 51 |
| Thymelicus | 3 | 106 | 32 |
| Vanessa | 2 | 105 | 24 |
| Coenonympha | 2 | 170 | 11 |
| Lycaena | 3 | 70 | 4 |
| Pyrgus | 2 | 11 | 0 |

**環境指標は栄養NUTRIではない。** 同地点に植物種数、草・花被度、裸地、撹乱管理、湿潤区分、pH、`site.cwm.pres.oek.f` はあるが、当該3 CSVに植物名×植被率、土壌窒素、Ellenberg **N**があるとは確認できない。`site.cwm.pres.oek.f` を勝手に N に読み替えない。属数は7（同地点での複数記録に限ると6）、HNS取得率80%以上でHNS差がある2種以上はAglaisの1属。植生種数に対する回帰なら可能でも、それは栄養応答の再現検証ではなく**別の問い**。

## 3. より適した候補のソース確認（未取得データを取得済み扱いしない）

**Guariento et al., Biodiversity and Conservation 2023**, [From diverse to simple: butterfly communities erode from extensive grasslands to intensively used farmland and urban areas](https://doi.org/10.1007/s10531-022-02498-3)：**南チロル93地点**で同じ年に蝶・維管束植物調査、植物のpresence/absenceから温度・光・水分・**栄養**・土壌反応のEllenberg指標を算出。デザイン上は今回の植物NUTRI×蝶属内HNSに合う可能性がある。ただし**原始補足XLSXの実ファイルはまだ取得・検証できておらず**、種別蝶記録やHNS対応種の実際の数も未確認。推定値は未報告。

**Curtis et al., Journal of Applied Ecology 2015**, [Butterfly abundance is determined by food availability and is mediated by species traits](https://doi.org/10.1111/1365-2664.12523)：英国UKBMS 54地点、27種、1,111蝶個体群、同地点165植物種の被度と寄主資源指数を解析。先行研究で既に**狭食性ほど寄主植物資源量と蝶個体数の対応が強い**ことを報告（公開論文の食性幅×寄主資源の交互作用係数−0.171, p=0.014）。植生データ DOI [10.5285/ca34e25b-6138-4b72-b32c-778d762118f1](https://doi.org/10.5285/ca34e25b-6138-4b72-b32c-778d762118f1) と導出データ [10.6084/m9.figshare.1494733](https://doi.org/10.6084/m9.figshare.1494733) が示されている。**データをこの作業で新たに取得し独立再計算したわけではない**。寄主**量**と栄養生態**嗜好**は別物であり、この既知の結果を自分の新しいメタ解析に数えたり再発見を主張しない。

## 4. 次の検証仕様と停止条件

最優先は**同じ場所・期間で**①植物種のEllenberg栄養指標（既報の実サイトN値、または真正植物種リストから元表に沿って計算）、②HNSありの**同属蝶種2種以上の比較**、③調査日・訪問労力あるいは0記録の定義、④空間ブロックあるいは地点反復、がそろった第三者データのソース監査。ソースが足りなければ適格性を明示してSTOPし、ドイツの植物種数・デンマークの処理施設距離・植生湿潤指標を任意の栄養効果に代用しない。

**論文構成に関する含意**：主解析の約1.41倍はコントロールで大きく変わり、親流域外予測も失敗した。実際のモニタリングを使った独立再現の条件も未成立。したがって現時点で「栄養増加が蝶の同属種選別を普遍的に引き起こす」「下水を介する群集間相互作用の機構が同定された」とは結論できない。

数字の再計算・元ファイルのGit blob・全属支持範囲は [external_standardized_monitoring_source_gate.json](../results/20261122_butterfly_within_genus/external_standardized_monitoring_source_gate.json)。元の実データは重複配布せず原研究を必ず引用。Draft PR #138は未マージのまま。
