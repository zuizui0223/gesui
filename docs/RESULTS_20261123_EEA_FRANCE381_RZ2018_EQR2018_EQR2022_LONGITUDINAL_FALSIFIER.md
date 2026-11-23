# 2018年の河畔土地被覆は2022年の底生動物EQRを追加予測するか――実同一河川381地点の縦断的反証

**2026-11-23｜フランス国内・実同一地点・同じBQE評価方法・元河道同一の実データ解析。探索的(posthoc)であり、河畔植生操作・実下水負荷変化の因果推定ではない。**

## 問題の再定義とこの解析の新しい点

もともと2018年の欧州13か国976地点で、**同じ河川観測地点の250–500m圏の公式Riparian Zones LCLU2018土地被覆**を2018年の底生無脊椎動物EQRに加えると、地理的別地域への予測MSEが地点均等 **+3.423%** 改善した([原PR #123](https://github.com/zuizui0223/gesui/pull/123))。しかし国別外部予測には不確実性があり、上流WWTP P×近岸植生の交互作用は追加予測 **+0.250%**、地域CIはゼロを含む。土地被覆はバイオマスや水→陸栄養流入ではない。

**新しい問いは時間的なもの**：「2018年の河畔景観は、同一河川における**2018年EQRで把握済みの生態状態を超えて**、2022年EQRを予測するか」。同年の地理的群集差を4年後の動態や富栄養化伝播と混同しないための反証的分析。

## 本当の同一地点・時間Sourceゲート（数値生態状態は未使用）

- 原2018 EEA Waterbase Biology benthos BQE `EEA_13-01-4` 真正2,634モニタリング地点、原数値表SHA256 `2a1a72d43656aca588a832dfae13860054c01a651e76bf25f984c5a374709b02`。[PR #17](https://github.com/zuizui0223/gesui/pull/17) source-reporting/Q/HydroRIVERS付き原Git.
- 2022 EEA真BQEは[PR #14](https://github.com/zuizui0223/gesui/pull/14)のoriginal 2,475地点、原Git blob `325f199e985d9862ffee5b0edb205c61a7be3833`。新しいActionsは同一リポのそのブランチから `git show FETCH_HEAD` で当該blobを取得し、`git hash-object` **厳密一致**。元生物EQRはSOURCE確認前に開かず、サイトID/国/評価方法/HYRIV_ID/座標だけでSOURCE適格を判断。
- original official Copernicus RZ2018 Level 1 56-class land cover、河川地点から **0–250mと250–500m**円環、元source近遠各>=70% classified polygon、実非水域地表>=10%の指定を固定。真2018 `RZ` site SOURCE **978**のうち、原2018数値EQRを持つのは**976**（残り2は原EQRなし）。原source-only [Actions #44004357868](https://github.com/zuizui0223/gesui/actions/runs/44004357868)、source ZIP artifact ID **14293433193**。
- 元の2018と2022両年に同じ`site_key`があるのは **1,167**地点。原2018/2022 `HYRIV_ID`同一 **1,167**。評価方法の旧SOURCE品質表「同一1,148」には**両年ともmethod空欄3サイト**が含まれていたので、それらを安易に同一手法と認めず、両年**非空欄かつ同じ方法なのは1,145**。元source上「同名評価法1,148」は文字列の空欄同士の一致を含む旧集計として併記する。新たに有利な群を探さない。
- 2018 mapped habitat sourceと2022再訪を重ねると **418**原同一site、そのうちsource年・実同一BQE法・HydroRIVERS河道ID・0.01度以内の元座標一致で **409**。国別FR **381**、ほか計28（HR8, NL9, NO8, BE2, IE1）。**多国の代表的な縦断群を作る条件はSOURCE STOP**。しかし国別のFrance sourceは事前指定の最低150site、5原3度ブロック、3つのsource-fold各>=20を達成: **381サイト、14独立3度ブロック、5元地理fold =122,58,48,33,120**。FR source-only fallbackだけGO。
- 2018上流報告されたP>0 **15/409**、2022報告P>0 **19/409**、**両年>0は14/409**。処理場のN/P負荷変動・排水の変化とEQRの因果反応を測る標本としては全く不十分。報告ゼロは実際のsource未曝露を意味しない。
- [2018→22独立SOURCE GO成功 #45095545359](https://github.com/zuizui0223/gesui/actions/runs/45095545359)、集約source artifact **14783096791**。[SOURCE-only事前固定仕様](https://github.com/zuizui0223/gesui/blob/analysis/eea2018-22-riparian-eqr-longitudinal-source-20261123/docs/PROTOCOL_20261123_EEA2018_RZ_LAND_TO_2022_BENTHIC_EQR_TEMPORAL_HOLDOUT.md)。

## 生態状態の実測時間変化と元予測仕様

- 実同一French site n381 の元**2018 EQR平均 0.68938**、元**2022 EQR平均 0.66809**、`ΔEQR2022-2018`平均 **−0.02129**。これは共通siteの記述差であり、河畔植生や下水による因果的低下・有意なトレンドと解釈しない。
- 同一381地点の2018/2022 EQRのPearson r **+0.81516**。異なる国ごとに方法が違うEQRを混ぜない。
- 事前指定の**元2018地理source 5fold**、地点ID重複なし、同じ3度ブロックのsiteが別foldに分かれないsource-leak検査PASS。train foldだけ標準化、Ridge alpha10（データサイズ381、3国重み合計3の外部4か国α10とは尺度が違う）。
- 生物学的応答 `y2022`、モデル `B=2018同site生態EQR+2018 lat/lon/logQ/logArea`、`L=B+2018 250–500m natural-like forest+managed grass+heath+wetland/urban/cropland fraction`、`D=L+2018 0–250 minus 250–500m natural-like class fraction`。元Copernicusは実土地被覆地図だが**真の葉栄養・バイオマス・2018→22土地被覆変化は測っていない**。

### 主結果：2018 EQRを知ったうえで4年後を予測する（時間的な追加予測）

| Future EQR2022 heldout: 増やした説明変数 | 381サイト均等 MSE改善 | 14地理ブロック均等 MSE改善 | 14ブロック固定OOF error bootstrap95（記述区間） |
|---|---:|---:|---|
| `B + 2018の広いLandContext L` | **−1.540%** | **−2.640%** | **[−6.440,+0.345]%** |
| `L + 河川から近い植生割合の差 D` | **+0.025%** | **−0.194%** | **[−1.273,+0.819]%** |

5foldで元の2018 EQRを入れたbaseline BのサイトMSE **0.01752151**、さらにLを入れると **0.01779126**。局所近岸Dまで足しても **0.01778674**。Lの残差予測信号は支持されず、Dは実質ゼロ・地域差の大きい不確実な値。14地域のうちL追加でMSEが減ったのは**3地域**のみ。水生生物の年内空間差を説明する景観と、年をまたぐsite単位の生態状態変化を説明する景観は同じ主張ではない。

**固定のΔEQR endpoint感度**：2018 EQRをbaselineに入れて `ΔEQR=EQR2022−EQR2018` を直接予測する場合も、L追加 **−1.934% site／−2.770% region [−6.775,+0.495]**、D追加 −0.011% site／−0.220% region [−1.194,+0.622]。方向・解釈は主結果と一致する。ΔEQRと将来EQRは代数的に重なった同一データの別表現で、独立した外部再現ではない。

[真正原EEA source同一地点の生物数値/空間CV Actions #45096054949 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45096054949)、集約original result artifact ID **14783103096**。

## 重要なposthoc：すべて「同じFrance381地点」で空間差と将来予測を比較

最初の欧州976siteとの比較では「標本集団が違う」ため原因が曖昧だったので、**上の負の将来予測結果を見た後で**新たに[fixed same-source診断仕様](https://github.com/zuizui0223/gesui/blob/analysis/eea2018-22-riparian-eqr-longitudinal-source-20261123/docs/PROTOCOL_20261123_FRANCE381_SAME_COHORT_SPATIAL_VS_FUTURE_EQR_FALSIFIER.md)を作成。source fixed381、2018実同じ土地被覆/5地理fold/元河道/同site EQRに対し、baseline B0=2018 lat/lon/logQ/logAreaのみ（**2018 EQRを入れない**）として、Lを追加した場合を比較。モデル選択やalphaを結果によって変更しない。以下の3つは**違うbaseline**であり、「同じ量の効果」が動いたのではない。

| 本当に同じ381site / 同じ2018 RZ /同じfold | Land Lを追加したsite均等MSE改善 | 14region均等MSE改善／95記述区間 |
|---|---:|---|
| **2018 EQR横断地理モデル**（過去の生態状態を説明変数に含めない） | **+2.469%** | +1.959% [−12.349,+12.703] |
| **2022 EQR横断地理モデル**（2018 EQRを入力しない） | **+5.135%** | +4.013% [−5.123,+14.050] |
| **2022 EQR将来モデル**（**2018 EQRを入力して追加予測を問う**、一次） | **−1.540%** | **−2.640% [−6.440,+0.345]** |

追加の河川近接0–250m vs遠方250–500m natural-like **D**：
- 2018 no-lag EQR site改善 **+0.168%**、region −0.050% [−2.426,+1.871]。
- 2022 no-lag EQR site改善 **+0.397%**、region −0.100% [−2.760,+2.175]。
- 2022 lagged EQR site改善 +0.025%、region −0.194% [−1.273,+0.819]。

**結果の意味**：同一フランスsiteでも2018年の土地被覆は、2018と2022両方のsite間の生態状態の違いに**ある程度の地理的予測情報**を持つ。一方で2018のEQRを知ると、2018年土地被覆は2022状態に**独立した追加予測をもたらさない**（事前に固定したモデルにおける結果）。この整合的な3-way比較は、景観–水生群集の横断相関をそのまま4年間の生物動態の機構へ拡張することの危うさを、同一site/sourceで数値的に示す。

しかしregion95は3つとも0を含む、過去の生物状態と長期の土地被覆が共通の環境特性を表すだけかもしれず、**「河畔林の生態機能がない」「植生が生態変化を防げない」とは言えない**。水生から陸生への羽化、陸上バイオマス、現地植生時系列はデータにない。

[同じsiteに制限した事後的な2018/22空間source元数値 #45096506660 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45096506660)、source-signed aggregate artifact **14783138192**。

## 独立コードによる数値監査

- 2022 lagged本解析を元コードとは別の**手計算NumPy 5fold train-only標準化+正規方程式Ridge alpha10**で再実装し、14原地理ブロック固定予測誤差+2,000block bootstrap、2 outcome（実2022 EQRとΔEQR）×2比較 **計4項目を全て再現**：[Actions #45096295649 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45096295649)、signed独立QA artifact **14783486858**。
- 381site同一サンプルの断面2018/22 B0/L0/D0の追加数値についても、別実装の**NumPy手計算Ridge10＋14地域誤差＋同じbootstrap**で4つの地理的MSE比較（2018/2022それぞれL0・D0）を完全再現：[Actions #45096643734 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45096643734)、signed独立QA artifact **14783028962**。
- この「独立監査」は本当の新しいフィールド調査ではなく、同一のoriginal EEA Biology/Copernicus元データに対する**異なる計算コードの再現**である。元のRZ2018 source、2018 EQR、2022 original Git blobのSHAを固定し、出所が異なる1975/2018/etcの研究を無理に混在させない。

## 下水研究への統合範囲と生態学的な一般性

**GO**：同一物理河川地点をつかった水域（底生生物の生態状態）–陸域（河畔の土地被覆）の**時系列的な環境–群集予測対応**。環境条件の空間的な効果は、独立した将来変化の予測には直接つながらない、という*methodological ecological caution*を支える探索的な外部年比較。

**SOURCE STOP**：
1. WWTP報告のN/P source陽性は409地点中両年14地点のみ、下水源の時間的変化の生物効果・nutrient-specific mediationは同定不可。2018/22の「非報告」を未曝露・真のcontrolに置き換えない。
2. 2018のRZ地図は実植物個体数・土壌Nや草地 biomassの測定ではなく、2022の同地点土地被覆変化も実測していない。河畔植生によるストレス緩衝／自立的な生物群集動態を因果的に主張できない。
3. 本解析の381地点はフランスに偏ったSOURCE適格空間。13か国976地点や2018/2022全標本への人口代表性を主張しない。2018 RZ源の分類実施とBQE評価法は選別された河川が対象で、5fold地理検証も完全な独立生態学的試験ではない。
4. **生態系サービス（浄化・送粉・昆虫バイオマス）、水→陸の食物網（羽化–捕食・蝶寄主）と下水処理場の物質移動は何も測定されていない**。US/EUの前の否定・STOPを削除しない。

## 投稿向け結論（限定する）

> Same-site riparian landscape composition retained cross-sectional predictive information about river macroinvertebrate ecological status in 2018 and 2022, but did not improve the spatially held-out prediction of 2022 ecological status once the observed 2018 status was included. This contrast was confined to 381 French monitoring sites and does not establish a causal influence, absence of riparian functions, or a wastewater-specific ecological pathway.

主張の軸を「陸から水への物質伝播が証明された」ではなく、「**同時期の空間相関を時間的な群集応答だと解釈してよいのか**」に置くことが望ましい。しかしこの単一の探索的フランス解析をもって普遍的な時空スケールの原理として単独論文化するにはまだ弱い。独立した反復期間／土地被覆の真の変化・同じsite群集の検出率と実源の量を追加できる場合だけ次段階へ進める。

**既存原稿・`main`・WWTP/source停止判定を改変しない。すべてreview-only Draft。**
