# gesui — Wastewater ecology / 下水による生態系の再編成

**中心的な問い：下水・処理水・栄養塩の人為的な再配分は、植物群集と水域・陸域の生態学的連結性をどう変えるか？**

蝶中心の [chocho](https://github.com/zuizui0223/chocho) から独立。第1の応答は**全植物相**。水生植物、河岸植生、淡水無脊椎動物、昆虫、河口・沿岸の生物へ段階的に拡張する。島や港は比較する地域の一部であり、研究の必須条件ではない。

## 引き継いだ探索的な結果

- HydroWASTEの下水処理場数25kmと外来／在来「候補寄主植物のGBIF記録種数比」は、国・大まかな地理の基本モデルで β **+0.1285**、95% CI [+0.0077,+0.2494]（1,368港、52か国）。記録量の代理で補正すると β +0.097、区間はゼロを含む。
- 2015年の**モデル下水発生量 WWp**は都市・農地・WWTP数などを調整して β **+0.05870**、CI [+0.01310,+0.10429]（1,055港、42か国）。R²の標本内改善は+0.00171。WWpは窒素負荷ではない。
- Ellenberg Nと植物の記録港頻度は偏順位相関 **ρ +0.210**。欧州6島の同一27沿岸ペアでは、候補寄主とその他の植物の相対N記録組成差が期間間で **+0.1494**（記録重み付けに敏感）。
- NIESの国別衛生経路N・CIRAD人由来Nポテンシャルでは、**下水由来Nが植物を増やしたという独立した機構は確認されていない**。

ここに記載した値は元の [chocho Draft PR #122](https://github.com/zuizui0223/chocho/pull/122) の監査結果。今回新しい回帰を計算した値ではない。

## ナビゲーション

- [研究計画と事前の検定](docs/RESEARCH_PROGRAM.md)
- [肯定結果・対照・限界の証拠台帳](docs/EVIDENCE_LEDGER.md)
- [データの意味と再現条件](docs/DATA_AND_REPRODUCIBILITY.md)
- [chochoとの分業・移行記録](docs/MIGRATION.md)
- [80ファイルの原本保存](archive/chocho_pr122)
- [完全一致したGit blobの照合台帳](provenance/CHOCHO_PR122_FILE_MANIFEST.json)

元PR122から **80/80ファイルを内容変更なしにコピー**。内訳：25旧Actions YAML、21解析スクリプト、18テスト、16監査・結果。旧Actions設定は archive 下にあり自動実行されない。**外部原データ・GitHub Actions artifact・旧chocho共通入力は複製していない**ため、移植先での独立実行はまだ未検証。

## Offline verification and new all-plant table contract

- `python scripts/verify_migration_integrity.py` verifies all original 80 Git blob SHA values **without network access**.
- `python -m unittest discover -s tests -p 'test_gesui_*.py' -v` runs the archival check and **synthetic-only** schema quality tests.
- `python scripts/validate_plant_site_panel.py data/site_period_template.csv` validates a site × period × plant taxon CSV schema. The supplied template is header-only, **not empirical data**. Measured/report-level/modelled nutrient mass is distinguished from modeled wastewater volume.
- `.github/workflows/gesui-offline-integrity.yml` runs these checks without invoking any archival HydroWASTE/GBIF/PANGAEA re-download workflow.

These tests establish preservation and input hygiene, **not** independent re-analysis of wastewater–flora relationships.

## Competing biological mechanisms

The new [water-subsidy vs nutrient-enrichment vs contaminant-stress framework](docs/MECHANISM_CONTRASTS.md) separates effluent-supported dry-season habitats from N-driven compositional filtering and water-quality stress. This is a **pre-analysis hypothesis set**, not an empirical result from chocho.


## 実際に完了した新規解析（2026年11月7日・探索的Draft）

- [河川ネットワークの検証・PR #3](https://github.com/zuizui0223/gesui/pull/3) — HydroWASTE 58,502施設の原本確認、欧州の条件適合4,713施設の **HYRIV_ID が第三者版HydroRIVERSの河段に4,713/4,713照合**。4,497件の上流2→放流河段→下流1の連続河段、3,569非重複候補から植物を見る前に16放流口を選んだ。推定放流口と対応河道のずれは中央値約25m、最大約47m。公式版の同一バイト／全トポロジーは未照合。[実データGIS CI成功](https://github.com/zuizui0223/gesui/actions/runs/42364866919)。
- [河川沿いの植物の初回観測・PR #4](https://github.com/zuizui0223/gesui/pull/4) — 16放流口×上流2地点／下流1地点=48地点の**実際のGBIF Plantae記録**、**1,402別speciesKey、5,394地点×植物種記録行**を取得。主解析は希薄化30記録＋観測量の均衡を満たす12放流口。モデル処理水割合と「下流変化−上流背景変化」の関連 **ρ=+0.1259（n=12、探索的p=0.6967）**、10記録の感度分析 **ρ=+0.1736（n=14、p=0.5528）**。明確な下水固有の植物群集効果は未確認。[元GBIF実行成功](https://github.com/zuizui0223/gesui/actions/runs/42365305227)／[独立5,394行の監査成功](https://github.com/zuizui0223/gesui/actions/runs/42365622805)。
- [先行20地点・PR #2](https://github.com/zuizui0223/gesui/pull/2) — 推定放流口対施設位置の植物記録（河川上下流比較ではない）。正の一貫した関係なし、ρ=−0.150、n=9、p=0.700。
- **実測栄養塩の新規ソース検証**：[Draft PR #6](https://github.com/zuizui0223/gesui/pull/6)／[Issue #5](https://github.com/zuizui0223/gesui/issues/5)。[全球河川栄養塩データPeters et al.2026 v5](https://doi.org/10.25400/lincolnuninz.30546371.v5)で実測TN/TP観測局と16放流口の上下流を照合。**1放流口のみ**TN/TP双方の上下流観測局メタデータが対応（河道距離25km以内）。ただし実際のTN原本は2局とも公開不可で**TN/TP同時主検証は実行不能**。TP単独は2012–13年24か月、下流／上流中央値1.167倍（2008–22年179か月、1.406倍）という河川濃度勾配が実測できたが、下水由来とは特定できない。
- **新しい独立6河川の実測TP×植物記録**：[Draft PR #7](https://github.com/zuizui0223/gesui/pull/7)。植物応答もTP濃度も見ずに元3,569河道候補から別の3,553候補を組み、**12個の独立水質観測局／異なる6流域／最短57km**のドイツ6放流口を固定。2011–22年**2,500件の実測TP、836か月分の上・下流ペア、GBIF植物18地点**を新たに取得。モデル処理水割合と実測TP勾配は**ρ=+0.829、n=6**（通常p=0.042、720順位置換の**正確な両側p=0.058**、LOO ρ+0.70〜+1.00）。植物種数勾配との対応は厳格30記録では3/6のみ、10記録感度で**ρ=−0.029、n=6、p=0.957**。**TPの探索的な正の空間関連はあるが、下水からのリン供給や植物群集への効果は立証していない。** [原本取得](https://github.com/zuizui0223/gesui/actions/runs/42371970951)／[独立再計算](https://github.com/zuizui0223/gesui/actions/runs/42372402396)。

**重要**：新しい実測結果はDraft解析ブランチにあり、mainへの結果マージはレビュー前に行っていない。都市化、土地利用、他の放流源、植物記録努力と河川距離の交絡は残る。小標本での探索的相関を、下水由来リンの因果効果や一般法則とはしない。元chochoの既存肯定結果は維持する。

## 欧州河川の下水インフラと実測栄養塩のマクロ検証（2026-11-07／Draft PR #8）

独立した新規の[河川マクロ解析 PR #8](https://github.com/zuizui0223/gesui/pull/8)で、Peters 2026 v5 の**実測**TP 1,418観測局／TN 236観測局を、HydroRIVERSの流下方向とHydroWASTEの上流50km内の**モデル処理水量**に接続した。原本のsource SHA/MD5検証、空間3°ブロック5分割の予測、上流と下流の負の対照、同一標本の測定値を固定し、外部API再照会なしの[独立数値監査SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42372545545)を完了。

- **TP全球型パネル n=1,418**（上流50kmに登録処理場があるのは85局）：実測TP濃度×上流処理水量の順位相関 ρ−0.00192。地理・河川Q/流域面積ベースラインより、上流処理水量追加で空間外MSEは**4.00%悪化**。下流の処理場量という負の対照を追加すると2.47%改善するため、人為土地利用やインフラ配置の交絡が残る。
- **独立の農地・人口密度情報があるTP n=510**（上流処理場あり73局）：同じ空間ブロック外予測のMSEは上流処理水量を足しても**1.09%悪化**、下流負の対照の改善は0.11%。TP元ソースにはN_deposition列がないため、最初の実行で510地点をゼロとして誤除外したコードを、既に固定していた共変量仕様に従って修正し、同じ原本で再計算した。
- **TN n=236**（上流処理場あり13局）では上流水量追加でMSEが0.66%悪化。農地・人口・窒素沈着がそろう局は20、そのうち上流処理場あり1局しかないため、TNの土地利用補正下の一般化は未検証。

**解釈**：上流50kmの登録処理場からの水量だけでは、実測TN・TPの地理的違いを独立に予測できない。これは実測した下水のN/P放流量の効果を否定する結果ではなく、**水量≠栄養塩の質量負荷**、実際の処理段階・流量時期・他施設・農業／都市化などの分離が必要という結論。上記マクロ負の結果と、[別標本のドイツ独立6河川の探索的な正のTP勾配（PR #7）](https://github.com/zuizui0223/gesui/pull/7)は、対象・曝露尺度・観測単位が異なるためそのまま置き換えない。いずれも生態学的因果証明ではない。


## 2022年の実測下水N・P排出負荷と新しい植物検証（2026-11-07、Draft）

- [**EEA 2022原本の実測N・P排出負荷×河川水質 PR #9**](https://github.com/zuizui0223/gesui/pull/9)：EEA公式UWWTDの処理場と一意の淡水（FW）放流口を照合し、**2,969施設**の報告負荷を河川網に結合。Peters 2026 v5の**同じ2022年**の河川TP実測 **1,291観測局**／TN **159観測局**と接続。上流50kmに対応する2022実測負荷があるのはTP **53局**／TN **26局**。固定3°地域ブロック予測では排出負荷追加によりTPのMSEが**1.97%悪化**、TNが**3.33%悪化**。TP農地・人口補正446局でも**1.33%悪化**。報告排出負荷に強い独立予測寄与は確認できない。**排水そのものの効果がないという主張ではない**。[元解析成功](https://github.com/zuizui0223/gesui/actions/runs/42375542151)・[地域分割を固定した独立監査成功](https://github.com/zuizui0223/gesui/actions/runs/42376519412)。
- **年間平均希釈の物理的な診断**：上流に報告負荷がある局に限り、処理場からの全報告量が年間を通じて河川平均流量に完全混合したと仮定すると、実測河川濃度に対する仮想濃度比の中央値は**TP 4.07%（53局）／TN 2.39%（26局）**。実測した下水由来割合でも季節的影響の上限でもない。**流量・希釈は下水の生態学的効果の重要な条件候補**。[原本・2022観測局CSV・監査・固定5分割を保存](https://github.com/zuizui0223/gesui/tree/analysis/eea-wwtp-nutrient-load-source-20261107/results/eea2022_measured_nutrient_load_river_20261107)。
- [**EEA実測P排出源あり／登録源ゼロの河川×全植物 PR #10**](https://github.com/zuizui0223/gesui/pull/10)：水質・植物の結果を選別に使わず、流量・流域面積・農地割合・人口密度と3°地域区分で**18組・36地点**を固定。新しいGBIF 2022–2024 Plantae記録から**1,006 speciesKey**を取得。厳格な30記録希薄化・記録均衡を満たすのは**5組のみ**で主解析の推測力不足。10記録の感度解析は11組、源あり側の記録種数が多いのは7組、**両側の正確符号検定 p=0.549**。植物個体数、植物の栄養選好性、土壌N/Pや真の曝露は未測定。[植物原本取得成功](https://github.com/zuizui0223/gesui/actions/runs/42377454765)・[独立再計算成功](https://github.com/zuizui0223/gesui/actions/runs/42377732661)。

**まとめ**：ソース品質・水文接続・同一年の実測河川水質・独立した植物記録までは統合できたが、下水特有の栄養塩による植物増殖／群集改変は未実証。次は種数だけでなく、**植物の栄養塩選好性、種組成の入れ替わり、低水量時の放流水寄与、他栄養源**に分ける。港・島・風・AIS/航路・海上Circuitscapeは引き続きchochoに残し、既存の肯定結果を保持する。

## 2026-11-08：植物NUTRI形質の監査と水生生物モニタリングへの拡張

- [**欧州Ellenberg栄養選好性の照合・Draft PR #12**](https://github.com/zuizui0223/gesui/pull/12)：18組・36地点の原本GBIF **1,006 speciesKey**に対し、受理済み維管束植物911種キー、欧州NUTRIとの完全一致624種キー。形質品質を満たす地点は25/36、比較ペアは**9/18**で事前の最低10組に届かず、下水効果の主検定は未成立。9組の記述値は源あり−対照 **−0.255 NUTRI**（3組正・6組負）で、植物の栄養選好性が低下したとは結論しない。
- **重要な観測・分類データ品質の結果**：源あり側のEIV照合率は**18組中15組で低い**（事後的品質診断の両側符号p=0.00754）。源あり側のGBIF記録数は13/18で多い。栄養指標の欠測と観察努力が非対称なので、標準化した植生資料や原産性情報を補わずに下水による植物群集の選択を主張しない。[原本と独立再計算成功](https://github.com/zuizui0223/gesui/actions/runs/42430531522)。元の欧州形質表全体や種別の数値を転載せず、出典・SHA・分類照合フラグ・地点要約のみ保存。
- [**EEA Waterbase Biology 2026**](https://www.eea.europa.eu/en/datahub/datahubitem-view/7bd881f0-59ca-4f0a-a43c-7e0cd0c66d6c)：2012〜2024年の河川底生無脊椎動物・魚類・水生植物・付着藻類などの**生物学的品質比（EQR）**を扱う公的データを次の独立応答に採用。2026版公式ZIP 132.2MBと内部EQR・位置情報の原本構造を確認。[ソース監査](https://github.com/zuizui0223/gesui/actions/runs/42430869887)。2022年の地点・分類群と上流のEEA実測排出負荷の照合を完了し、次の節の独立した水生生物EQR解析へ進んだ。


## 2026-11-08：実測下水排出負荷と河川水生生物の初回マクロ解析（Draft PR #14）

[**独立した水生生物EQR解析 Draft PR #14**](https://github.com/zuizui0223/gesui/pull/14)で、EEAの**2022年処理場別実測N・P排出負荷（年間t）**と、同年の河川の**底生無脊椎動物・付着藻類・魚類の生物学的品質指標EQR**を河川ネットワーク上で接続。生物結果を見ずに、2022年の調査地点・国・評価手法・50km上流／下流排出負荷・河川流量・流域面積・3°地域ブロックの5分割を固定。[物理ソース照合成功](https://github.com/zuizui0223/gesui/actions/runs/42432030617)／[実測EQR取得成功](https://github.com/zuizui0223/gesui/actions/runs/42432384196)。

- **底生無脊椎動物EQR（2,475河川地点、13か国、上流排出源あり147地点）**：国・生物評価法・地理・河川Q・集水面積を含む基準モデルに上流の実測P排出質量を追加すると、地域外予測MSE **+1.10%改善**、N排出質量では**+1.40%改善**。物理的に**下流の**処理場を入れた負の対照ではそれぞれ**0.35%／0.43%悪化**。予測改善した3°地域ブロックはP **32/45**、N **31/45**。排出負荷とEQRの単純順位相関は約ρ−0.072であり、下水が生物に良いという証明ではない。
- **付着藻類EQR（2,784地点、12か国、上流源あり110）**：予測改善はP+0.16%、N+0.18%と小さく、国単位の予測・非線形地理の事後感度分析で再現せず。**魚類（1,243有効地点）は4か国のみ**で事前の最低5か国を満たさず、予測モデルを適用しない。
- **再現性**：元データのファイルSHAを固定した別コードの独立再計算が[成功](https://github.com/zuizui0223/gesui/actions/runs/42436814278)。原本と[ソース・実測EQR・地域ブロック・監査一式](https://github.com/zuizui0223/gesui/tree/analysis/eea2026-waterbase-biology-bqe-20261108/results/eea2022_biology_bqe_20261108)を保存。国別除外／非線形地理の**事後感度分析**も[成功](https://github.com/zuizui0223/gesui/actions/runs/42437061882)：底生無脊椎動物では国を丸ごと除いてP+2.47%、N+3.00%、非線形地理を入れてもP+0.79%、N+1.00%。排出源あり147地点だけの元予測の条件付き誤差改善はP+14.4%、N+18.8%だが、これは事後的な下位群診断である。
- **判定と重要な未解決点**：上流の排出負荷は底生無脊椎動物EQRに追加の**予測情報**を持つ可能性がある。しかし源あり147地点内でPとN排出負荷の順位相関は**ρ=0.882**と高く、作用する栄養塩は分離できない。生物評価法の国差、ポーランド等でのEQR欠測／照合脱落、農業・都市化・他放流源・受水河川での実際の移流／希釈は未制御。**下水由来の因果効果や種多様性の増減と同義にしない。**

[次の独立検証・Issue #13](https://github.com/zuizui0223/gesui/issues/13)は国別生物評価法、モニタリング／負荷の欠測、河川栄養源の交絡を明示して追跡する。元のchochoの港・島・蝶の資源地理とgesuiの肯定・否定の各結果は、対象・推定対象を混同せず並行して保持する。

**追加のソース機構診断（2026-11-08）**：[再現テスト成功](https://github.com/zuizui0223/gesui/actions/runs/42437488249)。底生無脊椎動物2,475地点の同一地域分割予測で、単なる「上流に対象処理場あり」のダミー変数だけならMSE改善**0.63%**、施設数だけなら**0.42%**。実測報告P排出質量を使えば**1.10%**、N質量なら**1.40%**改善。処理場存在ダミーに加えて報告排出質量を入れてもさらにPで**0.34%**、Nで**0.73%**改善する。これは**結果確認後に行った追加診断**であり、報告された排出の「量」が施設存在のみより予測に有用である可能性を示すが、P/Nの独立した化学機構や因果的な生物影響は実証しない。付着藻類では同じ診断の改善がごく小さい。

## 2026-11-08：反証検証・EQR欠測・農業負荷の観測可能性

[**Draft PR #16：水生生物の下水排出負荷予測に対する反証と品質監査**](https://github.com/zuizui0223/gesui/pull/16)を別解析ブランチで完了し、[公式原本・実行監査・データを永続保存](https://github.com/zuizui0223/gesui/tree/analysis/bqe2022-missingness-wfd-pressure-20261108/results/eea2022_bqe_falsification_20261108)した。元のPR #14の底生無脊椎動物EQR地域外予測MSE改善（上流の実測報告P負荷+1.10%、N負荷+1.40%）は削除せず保持する。ただし**一般化を制限する欠測と交絡**を明確にした。

- **欠測の原因を公式EQR原本で特定**：[検証成功](https://github.com/zuizui0223/gesui/actions/runs/42438230754)。底生無脊椎動物の河川接続候補3,282地点から807地点（すべてポーランド）が元の数値EQRを欠き、源あり候補394地点中247地点（62.7%）が脱落。付着藻類では3,663候補から879（ポーランド802・オーストリア77）が欠測、源あり403地点中293地点（72.7%）が脱落。元のモデルは正しく除外したが、pandas nullable NAによる欠測集計が誤って0件になっていたため、**品質集計のみ修正**し、推定対象と予測MSEは変えない。
- **農業由来の独立負荷補正は地域カバー不足**：[EEA Waterbase Emissions 2026原本の集計単位・報告方法を検証](https://github.com/zuizui0223/gesui/actions/runs/42438695309)し、[農業面源NP1を固定生物地点へ結合](https://github.com/zuizui0223/gesui/actions/runs/42439118473)。2022年の農業由来P負荷が報告された底生地点は142（リトアニア・オランダの2か国）で、下水源あり20地点のみ。付着藻類では184（リトアニアのみ）、源あり20。元の農業P値は**河川流域地区RBD／下位地区単位のモデル推定値**であり、各地点の実測P流入ではない。未報告地点を農業負荷0にはしない。
- **排出負荷の「量」の探索的追加検証**：[元の地点・源あり／なし・国・地理的検証分割を固定し、源あり地点のNまたはP質量を199回入れ替えた実行](https://github.com/zuizui0223/gesui/actions/runs/42439198765)。底生無脊椎動物では国別および国×流量層別の両置換条件で、元の報告配置より予測誤差が小さい置換は**P・Nとも0/199**。有限MC指標0.005は**事後・交換可能性未確認の記述値**であり、確証的なp値でも下水由来の因果効果でもない。付着藻類の同じ層別検証はP42/199、N2/199で、元の予測改善が微小なことも保持。

次の独立検証は、栄養塩の実際の受水河川・処理方式・農業／都市化の共変量を同じ流域・年にそろえ、別の未開封年・地域で同じモデルを再現すること。[継続課題 Issue #15](https://github.com/zuizui0223/gesui/issues/15)。PR #16はDraftのままで、`main`や`chocho`に未査読解析をマージしない。
