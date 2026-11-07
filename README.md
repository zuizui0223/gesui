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
