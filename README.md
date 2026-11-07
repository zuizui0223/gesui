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
- **次の因果識別の入口**： [実測TN/TPと放流口を結ぶIssue #5](https://github.com/zuizui0223/gesui/issues/5)。2026年公開の[全球河川栄養塩測定データ](https://doi.org/10.1038/s41597-026-07625-1)と[EEA Waterbase 2026](https://www.eea.europa.eu/en/datahub/datahubitem-view/fbf3717c-cd7b-4785-933a-d0cf510542e1)を地理・水文・年月で照合する計画。**照合と栄養塩の生態効果解析は未実行**。

**重要**：この新しい結果はDraftブランチに存在し、mainの研究用データを無断で合併していない。河道上の地点の種記録差は、直接の下水N/P暴露、植物個体数、定着・分散、因果効果を意味しない。元chochoの肯定結果と矛盾する分析ではない。
