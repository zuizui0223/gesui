# Supplementary notes: exploratory interfaces, external studies and source-gate STOPs

日本語草稿の補足記録。主文 [MANUSCRIPT_FULL_20261119.md](MANUSCRIPT_FULL_20261119.md) と図1の説明 [FIGURE1_LEGEND_AND_BOUNDARIES_20261119.md](FIGURE1_LEGEND_AND_BOUNDARIES_20261119.md) に対応する。本補足中の数値を本文の独立した確認的生態効果として二重計上しない。

## S1. フランス同一河川地点の2018→2022年TP・NDVI変化

同じ河川監視地点で2018年・2022年の総リンTPと底生EQRを持つ263元候補から、河道位置・夏季同月の実測の条件により178地点（58 H7・23 H6）を先に選定。河道に沿った近距離20–100 mと遠距離500–1,000 mのSentinel-2 NDVI差を比較するため、各帯の有効自然被覆画素100以上・両年同月2か月以上を事前の品質基準とした。

元の単一衛星シーンでは116/178地点にしか満たず、事前の必要120地点より少ないため主検定はSTOP。事後的な同月最大2シーン補完により得た探索的122地点・44 H7・22 H6では、ΔTPとΔ近接−遠方NDVIのSpearman ρは−0.0276、親H6再標本化記述95%区間［−0.436,+0.316］。背景モデルにTP変化を追加した未知親H6予測改善は−3.73%［−9.19,+2.25］で悪化方向。同じ標本のΔ底生EQRとΔNDVIのSpearman ρは+0.006で、底生EQRの時点は必ずしも同月ではない。

**2018年上流WWTP報告リン陽性は元263地点のうち1地点、位置・季節を揃えた178地点で1地点、最終122地点で0地点。** 未報告は下水が実際に無いことと同義ではない。この事後比較を下水由来栄養塩の陸域移行に対する確証的「効果ゼロ」としない。

出所： [草稿PR #124の原source票・監査履歴](https://github.com/zuizui0223/gesui/pull/124)、[フランス122地点の元Actions #44020568960](https://github.com/zuizui0223/gesui/actions/runs/44020568960)、[元source負荷陽性台帳](https://github.com/zuizui0223/gesui/blob/docs/wastewater-riparian-interactions-outline/docs/manuscript_riparian_interactions/evidence/fr122_2018_upstream_WWTP_report_support_audit.csv)。

## S2. ドイツ6放流口の水質・周辺植物

2011–2022年の同月・上流下流TPと近傍植物群集の栄養立地指数を対比できる6放流口を探索的に集計。植物側に十分な出現記録がある4対では全て下流側の植物NUTRIが低く、そのうち3対で下流TPも低かった。TPが下流で高い1対でも植物NUTRIは低かった。河川監視地点は放流口から約7–25 km離れるため、同じ放流水を受けた植物個体を直接追っていない。今回の主張「下水が植物の栄養好みを高める」には合わず、上下流の観測的対比として補足のみに保存する。

出所：[ドイツ6地点の元source比較](https://github.com/zuizui0223/gesui/blob/docs/wastewater-riparian-interactions-outline/docs/manuscript_riparian_interactions/evidence/DE6_interface_contrasts.csv)、[PR #7](https://github.com/zuizui0223/gesui/pull/7)。

## S3. SOURCE-only点源 vs 陸へのN投入：行政区269→共通25地域でSTOP

JRC GREENの2014–18年平均水へのモデル点源N（家庭＋産業、トン/地域）と、陸へのEMEP 2015大気N沈着（mg/m²/年）・DAYCENT2010–19鉱質/有機N施肥・N固定（kg/ha/年）をGISCO2016 NUTS2で照合した元SOURCE-only6変数完全ケースは269地域・28か国。点源総量の対数と大気N、鉱質、肥料有機、固定の対数値の相関は順に+0.051/+0.069/+0.011/+0.166。国・地域面積を考慮した国内国等重み相関は+0.084/+0.025/−0.268/−0.021。有機Nの記述95%区間［−0.471,−0.029］、点源と陸3源の地理条件付きpartial R²は約0.111。低い単純相関を栄養源の地理的・因果的独立性とは読めない。

同じ2018河川底生・2017–19植物出現・全点面源を正式GISCO2016 1:100万境界で合わせると25 NUTS2地域・6か国・43 H6。事前の40地域・8か国などのSOURCE下限を満たさず**生物の水×陸2×2モデルは未実行**。2016底生動物だけを対象にした別の事後探索では、元EQR有効67 NUTS2・11か国で点源Nの背景追加予測改善+4.86%、陸源3変数の追加は−30.72%。点源Nの固有改善+4.06%の国ブート95%区間は［−3.21,+7.66］でゼロを含み、強い点源優位の基準は不合格。新しい独立実験や2018/2022 WWTP実測の再現として数えない。

出所：[行政区と生物共通標本のSOURCE-STOP PR #125](https://github.com/zuizui0223/gesui/pull/125)、[国内・面積・水域単独事後探索のPR #126](https://github.com/zuizui0223/gesui/pull/126)、[水域2016の元Actions #44436814968](https://github.com/zuizui0223/gesui/actions/runs/44436814968)。

## S4. JRC格子由来の元流域負荷モデル：全226 H7でも独立源不足

元EIV品質適格の2017–19植物H7は226。JRCの0.1° EMEP大気N・1 km DAYCENT肥料Nから作った元native river basin averageとGREEN水点源Nモデルの原ポリゴンを、元GBIFで使った真のH7 WKTへ重ねた。元226 H7×5モデル源＝1,130セルのsource取得は全件完走した。

複数元sourceポリゴンの和集合がH7面積90%以上を覆うH7は221。**それぞれの主要sourceについて単一元JRCポリゴンがH7面積90%以上を覆う**のは204 H7。元EsriJSONリングの大きな修復を別途厳格に考慮した保守的適格は19 H7、独立土壌側native source流域12。**全ての幾何修復を正しいと認める最も楽観的な境界でも、204 H7に含まれる独立土壌側元モデル流域は35であり、先に定めた最低100流域に届かない**。したがって新しい陸域生物モデルはSOURCE-STOPのまま。モデル元ラスターのH7直接ピクセル集計を実施したという主張はしない。

出所：[全226 H7 official source-only実行 #44449479977](https://github.com/zuizui0223/gesui/actions/runs/44449479977)、[全修復を許した独立sourceの上限 #44451436297](https://github.com/zuizui0223/gesui/actions/runs/44451436297)、[PR #128](https://github.com/zuizui0223/gesui/pull/128)。

## S5. 先行研究の公開原本再計算の扱い

カナダGrand Riverの8処理場における羽化昆虫と河岸クモのマイクロプラスチック、Bow Riverの安定同位体・微生物叢、および米国産業排水・PFAS研究の元表は、原著者による**公表済み研究の公開データを再集計した二次資料**である。時期・物質・追跡対象が異なり、本研究の独立な欧州栄養塩–植物–蝶の反復流域として合算してはいけない。これらの経路の存在は原著の結果として序論で引用し、独立再計算の監査元は[草稿PR #124の履歴](https://github.com/zuizui0223/gesui/pull/124)に残す。

参照する原著：[Environmental Pollution 2025](https://doi.org/10.1016/j.envpol.2025.125948)、[Environmental Pollution 2026](https://doi.org/10.1016/j.envpol.2026.129116)、[Global Change Biology 2023](https://doi.org/10.1111/gcb.16934)、[Ecology Letters 2024](https://doi.org/10.1111/ele.14401)。先行する物質フラックスの現地証拠を「未発見の現象」と書かない。

## S6. 解析の透明性

元の観測コホート、固定閾値、失敗ジョブ、元資料SHA、欠測原因・国別支持が検証可能なActions artifactおよびPRに残されている。本文の各数値の標本単位・観測形質・検定の結論については [主張と原本根拠の監査表](CLAIM_EVIDENCE_LIMITS_20261119.md) を優先する。SOURCE-STOPは効果ゼロではなく、負荷源と生物結果を同じ独立地理単位で比較できなかったことを意味する。
