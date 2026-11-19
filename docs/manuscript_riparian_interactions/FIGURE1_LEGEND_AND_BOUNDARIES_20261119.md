# 図1 A–C：投稿図版用キャプションと原図の区別（2026年11月19日）

**主図は1枚（A、B、Cの3パネル）**。本文 [MANUSCRIPT_FULL_20261119.md](MANUSCRIPT_FULL_20261119.md) に対応。[元のMermaid構成](README.md) を保持したまま、数値・観測単位・実測／モデル／仮説の区別が消えないキャプションを固定する。

## Figure 1. Within-system nutrient-related associations and the unresolved aquatic–terrestrial pathway

**A | Distinct within-system associations measured at different observational supports.** In Austria (2022), reported upstream wastewater-treatment-plant phosphorus load normalized by discharge provided additional predictive information for measured total phosphorus and benthic macroinvertebrate ecological quality ratio (EQR). In France (72 distinct H7 subcatchments, 139 H7-by-period observations), observed stream total phosphorus was negatively associated with benthic EQR (standardized β = −0.434, 95% interval −0.654 to −0.279), while the nutrient-affinity score of recorded terrestrial plant assemblages was positively associated with the nutrient affinity of known larval host plants of recorded lycaenid butterfly species (β = +0.634, interval +0.116 to +0.880). The stream phosphorus–plant affinity association was uncertain (β = +0.027, interval −0.170 to +0.229). Solid arrows represent estimated statistical associations or predictive contributions, not causal fluxes. Austrian and French estimates arise from separate samples and have different scales; their coefficients are not multiplied to infer mediation.

**B | The limit of a same-site, matched-month boundary test.** A planned comparison of change in measured stream total phosphorus between 2018 and 2022 with concurrent changes in the near-bank versus distant riparian vegetation index did not pass its pre-specified single-scene satellite coverage threshold (116 of 178 eligible monitoring sites; minimum 120). A two-scene post-hoc 122-site analysis is reported only in Supplementary material, not presented as a confirmatory null result. No site in that post-hoc sample had positive 2018 reported upstream WWTP phosphorus load; the analysis therefore could not assess the absence of a wastewater-specific effect on terrestrial vegetation. A dotted arrow denotes an attempted but unsupported test, not proven ecological disconnection.

**C | Three hypothesized, unmeasured routes from wastewater treatment to butterflies.** Route A: treated effluent enters a river and may reach riparian soils and host plants through inundation, deposition or infiltration, potentially entering butterfly larvae through their hosts. Route B: adult butterflies may contact or ingest water, salts or organic substances from damp substrates; observed puddling does not itself demonstrate nitrogen assimilation. Route C: sewage sludge is processed into a fertilizer or soil amendment and transported to documented application fields, where sludge-derived nitrogen may be incorporated into host plants and larvae **without passing through the river**. The source treatment works, receiving river and application field need not be in the same catchment. Aquatic insect emergence to riparian spiders is a separate contextual route. All arrows in C are hypotheses, not measured fluxes in the present data. Natural-abundance δ¹⁵N alone cannot identify effluent versus sludge because their source signatures can overlap.

## 日本語キャプション（投稿前確認用）

**図1｜水域内と陸域内の栄養関連の生物指標の対応と、未解決の水陸境界。** **A** オーストリア2022年の処理場報告リン負荷／流量が河川実測TPおよび底生EQRの地域外予測に加えた情報と、フランス共通139 H7×期間（72 H7）における河川TP–底生EQR（標準化β−0.434）、植物群集NUTRI–シジミチョウ科の既知寄主HNS（β+0.634）の正負の関連を示す。河川TP–植物NUTRIは不確実（β+0.027）。オーストリアとフランスは別標本であり、係数や予測改善を同じ効果量として掛け合わせない。**B** 2018→2022年の同一河川地点・同一夏季観測月におけるTP変化と河岸近接−遠方NDVI変化の比較は、単一シーンの元事前品質基準120地点に対し116地点しか残らずSTOP。事後的な122地点の数値は補足にのみ記録し、この標本では2018年の上流WWTP報告P陽性が0であるため、下水特異的な水陸影響なしとは結論できない。**C** 下水処理から蝶への仮説的3経路：**A** 処理水→河川・河岸寄主→幼虫、**B** 水際の湿った基質→成虫の口吻接触・摂取、**C** 汚泥由来肥料→施用農地→寄主→幼虫を示す。処理場と汚泥施用農地の位置を同一流域と仮定しない。羽化昆虫→河岸捕食者は別の既知候補経路。実線は観測的な予測関連、点線・破線は品質不足の検定や未観測の機構であり、物質の実測フラックスではない。

## 図版のデザイン仕様（実測・推定・仮説を誤認させない）

Aの水域と陸域は同じ色で統一しない。数字の単位はフランスの標準化係数、オーストリアのMSE改善率として別に記載。TP→植物NUTRIは実線ではなく不確実な推定区間が0を跨ぐことを注記する。Wは静的なモデル下水関連地理であり、実測処理場の到達量ではない。Bは数値効果ゼロではなくSOURCE/QUALITY STOPの位置として表示する。Cは本文に実際の安定同位体フラックス測定があるかのような矢印幅・太さを使わない。

**配色の意味だけで結果を伝えない。** 凡例に「直接の観測関連」「モデル・予測」「未実施または仮説」を必ず文字で併記する。H7のGBIF陸域記録は河岸距離内に限定されていないため、図中で河岸植物の実測個体群と同義に描かない。

## 対応元データ
- [フランス139共通標本の元係数](evidence/fr_interface_paths.csv)
- [オーストリアの原本・地点重複QA](https://github.com/zuizui0223/gesui/pull/118)
- [フランスの衛星SOURCE STOP・事後解析](SUPPLEMENTARY_ANALYSES_AND_SOURCE_STOPS_20261119.md)
- [主張・元source・不確実性の監査](CLAIM_EVIDENCE_LIMITS_20261119.md)

## 新しい陸域汚泥経路Cの資料

[EU・日本の汚泥肥料利用、¹⁵N標識・自然同位体・既知の昆虫影響文献と野外設計](SLUDGE_TO_HOST_BUTTERFLY_ROUTE_C_LITERATURE_20261119.md)。既存の陸域DAYCENT有機肥料モデルを汚泥施用記録と読み替えない。
