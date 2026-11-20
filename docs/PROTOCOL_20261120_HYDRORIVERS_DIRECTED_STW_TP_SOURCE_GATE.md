# 2026-11-20 実際の下水処理水TP採水点とChemPopの有向HydroRIVERS河道を照合するsource-only仕様

**SOURCE FIRSTの事後探索。処理施設の変更や下水由来Pの因果的な群集選別を実証する検定ではない。**
対象はEA公式Water Quality Explorer現行API 2026のWater Company「SA＝最終処理水」採水地点で、**2017/18当年に実際に0348総リン（TP as P、mg/L）が測定された地点**。元SOURCE起点の原資料は2017年67・2018年71 native SA source points（共通67）と、ChemPop v2の2017年389・2018年366の形質得点可BIO site。前処理のsource frameは2017または2018の生物調査源適格463siteから2km直線以内にある現在登録SA採水点367件であり、**全国の全処理場には当たらない**。この選択法は実際のBMWP指標を見ずに事前固定した。

## 河川グラフの公式原ソース

[HydroSHEDS HydroRIVERS v1 European official ZIP](https://data.hydrosheds.org/file/HydroRIVERS/HydroRIVERS_v10_eu_shp.zip)をpublisherから取得し、元ZIP SHA256 **500da7d36ceee0aa4c82dd625f5c53d80ba2f7ac14ad022cd3a6601bf55068e5**、938,544 original reachesと元CRS EPSG4326、正式フィールド HYRIV_ID/NEXT_DOWN/LENGTH_KM/MAIN_RIV/HYBAS_L12をsource-onlyで事前確認済み。[HydroRIVERS SOURCE #44679393760](https://github.com/zuizui0223/gesui/actions/runs/44679393760)。図郭はWGS84 lon −9〜4、lat49〜58（英国に沿う有向ルート最大50kmの余裕域）。元全生物・source資料を英国BNG(EPSG27700)で統一。元河道グラフはWGS84→BNGへ投影し、出現可能/不可能な小流路を推測補完しない。

HydroRIVERSは原対象河道の条件が集水面積≥10km²または平年流量≥0.1m³/sで、都市側溝/小河川/処理場の原排水管を含まない。**原登録SA採水点が本当の放流水が混合する位置とは限らない**。したがってsource採水点の河道スナップができても「真正施設アウトレット位置を再現した」とは決して記載しない。

## 元ソースと同年の河道SOURCE判定の凍結

1. **源試料**：原 WQE native pointNotation、TP=0348、unit mg/L、sampleMaterial FINAL SEWAGE EFFLUENT、samplingPointType SA、phenomenonTimeが2017または2018に一致する原試料だけ。水道会社以外UA、storm overflow SB、非水域試料・計測N/Pでない行は除外し、源未報告は0と扱わない。
2. **生物現地試料**：original EIDC ChemPop v2の７地域の実採集日2017/18とBWMP科値付与率70%以上、得点3科以上、元source環境欠測なしという原SOURCE品質ゲートのみ。**まだ生物BMWPの値は開かない**。
3. **河道にのる位置**：元source SA採水点、元BIO両方で、正式HydroRIVERS v1最近傍lineから**500 m以内**を要求し、複数異なるreachが最近接距離±1mの河道分岐で同等ならSOURCE AMBIGUOUSとして除外。小流路欠落のための非snapはSOURCE UNROUTABLEに残し、処理場の影響なしと判定しない。
4. **方向と流路**：SOURCEのHYRIV_IDから正式NEXT_DOWNを辿る。50km以内で元BIOのHYRIV_IDに達する場合だけ**潜在的「source側が上流にある河道候補」**と数える。逆向きBIO→SOURCEが辿れるものは「sourceが生物地点の下流にある対照」、同じ reach はsegment内の前後位置が不明として別の不確実群（TRUE順序不明）、異なるbranch/50km以上/グラフ外は「検証対象外」を区別。どれも元source生物個体数効果ではない。
5. 真正HYBAS_L12、MAIN_RIVは対象河道の独立性の目安として源別に数えるが、同一大河川内の複数HYBAS_L12を独立5実験と無条件に数えない。元施設のoutfallやsource→受水到達は座標/管路や放流水量の元台帳で確認されていないため因果モデルSOURCE STOPは別途維持。

## 事前固定のSTOP

原研究での目標ソース規模は**実際に有向接続した異なるBIO site≥50、BIO site×year≥100、独立水文block≥5**。今回のgraphは対象STWを最初から2km以内の地理点でsource selectionしているため、ここから外側の上流STWの寄与が落ちている。原SOURCE結果がこの下限に満たなくても**0から好結果が出るまでスナップ500m→2kmやバッファ50kmを後付けで緩めない**。真のアウトレット座標と実際のsource-year N/P質量を持たないため、閾値を通過しても「下水TPによる科BMWPの因果選別」主モデルはSTOPを保持する。

流下距離に不確実性があり、HydroRIVERS v1 NEXT_DOWNは大河川の流向の近似であって放流水の管路、河川横断水質フラックスや下水の栄養起源を分離しない。源のTP濃度はmg/L、元処理施設質量kg/yrもTNも未取得。全成果は原zip/ソースSHA・無名件数・サンプル生物source QAのみ公表。元処理場・EA地点のnative codes、座標、化学原数値、元BWMP科scoreと生物原site行を公開Gitに転載しない。

結果がSOURCE STOPでも科学的には「確認に必要な有向流域の支持条件が公開データでは不足」として報告し、元の英国BMWPにおけるEDF⇢栄養WQ予測性という確実な観測結果を捨てない。
