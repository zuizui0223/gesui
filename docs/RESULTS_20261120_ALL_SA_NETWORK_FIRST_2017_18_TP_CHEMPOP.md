# 2026-11-20｜全国の実測TP処理水源を河川網で先に選別：148異なるBIO地点・237地点×年

**SOURCE-only・成功：[GitHub Actions #44681039804](https://github.com/zuizui0223/gesui/actions/runs/44681039804)、集約JSON artifact ID 14604522034。** 英国の処理水源を対象とする後付け探索的なデータ支持監査であり、生物BMWPの新しい効果係数、処理場の年間N/P質量、源→受水の実際の化学的到達や生態系サービスを直接推定したものではない。

## 2km前抽出の選択上の問題

以前の[河川網解析 #44679903035](https://github.com/zuizui0223/gesui/actions/runs/44679903035)は、同年BIO地点から直線2km以内にある現行SA採水点367を**先に**抽出。そのうちTPを実測した2017年67、2018年71点でHydroRIVERS50km流向を調べ、上流reach候補BIOサイトは2017年36/2018年43、和集合51地点・79 site-yearだった。この実装では**直線2〜50km上流の放流水源を母集団から除外**し、流下50kmという曝露定義を正しく覆わない。

今回の追加探索は[元事後source protocol](PROTOCOL_20261120_ALL_SA_NETWORK_FIRST_HISTORICAL_TP_GATE.md)を使い、公式全国の5,722 SAを**先に**HydroRIVERS上でBIO地点との上流・逆方向placebo・同reach順序不明候補へ分け、その後初めて原2017/18の実際の処理水TP観測で選ぶ。BMWP形質値や好ましいモデル係数を見て源の数や年を決めていない。

## 原SOURCEと元のデータ品質

- [EA Water Quality Explorer v1](https://environment.data.gov.uk/water-quality/view/landing)、2026年現在登録のWater Company FINAL/TREATED EFFLUENT native SA採水地点：**5,722**。全23公式ページを取得、原BNG座標の有効率100%。現在登録されているだけで2017/18年の実測や稼働を意味しない。
- [HydroRIVERS v1 Europe ZIP](https://data.hydrosheds.org/file/HydroRIVERS/HydroRIVERS_v10_eu_shp.zip) 正式SHA256 **500da7d36ceee0aa4c82dd625f5c53d80ba2f7ac14ad022cd3a6601bf55068e5**、元欧州938,544河道、英国周辺bbox28,830 reach。原NEXT_DOWNで50km以内の有向関係を確認、BNGへ座標変換。公式SA5,722中500m以内に原河道へ一意スナップしたのは**2,992**、2,729に500m以内の河道がなく、1点は同距離で複数河道が競合。BIO2017/18の和集合463地点では442スナップ、21は未対応。後者は河川流下不明であり下水影響ゼロではない。
- ChemPop v2原2010–18 source score-eligible生物7,836原analysis／707site。2017年389、2018年366、年和集合463異なるBIO地点。元生物のBMWPスコア・調査結果を今回のデータ選択には使用しない。
- 登録SAから上流・逆方向・同reachいずれかの河川ネットワーク幾何関係にある源候補は**1,693点**。正式EA APIの観測上限は**250件/ページ**。[公式API移行説明](https://dsp.environment.data.gov.uk/support/announcements/275811447/1156939777)に合わせ、元native pointNotationをコンマ連結40点/回、43 batch、90 pagesで年別TP code0348を完全照合。原元観測8,810件のうちwater company SA・FINAL SEWAGE EFFLUENT・TP as P・mg/L・非負数値を満たす**8,753件**。54件は非数値など、3件は元material不一致で除外。
- 真の2017年TP実測があるSA採水地点 **281**、2018年 **293**、**275点は両年で実際のTP測定あり**。これらはsampling-point notationsであり、処理場本体の施設数や原permit/outfall数と同値とは確認できない。2017測定281点のうち41、2018測定293点のうち29は「別年のBIOとネットワーク候補であるためsourceデータ取得はしたが、**その年のBIOとの有向関係はない**」ものとして保管し、その年の曝露陽性から除外した。

## 同年TP実測source・50km有向河川網の地理スクリーン

| 元source品質指標 | 2017年 | 2018年 |
|---|---:|---:|
| 水道会社最終処理水TPが実際に測定された元SA地点 | **281** | **293** |
| 同年BIOに有向／逆向／同reachの関係がある源SA地点 | 240 | 264 |
| **NEXT_DOWN上流reach内50kmのTP実測SA候補を持つ異なるBIO地点** | **118** | **119** |
| その上流候補に対応する原SA採水地点 | 171 | 190 |
| 上流SA–BIOの有向reachペア数 | 282 | 314 |
| **逆方向：SAがBIOの下流にあるplacebo候補BIO地点** | **170** | **176** |
| 同reachで採水点とBIOの順序不明の候補BIO地点 | 61 | 71 |
| 少なくとも2SA上流幾何候補があるBIO地点 | 63 | 72 |

**元の上流候補の和集合は148異なるBIO地点・延べ237 site×year。** 元固定の50 distinct site・100 site×yearという数値source品質基準は通過。BIOの対応元HYBAS_L12ラベルは113、MAIN_RIVラベル33だが、必ずしも互いに独立な5個以上のhydrological blockを意味せず、旧「真の独立5流域」条件や生物学的因果試験の基準がすべて通過したわけではない。

逆方向placeboのBIO地点170/176は上流候補118/119より多い。方向を見ずに近い施設をsourceとして扱うモデルは偽陽性を生みやすい。反対に、逆方向であることだけで生物に別の上流源がないとは言えない。同reachの前後不明は上流へ勝手に割り当てない。

## 識別できなかった重要事項とSOURCE STOP

1. **採水点は物理放流口と必ずしも一致しない。** 元WQE SAは処理水の成分を測るsampling point、許可されたend-of-pipe outletの位置とは別の場合がある。[GOV.UK WWTW採水指針](https://www.gov.uk/government/publications/waste-water-treatment-works-treatment-monitoring-and-compliance-limits/waste-water-treatment-works-treatment-monitoring-and-compliance-limits)では複数outletや統合sampling pointを想定。したがって今回のNEXT_DOWNはSA sampling point reachの幾何学的方向を示すだけで、処理水がその河川経路を流れた実証ではない。
2. **TP mg/Lは下水年度P質量kg/yrでない。** 正式な元施設排水流量・原施設別N/P報告量・実際の受水化学の源別濃度到達が未取得。英国EEA2018の両実測NP正値FW適格施設には英国報告元なし、Wessex原CSVはActionsからHTTP403、Wessex元日別流量ファイル公開停止。0補完不可。
3. **HydroRIVERSの小流路欠落とsource snap誤差。** HydroRIVERS v1は河川規模の閾値により小川・側溝・施設内管渠を網羅しない。処理水原採水点を500mで元河道へスナップしているので、近くの別支流へ誤投影される可能性。正式outfall位置や同一reach内の順序を確認できないSOURCEは不確実として残す。
4. **下水による生物群集変化を今回のsource解析だけで主張しない。** 同じ英国ChemPopでは元モデルEDFと科BMWPの負関連、河川の3栄養化学種を調整するとEDF追加地理予測−0.519%という陰性、目内形質ラベルnull p=.206という未支持が存在する。SOURCE候補の増加により元結果を因果的に再解釈してはいけない。

## 次段階（方向ロバストネス）

元SOURCEの地理判定がスナップ上限100/250/500m、流向10/25/50kmの選択にどう依存するかを検定し、元固定50kmの主判定を後付けで変えない。上流reachと逆方向placebo・同reach順序不明の全年原件数を併記する。公式permit/outfall位置が許可された形で確認されない場合、実測処理場別N/P負荷による生物効果モデルは引き続きSOURCE STOP。必要なら結果を明確に事後探索の「候補処理水source監視地点との非因果的関係」と定義し、元river chemistry3種・別流域予測・placeboを守る。

本結果と旧2kmパイロット、Wessex403、EEA UK欠測、ChemPop BMWP/EDFの陰性を削除せず、源の原監査だけをPR #136のDraftに保存する。元Source測定TP行、原地理、source site IDを公Gitへ再配布しない。
