# 2021年の実測下水放流水質まで接続：米国128湿地・160許可・344登録放流口の全国マクロ解析

**2026-11-11／Draft／実測排水データのSOURCE-ONLY監査完了。植物多様性への下水因果効果は未検証。**

## 研究の主線

世界1,055港・42か国ではモデル下水発生量と植物・動物・菌類のGBIF記録数に再現的な正の相関がある。ただし植物記録数を50件に希薄化すると植物種数の正の効果は小さく、不確実性区間がゼロをまたぐ。米国の標準化NWCA2021湿地植生では、登録下水放流口から10km以内（103湿地）と20km以上（137湿地）を比べると、全植物Shannon多様性 H_ALL の距離係数が**地点等重みでは−0.1204 [−0.2506,+0.0097]、元のNWCA調査面積重みでは+0.1806 [+0.0104,+0.3509]**に逆転する。重み上位5%を上限制限すると加重CIもゼロをまたぐ。元の調査重みのKish有効標本数は240実地点に対し約61地点相当である。つまり正の関連があっても強い重み依存性を持ち、実測排水影響の証拠ではない。[旧Draft PR #49](https://github.com/zuizui0223/gesui/pull/49)。

今回の目的は**2021年にN/P排水実測が行われた登録放流口そのもの**が、この全国湿地標本の近くにどれだけ存在するのか、植物の値を見ずに数えることである。水文的に下流の湿地へ放流水が届くかはさらに独立検証を要する。

## 1. 植物測定を見ずに対象放流源を固定

- 2021年NWCAの初回確率標本・植物と水質の両方が元IDで対応する**568湿地**。2026年11月8日EPA ECHO全国ICIS-NPDES登録放流口（原本SHA-256 `f0d1310b1c07085c42527e2ede5cf97ab83d517d50548b7dadb7621a44aa6160`）。下水関連施設の**NAICS 221320 または SIC 4952**、LATLONG_TYPE=Permitted Feature、2021年許可日付としてあり得るものだけを対象とした。
- 元の568湿地で登録放流口から**10km以内は128地点（32州）**。元の2021植生比較の同州×生態地域内のNear103とFar137は変更しない。近接側の登録**160行政許可**をSHAで事前固定し、許可単位ではなく**実際の登録放流口番号**まで追加で照合すると**344異なる放流口**、湿地×放流口の候補組み合わせは**392**であった。344放流口は2026年登録座標であり、現地の2021放流口を現場確認したものではない。
- 160許可・128湿地・344放流口は、2021実測N/Pの有無を見る前にファイルSHAとともにGitHub上の独立SOURCEファイルに凍結した。

## 2. EPAの2021年の実測排水報告を取得

米国EPA公式の**ICIS-NPDES Discharge Monitoring Reports (DMR)** 全国年度原本、FY2021およびFY2022を使用。年度ファイルは会計年度のため、2021年**暦年**の測定を揃えるには二つの原本が必要。公式の`EXTERNAL_PERMIT_NMBR`と`PERM_FEATURE_NMBR`を保持し、事前固定160許可の2021年`MONITORING_PERIOD_END_DATE`に限定。実測対象は**総窒素 as N = 00600、総リン as P = 00665、流量 = 50050**で、最終外部放流（PERM_FEATURE_TYPE_CODE=EXO、測定位置=1/2/EG）かつ数値・非NODI・非検出限界表記だけを残した。

EPA公式はDMR測定値を濃度・数量など**異なる単位と統計基準**で報告する。今回のSOURCE版では原数値の勝手な合計、N/P質量の推定、未測定をゼロとする処理は行っていない。

| 実測報告の存在（2021暦年、全国160候補許可） | 確認できた許可数 |
|---|---:|
| **総窒素TNの数値** | **34** |
| **総リンTPの数値** | **79** |
| **TNとTPの両方** | **31** |
| **放流水量の数値** | **117** |

許可だけで集計すると、128近接湿地の**76地点**に、2021年の総Nまたは総Pが報告された許可が近くにあった。しかし同一許可の**別の放流口**で測定された可能性があるため、これをそのまま湿地の栄養塩曝露と見なしてはならない。

## 3. 厳密な放流口番号一致で、さらに絞り込み

登録放流口の**許可番号とFEATURE番号が完全一致**する2021年DMRだけを照合。元の344放流口のうち、数値TNあるいはTPが実際に報告されていたのは**86放流口**であった。**そのうち76放流口では同じ放流口番号の水量記録も確認**できた。ただしN/Pは別月または別統計量でもあり得るため、濃度×流量の無条件の負荷変換はできない。

| 2021年の厳密な同一登録放流口での報告条件 | 登録放流口数 | 128元湿地中の近接候補数 |
|---|---:|---:|
| 実測総窒素 TN | **34** | 31 |
| 実測総リン TP | **83** | 66 |
| **TNまたはTP** | **86** | **71** |
| TNとTPの両方 | **31** | 26 |
| 2021年の水量記録 | **124** | 97 |

この**71湿地は26州、81行政許可**に分布する。元の240地点・27州の標準植生Near103のうち、登録放流口番号が2021年の実測総Nまたは総P報告と一致した地点は**51地点**だった。遠隔137地点は元の距離対照であり、未登録下水源や流域外からの流入がないことまで保証するものではない。2021年の測定報告は1〜4四半期に及ぶ。

事前に固定した**実測報告のSOURCE可用性条件（20近接湿地以上、10行政許可以上、6州以上、2四半期以上、実測N/Pを持つ放流口の半分以上に水量記録）を満たした**。これによって**次のNHDPlus水文流路確認へ進むための規模はある**が、現時点で植物の下水反応を推定する許可ではない。

## 4. 科学的な現在地

**確実に成立した肯定的事実：** 米国の標準化された2021年植物調査を持つ湿地の周辺に、同年に総N/総Pを実際に測定報告した下水放流口が、26州規模で存在する。これは『下水処理場の位置が近い』だけよりも一段強い曝露**候補**である。

**いまだ未確認の本丸：** 登録2026座標＝2021実際の流入口とは限らない。湿地植生が存在した場所への**河川流下経路・到達・希釈・塩水・地下水**が確認できていない。現地のTN/TP観測値が上流の登録WWTP由来と同定されたわけでもない。また原DMRは一部の窒素・リンが未報告の許可もあり、欠測は排水ゼロではない。

この条件で今すぐ2021年の植物Shannon H_ALLを『実測下水栄養塩負荷』に回帰するのは誤り。**次の独立したSOURCEステップはUSGS NHDPlusV2の有向流路によって、N/P測定のある86登録放流口から71湿地の流路を区別すること**。同じHUC12内・直線距離10km以内のだけでは実際の下流曝露を立証できない。2021DMRの数値を負荷に換算するには、報告期間・統計量と濃度／数量単位の整合を別途確認する必要がある。

## GitHub原本・独立監査

- [元160許可のSOURCE凍結](https://github.com/zuizui0223/gesui/blob/analysis/nwca2021-POTW-DMR-2021-source-20261111/results/us_NWCA2021_plant_TNTP_20261110/EPA_NWCA2021_160_SEWAGE_PERMITS_BEFORE_DMR2021_SOURCE.json) — [原本SOURCE固定SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42960834382)。
- [元344放流口の許可×FEATURE番号と568湿地の距離SOURCE凍結](https://github.com/zuizui0223/gesui/blob/analysis/nwca2021-POTW-DMR-2021-source-20261111/results/us_NWCA2021_plant_TNTP_20261110/EPA_NWCA2021_128WETLAND_160POTW_EXACT_NEAR_OUTFALL_FEATURE_SOURCE.json) — [公式ECHO元ファイルSHA照合SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42961350082)。
- [公式FY2021・FY2022 DMR全国原本を160許可へ絞った2021暦年レポート](https://github.com/zuizui0223/gesui/blob/analysis/nwca2021-POTW-DMR-2021-source-20261111/results/us_NWCA2021_plant_TNTP_20261110/EPA_NWCA2021_ORIGINAL_160POTW_2021_DMR_TNTP_SOURCE_FEASIBILITY.json) — [FY成功成果物の本来2021ソース統合SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42961956869)。
- [元344放流口×実際の2021 N/P/水量SOURCE一致、71湿地・26州](https://github.com/zuizui0223/gesui/blob/analysis/nwca2021-POTW-DMR-2021-source-20261111/results/us_NWCA2021_plant_TNTP_20261110/EPA_NWCA2021_EXACT_344POTW_OUTFALL_2021_TNTP_DMR_NHDPLUS_SOURCE_GATE.json) — [厳格な2021 N/P放流口照合SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42962030306)。
- [別実装で元の128/160/344/392と2021N/P/流量報告件数を再確認](https://github.com/zuizui0223/gesui/blob/analysis/nwca2021-POTW-DMR-2021-source-20261111/results/us_NWCA2021_plant_TNTP_20261110/EPA_NWCA2021_344OUTFALL_2021_DMR_INDEPENDENT_SOURCE_QA.json) — [独立SOURCE監査SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42962369321)。

出典：[米国EPA ECHO公式DMRデータ](https://echo.epa.gov/tools/data-downloads)、[DMR要素辞書](https://echo.epa.gov/tools/data-downloads/icis-npdes-dmr-summary)、[USGS NLDI有向流路API](https://api.water.usgs.gov/docs/nldi/navigation/)。

**研究はDraft。`gesui/main`、`chocho`、既存のPR #49の元の植生結果にはマージ・上書きを行わない。**