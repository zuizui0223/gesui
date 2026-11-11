# 全米下水×標準化湿地植物：2021年実測N/P放流口からのUSGS流下接続の原本判定

2026-11-11／Draft／**流下経路のsource-only判定、追加の植物の正負の回帰は行っていない**。既存の全球42か国の正のGBIF記録関連と、米国2021年NWCA標準化植物の重み依存の近接関連を変えない。

## 1. 独立した実測下水排水源の確保（既存PR #54から固定）

- 2021年の米国EPA NWCA確率標本568湿地（植物と現地水質の両方に原本SITE_IDが存在）、EPA2026年11月8日のNPDES全国放流口報告書のうち下水NAICS221320/SIC4952・Permitted Feature実座標・2021許可可能日付を満たす登録放流口。
- 2021年の実測DMRを暦年としてFY2021+FY2022から取得。正確に同じNPDES許可番号+放流口FEATURE番号で **実測総窒素00600または総リン00665が数字として報告されていた86登録放流口**と、その登録点から10km以内の **71湿地・26州**、100個の2021実測放流口×湿地ペアに限定。この71湿地は2021年排水の現地受水確認を意味しない。
- 元のNWCA同州×生態地域の240地点（Near103/Far137）では、Near103のうち51地点に同一放流口の実測DMRが対応。**植物H_ALL/Shannon、FQAI、NWCA現地TN/TPを新たに見ていないSOURCE集団**。

## 2. 新しいUSGS NLDI NHDPlusV2原本を本当に取得

公式USGS [NLDIのキャッチメントCOMID](https://api.water.usgs.gov/docs/nldi/feature-sources/)と[hydrolocation](https://api.water.usgs.gov/docs/nldi/)を、元の71湿地・86登録放流口の座標でそれぞれ照合。NLDIの`comid/position`は元の点がどの集水域にあるかを判定し、必ずしも川そのものに接することを示さない。別の`hydrolocation`の河道上の位置と元座標との距離を計算し、両者のCOMID一致、距離**200m以内**を厳格な河川近接候補とする（500mは事前に定めた感度表示のみ）。

| 厳密なSOURCE対応段階 | NWCA湿地 | POTW登録放流口 | 備考 |
|---|---:|---:|---|
| 元の2021TN/TP数値DMRと同じ登録放流口で10km以内 | **71** | **86** | 2026登録位置・2021実測報告、流路不明 |
| USGS NHDPlusのcatchment COMID取得 | **71** | **83** | 流域インデックス、河道近接とは異なる |
| USGS hydrolocationで200m以内＋元COMID一致 | **20** | **35** | API成功例での河道近接候補 |
| **両端が200m以内**＋別の河道COMID＋元の100登録ペアを維持 | **10** | **11** | 11候補ペア／11行政許可／8州。まだ下流か未判定 |

- 元の100湿地×実測放流口ペアのうち**11ペア**が双方の河道200m条件＋別COMIDを満たす。他に**1ペアは同じCOMID**で、どちらが上流かNLDI whole-reach表現からは決められない。元の200m strictと異なる**500m modelled location**の組み合わせは16ペアあったが、事前の厳格ゲートへの代入は不可。
- **API欠測が大きい**。NLDI hydrolocation **157地点中70地点でHTTP成功応答が得られず**（湿地28／放流口42）。対してcatchment positionは湿地71/71、放流口83/86で成功。こうした技術的未評価地点は、河道につながらない地点でも、TN/TP曝露がゼロの地点でもない。したがって**10湿地という値は現在のAPIと200m設定で証明できる候補数**であって、物理的な全接続湿地の上限ではない。

## 3. 前もって決めた流路・生物解析のSTOP判定

| NHDPlusV2下流流路の探索へ進むためのSOURCE要件 | 最低数 | 実際に確認できた数 |
|---|---:|---:|
| 別COMIDの河道200m条件が両端で成立する湿地 | **20** | **10** |
| 対応する実測TN/TP報告POTW放流口FEATURE | **15** | **11** |
| 異なる行政NPDES許可 | **10** | **11** |
| 州 | **8** | **8** |
| 原の比較で同州×生態地域の遠隔対照湿地 | **30** | **137** |

**最初の二要件が失敗したため、USGSの下流DM／分流DDの50ネットワークkm内流下検索は実施していない。** ソースJSON内では下流接続件数を**null（未判定）**とし、「ゼロの下流河川接続が確認された」という誤った表示を避けた。したがって今の科学的な結論は「下流接続が存在しない」ではなく、**『2021年測定POTWと2021年標準化湿地で、先に決めた厳格条件のもとでは下水曝露の生物解析に足る独立SOURCEを確定できなかった』**である。

USGS公式の[NLDI navigation仕様](https://api.water.usgs.gov/docs/nldi/navigation/)は河道全体の上流/下流の流路を返すが、湿地生態系への水の流入・地下水供給・潮汐の逆流・河口の塩水混合などは別問題。かりに同じ下流流路を検出しても、2021年WWTPの実排水が2021年湿地の植物の根に到達したとは言えない。

## 4. それでも残る肯定結果／次の研究の優先順位

- 元の全球1,055港・42か国の2015年モデル下水発生量×2011〜2024年GBIF全Plantae記録件数：β**+0.3156** [95%CI +0.1816,+0.4495]、未希薄化植物記録種数β**+0.2288** [+0.1094,+0.3483]は維持する。ただし50記録希薄化種数β+0.0126（CIゼロをまたぐ）。動物界・菌類の記録数も正であり、観測機会の集中は主要な代替説明。
- NWCA240湿地・27州で元の登録POTW10km以内対20km以上の植物Shannon差は、等重み**−0.120** [−0.251,+0.010]／NWCA面積重み**+0.181** [+0.010,+0.351]に符号反転し、重み上位5%制限で正の区間はゼロをまたいだ。これも**行政放流口の近接記録との観測関連**の範囲を出ない。
- この71湿地×86放流口では新しい下水に対する植物Shannon係数をフィットしない。次に本当に下水効果を研究するなら、広域の**source-verified水文流域標本**を最初から対象集団として選び、NWCA内の少数の河道近接点を救済するために半径・NHDソースの閾値を緩めない。原の2018/2022年欧州河川底生動物EQR×公式下水N/P放流量の複数年マクロ結果は、別の大規模証拠として優先するのが妥当。

## 5. 完全な原本・再現性

- [原の86POTW×71湿地の2021年TN/TP計測放流口SOURCE成功・独立監査（PR #54）](https://github.com/zuizui0223/gesui/pull/54)。
- [新しい公式USGS NLDI API（公開サンプル）位置＋下流flowlineレスポンス検証 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42963556338)。
- [71wetland/86outfallの公式USGS API・4分割原本と厳格200mSOURCE判定成功](https://github.com/zuizui0223/gesui/actions/runs/42964085386)、[結果JSON（下流未実施は null に修正）](https://github.com/zuizui0223/gesui/blob/analysis/nwca2021-POTW-2021-USGS-NLDI-direction-source-20261111/results/us_NWCA2021_plant_TNTP_20261110/EPA_NWCA2021_POTW_2021_NP71WET_86OUTFALL_USGS_NLDI_STRICT200M_DIRECTION_SOURCE_GATE.json)。
- [原のUSGS成功4シャードを再利用しDOWNSTREAM未実施をnullにするオフライン再集計 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42965476361)。
- [別スクリプトで71/86/157地点のソース・COMID・200m/500m基準・API HTTPステータス・NHDPlus適格性を独立監査 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42965721292)、[監査JSON](https://github.com/zuizui0223/gesui/blob/analysis/nwca2021-POTW-2021-USGS-NLDI-direction-source-20261111/results/us_NWCA2021_plant_TNTP_20261110/EPA_NWCA2021_NLDI_POTW71WET86_SOURCE_INDEPENDENT_QA.json)。

**全過程はSOURCE-ONLY。未判定サイトを汚染なし・排水なし・接続なしには置換しない。元のNWCA植物Shannon/H_ALL/FQAIや湿地現地N/Pの生物応答値を開かずに停止し、Draft PR #55以外のファイルは変更していない。**