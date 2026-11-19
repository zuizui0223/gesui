# 2026-11-19｜JRC-native流域の点源N×陸域N投入：元データ全件・ジオメトリ・源のみの地理対応

**結論の範囲：SOURCE-only 9,739 native-ID面積照合候補。30/30で元ジオメトリの相互一致、しかし9,739流域全ポリゴンの一意性・空間的独立性・生物応答の変化は未検証。**
従来2018河川×2017–19植物・蝶のNUTS2共通2×2は正式SOURCE-STOP（25地域・6国）であり、この新たな流域SOURCE結果によってそのSTOPを取り消さない。2016河川EQRの既存の点源追加予測改善も国等重み95%区間ゼロ横断のまま。

## 原データと正式な再現性

JRCの [Inputs to soil](https://water.jrc.ec.europa.eu/arcgis/rest/services/inmap/MIS_Inputs_to_soil/MapServer) の河川流域レイヤー5/10/15/20は元地理コード Basin_3035_HydroID（MEAN_mgsqmy または MEAN_kgha）。[Losses to water](https://water.jrc.ec.europa.eu/arcgis/rest/services/inmap/Losses_to_water_MIS/MapServer) layer5 は GREEN 2014–18の点源 PsN（BasinID, AreaSqKm、家庭＋産業を含む）。異なるフィールド名を**無検証でID一致とみなさない**ため、全ID/Areaの照合と元ポリゴン検証を別々に実行した。

| SOURCE レイヤー | 元features | 有効モデル値 | 欠測（真ゼロではない） |
|---|---:|---:|---:|
| EMEP 2015 大気N沈着 native HydroID | 74,016 | 36,483 | 37,533 |
| DAYCENT 2010–19 鉱質N施肥 HydroID | 12,338 | 9,823 | 2,515 |
| DAYCENT 2010–19 有機N施肥 HydroID | 12,338 | 9,823 | 2,515 |
| DAYCENT 2010–19 生物N固定 HydroID | 12,338 | 9,823 | 2,515 |
| GREEN 2014–18 水への点源N BasinID | 71,779 | 71,779 | 0 |

元HydroIDを別4土壌レイヤーで一意に照合したところ、主要3陸源（大気沈着、鉱質・有機施肥）の数値がそろうのは**9,740流域ID**。そのうち点源BasinIDの整数IDおよび元公式AreaSqKmが**全一致の候補**は**9,739流域ID**。これは元の9,740の流域のうち1つが不一致/未対応であるという意味だが、その1流域の原因分類は本段階では未確認。元公式属性サンプルのJRC HydroID6/BasinID6は共通面積143,655.77 km²で、別々に取得した元ポリゴンも相互の面積重複100%。

- [元5層のレコード/ID全件確認 #44439326535](https://github.com/zuizui0223/gesui/actions/runs/44439326535)
- [全元属性・欠測・HydroID/BasinID面積照合 #44440277477](https://github.com/zuizui0223/gesui/actions/runs/44440277477)／[元属性候補CSV＋原本API SHA artifact](https://github.com/zuizui0223/gesui/actions/runs/44440277477/artifacts/14491365519)

## 固定30IDの地理証明（元形状をそのまま別途取得）

全native候補9,739のIDを文字列ハッシュで固定順序に並べた先頭**30ID**をあらかじめ選び、JRC EMEP土壌HydroIDとGREEN点源BasinIDの元PolygonをEPSG3035の同一座標系で取得した。両方向のポリゴン面積重複率が**全30/30流域で1.000**、元areaとGIS面積の比も全30/30で1.000。先に固定した「少なくとも27/30が両方向95%以上」とする地理支持ゲートをPASSした。サンプル30流域の異なるID同士の相互面積重複が5%を超えるペアは0件。

- [単一ID6ジオメトリの独立pilot #44440914367](https://github.com/zuizui0223/gesui/actions/runs/44440914367)
- [30/30元Polygon厳格一致 #44441131398](https://github.com/zuizui0223/gesui/actions/runs/44441131398)／[原形状比較の元SHAと対象ID artifact](https://github.com/zuizui0223/gesui/actions/runs/44441131398/artifacts/14490309516)

30流域のgeometry一致は「元JRC HydroIDとBasinIDが本研究の共通流域キーである可能性」を強めるが、**全9,739流域の形状・親子構造の完全一致**を証明するものではない。異なるJRC流域IDも集水域の階層で依存し得る。したがって9,739件を独立サンプルとしてのp値・CIや効果量のメタ解析にしない。

## 広域source地理の記述的比較（生物値は一切未読）

上の全属性・30地理サンプルのSHAチェックが通った場合だけ、同じ9,739候補に対し、log1p(水への点源モデルN t/元流域km²)と各log1p(陸へのモデルN平均)のPearsonおよびSpearman相関を全源で算出した。

| 陸域N投入 | Pearson(log点源N密度,log陸N平均) | Spearman（順位） | 点源N>0の候補だけのPearson |
|---|---:|---:|---:|
| 大気沈着 | **+0.0890** | +0.2646 | +0.0631 |
| 鉱質肥料 | **−0.0035** | +0.0795 | −0.0592 |
| 有機肥料 | **−0.0168** | +0.1599 | −0.0678 |
| 生物N固定（副次） | **+0.0438** | +0.1058 | +0.0476 |

- 共通候補9,739：点源モデルが正は**7,761流域**、元モデル真ゼロは**1,978流域**。元流域面積分布5%/中央値/95%は**2.53/17.94/703.15 km²**。これらは独立な9739地域のサンプリング精度ではない。
- 大気沈着に対する点源N密度のPearsonは小さいが、Spearman+0.265で非線形や順位の対応が見える。さらに元流域面積四分位の最大面積群ではPearson **+0.191**、最小面積群では **+0.058**。各出所の観測支援と流域面積に依存する。各負荷源は2015/2010–19/2014–18に測定年と元モデルが異なる。
- この規模での単純相関の小ささは、源モデルの地域的非冗長性と整合的だが、点源と面源の物理的独立、流域植物・蝶の生態選択、河川内EQRの因果メカニズム、放流水の陸域到達を証明しない。国・気候・都市化・親流域による依存を調整していない。p値は計算しない。

[全4源9,739候補のsource-only地理検証 #44441591484](https://github.com/zuizui0223/gesui/actions/runs/44441591484)／[全数値JSON artifact](https://github.com/zuizui0223/gesui/actions/runs/44441591484/artifacts/14491521088)。

## 次の真のボトルネック

**元H7（2017–19のEIV品質適格226流域）との、厳格な原Polygon⇄JRC native polygon重複**が必要。元H7のGBIF genuine query WKTはsource-frozenされているため、JRC flow-networkと面積交差を**生物NUTRI数値を開く前**に確認できる。ただしJRC大流域はH7と1対1ではなく集水域の親子の重複や上流点源負荷の二重計上に注意。下流の葉N濃度や蝶の摂食はこの解析では測定していない。

このsource-onlyゲートが通るまでNUTS2時点の2×2を新地理単位で再実行せず、主論文の要旨・結論を変更しない。mainや草稿PR #124、および旧SOURCE PR #125/#126は未変更。新しいnative source解析は独立Draft PR #127に保存する。
