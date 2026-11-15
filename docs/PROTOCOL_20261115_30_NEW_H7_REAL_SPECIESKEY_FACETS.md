# Source-only spec: 30 genuine previously count-eligible external H7 × two periods × six taxon speciesKey facets

**作成：2026年11月15日。研究の植物NUTRI・蝶HNS回帰効果を一切見ずに元source-count監査 #43723994932 成功後の追加段階。元の30/40 H6 source STOPは新しい完全適格コホートができるまで不変。**

## 決定された地理・分類群

- 元EEA18/22元H7監視局6,330行、公式欧州HydroBASINS L6/L7、元454 H7×periodは既存の原本SHAで照合済。
- 一期のみEEA記録の267 H7等から、前段の本来のsingleton真正H6=56中12 H6に40 source-geography candidate H7がある。そのうち**30 H7／10 H6**は真正の旧GBIF原本1,219流域×2期間×3 taxa（Plantae≥50, Lycaenidae≥10, Animalia≥50）を両期間で満たす。SOURCE: [Actions #43723994932 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43723994932)。
- **30候補全件を対象**に5蝶科＋全Plantaeの元GBIF speciesKey facetを追加取得。結果を見て10 H6×1候補へ勝手に絞り込まない。元source rank（EEA dual-era優先、監視局数多い順、H7数字昇順）は前GISで固定されている。
- 追加H7は過去の元232 H7でないこと、真正H6 parentが既知元singleton H6であることをゲートで確認。決して同じ流域の2期を2つの独立地理反復にしない。

## 原本取得方法

- 追加地理: [Actions #43197894183](https://github.com/zuizui0223/gesui/actions/runs/43197894183) の元Level7、1,219河川流域WKT source SHA256 `47f597e0a798cd6155ee9fb69ce8916ed5629951d046d114feb9f667963310e4`。各source polygon part WKTの個別 SHAを照合する。港radiusやH3へ置き換えない。
- Butterfly family original native GBIF taxonKey は [Actions #43204064427](https://github.com/zuizui0223/gesui/actions/runs/43204064427) の元 `verified_butterfly_family_taxa` をSHA256 `afc524ebfa353e4d1c6ea51d208ba529a85fedb57c4d0f28d83545010bb76189` で保護し、Lycaenidae=5473等5科固定。
- 過去の元234 H7で実際に使用した `scripts/gesui_H7_original_fivefamily_plants_speciesKey_retention.py` の `task` 関数をそのまま再利用。GBIF `limit=0`, `facet=speciesKey`, `facetLimit=1000`, `facetOffset` pagination、`hasCoordinate=true`, `hasGeospatialIssue=false`, `occurrenceStatus=PRESENT`, `year=2017,2019` or `2020,2022`。元各WKT partsの単位記録から speciesKey別件数を保存、SHA付きoriginal source URL/レスポンスreceiptを保持。
- 6 groups = 全Plantae（taxonKey6）＋Hesperiidae, Lycaenidae, Nymphalidae, Papilionidae, Pieridae（family taxonomy original）。30 H7×2期×6 groups = **360 source facets**。3 deterministic H7 hash shards×2 concurrent workers/shard、最大2 shards並列。最初の2 H7のsource-only pilotで分類学的パラメータとAPI可否を確認してからfullを実行する。
- 元HTTP429またはAPI source failuresは未取得のままsource_status=`missing`。0記録とは混同しない。以前の成功レスポンスの再取得を強行したり、生物応答が望ましくないsource shardを再実行しない。

## 次の段階で評価する予定の *SOURCE-ONLY* 適格条件

- 各H7×両期間、Plantae源speciesKey facetsの元counts、真正種ID、記録数≥50は原本で検証済。
- 既存のEIV数値NUTRI 4,011 scored keysと受理名・accepted vascular exact matchを照合し、**元全Plantae種の≥40%がNUTRI対応し、10個以上の数値speciesKeys**（未同定新キーを成功扱いしない）。
- 蝶5科の真正 speciesKey源から元303種固定HNS数値あり186種と比較。元索引記録数≥10、HNS数値つき種≥2、5科総元記録の≥70%がHNS数値対応を仮のsource品質検査とする。研究後の係数や実際の蝶個体数は使わない。5科外や分類不能はsource missing/unknownとする。
- 元NUTRI対応辞書に**新規追加H7でのみ出現した種**がある場合、現時点で未採点（unresolved）にして下限品質を報告。これが足りなければexact GBIF accepted taxon再照合＋同じ元EIV辞書の原本SHA検証によるsource-only追加取得を独立ステージとする。現地葉Nを推定しない。
- これらを2期とも満たす新規H7が少なくとも10の**異なる元singleton H6**をカバーすることが、次の真のH6地域内モデルの最低条件。ただし単期EEA源の新H7を含むため、元のW/EEA/EQRと同一モデルと呼ばない。**元451 H7×期間との全source共変量再構成が必要**。
- 1 H6内に品質源を満たす複数候補がある場合は、前のsource-only監視局数・年・H7 IDで決まった順位最上位を採用。効果係数の大小で差し替えない。

## 解釈

この段階は**source-only taxa identities/data availability audit**であり、新しいNUTRI/HNS傾き、植物と蝶の因果媒介、排水窒素の伝達、生態系サービスは推定しない。GBIF indexは観測記録で実個体数・実際の植食関係ではない。旧PR #103の30/40 source STOP、6 EIV trait joint STOP、元219 H7時間的非同期、元水EQR×蝶残差不確実性は保存する。main/chochoを変更しない。
