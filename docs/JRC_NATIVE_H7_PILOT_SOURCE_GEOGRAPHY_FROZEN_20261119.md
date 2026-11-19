# 2026-11-19｜元GBIF H7植物流域とJRC native点源N/土壌N流域の独立SOURCE-only空間照合pilot

**事後探索のための新しいデータ取得設計。EIV NUTRI・蝶HNS・EEA EQRの数値は一切使用しない。過去のNUTS2 2×2 SOURCE-STOPを解除しない。**

JRC-nativeの原土壌4層、元点源（家庭＋産業）の全属性と公式流域Areaが一致する候補は9,739 native IDで、予めハッシュ選択した30の元Polygonでは土壌HydroID＝点源BasinIDの相互位置が100%一致した。ただし別の元H7（HYBAS_ID_L07）への直接ID対応や同じ受水河川の曝露は未検証。

## pilot H7の固定

元454 H7×2期間のうち2017–19を使用し、元EIV source-only品質フラグ451/454から適格な**226 H7**を固定。元2,808 GBIFリクエストfacetの2017–19 Plantaeの source_request_receipts に含まれる genuine WKT PolygonとSHAを使用し、EPSG3035で面積を計算。**EIV数値・蝶HNS・水域EQRを開かない**。

全226 H7のID文字列をsha256("gesui-20261119-JRC-H7-pilot:"+H7)で並べた**先頭12 H7**をpilot固定標本とする。サンプルの地理、元源負荷、結果の符号で選び直さない。JRC点源N river-basin layer #5 に対して、元H7のEPSG3035 bounding envelopeで公式ArcGIS Polygonの候補を取得し、元JRC候補9739 IDのうち実ジオメトリがH7と重なるものだけを使用。native JRC POINT BasinID/HydroIDや面積は照合済み元SOURCEの表を必須とする。

## pilotの source-only 合否（取得前に固定）

- JRC公式 ArcGIS native Polygon の元IDがAPIからID欠落/重複なく取得可能。1つのH7のbbox検索が3,000件を超える・ArcGISページ欠落・Polygon変換不能ならpilot SOURCE STOP。API 403は点源N=0ではなく取得失敗。
- 元JRC-native候補Polygonを元H7で切り、そのunion面積/H7面積が**80%以上**かつ、JRC候補同士の重複による「クリップ面積和/union面積」が**1.05以下**であるH7が、pilot12 H7中**8 H7以上**であれば次段の全226 H7 source地理取得を検討できる。元pilot有効H6も少なくとも4親流域。
- いずれか未達の場合にはこのpilotはSOURCE STOP。H7を任意の点源nativeポリゴンへ最近傍で割当てないし、地理重複のあるJRC流域を独立サンプルとして重複集計しない。
- このpilotの合格は**source geometry取得技術**についてのみであり、植物/蝶・底生の効果の肯定ではない。NUTS2実験の水陸source閾値を後から緩める形にはならず、元H7で新しい独立sourceゲートが必要。JRC点源モデルは河川に排出された（実測到達でない）N、陸投入は大気/農地源の別年モデルである。

確認結果・失敗行・元REST response SHAを全て保存し、source候補を使って相関の高いH7を選んで再モデルすることを禁止する。


## 2026-11-19：初回12H7で被覆100%が出た後の追加・厳格な解像度監査

先に固定していたsource-onlyの「元H7の80%以上がJRC SOURCEポリゴンunionで覆われ、重複1.05以下」の技術的な第1段は**12/12 PASS**（[元Action #44442417933](https://github.com/zuizui0223/gesui/actions/runs/44442417933)）。ただし元ポリゴンの中に大きなJRCネイティブ流域があり、H7がその内部に丸ごと含まれる場合もcoverage=100%となる。このため、元被覆ゲートだけでは小流域に局所的な独立曝露を割り当てられないことが**結果の確認後に明らかになった**。

そこで別の**post-pilot stricter diagnostic**を追加する。これは元の12H7を再選定せず、SOURCE-onlyであっても前向き事前登録ではない。JRCの元 native Basin AreaSqKm と各H7自身の genuine Polygon面積の比をすべての対応候補で算出し、H7を覆うnative Polygonのうち native面積がH7の**1×、2×、4×、10×以内**のものだけから計算した被覆率を全段表示。事後の「より細かくできる」主張は**4×以内でH7面積80%以上**を覆うpilot12H7のうち**8 H7以上**、さらに主被覆するJRC点源/土壌共通の優占ネイティブBasinIDが**8種類以上**のときに限る。

この追加閾値を満たさない場合は、元技術ゲートPASSは保持しつつ**H7で独立曝露を作れるかについてはSOURCE-STOPまたは要新しい地理資料**と判定する。モデル結果の符号を見てないため結果選択はないが、SOURCE pilot後に追加したことを公開記録する。JRC本流域の重複・親子関係、複数H7が同じnative IDに帰属する問題を無視しない。新しい植物NUTRI/蝶HNSの回帰は、独立単位と欠測が完備するまでしない。
