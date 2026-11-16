# 2026-11-16 次のsource-only反証：別の水質観測局を同じH7で使えるか

**元の同一EEA Biology monitoringSiteIdentifier+Scheme+2018/2022 yearにおける真のTN+TP source availabilityが、公式Part1+Part2合算でも上流報告負荷と9 H7×period/8 H7/8 H6/3国でSTOPした後の追加の地理的source検討。元の厳格同一局ゲートはSTOPのまま維持する。**

公式EEA Waterbase ICM 2026には公開`WISE6_SpatialObjects_DerivedData-csv.zip`（WebDAV original name、約2.8MB）があり、元の年別・試料別水質測定の監視局座標を別の局に対して得られる可能性がある。

## Source-only段階と厳格な境界
1. **第1段階（本作業）**：真正公式公開ZIPの元バイト数・SHAを検証し、内部CSV名・全ヘッダ・レコード数と座標／監視局ID/識別子schemeの存在・非欠測件数だけを監査。**観測局の実値や水質N/P数値はまだ評価しない**。この段階で生態係数をfitしない。
2. 第1段階の結果、元同一監視局を外れる空間結合が可能でも、これは**完全に別の推定対象**である。元の同一EEA Biology site×yearの水質測定を示したことにはならない。新しい水質局が元H7の区域内に入るという事実だけでは、下水上流源の流下経路・同一河川の観測値・河畔植物の生息地点への投下量は特定できない。
3. 第2段階へ進むなら、元H7流域の公式HydroBASINS原本多角形、WISE6実際のmonitoringSiteIdentifier+Scheme＋WQ2026 Part1/Part2の年と数値を**生物応答より先に**結合し、少なくとも同じ真H7内にあり、受水側の同じ河道ネットワーク上／地理的な上流下流整合ができるかsource onlyで検証。H7外の近傍局を直線距離だけで同じ流域に移すことは禁止。
4. 同じH7でも異なる水質局によるTN/TPは**流域内環境濃度proxy**であり、元生物局と同局年に測定した値とは異なる。上流年別WWTP報告源なしはnon-exposedではない。真のTP/TN条件付き比較の閾値や国・H6の独立反復を目的変数で変更しない。
5. アクセス成否と取れる監視局座標形式だけを結果として保存し、元の下水生態系のメタ解析に必要な真の処理水N/Pの寄与割合・植物葉N・幼虫寄主利用を捏造しない。

**状態：** 主研究の正負を含むPR #106はmainにマージ済み。追加source検証はDraft PR #107に記録、main/chochoを改変しない。
