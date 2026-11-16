# 2026-11-16 元EEA Biology観測局と公式2026 ICM水質観測局の空間候補監査（事後source-only）

**起点：検証済みmainの肯定/否定はPR #106で統合済み。PR #107内の公式2026年水質原本Part1+Part2で、生物調査局同一ID/2018・2022年の有効TNTPが101 station-year、うち上流WWTP報告正は12 station-year・9 H7×期間・8真正H6であり、固定30 H7×期間の曝露source gate STOPを確認した後の探索である。前段の同一観測局ゲートを緩めない。**

## 独立検証対象と前段の固定条件

- 生物監視局：元EEA2018/2022＋EPA2018/19地理原本6,330行のうち、元454 H7×期間と一致する真正EEA station-year **3,525**、真正distinct物理site **2,485**。元 source SHA256 `a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89`。source下流河道座標が得られない場合は生物サイトの近傍局検索を実行しない。
- 水質の位置原本：EEA公式2026年 Waterbase ICM `WISE6_SpatialObjects_DerivedData-csv.zip`、原本ZIP **2,797,036 bytes / SHA256 aa505be6f85d06c675e0d0fcb148a9148f80162f8e9590740a887cc95a3e349d**。内部 `Waterbase_v2025_1_S_WISE6_SpatialObject_DerivedData.csv` **90,859行**。位置: lat/lon、識別子: countryCode、monitoringSiteIdentifier、monitoringSiteIdentifierScheme。非欠測座標71,592件。原zipはGitHubに再配布せず、derived source-only counts。
- 同一観測局のTN/TPを取得済み：公式2026 Part1 original SHA256 `ff779521749b40b952cce35a4671240c3b6c47d87f32a8fef9240a4ac503a082`、Part2 `315396800d76e46db415164851fb4f99765b7ec2654828fc42782d902e679a7a`。その真正年別測定源と報告源のjoinは**別source gate**。この段階で新しい水質数値や蝶・植物・EQR値を開かない。
- 監視局の同一キーは `(countryCode, monitoringSiteIdentifier, monitoringSiteIdentifierScheme)` を必須とし、文字列完全一致の前に表記正規化・未知コード補完は禁止。同じ水質WQ ID＋schemeで複数の異なる座標がある場合はambiguousとして除外、平均座標を作らない。複数recordsが同一座標ならde-duplicateのみ。
- 元の生物局の座標が取得できる場合に限り、別WQ局までの**測地線直線距離**の閾値を結果前に固定：≤250m、≤1km、≤5km、実際の候補ID数と元生物station年/H7/H6との重複数を全件記録。≤1kmのケースを次段のsource-only候補とするが、**これだけでは同一河道・同じH7流域・放流水到達を意味しない。**
- 次段で公式HydroBASINS本物L7多角形を用いて両地点が同じH7か検証、河川線の原本HYRIV_IDと流下方向も照合。未確認段階の距離だけの結果を水生-陸域の近傍生態効果やTN/TP同局実測と呼ばない。
- 本人が求める下水固有の生態系検定には、**実際の受容河川TN/TPが下水由来と特定可能**であることに加え、河岸植物の局所出現/被度、蝶の現地幼虫寄主利用/個体群量・共通調査日が必要。元NUTRI/HNSはGBIF生息地嗜好であり直接実測ではない。

## ゲートと次段階

1. まず原本schemaと座標の対応率を出す。EEA Biology native siteキー2,485のうち元水質SpatialObjectsに同一キーがある数、座標のある数、生物側の座標の有無と品質。生物サイトの座標が不足すれば **STOP_NEAREST_WQ_SITE_GEOGRAPHY** とし、対象H7の公式EEA biology site座標ソースが必要。
2. 座標がある場合でも、距離のみの候補は **SOURCE_OPPORTUNITY_ONLY_NO_HYDROLOGICAL_VALIDITY** として、≥30報告正H7期間、≥20 H7、≥10真H6、≥3国（元の厳格same-site-year source gate）達成したとは扱わない。**同一局で101/12、source報告正9 H7×期間という元判定は不変**。
3. 水質2026 Part1/Part2を再読込するときは、取得元ZIPと本当のsample year・total N EEA_31615-01-7／total P CAS_7723-14-0・RW・W・標準単位・LOQ・≥3 distinct water samples・同じschemeを厳格に保持。観測元sourceを別途照合するまではゼロや欠測と記さない。
4. **生物応答の係数は本検定で一切計算しない。** main・chocho・既存Draftを改変しない。PR #107だけで記録し、否定結果/STOPを消さない。
