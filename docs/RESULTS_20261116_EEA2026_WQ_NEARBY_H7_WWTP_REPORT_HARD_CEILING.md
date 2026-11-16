# 2026-11-16：公式水質観測局を近傍同H7で追加できるか——原本監査の結果とSTOP

**対象：EEA Waterbase Water Quality ICM公式2026年版・生物監視局の2018/22年、下水由来のN/P到達と植物―蝶の生態媒介を実証できるかについての事後source-only検証。生物結果は今回開いておらず、実測負荷による生態効果・因果SEM・個体数を推定していない。**

## 1. 位置データの品質と同一ID再現

- 真正EEA Biology EQR source：原2018/2022の物理生物監視局2,485・3,525 station×year、原HydroBASINS 454 H7×period/232 H7/90 true H6。元EEA source SHA256 `a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89`。**2,485局すべて原 source lat/lonあり**。
- 真正公式水質位置：EEA 2026 [WISE6_SpatialObjects_DerivedData-csv.zip](https://sdi.eea.europa.eu/datashare/public.php/dav/files/ZcJ37WK6oaWonFJ/WISE6_SpatialObjects_DerivedData-csv.zip)、2,797,036B、原本SHA256 `aa505be6f85d06c675e0d0fcb148a9148f80162f8e9590740a887cc95a3e349d`。90,859 source rows；緯度・経度両列が非空71,592行；ID+scheme+国まであり座標非空69,821行；正規WGS84有効座標69,820行。空欄でないが原本完全なIDを持たない座標が1,771行、実座標として不適切な(0,0) placeholderが1行。IDなしを誤って座標欠損と解釈しない。2,485旧生物siteキーの同一`country|monitoringSiteIdentifier|Scheme`は公式WQ空間sourceにも2,485すべて存在した。
- 原本WQ座標（約7万station）はGitHubに転載しない。原本site/年/座標と真のH7は原 source Action SHAで再現可能。
- 検証：[source coordinate QA #43804308185 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43804308185)。

## 2. 元H7公式多角形の厳格照合

欧州 [HydroBASINS v1c Level07公式ZIP](https://data.hydrosheds.org/file/hydrobasins/standard/hybas_eu_lev07_v1c.zip)、20,992,584B、SHA256 `b9b810f4db271f0fa6b39fecf01784041b948944b3c3f296a393766388e86b44`。対象の232 H7 polygonを全て検証、invalid0、元EEA Biology3,525 source site-yearの元H7内に位置しない局**0**。別IDの公式WQ局が、元生物局と同じH7多角形内に位置する場合だけを機会の候補として認める。地理distanceは大円距離、元station keyとID/scheme一致を保持し、自己と別IDを混同しない。

| 元source strict geo site opportunity | 250m | **1km** | 5km（感度） |
|---|---:|---:|---:|
| 別WQ局が同一H7にある元生物サイト | 86 | **326** | 1,738 |
| 元生物サイト×年 | 110 | **402** | 2,343 |
| 対応する異なるH7 | 40 | **99** | 205 |
| 真正H6 | 31 | **59** | 83 |
| 区分国 | 8 | **10** | 15 |

1kmのsource candidateは429個の異なるWQ実station identifier（完全`country|site|scheme`）、生物station-year×別WQ stationの組575。≤25mの同一位置・別IDはalias候補としてフラグ（≤1km候補のsite-year×WQうち42組）を保持。別WQ局が生物stationと全く同じ測定所・同じHydroRIVERS reach・上流WWTP sourceに接続されるとはまだ判定していない。

[source source verified #43804933915 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43804933915)。

## 3. 近傍水質局の2018/2022年原本Part2での水質観測機会

生物応答・TNTP実測値を見ずに上述1km同H7の575組を凍結。その後、公式2026年版**Part2** `WISE6_AggregatedData-csv.zip` (170,422,309B、SHA256 `315396800d76e46db415164851fb4f99765b7ec2654828fc42782d902e679a7a`) の原本6,450,192行を再取得・全行走査。575候補に出る**531 WQ station×yearキー**の2018/2022年を検査。RW河川、水matrix W、総N`EEA_31615-01-7`/総P`CAS_7723-14-0`、元素別標準mg/L、平均値・LOQ品質、≥3試料、各 nutrient×station×yearにちょうど1本のvalid年値の条件。

- このPart2で、有効な**TN/TP両方が同一年・別WQ局**にある近傍sourceは **1生物station×year**（1 H7×period）。その上流WWTP報告正は **0**。
- 生物stationと近傍WQの**同じH7内・≤1kmの地理近接**は確認したが、河川流下線・源距離・流域中の栄養流量を確認できていない。公式Part1の個別水質試料をここではまだ重複ダウンロードしていないため、Part2の1件を「実在する化学測定すべてが1件」とは呼ばない。
- [source source verified #43805349153 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43805349153)。

## 4. より強い結論：上流WWTP報告正のsource上限で既に通過不能

近傍WQ局にTNTPがあるかどうかを**見ずに**、元EEA Biology3,525 station×yearのsource positive=132を、固定済み1km・同H7の402 source Bio station-yearに照合。近傍WQのTNTPが**全件揃っていたとしても**：

| 元のWWTP報告正に関する機会上限 | 250m・同H7 | **1km・同H7** | 固定閾値 |
|---|---:|---:|---:|
| 上流WWTP報告正の生物station×year | 19 | **47** | — |
| source positive H7×period | 10 | **16** | **30** |
| source positive distinct H7 | 8 | **13** | **20** |
| source positive native H6 | 8 | 12 | 10 |
| source positive country | 6 | 7 | 3 |

同位置aliasの疑いがある**≤25mを除いた1km**ではreport-positive source H7periodは14、H7は12。よって、元の固定条件30 H7×期間・20 H7という二つの閾値を満たすことは**数学的に不可能**。源に一致する施設報告なしは無放流ではないし、正の報告源の施設から近傍WQ stationへN/Pが流れたことの証明でもない。

このsource hard upper bound は**公式Part1を追加しても破れない**。Part1の約96,597,294行/1.639GB original ZIPを再び重複取得しても、報告正site×yearの数は増えないため、この元近傍1km source cohortでは追加大容量試料の処理は行わずSTOPとする。[hard ceiling #43805542822 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43805542822)。

**結論：`STOP_1KM_SAME_NATIVE_H7_REPORTED_POSITIVE_SOURCE_HARD_CEILING`**。元の2018・2022年公式Part1+Part2、同一EEA Biology局の厳格両TNTP12報告正station×year/9 H7×period、そして今回の「近傍別局でも上限16 H7×period」は別設計・別推定対象であり、統合して解釈や効果量を作らない。

## 5. 生態学的な評価と具体的な次段

- **すでに確認できた肯定：** `main` 統合済みPR #106の欧州同一H7での植物NUTRI―蝶HNSの正の地理的対応、真正H6内部183 H7・40 H6の栄養生息地嗜好・属レベル形質構成の正、米国EPT physical trait filtering。これらの結論はそのまま残る。
- **同時に残す否定：** 同一H7の時間的同期は出ず、同属内EIVラベル入替p=0.6864、W→植物NUTRI/蝶HNSは真正H6内部で不確実、新しい10 H6の外部主指標予測も不確実。実測の受水N/Pが下水のものだと識別できない。
- **合理的な次段は、既存の454 H7×periodにこだわってsource-positive地点を増やすのをやめること。** 公式EEAのもっと広い2012–2024 Biology原本にある異なる水生種群・測定局を探す／年別UWWTD施設source reported loadの対応年と重なる生物観測時期を選ぶ、あるいは施設レベルで実測流入→下流近接河道濃度→植物/蝶の真の現地データが同時に得られる独立監視系を事前定義する。元の454/183H7の肯定を実測下水の因果効果と解釈しない。
- 元の源に基づく方向特異性が弱いこと（PR #92・#93）も残る。H7全域におけるWWTP報告年との関連を生物間相互作用と断定しない。
- 5kmで別局を見つけても、**同じHydroRIVERS river segment・流下方向の保証なしに生物試料・環境負荷を結ぶことをしない**。単なる距離緩和で30 H7periodを作らない。

## 再現性と未マージ範囲

- 原GitHub Actions：#43804308185, #43804933915, #43805349153, #43805542822（source-original SHA verified, successful）。
- コード：`analysis/20261116_reported_np_source/audit_biology_wq_spatial_source_v1.py`、`audit_official_H7_biology_wq_geometry_v2.py`、`scan_eea2026_part2_nearby_sameH7_TNTP.py`、`nearby_WQ_reported_positive_upperbound_sourceonly.py`。
- 元CSVと原本source ZIPは再配布しない。各runは元IDのsource-only候補台帳、repoにはsource由来の結果summary JSONだけ保存：
  - `results/20261116_reported_np_source/EEA2026_BIO_WQ_SPATIAL_SITE_QA.json`
  - `results/20261116_reported_np_source/EEA2026_SAME_NATIVE_H7_NEAR_WQ_CANDIDATES_BRIEF.json`
  - `results/20261116_reported_np_source/EEA2026_NEARBY_SAMEH7_PART2_TNTP_SOURCE_BRIEF.json`
  - `results/20261116_reported_np_source/EEA2026_NEARBY_1KM_REPORTED_SOURCE_HARD_UPPERBOUND.json`
- この追加ゲートは[Draft PR #107](https://github.com/zuizui0223/gesui/pull/107)。**main・chochoは更新せず**、元PR #106の肯定・否定と以前のsource STOPは不変。
