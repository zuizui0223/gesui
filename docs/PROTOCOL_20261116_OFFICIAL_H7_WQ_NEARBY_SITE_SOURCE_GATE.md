# 2026-11-16 追加固定：250m/1km近傍のWQ局が本当に同じHydroBASINS H7内か

**公開結果確認後のsource-only対抗検定。元の同一EEA Biology監視局×年でTN/TPとWWTP報告正をともに得られた9 H7×期間のSTOPは不変。**

### 元資料と識別

- 2018/2022年EEA Biology 2,485 distinct native観測局・3,525局×年、元H7×期間454行/232真正H7/90 H6。観測局の元緯度経度と国/原本ID/識別子scheme。原station source SHA256 `a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89`。
- EEA Waterbase ICM 2026公式 `WISE6_SpatialObjects_DerivedData-csv.zip` の2,797,036B、source SHA256 `aa505be6f85d06c675e0d0fcb148a9148f80162f8e9590740a887cc95a3e349d`。全90,859行、identifiers原本完備88,914行、source lat・lon両方非空71,592行、そのうち識別子完備で両方非空69,821行、正規WGS84 geocoded 69,820行（0,0の無効プレースホルダー1件）。元WQ局実値はリポに再配布しない。
- 真正Europe HydroBASINS native Level7 official `hybas_eu_lev07_v1c.zip` SHA256 `b9b810f4db271f0fa6b39fecf01784041b948944b3c3f296a393766388e86b44`。元生物対象232 H7ポリゴンだけを読み取り、各WQ候補は**元生物局と同じ本物のH7ポリゴンに両地点座標が含まれる場合のみ**厳格な空間候補とする。PFAF、近傍H7代理、曖昧な行政範囲で補完しない。
- H7源の親H6は元`H7_454_HNS_panel.csv`の genuine repaired polygon parent source に限定。新しい2026WQ座標で新たなH6親を推測しない。

### 候補規則と固定感度

1. 生物局元2,485 site-IDは元2期の`(countryCode, site, scheme)`に厳密に対応し、別のWQ局も同国・IDキー重複なし・WGS84座標有効。元生物局から測地線直線距離 **250m、1km、5km** の3段階を事前に固定。**主要の機会指標は1km**、補助は250mと5km。
2. 元生物局そのもののH7と公式多角形を照合。元H7に入らない元生物局は品質異常として除外。H7内に別のWQ局があっても、**同じHydroRIVERS河段・流下方向の一致は未確認**として区別。
3. 同一座標(≤25m)のID違いは観測局aliasの疑いがあり、別の生態学的調査地点としない。各閾値で全候補、>25m距離を持つ候補、元生物局×年/H7/H6/国数を比較し、肯定だけを選ばない。
4. 次段の水質検査で2018/2022年の総N/総Pを判定する際は、**最初に全ての同一H7/≤1km候補ID**を固定し、その後で元EEA2026 Part2 6,450,192行／Part1 96,597,294行を本物のZIP SHAで照合する。生物応答値やNUTRI/HNS、元EQRは検査前にも検査後にもsource選定には用いない。
5. 近いWQ別局の実測TNTPは同一EEA Biology局での測定とみなさず**別サイトの受水環境proxy**と呼ぶ。上流正報告の少なくとも30 H7×期間/20 H7/10 H6/3国という元のsupport gateを事後変更しない。たとえ近傍WQ局で数を確保しても、河川の同じHydroRIVERS reach、流下方向・施設位置が同一でない場合は生態学的な因果パス推定をSTOP。
6. N/P記録ありでもそのN/Pに下水が実際に寄与した割合は不明、報告源「未一致」は無曝露ではない。既存の水域EQRと陸域植物/蝶の正負を、新しい水質局による別地点代理と混同しない。

### source-only停止規則

- 元source SHA、station国/ID/scheme、H7本物の多角形ID、WGS84座標の確認に失敗したら **STOP_SOURCE_GEOGRAPHY**。候補局の水質値にアクセスしない。
- ≥1km近接候補があっても同一H7が少ない場合、**近傍で生物観測や処理水N/P影響を測ったことにはならない**。同じH7内の実座標と年別WQ観測を照合してから再評価。
- これは先行main肯定群の再解析・新メタ解析ではない。作業記録はDraft PR #107に保存、main/chochoを改変しない。

## 2026-11-16 source-only Part2実測値の取得規則（化学値の再取得前に固定）

公式H7多角形上で**同じ真正H7内・別WQ局・1km以内**という事前条件を通過した原sourceは、326 Biology実局・402 Biology局×年・429 WQ別局・99 H7・59 H6、station-year候補キー組575（元GIS Action #43804933915）。対象情報は既に取得した候補台帳の生物局ID×scheme×year・WQ局ID×scheme×country・距離とH7/H6だけ。

最初に **publisher 2026年版Part2** の原本ZIP SHA256 `315396800d76e46db415164851fb4f99765b7ec2654828fc42782d902e679a7a` を再照合する。2026年版Part2の真のannual source 6,450,192行を読み、凍結した575組のWQ局×2018/2022年の水質測定があるかだけチェック。真のTN `EEA_31615-01-7`／TP `CAS_7723-14-0`、RW河川・W水、mg{N}/L・mg{P}/L、original sample count≥3、LOQ未満を無効とし、site+scheme+country+year+nuntientで採用できるannual source行が**ちょうど1**でなければ不確実とする。元の数値TNTPの大小で候補を選ぶことは禁止。

同じ生物局×年に複数の近傍WQ候補がある場合、**主要感度**は「水質測定の有無を見る前に距離が最小のWQ局（同距離なら原ID）」を採用し、この近い局に真の同年TNTPがあった場合のみeligible。**副次の源取得機会上限**として、「候補全体のうち最低1別WQ局に両TNTPがある」場合も報告する。両者を混同せず、同一観測局研究の元ゲートSTOPを残す。25m以内の局ID別名可能性は別集計し、近傍別IDを独立現地反復と誤認しない。

生物局の上流WWTP正の報告があるかは、既存の**原本source-only**公式Part1/Part2合算Action #43800794527が保存した「元3525 EEA station×year」のsource報告フラグだけで結合し、元source positive station-year=132と元report-positive H7period=56の母集団を変更しない。元生物応答EQR、PlantN、蝶HNSは開かない。

新しい近傍別局源の30 H7×期間/20 H7/10 H6/3国は**補助的なSOURCE OPPORTUNITYの最低規模**であり、同じ生物観測局での実測の閾値を達成した意味ではない。近傍別WQ局の数を満たしても、河川ネットワークの同じ河段・流下方向／流域内源の具体的到達は次に検証する。Part2源が疎なら、Part1オリジナル96,597,294 sample rowsとのunionをsource-onlyで追加するが、同じ観測局原本Part1+Part2を独立研究として重複カウントしない。Part1の成功済み元の strict same-site-year STOP を変更しない。
