# 2026-11-16：EEA公式2026年版TN/TPの同一観測局・同一年の原本完走と最終Source gate

**Draft PR #107／元main統合PR #106の肯定・否定結果は維持／実証したのは原本到達・監視局化学データの有無。新しい生態効果モデルは一切実行していない。**

## 1. 技術的な取得問題は解消した

EEA公式[Waterbase Water Quality ICM 2026](https://www.eea.europa.eu/en/datahub/datahubitem-view/fbf3717c-cd7b-4785-933a-d0cf510542e1) の原本は公開Nextcloud（`/datashare/public.php/dav/files/ZcJ37WK6oaWonFJ/`）で直接取得できることを[Action #43797712726](https://github.com/zuizui0223/gesui/actions/runs/43797712726)で確認した。これは [WISE SoE](https://water.europa.eu/freshwater/resources/wise-soe-data-collection) の更新が保証されない旧Discodataと明確に区別する。

- **Part 2、年別site mean：** 正式ファイル `WISE6_AggregatedData-csv.zip` 170,422,309 B、公式原本SHA256 `315396800d76e46db415164851fb4f99765b7ec2654828fc42782d902e679a7a`、CSV全6,450,192行。源SHA・全行走査[Action #43798318436 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43798318436)。
- **Part 1、試料別：** 正式ファイル `WISE6_DisaggregatedData-csv.zip` 1,639,231,716 B、公式原本SHA256 `ff779521749b40b952cce35a4671240c3b6c47d87f32a8fef9240a4ac503a082`、ZIP method9 **DEFLATE64**（Python`zipfile`不可）。事前の元ヘッダsource検査[Action #43799308056](https://github.com/zuizui0223/gesui/actions/runs/43799308056)のあと、`7z x -so`で**28,458,399,518 BのCSV**を解凍内容を全保存せずに走査。全**96,597,294行**を検査し、元日付は`phenomenonTimeSamplingDate`（Part 1には`phenomenonTimeReferenceYear`列がない）から2018/2022年を抽出。[Action #43799818906 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43799818906)。

いずれも生物調査元H7×period=454、232 H7、90真正H6から、EEA Biology 2018・2022の**元3,525 native physical station×year（2,485固有station）**を先に固定。原station SHAs、元2018/2022 WWTP正報告site-year=132、正報告H7×period=56を固定した上で、水質を読んだ。

## 2. 元EEA EQR局と同一年の総N・総Pだけを受理

- `monitoringSiteIdentifier`と`monitoringSiteIdentifierScheme`、年2018/2022、country code、`RW`河川、`W`水試料を**完全一致**。
- 総窒素`EEA_31615-01-7`および総リン`CAS_7723-14-0`、単位それぞれ`mg{N}/L`と`mg{P}/L`。硝酸塩・リン酸塩を総N/Pとして代用しない。
- Part 2はその局年で元報告の試料数≥3、有限値、品質フラグ（below LOQ）受理。Part 1は**元の`sampleIdentifier`（ない場合厳格日付/値）で重複を除いた異なる水試料3件以上**。Part 1からの年meanは個別値の算術平均として記述し、同日性を他の生物観測に仮定しない。
- Part 1とPart 2を合わせるときは**局×年×栄養塩**を1単位とし、両方有効なら Part 1優先、片方だけならその観測値を残す。同じ実測を独立研究・2反復に二重計数せず、両表間の元値不一致があっても符号に応じてレコードを除外しない。

## 3. 原本 source gate 実測結果

| 指標／同じ3,525 EEA Biology元station×year | Part 2 年別 | Part 1 個別から算出 | Part 1＋2 source union |
|---|---:|---:|---:|
| strict TN対応（site×year） | 20 | **166** | **186** |
| strict TP対応（site×year） | 20 | **242** | **262** |
| **同じsite×yearでTNとTP両方** | **20** | **81** | **101** |
| **そのうち上流の年別WWTP N/P報告に正の値** | **1** | **11** | **12** |
| WWTP報告正・TNTPとも有効な H7×period | 1 | 8 | **9** |
| WWTP報告正・TNTPとも有効な distinct H7 | 1 | 7 | **8** |
| WWTP報告正・TNTPとも有効な真正H6 | 1 | 7 | **8** |
| WWTP報告正・TNTPとも有効な国 | 1 | 2 | **3** |

- Part2：2018年局年8件、2022年12件で**ベルギーのみ**、内上流報告正1件。元soure quality breakdown [Action #43798670501 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43798670501)。
- Part1：元個別試料CSV2018年5,800,109行、2022年9,130,406行が存在し、源station+国の一致で約183.5万行。元観測局×年で検証済みTN 166・TP 242・両方81、上流報告正11。元source抽出 **完全一致**かつ3件以上の実測試料であり、全EEA集団の無測定を示すものではない。
- **Union source gate [Action #43800794527 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43800794527)**：真正site×year単位で Part1優先・Part2補完によりTNTP 101、上流報告正12、**9 H7×period／8 H7／8 H6／3国**。源全部をGitに再配布せず、元SHAとsource-only抽出結果・原run IDsで再現可能にした。
- 2026-11-16に**Part1の完全原本取得と集計式を初めて実現**したことと、実際の監視源標本が不足することは別。前者は技術的成功、後者は科学的STOP。**旧Discodataが古い→水質が全くない**という誤った解釈は棄却した。

## 4. 固定済みソース支持条件と科学的判定

| 事前に固定した最低支持条件 | 必要 | 原本 union | 判定 |
|---|---:|---:|---|
| WWTP報告正＋両栄養塩のH7×period | 30 | **9** | 未達 |
| distinct H7 | 20 | **8** | 未達 |
| genuine H6 | 10 | **8** | 未達 |
| 国・測定報告区分 | 3 | **3** | 達成 |

**`STOP_SAME_SITE_YEAR_REPORT_POSITIVE_MEASURED_TNTP_SUPPORT`。** 上流報告源の真の強さと、同じ観測局×同年の実測TN/TPのパスを、他地域へ外挿できるだけのsource supportはまだない。少標本の全試行で肯定的な方向を恣意的に選ぶため閾値を下げるべきではない。

また、得られたのは**受容河川の環境水質**であって、上流処理場からの**実際の投下N/Pの寄与率を同位体等で識別**したわけではない。河岸植物GBIFと蝶GBIFは流域スケールであり、EEA生物測定局と同じスポット・訪花日・個体群サンプルではない。「水→植物→蝶」の因果媒介や生態系サービスは実証していない。

## 5. 次に進むための合理的な選択肢

1. **固定済みsource gateを無理に通した扱いにしない。** TNとTPの同時site-yearに限定した機構候補は今の原本では不成立。TPだけ、TNだけの対象拡張は必ず**独立の事後探索的テーマ**として別のフリーズ仕様とソース品質チェックを要し、同一site-yearのTNTP共観測と同じ推定量とは呼ばない。
2. **別の実測源**（国別水質監視やPeters 2026 true TN/TPなど）を使うなら、元EEA EQR物理stationと報告施設の上流系をソースだけで照合できる独立データ集を探す。河川全体の広域TNTP平均を単一EEA局の水質として置かない。
3. **本来の下水の生態学的帰結**は、同じ受水地点の処理水由来割合・水温・流量・栄養濃度・植物（被度・寄主）・蝶（種数、観察努力、幼虫摂食）が得られてから検証。主結果は「属レベルに保存された植物栄養環境嗜好と蝶寄主HNSの広域・H6内での構成対応」という限定を維持。

## 再現性とマージ状態

- 元統合済み [gesui/main PR #106](https://github.com/zuizui0223/gesui/pull/106) の水域EPT・植物―蝶の**肯定・否定・STOP**結果は変更していない。
- 新ソース検査は別の [Draft PR #107](https://github.com/zuizui0223/gesui/pull/107)。元原本一致のscriptsとSOURCE receiptsを積んだが、**mainには未マージ**。
- Part1 96.6百万行 [#43799818906](https://github.com/zuizui0223/gesui/actions/runs/43799818906)／Part2 645万行 [#43798318436](https://github.com/zuizui0223/gesui/actions/runs/43798318436)／Part1+2厳格union [#43800794527](https://github.com/zuizui0223/gesui/actions/runs/43800794527)。生物応答係数は新規に推定していない。
- 原本公式ZIP2本（個別は1.64GB、年別は170MB）や完全な28GB解凍CSVはGitHubへ再配布していない。アップロードした原本一致の3,525局年分の限定抽出結果も、処理水施設の報告書の「ゼロ汚染」という誤解を誘う用語を使わない。


## 6. TNだけ／TPだけの別の探索的source機会

メインのTNTP同時site-year支持が未達と確定した**後**、別の「TN単独」「TP単独」の source-only gate を新規に固定した。先に固定した`docs/PROTOCOL_20261116_EEA2026_SINGLE_TN_TP_SOURCE_REPLICATION_POSTHOC.md`に従い、同じ元3,525 site×year、元Part1/2 source unionを、栄養塩の測定値の**有無のみ**で3通りすべて集計。生物応答/EQR/HNS/NUTRI/Wの数値を読まない。試料品質・国・報告源・H7/H6等は同一。

| 厳格な条件 | TN単独 | TP単独 | TNTPとも |
|---|---:|---:|---:|
| 全EEA元site×yearで有効な観測 | 186 | 262 | 101 |
| 上流WWTP N/P報告正と一致したsite×year | 12 | 12 | 12 |
| 報告正のdistinct H7×period | 9 | 9 | 9 |
| 報告正のdistinct H7 | 8 | 8 | 8 |
| 報告正のgenuine H6 | 8 | 8 | 8 |
| 報告正の国 | 3 | 3 | 3 |
| **同じH6内で報告正かつ当該栄養実測のある異なる2 H7が共存** | **0 H6** | **0 H6** | **0 H6** |

**結論：** TNかTPに別々にしても、実際に水質値と上流処理場報告が重なったsource populationは**12 site×yearから増えない**。メインゲート（30 H7×period/20 H7/10 H6/3か国）だけでなく、より厳密な同一H6内の正の報告源強度比較（2 H7以上のH6を10）も **0/10** でSTOP。TN-only/TP-onlyの効果推定を試みない。報告なしの地点を真の無放流対照にはできない。

真正元source-only [Actions #43801308202 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43801308202)。コード`analysis/20261116_reported_np_source/audit_totalN_totalP_individually.py`、結果原本ログを`results/20261116_reported_np_source/EEA2026_SEPARATE_TOTALN_TOTALP_SOURCE_GATE.json`に保存。メインのTNTP同時のSTOPとPR #106のmain統合結果は変更していない。
