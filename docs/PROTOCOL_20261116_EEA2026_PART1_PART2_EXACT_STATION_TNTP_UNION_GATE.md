# Source-first: 2026年版EEA Waterbase ICMの個別試料と年集計を用いるTN/TP検証仕様

**2026-11-16、元研究の同一流域W・植物NUTRI・蝶HNS関連と2026 Part2の現地source欠測を確認した後の事後設計。新しいPart1の試料結果を開く前に固定。独立前向き事前登録ではない。**

## 原本とコホート

- **2つの公式2026出版ファイルは区別**する：[EEA公式 Water Quality ICM 2026](https://www.eea.europa.eu/en/datahub/datahubitem-view/fbf3717c-cd7b-4785-933a-d0cf510542e1) の`WISE6_AggregatedData-csv.zip`（Part2、170,422,309 B、SHA256 `315396800d76e46db415164851fb4f99765b7ec2654828fc42782d902e679a7a`）と`WISE6_DisaggregatedData-csv.zip`（Part1、1,639,231,716 B、SHA256 `ff779521749b40b952cce35a4671240c3b6c47d87f32a8fef9240a4ac503a082`）。Discodata SQL latestは旧源であり、2026最新の母集団とみなさない。
- 元EEA Biology EQR source（EEA18／EEA22）の**真正3,525 monitoringSiteIdentifier+monitoringSiteIdentifierScheme+2018/2022 year**／2,485固有physical site（H7×期間454、H7=232、H6=90）。元station割当JSON SHA256 `a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89`、元H7パネルSHA256 `1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca`。
- 元2018・2022年EEA上流報告下水N/P mass元station blob SHA1 `0f30bfbc3dac13f6a182c8a910e9899c04778322`および`325f199e985d9862ffee5b0edb205c61a7be3833`。元3,525 site×yearのうち**132 site×year・56 H7×期間**にのみ報告された正の上流N/P負荷を持つ。無照合はゼロ曝露・非放流の対照ではない。
- 元生物EQR結果の大きさや植物NUTRI・蝶HNSの係数を、**source cohortを決める前に見ない**。今後の解析は新しいsource-supported site-year/H7×periodとして成立したときのみ。

## 元水質の厳格選別

- 生物学的調査と同じ`monitoringSiteIdentifier`＋`monitoringSiteIdentifierScheme`＋calendar year 2018／2022、country code、河川`RW`、水試料`W`を完全一致。
- 元TN determinand `EEA_31615-01-7`＋単位`mg{N}/L`、元TP determinand `CAS_7723-14-0`＋`mg{P}/L`。総Nと硝酸N、総Pとリン酸Pを代替しない。
- Part1個別試料は2018/22その年に測定され、`resultQualityObservedValueBelowLOQ`がfalse／empty、非負で有限、元`sampleIdentifier`により重複を除外して**各栄養塩3件以上の異なる水試料**がある局だけを採用。試料間の実測値を算術平均した年別局値を算出。
- Part2年集計は測定原本`resultMeanValue`、`resultNumberOfSamples>=3`、LOQ flag false、年/単位の確認が成立した局のみ。元両表の年別局が重複しているなら**Part1優先**とし、原本値の一致・違いを記述する。Part1 TN/Part2 TPを単純に組み合わせて「同日測定」とは記述せず、同じ年・局に独立して値があることだけを記録。元site×year×TN/TPに複数の互いに不整合な測定集計がある場合はそのままflag、恣意的平均値を作らない。
- 元Part2は正式な6,450,192行を全走査済み。厳格同局TNTP 20 site-year（すべてBE）、そのうち元上流報告正は1 site-year。**Part2はsource gate STOP**（先行レポート参照）。この結果をPart1の先入観にしない。

## Part1および源の合併（結果の方向とは独立）

- Part1のみで正の上流報告がある同局同年TN/TPの有効site×year、H7×期間、distinct H7、H6、国数を記述。厳格母集団の事前最低**30 H7×期間、20 H7、10 genuine H6、3か国**をすべて達成した場合のみ次の生態関係を見る。
- Part1とPart2の和集合は、site×year×nutrientで重複処理・ロジックを**原本の一致条件だけで固定した別gate**として記録。個別試料と集計表を独立野外標本として二重計数しない。Part1で解けた正例が増えても、EEA年度の水質測定が実際の下水由来物質フラックスであるとは断定できない。
- 正の上流報告源あり+同局同年水質が揃うH7が充分でも、観測N/Pと放流源N/Pの相関については農地源・処理段階・混合希釈・H7内の別流域位置・報告不能の施設の交絡が残る。

## 作用検定へ進む場合（厳格source gate合格時のみ）

- EEAの同一生物調査点における水質N/Pは、**真の河川中のTN/TP濃度**。静的WやEEA上流の報告年別t/yrは異なる変数。水質濃度との関連を因果的な下水効果に読み替えない。
- 推定単位・サンプル数と、国×生物EQR評価方法・土地利用・河川流量・GBIF観測努力を元mainと完全に揃え、元H6 block out-of-fold予測を使用する。完全に同じH7×期間に植物NUTRI・蝶HNSを存在し得るが、**同じ河岸スポットではない**ことを必ず記す。
- 本段階で扱うのはsource availabilityであり、新しい正の β、SEM/媒介/生態系サービスを作らない。元PR #106で統合した肯定・否定（時間的非同期、属内NUTRI不確実、報告源WのH6内不確実）は変更しない。

**研究の優先：** 作用パスに必要な局所曝露実測がなければ次の回帰へ進まず、年度一致のサンプル源と河川ネットワーク位置が不足することを明示する。
