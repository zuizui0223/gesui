# 2026-11-16 事後追加source gate：EEA2026 TN-only／TP-onlyの独立流域反復は足りるか

**Part1/Part2同一局・同年TNTP unionが101 site×year、正の報告WWTP sourceは12 site×year（9 H7×period／8 H7／8 H6／3国）で主Source gate STOPとなった後の別設計。事前登録や陰性結果の上書きではない。**

## 母集団を維持
元454 H7×period/232 H7/90 H6に位置するEEA benthic source計3,525 exact original native station+scheme+calendar year（2018/2022）のみ。2026年EEA公式Part1（original SHA256 `ff779521749b40b952cce35a4671240c3b6c47d87f32a8fef9240a4ac503a082`）とPart2（`315396800d76e46db415164851fb4f99765b7ec2654828fc42782d902e679a7a`）を同じ局年・各TN／TPでPart1優先・Part2が有効な場合のみ補完して作った、source-verified 3,525行unionのみ使用。局年の重複を絶対に許さない。

## 対象と結果を開く前のゲート
- 主ソースは総TN`EEA_31615-01-7`と総TP`CAS_7723-14-0`、厳格water RW/W、同局・同年、LOQ条件及びそれぞれ元3試料以上。硝酸・リン酸等にすり替えない。
- **TN単独**、**TP単独**、元の**TNTP両方**をすべて集計して開示。各々、原本上流WWTP源報告が正のsite×year数、異なるH7×period／distinct H7／真のH6／country／2018/22期間を計数。少なくとも`30 H7×period、20 H7、10 H6、3 country`という元主Source gate条件は維持。TNとTPどちらでも合格基準を下げない。
- 追加の厳しい親H6内の強度比較：各真正H6で`別々のH7二つ以上`に、源報告正かつ当該実測TNまたはTPが存在すること。≥10独立H6を別ゲートとする。たとえ国数・H7数を満たしても、2期間を2空間反復として数えない。
- 生物目的変数W/NUTRI/HNS/EQRは読まない。元源報告positiveがないサイトも施設なしとは言えず、負の真の対照に使わない。
- 元のTNTP union主ゲートSTOPを維持し、この探索的別課題の通過があっても「元の同時TNTPでの下水媒介合格」とは呼ばない。排水の報告質量と受容河川の実測TN/TPの関係は源独立性が未確定。

## 判定
- 単独TN/TPが成功→次の**別の**非因果的な水質構成モデルのsource検証が可能かを考える。生物効果をここで推定しない。
- 不成立→物理局を増やす／国別の観測局横断照合／真の処理水混合率の独立測定源が必要。モデルや閾値を変えて元の肯定的な蝶HNS関連を救済しない。

**元main PR #106は既にマージ済み。source実験はDraft PR #107だけに保存し、main/chochoの既存研究結果は編集しない。**
