# EU H7 454 × period: 全Plantae speciesKey → 欧州NUTRI完全照合と蝶HNSの条件付き共変動

**2026-11-15 凍結仕様、元データ確認後・結果を見る前の追加実装／既存PR #101後の探索的研究。独立した前向き事前登録ではない。**

## 第一段階：元植物種の完全照合

- 元パネル: 真正454 H7×period、232 H7、90 H6。SHA256 `1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca`。2017–2019／2020–2022のみ。H7-periodの重複・脱落禁止。元の `country` は15 source labels、`country_holdout` は11 H7多数派ラベルで別物。
- 元Plantae speciesKey 234 H7×2期間×6 facet（original 2808 facets）のうち、元454 H7×期間のみ。真正JSON SHA256 `251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4`。元Union speciesKey **11,511**種（先行の源資料ゲートで確認）。同じGBIF記録を独立反復と数えない。
- 元の独立したPR #12植物1006-key speciesKey crosswalkを再利用（EIV原本 SHA256 `ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba`）。既知624種は既存辞書との照合一致を用いる。他の元speciesKeyはGBIF公式API `/v1/species/{speciesKey}` でIDごとに最新タクソノミーを取得、エラー／見つからないキーは**未知**として保存、植物欠如としない。
- 欧州植物N形質の原本: Tichý et al. (2023), *Journal of Vegetation Science* DOI:10.1111/jvs.13168 の公式 `EIV_Europe.csv`、ソースSHAを固定。元行数8,908、数値NUTRIは1–12、種名は原本の `Unnamed: 0`。原本の曖昧な名称重複を捨てる。正確な二名法のみ（完全一致、accepted `SPECIES`、`Tracheophyta` phylum）、属フォールバック／fuzzy／補完なし。
- 1006-key外部辞書との共通キーは古いGBIF original nameを維持し、今日の再照合との差異は別報告。年代によるタクソノミーのdrift、欠測は原因が不明な限りサイト0としない。API失敗を `NO_MATCH` と同一視しない。
- 実行は決定論的な8 shard（GBIFキーのhash）、最大同時3 shard・1 shardあたり最大6並列HTTP。source errorsは各keyでHTTP statusと理由を記録し、API再試行後も欠測ならそのまま残す。レスポンス失敗を正の値に代入しない。

## 第二段階：**形質アウトカムを開く前**の対象支持確認

元454 H7×期間それぞれに次を算定する（HNSやWの回帰はまだ行わない）。

- 元GBIF植物記録数≥50（元の454行はすべて既に満たす）
- Accepted vascular speciesKeyの数を保存。
- **全Plantae source speciesKey**に対する実測NUTRI scoreカバー率≥40%、EIV得点種10以上。
- 生態学的結果を見る前に、10数値種・40%を満たすH7×期間数≥100、異なるH6≥40、元 `country_holdout` ≥5、かつ全キーに対するGBIF taxonomy API status OK ≥99%という診断ゲートを固定。未完取得が1%以上の場合は環境モデルへ進まない。各サイトで未知API種が多い場合は形質平均の解釈とその割合を別記。

前段の624得点種だけで112 H7×期間、39 H6、6 country_holdout。**この既存辞書は下限**であり、新規約10,500キーの照合後に可観測性が改善するかは未確定。

## 第三段階：モデル（第二段階が成立したときのみ）

### 固定アウトカム

主 `Plant_N_incidence`：同じH7×期間に出現した**EIV得点確認済み種のみ**の種均等NUTRI平均。
副 `Plant_N_record`：得点確認済み植物speciesKeyの元GBIF記録数で重み付けした平均。
記録数は観測インデックスで植物被度・実個体数ではない。

### 明示的なモデル比較

**A: W→PlantN 条件付き関連（未検証）**

`Plant_N_incidence ~ standardized W + source geography/effort controls`。Wは静的モデル排水関連地理、河岸の実測放流Nではない。

**B: PlantN×Butterfly HNS 追加説明力（未検証）**

`Butterfly_HNS ~ standardized W + standardized Plant_N_incidence + source controls`。蝶HNSは同じ454行から取った元の記録重み付き平均、補助は種均等平均。前段で用いた source-locked controls と国×method・期間を保つ。植物Nを**媒介変数と呼ばない**。共有Wや都市化・土地被覆・記録努力による交絡はまだ残り得る。

**検証原則：** H6を跨がない決定論的5分割、train-only標準化、固定Ridge alpha=10、H6均等 heldout MSEと元H6単位の2000回固定誤差再標本化による記述95%区間（モデル再学習CIではない）。観測回帰はH6クラスタWLSの記述95%区間。国丸ごとのleave-country-out予測は補助。係数と予測指標を混同しない。標本閾値、調整変数、ランダムseedは事後に変えない。

### 判定／中止

- 主要な植物NUTRI支持基準が未達、原本checksum不一致、taxonomy lookup失敗多発、またはモデルの交絡・線形従属性が大きい場合、原本品質監査だけを公開し、推定・因果パスは**STOP**。
- Wの植物N係数が正でH6-cluster95 >0でも、H6留保MSE改善95%が0を含むなら「地理外挿未確認」。逆も同様。
- W→PlantNとPlantN↔HNSがともに正でも、それは同一流域での関連であり、`W→Plant→Butterfly` の因果媒介ではない。
- PR #96は同じH7で水EQR―蝶HNS残差相関がゼロを含んだ。結果が不利ならそのまま残す。米国EPT XEMBEDと欧州H7 Wは曝露が異なる、単一効果量メタ解析をしない。

**科学的責任：** 局所実測N/P、葉N、幼虫利用、成虫吸塩、水生昆虫羽化量はこの解析に存在しない。植物の種組成と潜在的寄主形質に限定して結論を書く。ソース未承認の原生データや写真は再配布しない。
