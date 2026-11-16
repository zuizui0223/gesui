# 2026-11-16 — 既存原本の「上流報告N/P」対「下流報告N/P」で植物NUTRIと蝶HNSを反証

**前の正結果と元PR #90〜#93の負の対照を読んだ後の、事後探索的解析。独立事前登録ではない。2026-11-16の検定結果を見る前にコードと閾値を固定する。**

## 中心的な疑問

元454 H7×期間のうちEIV形質が揃う451行では、静的なモデル排水関連地理Wと植物NUTRI、植物NUTRIと蝶HNSが正の地域共変動を示した。しかし現地の下水由来Nが河岸へ届いたかは未測定。元の別実装PR #92・#93では、**実際に行政に報告された2018/2022年の処理場年間N/P排出質量**を上流・下流方向で地理的に対応づけると、EQR・植物S50・蝶記録重みHNSで上流だけの安定した予測効果は確認できなかった。今回、**新たに得られた植物種均等NUTRIと蝶種均等HNSに、同じ地理的負の対照を適用する。**

## 真正原本と接続単位

- 基本母集団：`results/20261115_sync/hns_macro/H7_454_HNS_panel.csv` SHA256 `1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca`、454 H7×period、232 H7、90真正H6。源のW、crop、urban、riverflow、自然流出、観測努力、country×method、期間、蝶数値scoreカバレッジを保持。
- source方向：元PR #92 `results/20261115_np_direction/H7_454_DIRECTIONAL_SOURCE_AUDIT_PANEL.csv`、原Git blob SHA1 `702131afc7d1a7bfb15e80357a506754cae7a70a`。真の上流・下流EEA処理場報告源を原元報告施設まで辿り、2018/2022年の年次N/P排出**報告質量**と流量モデルから作った仮想完全混合濃度単位、及び報告源有無（0=一致した正の報告がない、真の無下水ではない）を元ソースのまま保持。元の`upstream/downstream_report_positive`は方向づけされた**報告の登録メタデータ**であり介入・実際の栄養曝露ではない。
- 植物：`EU454_H7PERIOD_NUTRI_TRAIT_MEANS_SOURCEONLY.csv`を原source v2 [Actions #43718069015](https://github.com/zuizui0223/gesui/actions/runs/43718069015)から取得。原欧州EIV 4,011 accepted vascular speciesKeyのみ数値、175 accepted雑種に補完なし。原451 H7×period・90 H6・11 majority-country groups。
- 蝶：元の`primary_incidence_mean`（観測種均等HNS）を主、元の`primary_record_mean`を副。植物種均等NUTRIを主、記録重みNUTRIを副とし、一つの統合効果に平均しない。
- すべて`h, period`のユニークな元種IDを保持。source direction panelの`h6, W, primary_record_mean`を同一キーの元パネルと数値チェック、最大誤差1e−10、すべて整合しなければ停止。国の定義が合わない行数は元source ID差として報告し、勝手に置換・削除しない。

## 対象ゲート（生物アウトカムを開く前）

1. 原方向パネル454/454、前ソースの4分類上流/下流報告源(UP-only=21、DOWN-only=26、both=35、neither-positive-report=372)の元同定を検証。
2. 原451 EIV scored H7×periodを重複なしにjoin。source-positive upstream ≥30 H7×period、≥20真正H6、≥3独立国群、下流positive ≥30、2×2の3正報告セル各≥10を満たすこと。境界の報告0を真の非曝露と解釈しない。
3. 2期間は同じH7の擬似独立地理サンプルではない。Wと reported annual load index は3年間の生態GBIF記録と時点が完全一致せず、年次の排出変化・因果的な時差応答は推定しない。
4. 報告上流・下流が流域/施設共有で空間的に強く相関する可能性を想定し、完全分離や高いVIF・デザイン列ランク不足はすべて報告。

## 固定した解析と反証

- **一次検定（植物Nの応答）**：元`P_GEOG` + W + 国×method/期間 + DOWN positive-report indicator を固定基準とし、UP positive-report indicator を1つ追加したときのH6クラスタWLS係数（1原サンプルSDあたり）と真正H6 5fold外予測MSE改善（Ridge alpha10、train-onlyスケール、H6均等、H6固定誤差2000回bootstrap95）を計算。UP→DOWNの反対向き（UPを先に固定してDOWNを追加）も全て表示。
- **二次検定（蝶の種均等HNS）**：元`B_GEOG`（W含む）+ 元植物NUTRI種均等値 + DOWN indicator を基準に、UP indicatorを追加。これも逆方向を含む。植物NUTRIを足さない副次モデルでもUP/DOWNを調べる。
- **強度の感度**：上流/下流の原本報告N質量÷流量の仮想等価濃度`log1p(N_equiv)`、P同様`log1p(P_equiv)`を使う。UP positive indicator とDOWN positive indicatorを保持した基準にDOWNの強度を加えてからUPの強度を追加、逆の順番も監査。両方向/両成分を全結果報告。強度は**受水河川の実測TN/TPではない**。
- **科学的停止**：4セル・H6・国・原source整合未達、H6 5foldの列ランク不足、元比較不能なら**STOP**、都合よく源なしを清浄対照に戻さない。
- 4予測目的（植物種均等N、蝶種均等HNS植物調整、蝶種均等HNS植物非調整、植物記録重み）・2曝露方向・N/Pの多くの事後比較なので、一つの有意だけ抜き出して単独の「確認的発見」と呼ばない。係数95とH6外予測区間を区別する。十分な報告源あり／なしのsource同定は交絡の解消を保証しない。

## 言えること

UP負荷がDOWN負の対照より安定して予測するとしても、上下流の登録処理場N/Pが同一の受容植物区間に届くことを直接観測したわけではない。UPとDOWNが同程度・符号が不確実なら、元のWと植物・蝶の栄養嗜好の正の地域相関を`下水が栄養供給した`と読む根拠は一段弱まる。

河川の同一地点/同時点TN・TP、土壌/葉N・植物被度・幼虫寄主利用・蝶個体数、実際の水生羽化物を介した物質フローはどちらの方向結果でも未検証。元米国EPTの河床埋没度は独立の物理曝露であり、同一効果量にしてメタ解析しない。既存mainの正/負/STOP、原本成果物は改変せず、別Draft PRで新たな知見として監査する。
