# EPA2023 新規214流域の河床情報がなぜ移植できないか：粒径・流量・都市化のSource-only照合と厳密な反証

Date 2026-11-22 JST. **すべて事後的／探索的な観測解析。下水の曝露、因果機構、水→陸輸送の証明ではない。**

## 問いと元データ

2023–24年のEPA NRSA渡渉型（BERW/WADEABLE）1,165 UID・693 HUC8で、過去2008–09、2013–14、2018–19の3回の実測EPAから学習したEPT個体割合モデルは、現場の河床・河岸の物理指標6個を追加することで**2023全国HUC8均等の予測MSEを+7.057% [2.447,11.582]%改善**した。ただし元の旧期SOURCEに出現した479 HUC8（訓練からは全削除、2023 UID879）では+10.50%に対し、旧期のSOURCEに一度も出現しない214 HUC8（UID286）では−1.76%（区間は0を含む）。[母解析PR #141](https://github.com/zuizui0223/gesui/pull/141)、[旧期履歴の214ペア7環境マッチングPR #142](https://github.com/zuizui0223/gesui/pull/142)。

重要：**「新規HUC8」とは過去3期の解析に使ったEPA原本SOURCEに存在しなかったHUC8であり、米国の全調査履歴に新しい河川という意味ではない**。どちらの群もモデル学習から除外（old training 1,185 UIDs/551 HUC8、2023 target 1,165 UIDs/693 HUC8、target-train共通HUC8=0）。歴史的source membershipそのものが生態学的状態を変えるわけではない。

## 本質的なフィールド測定定義

[EPA NRSA2023/24 Field Operations Manual: Substrate Embeddedness](https://www.epa.gov/system/files/documents/2023-06/NRSA%202324_FOM_Wadeable%20V1.0_0.pdf) は、砂・シルト・粘土などを**定義上100%埋没**、岩盤を0%埋没と入力する規則を明示する。従ってXEMBEDは、河床の「既存の礫が砂に埋まっている現象」だけでなく、**そもそも砂や泥が占める割合**も含む複合指標であり、礫の独立した埋没現象を示す因果変数だと断定できない。

[原EPA2023 SiteInfo/PhysicalHabitat水質等Source-only schema gate #45021024628](https://github.com/zuizui0223/gesui/actions/runs/45021024628)で、新規に調べる候補は**2023原本1165 UID／693 HUC8の両群で成立**。元原表の変数名・測定単位・NARS flagsを保持した。

| Source-only original 2023 field | 1,165対象中source valid visits |
|---|---:|
| **PCT_FN** percent silt & finer | **1,165** |
| **PCT_SAFN** percent sand & finer | **1,165** |
| **LSUB_DMM** original log geometric mean bed particle size | **1,165** |
| EPA **MEAN_FLOW** NHDPlusHR modeled mean annual stream flow cfs | **1,157** |
| **WSAREASQKM** upstream catchment area | **1,165** |
| **URBN_NRS23** urban category | **1,165** |
| **XDEPTH** stream depth | **1,164** |
| **W1_HAG** riparian human/agriculture disturbance | **1,165** |
| **LRBS** original source relative bed stability | **not present in this specific EPA PHab metrics CSV** |

全693 HUC8に各流域の原本測定平均をとって、純環境統計だけを計算：
- XEMBED vs PCT_SAFN **r = +0.93965**
- XEMBED vs LSUB_DMM **r = −0.90009**
- XEMBED vs PCT_FN **r = +0.58985**
- PCT_FN+PCT_SAFN+LSUB_DMMによる**同一標本のOLS記述R² = 0.90278**。この値は同一流域内の組成関連を示すもので、未知流域予測や粒径の直接の作用の証拠ではない。

## 追加のSOURCE ONLYマッチング契約

[環境ブロック・追加sourceプロトコル](https://github.com/zuizui0223/gesui/blob/analysis/nrsa2023-novel-huc-substrate-domain-20261122/docs/PROTOCOL_20261122_NRSA2023_NOVEL_HUC_SUBSTRATE_HYDROLOGY_FALSIFIER.md)。以前から確定した元7変数の**214ペア/HUC8 exact SHA**を旧3期原本で再現し、底生数値を一切使わずマッチング用追加共変量を決定。

- **grain block**：従来の緯度経度、河川次数、実測河川logTN/logTP、XEMBED、log河幅＋原本 **PCT_FN, PCT_SAFN, LSUB_DMM**。新規214/旧SOURCE479のHUC8で欠損なし、**214対**、9 eco9すべてに対応。ペアSHA **`8c78c9060e0d512666062760f0467baaa54f6cfabeeaf76867ad609e3f1ea8e2`**。マッチ後10変数の標準化平均差は全てabs<0.10。
- **hydrology/landuse block**：元7変数＋**log10(MEAN_FLOW cfs), log10(WSAREASQKM km²), HUC8 urban割合**。旧SOURCE476 HUC8、新規213 HUC8の追加source成立、caliper3.5に固定した1対1 matchingで**209対**、9 eco9に対応。ペアSHA **`4d21361aa86707039d4afd62cfa11c70f629c50ea01e1a0663e4c00930f254b2`**。10変数のマッチ後abs SMDは全て<0.10。
- モデルに関係する**旧期トレーニング551HUCの分布**だけで元7環境変数を正規化。追加変数は旧3期統合データに存在しないので、**2023全693HUCのsource-only分布**で正規化（生物応答不使用）。経時的な値の互換性を根拠なく仮定しない。共通caliper3.5、新規source coverage≥60%（129流域）、70ペア、5 eco9、両群それぞれ80%新変数coverageという閾値は最初に固定。
- [追加source-only original EPA全old-era再構成 #45021173233 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45021173233)、aggregate source-only artifact **14753442780**。source cohort and paired SHA固定を、生物EPTの新しい誤差を開く前に完了。

## 実際の過去3期モデルを変えず、2023 EPTに対する214/209対の外部予測誤差を再集計

結果[実際のEPA原本全4調査期の再構築＋source-frozen後のモデル評価 #45021358057 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45021358057)、signed aggregate outcome artifact **14754236252**。**新規2023年のEPTは旧期学習・尺度・モデル選択に一切使っていない**。モデルは過去3調査期学習1,185 UID/551HUC8の元の同じ重み付きRidge10、基準B+実測河川N/Pと、source追加6PHab・XEMBEDのみ、他5PHabを調整後のXEMBEDを比較。対象HUC8均等に二乗誤差を集約。2,000組whole-pair fixed-prediction bootstrapを記述95%として評価。

| 元の2023年EPT予測における追加情報 | old-source出現HUCの予測改善 | new-source HUCの予測改善 | new−old差（paired 95%） |
|---|---:|---:|---:|
| **粒径組成追加調整・214対** ||||
| **XEMBED after N/P** | **+9.955%** | −1.778% | **−11.733 pp [−20.108,−3.202]** |
| XEMBED after N/P+他5PHab | **+8.624%** | +0.058% | **−8.566 pp [−16.146,−1.113]** |
| 6 PHab after N/P | +8.429% | −1.758% | −10.186 pp [−21.080,+1.459] |
| **推定流量・集水域・都市化追加調整・209対** ||||
| **XEMBED after N/P** | **+11.025%** | −1.489% | **−12.514 pp [−20.048,−4.493]** |
| XEMBED after N/P+他5PHab | **+8.602%** | +0.296% | **−8.306 pp [−14.571,−1.641]** |
| 6 PHab after N/P | **+10.863%** | −1.341% | **−12.204 pp [−22.997,−0.831]** |

数値は同一対象のbaselineとfull modelに対するHUC8均等**MSE百分率改善**。2種類のマッチングはほぼ同一の新規214流域を共有するので**独立したフィールド再現とはみなさない**。計算した計8件のposthoc comparisonsに多重性補正は掛けていない。**「既知群は有意、新規群は非有意」だけをもって2群の差を宣言していない**。固定HUC8ペアを丸ごとresampleした差そのものの記述区間を報告する。

## 独立モデル・統計数値検証

[全EPA原本→4調査期旧期再構成→213/214マッチSource SHA→別実装NumPy係数とペアbootstrapの独立検証 #45021574096 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45021574096)、signed original QA artifact **14753458308**：
- 新粒径214組、河川水理209組のsource-onlyマッチSHA厳密一致。
- 元のscikit-learnとは異なる**加重NumPy正規方程式 Ridge10**で過去3期のみから2023 EPTを独立予測。モデルごとにHUC8内損失と2,000ペアbootstrapの95%区間を再計算。
- 2ブロック各4比較＝8セットの**全paired loss contrasts/差/CIを1e−10以下の許容誤差で完全再現**。
- 独立コードは原本の再計算監査であり、追加の独立生態系現地観測を意味しない。原本raw site-level XEMBED/EPT/coordinateデータをリポにコミットしていない。

## 生態学的な意味と、言えないこと

1. **XEMBED指標は河床粒径組成の変動の代理にもなっている**。特にsand+fines割合だけで強い相関 r=.94、原本3粒径の記述R²=.903。このため「独立した堆積ストレスの機構が証明された」とは言えない。
2. しかし**粒径組成・流況規模・都市化の平均的環境差を新しく調整しても**、以前のsource記録のある流域と未記録の流域とのEPTの**予測可能性**の差が消えない。新規流域で河床情報の予測移植が無条件に成り立つという普遍性には強い留保が必要。
3. **観測履歴が生態応答の原因であるという因果推論は不可能**。そもそも過去3期sourceはEPA元の全調査ではなく、物理・水質・底生・時に付着藻類の測定が整ったsource-population。観測履歴ラベルは、未測定の河川流況、河床の質、種プール、地質、サンプル検出、過去の選択性などを代替している可能性がある。
4. **相対EPT割合のみ**。EPT order composition/family-specific反応、実際の種多様性、羽化乾燥質量、土砂輸送量、日照、水温、実流量は今回の損失表には含まない。次の独立仮説は2023年の同一sourceでEPT order/family構成を分類群ブラインドで固定し、XEMBEDと調査方式・群集組成の予測ギャップを判定すること。
5. 同一EPA調査プログラムの4時期という外挿であって、別大陸／別モニタリング制度への普遍性は証明できない。既知の河床粒径と底生動物の生態学に対する新規性は**環境指標のsource-domain transferabilityの限界を本物のsource-only matchingと独立誤差検証で定量化**した部分であり、新発見の河床機構ではない。
6. **下水処理場からのP量を直接再構成した分析ではない**。[gesui PR #139](https://github.com/zuizui0223/gesui/pull/139)のWWTP TP×樹冠source GO失敗とは異なる。水域→陸域の動物、河畔の植物・蝶、糞や汚泥、実際の栄養輸送、サービス連鎖とは繋いで因果図にしない。

## 投稿進捗・候補
探索的source-domain解析の統合としては **「Sediment-composition proxies transfer across monitoring eras but not uniformly across historically unsampled catchment domains」** を仮候補とできる。ただし「粒径に基づく限定」まで確証したことにはならず、**新規HUC群での予測欠如はHUCソース選択により条件付けられた記述結果**である。今すぐ総合環境学の普遍性論文にするには、独立したモニタリング源または群集組成・source sampling designの明確な反証がさらに必要。ここで新しい共変量を効果が出るまで追加する探索は停止するのが適切。
