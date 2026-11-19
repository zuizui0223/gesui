# 2026-11-19｜Feurs処理場：水域・陸域二出口の年次固定N/P SOURCE-only監査

**重要：同じ施設から二つの出口があるが、同一年の水と陸の年間実測N/P質量はそろっていない。水陸二出口の質量収支・生物への因果効果・媒介はSOURCE-STOP。** 下表は公的な濃度・施用量の記録からの参考推算であり、実際の年間排出負荷や植物の栄養吸収ではない。

## 2018年：農地で実施済みの汚泥施用

| 入力・推算 | 数値 | 属性 |
|---|---:|---|
| 元STEU | 0442094S0002 | 同処理場、Loire川に排水 |
| 施用面積 | 161.53 ha | 実施済み農地施用面積 |
| 汚泥乾物（石灰込み） | 481.31 t | 実施済み固体施用 |
| 平均総窒素N | 4.76% dry mass | 分析報告平均 |
| 平均P₂O₅ | 2.96% dry mass | **P元素ではなく酸化物換算** |
| 固体N量の概算 | **22.910356 t N** | 481.31t × 4.76% |
| 固体P₂O₅量の概算 | **14.246776 t P₂O₅** | 481.31t × 2.96% |
| 元素P相当の概算 | **6.217675 t P** | P₂O₅ × 0.436426817 |
| 独立した平均N/ha×施用面積 | **23.583380 t N** | 146 kg N/ha × 161.53ha |
| 上記N概算2方式の差 | **+2.94%** | 面積式と組成式の相違 |
| 同年水出口の実測TN/PT積算質量 | **NA** | 取得できていない |

注意：(1) 平均試料含有率を全施用乾物へ一律適用したスクリーニング推算であり、ロット加重集計した実測のN/P投入総量ではない。(2) 汚泥と別途施用された液体351m³のN/Pを含まない。(3) 農学報告の「リン92 unités/ha」がP₂O₅なのか元素Pなのかは元定義未検証のため上式に使用しない。(4) 4.76%/2.96%の乾物分母と石灰補正、PDF原本SHAとページ画像の独立再取得は残課題。(5) これは処理場からの土壌への投入可能性であり、作物の実同化や蝶への移行ではない。

**公的原資料**：[Feurs 2018年農地汚泥実施用 bilan（DREAL）](https://www.auvergne-rhone-alpes.developpement-durable.gouv.fr/IMG/pdf/10_step_feurs_cas_par-cas_annexe9_epandage_2018bilan.pdf)、本文p.5,11–12,14,16–19。2018年の処理水放流先はLoire川。元PDFの画像・原本SHAはアクセス制限で未確保（テキスト索引原文による暫定確認）。

## 2019年・2021年：処理水の年量と放流平均濃度

| 年 | 年間処理水量 m³ | NTK平均 mgN/L | PT平均 mgP/L | NTK×年水量 t（参考値） | PT×年水量 t（参考値） |
|---|---:|---:|---:|---:|---:|
| 2019 | 1,081,393 | 5.15 | 0.45 | **5.569174** | **0.486627** |
| 2021 | 1,358,574 | 6.60 | 0.52 | **8.966588** | **0.706458** |

元の平均濃度に年間体積を掛けた参考値であり、**流量加重平均と検証できていないので年間実測負荷とはみなせない**。とくに**NTK（Kjeldahl窒素）はTNではない**。2019年と2021年の処理水量は入口側バイパス244,427m³／329,358m³を含まず、未処理水等の質量収支も閉じない。報告書の固形汚泥乾物502.68t／663.69tは年生産量で、**実際に施用した圃場・年月・N/P量ではない**。異なる2018年の実施用と併せて同年の水陸輸送収支を主張してはならない。

**公的原資料**：[Feurs 2019年下水道サービス報告（SISPEA）](https://www.services.eaufrance.fr/sispea/referential/download-rpqs.action?collectivityId=79907&rpqsId=299841)、[Feurs 2021年下水道サービス報告（SISPEA）](https://www.services.eaufrance.fr/sispea/referential/download-rpqs.action?collectivityId=79907&rpqsId=568453)。原本PDFのSHAは未取得。

## 再現・SOURCE gate

- 元入力：analysis/20261119_fr_sludge_routeC/FR_FEURS_PRIMARY_SOURCE_YEARS_2018_2019_2021.json
- 推算：analysis/20261119_fr_sludge_routeC/fr_feurs_dual_outlet_source_mass_gate.py
- 正しさ・反例：tests/test_fr_feurs_dual_outlet_2018_2019_2021.py
- GitHub Actions：.github/workflows/gesui-fr-feurs-2018to21-dual-outlet-source-gate.yml
- 厳密なペア実測質量が存在する年：**0**。実際の生物的曝露の同定：**0**。生物EQR/NUTRI/HNSを本監査で読まない。
- EEA 2018の両N/P年間実測負荷が正かつ一意淡水口の元解析標本には**FR由来が含まれていない**（別[PR #17](https://github.com/zuizui0223/gesui/pull/17)のQC表）。これはFeursの全公的排水データが存在しない証拠ではない。
- 次段は**Feurs 2018同じ排水口の流量加重TN/TPまたは同年実測積算質量**、汚泥ロット別N/P₂O₅、同年他施設の実施用台帳を新規取得すること。旧実施用自治体≥30・独立親H6≥10・実施用支持H7≥20のSOURCE閾値を下げず、未記載を「未施用=0」としない。

**判定：source bifurcation verified for Feurs 2018; quantitative paired source mass STOP; biological pathway and causality STOP.** mainと論文草稿PR #124は未変更。


## 2018年の真正な水域への地理接続（県令で固定）

[Loire県官報特別号2018-019、2018年3月8日公表](https://www.loire.gouv.fr/index.php/contenu/telechargement/5025/45942/file/raa_special_du_8_mars_2018.pdf)に掲載された**Feurs下水処理施設の操業許可 42-2018-03-01-007**のArticle 2（原印刷p.135）には、EPSG:2154（Lambert 93）の別々の地点が明記されている。

| 対象 | X Lambert 93 (m) | Y Lambert 93 (m) | 説明 |
|---|---:|---:|---|
| 処理場 | 794449 | 6516744 | 施設所在地 |
| **河川放流口** | **794216** | **6516746** | Loire川右岸の正式な放流口 |

同県令にある**5,800 m³/日**は設備の計画日平均処理量、**33,300 EH**は設計人口当量で、2018年の実際の水量でも年間N/P負荷でもない。これによりSOURCE段階の**正しい放流口位置**は取得できたが、水側の2018年実測TN・PT排出量の欠測は解消しない。PDF検索可能本文で座標・Articleを確認し、画像ページの取得はゲートウェイの502で失敗したため元PDFバイナリSHAは未固定。

2018年の未検証植物・蝶・底生EQR応答を用いた因果推定は依然としてSTOP。
