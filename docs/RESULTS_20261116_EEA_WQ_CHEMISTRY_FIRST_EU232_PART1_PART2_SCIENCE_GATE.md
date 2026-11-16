# 2026-11-16 公式EEA実測TN/TP起点の水陸連結：原本源・標本支持判定

**科学的状態：公式原本取得・物理監測局×年の源照合PASS、植物・蝶の生態回帰は厳格 science gate STOP。** 無検出・不在・負の効果ではない。元の水域・陸域肯定結果は維持し、主リポ main/chocho は変更していない。Draft PR #111 (stacked on original source-only PR #107)。

## 1. なぜ旧454 H7×期間からの測定局一致をやめたか

PR #107 は、2018/22同一EEA底生EQR観測局と同一年のPart1+Part2厳密TN・TPがそろい、かつ上流処理場の正の年間N/P報告も確認できたものが**9 H7×期間/8本来H7/8真正H6**しかないと原本で確認した。近傍1kmに別WQ局を探しても、WWTP報告正で最大16 H7×期間、閾値30に数学的に達しない。

本研究は、固定した旧232 HydroBASINS H7（GBIF植物と蝶の真正454 H7×期間）内にある公式**全水質WQ局から始める**別の母集団に変更。EEA底生EQR物理stationの同一性、上流WWTP施設の報告有無は標本選択に使わない。これは水域EQR測定局の水質を代替するのではなく、**流域内の別のambient水質**を新たに問うのであり、元のPR107不通過を通過扱いにはしない。

## 2. official SpatialObjectsと原HydroBASINSの幾何 source gate

[原本source Actions #43816554985 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43816554985)。

- EEA WQ SpatialObjects 原本ZIP SHA256 aa505be6f85d06c675e0d0fcb148a9148f80162f8e9590740a887cc95a3e349d。90,859 original source行、69,820 source有効WGS84座標の完全な国・ID・schemeを持つ監測局。位置不明/0,0や位置の異なる同一IDを除外。
- HydroBASINS Europe Level7公式ZIP original SHA256 b9b810f4db271f0fa6b39fecf01784041b948944b3c3f296a393766388e86b44。旧232 H7全多角形有効。WQ局が厳密に1つのoriginal H7を含むときのみ採用。
- **12,139 official WQ station IDs** が元232 H7のいずれかにあり、2018/22×元454 H7×期間と対応する地理機会は23,717 station×year。元GBIF basin majority-countryとの完全一致を主にすると **11,540 WQ局・21,711 station×year候補**。原source候補SHA256 c5cde7bd8134b254c5e95b06849292fa260fd6f01832e6b6d2e1994a903c5910。
- 地理を源・化学値・生物の傾きから独立に先に固定。国不一致の地点を「生態効果が弱かったから」除外したのではなく事前定義のprimary。何れも同じH7 polygonがその生物観察点への実際の栄養輸送経路だとは保証しない。

## 3. 公式EEA Water Quality ICM 2026 Part2（年次原報告のみ）

[原本Part2 Actions #43816554985 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43816554985)。公開 ZIP 170,422,309 bytes、SHA256 315396800d76e46db415164851fb4f99765b7ec2654828fc42782d902e679a7a、原年次CSV 6,450,192 source行。

RW河川、水matrix W、TN EEA_31615-01-7 / TP CAS_7723-14-0、元元素mg/L単位、3つ以上の元試料・LOQ品質を要求し、同一WQ物理局×年で**両方の原年平均値が１行ずつ有効**ときだけ採点。

| 公式源 | 原TNTP同一局×年 | 主country一致 | 主H7×期間 | 主distinct H7 | 主真正H6 | 主国 |
|---|---:|---:|---:|---:|---:|---:|
| Part2 年次報告 | 54 | 32 | **5** | 3 | 2 | 1 |

結果：**STOP_SOURCE_CHEMISTRY_PANEL_SUPPORT**（事前30 H7period/20 H7/20 H6/3国）。Part2では全国のWQが測定していない、という意味ではない。Part1の個別観測も含めた補完前の源情報しか示さない。

## 4. 公式EEA WQ ICM 2026 Part1（別源の個別採水原本）により実質的拡張

既存Part2結果を見たあと、Part1用の別のsource-only仕様を文書で先に固定。**旧のPart2 STOPを改変せず**同一original geographic candidateを使用。公式原ZIPは1,639,231,716B、SHA256 ff779521749b40b952cce35a4671240c3b6c47d87f32a8fef9240a4ac503a082、DEFLATE64を7z stdoutにストリーム展開し、原個別試料**96,597,294行を全走査**した（[成功Actions #43817344055](https://github.com/zuizui0223/gesui/actions/runs/43817344055)）。

- sourceのgenuine original sampling-date 2018/2022、river RW/water W、厳密同一物理局・国・site identifier scheme、同じTN/TP決定因子の元値で元素単位を照合。各TNとTPは少なくとも3件の異なる妥当試料、重複試料とLOQ未満を除去して年平均を求めた。Part2年次と混ぜず**Part1単独**で照合、実測された原試料平均の分位値を見て地理選択は変えない。

| 公式源 | 原TNTP同一局×年 | 主country一致 | 主H7×期間 | 主distinct H7 | 主真正H6 | 主国 |
|---|---:|---:|---:|---:|---:|---:|
| **Part1 個別原試料** | **632** | **621** | **37** | **30** | **24** | **5** |

**PASS_PART1_SOURCE_SUPPORT_ONLY_NO_ECOLOGICAL_CLAIM**：当初固定した「源の可能性」最低閾値30 H7×period /20 H7 /20 H6 /3か国をすべて満たした。ただしPart1が少ない国もあるため、5か国を欧州全域の代表と一般化できない。生物効果が正/負かはここでは一切読んでいない。

## 5. それでも生態モデルの独立したサンプルゲートはSTOP

source gateより厳しく **生態主モデルを開始する前に**最低100 H7×期間、40独立H6、5か国、原EIV numeric 10種+coverage40%を事前固定。元GBIF454 H7パネル（sha256 1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca）と真正EIV4,011原受理種の公式trait v2に同じH7×期間のWQ Part1源をjoinし、元451trait-eligible sourceで確認する。

[科学source gate Actions #43818637326 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43818637326)：
- 本来のEIV支持植物・蝶元H7パネルへの同一H7×期間対応 **37行、30異なるH7、24 genuine H6、5か国**。
- source geographyとtrait coverageは一致し欠落追加0。しかし**37/100 H7×期間、24/40 H6**で事前基準未達。
- 結果コード **STOP_BIOLOGICAL_SAMPLE_SUPPORT**：事前の水質TN/TP→植物EIV NUTRI→蝶HNS解析で新しいβ・p・地域外予測改善を**１つも推定していない**。小さい標本で符号が正のモデルだけ探して選ばない。

**数学的な源上限**：Part1とPart2をどのように真の元物理station×year×nutrientで重複除去してunionしても、元のPart1単独37 H7period/24 H6とPart2単独5 H7period/2 H6から、最大でも **42 H7×期間、26 H6** を超えない（実際は重複で少なくなる）。したがって原232 H7の両源unionでも、事前の生態モデル100 H7×期間/40 H6は満たせない。非必要なPart2の170MB再取得と、結果を見た閾値の緩和は行わない。

## 6. 水域の肯定結果と陸域の肯定結果をつなぐ正しい解釈

- 既に支持：欧州上流処理場年報告N/Pが水生底生EQRをわずかに予測する2018/2022結果、米国河床埋没XEMBEDとEPT固着形質・耐性の3期反復、欧州の植物栄養嗜好NUTRIと蝶寄主HNSの同じ流域での空間的対応。
- 不確実：水域EQR―植物NUTRIは同一H7の負の条件付き空間関連（β−.0345）を示すがH6留保・時間差は不確実。実測TN/TPに対し原232 H7で十分な標本を確保できないため、**同じ実測水質の因果的連鎖**は依然検証されていない。
- 新しい肯定は「真正EEA水質原本と正しい流域地理の突合によりPart1で**37 H7periodの実測TNTP×元GBIF trait opportunity**を再現可能に確保した」技術・サンプル構築の成功であり、水質が植物・蝶を増やしたという生態効果ではない。
- 本来の次段は、water-quality-firstに元の232 H7枠外の水質測定局をさらに採用し、それらのH7で5蝶科＋植物数値EIVの真正GBIF source facetを新規取得し、まず生物値を開かず国・H6・trait coverageの源成立性だけを判定すること。あるいは測定局近接の河岸植物・幼虫記録の本当の局所比較が必要。旧EQR源454流域へ後付けで絞って「下水効果」と呼ぶのは避ける。

再現コード/原本source：analysis/20261116_chemistry_first_eu232 (station geography, Part2, Part1, ecology_if_source_pass)。各原大容量EEA ZIP/個別sourceレコードはリポジトリに転載せず、GitHub Actionsに有限の派生集計とSHAのみ保存。既存main/chochoへの自動マージなし。原標本の欠測は非排水・零濃度とはしない。
