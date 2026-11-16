# 下水―水域―河畔生物を一本の論文にするための真正source監査（2026-11-16）

**研究区分：原本検査済みのsource feasibility、効果量の推定なし。** 空間・時間・生物群の統一は、同じ流域IDを結合するだけでは正当化できない。今回の判定は「現行公開sourceで下水の水陸連鎖を因果的に検証可能か」。以前の水域・陸域の肯定・否定結果を変更しない。別Draft PR #115、main/chocho 未変更。

## 1. 既に支持される解析と、一本化に足りない観測

- 欧州施設2018/22年報告N/P上流負荷→本物の底生EQRの予測改善は約1%だが、方向に弱い推論・報告欠損/国の偏りがあり、生物多様性向上や局所曝露の因果効果ではない。
- フランス73 H7/24真正H6の実際の同じ河道TPとI2M2底生EQRは、真H6外の予測を+24.6%改善。一方、実測TP→同じH7植物NUTRI/蝶HNSの追加予測改善は不確実（負値）。植物NUTRI→蝶HNSは実測TP調整後+30.2%（正）、ただし両指標は共通の寄主生育地栄養嗜好情報と属構成を共有。水陸の因果栄養移動・現地幼虫摂食ではない。原本PR #114。
- 米国EPT clinger−/汚濁耐性PTV+は河床物理勾配XEMBEDの結果であり、欧州の下水の負荷と同一曝露ではない。
- つまり「下水が水域を変えて、その結果植物・蝶が変わる」を現在のsourceから計算したSEMにするのは不正確。

## 2. 今回のsource-onlyの新規検査

### A. 16既存処理場×48旧GBIF植生（2011–24）からの形質取得可能性

旧元PR #3/#4による処理場独立16（AUT/CZE/DEU/DNK各4）、U2/U1/D1の48サイト、公式原GBIF植物speciesKey×地点 **5,394 source行、1,402 source speciesKey**。EIV公式11,511元植物種のうち4,011厳格受理維管束植物に数値NUTRIがある。

真のsource original EIV scorer availabilityだけで各旧GBIF siteに >=10得点可能植物種、>=40% scored species, >=20 scoreable records要求：**44/48 site、13/16 outfallがU2 U1 D1全て通過（4国）**。[Actions #43845610859 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43845610859)。**ただし**元中間河道centre U1→D1は真の水流路で min4.096、median8.076、max16.479km、**2km以内0/16、5km以内2/16**、元2011–24の長期資料に対応する原WWTP年別実測N報告元項目は**0/16 source施設**（missingは無放流を意味しない）。単に植物EIVを再計算しても、局所下水影響を推定できない。

### B. 従来の河道中心点を捨てて、真正放流口±1.5kmの地点を地理的に作り直した

古い原16施設のsource真正放流口位置とHydroRIVERS U2→U1→outfall→D1四reachのNEXT_DOWN/line端点・snapを独立検証。真の河道polylineに沿って、元放流口から上流1,500mと下流1,500mの**新しい2地点**を構築した。両地点の次段GBIF/field radius=250m、源Outfall snap≤100m・河道接合≤100m・真のflow順序・250m buffer非重複をすべて検査。**16/16施設で地理的source gate PASS（4国）**。[Actions #43846527119 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43846527119)。

元source HydroRIVERSはSHA固定されたサードパーティ変換で、公式原HYRIV v1のbit単位topology独立確認ではない。座標生成の成功は同じ生態観察の成功でも年次化学N/P受水の検証でもない。

### C. 新しい局所250m径で、同じ2017–19と2020–22の真正GBIF録機会があるか

元32新地点×原2期間（2017–19,2020–22）×2群(Plantae, 非焦点Animalia)の**128 API source count**が全て成功、先行研究のGBIF record人数や生物個体数でなく公開original georeferenced occurrence totalのみ確認。標本gate=同じ原Outfall上下流に各20以上original plant GBIF観測があり、>=10独立WWTP/3国。

結果：**両期間とも0/16源処理場が両局で20件以上**。旧元2011–24 1km径/48河道centreの記録をそのまま短期250m径/真の放流口へ当てはめられない。[Actions #43847165094 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43847165094)。

### D. 独立して凍結した半径500m/1kmの再試験（250m gateは取り下げない）

流路上でUP/DOWN中心は各1.5km離し、原二地点の水平距離は**最小2.105km**。250m/500m/1kmの円バッファは原16全地点で重複しない。元3年窓・放流地点・GBIF filterは変更せず、半径500m/1kmだけを別途source機会として128回のoriginal official GBIF Plantae source countで比較。[Actions #43847809172 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43847809172)。

| 上流・下流両方、1地点20件以上original GBIF植物recordを得た処理場 | 2017–2019 | 2020–2022 |
|---|---:|---:|
| **250m primary（以前に凍結）** | **0/16** | **0/16** |
| 500m posthoc source sensitivity | 2/16（2国） | 3/16（2国） |
| 1km posthoc source sensitivity | 4/16（3国） | **10/16（3国）** |

**意味：** 年・物理距離を固定すると source植物観察が消え、GBIFを1km径に広げれば2020–22年について一部が成立するが、2017–19年は4/16しか残らない。両期間を一組の独立施設として>=10組確保する縦断解析は不成立。0件は個体群絶滅や植物種数=0ではなく観測機会不足。

これは「生態系の真の水陸境界」の実証ではなく、**観測位置・時間窓を変えるとデータ支持標本自体が大きく変わる**という元sourceの限界。1kmの植生は水際250mの河畔特異生物群集ではなく他の土地利用/庭などの記録を含み得る。2020–22の10 outletにEIV数値を計算しても、実際の曝露N/P源は未対応のまま。

## 3. 下水に注目する正当な理由と先行研究ギャップ

既出例：
- Behrens & Bernhardt 2025 Freshwater Science DOI 10.1086/735929：都市河川の放流影響域で羽化昆虫から河岸へ移動する生物量と金属量が増える例。水陸の資源補助と汚染物質の両面はすでに実証されている。
- Diesbourg et al. 2025 Environmental Pollution DOI 10.1016/j.envpol.2025.125948：下水放流水由来資源が水陸食物網へ移ると同時に昆虫群の微生物応答が分類群依存。
- Wardlaw et al. 2026 Environmental Pollution DOI 10.1016/j.envpol.2026.129116：8処理場上下流の幼虫・羽化昆虫・河畔クモにおける微細プラスチック移動も既出、明確な普遍的濃縮は支持されない。
- 2025 Ecological Monographs DOI 10.1002/ecm.70025：4欧州河川で水生昆虫の分散形質が河岸クモの食物網連結に関連。
- de Guzman et al. 2026 Hydrobiologia DOI 10.1007/s10750-025-06103-w：河川水系への排水追加が緑色食物網の資源利用を刺激しても食物連鎖長は変わらない。

したがって新規性を「下水が水陸をつなぐ初めての発見」に置かない。**下水という位置・日時・負荷を特定可能な人為流入源の資源補助と汚染ストレスが、どの機能的資源アクセス経路・生態系間連結性の条件で生物の反応へ変換されるか**が未解明な一般問題。

## 4. 一本化論文の本当のGO条件

一本の明確な一次仮説：
**Wastewater-associated changes propagate across aquatic–terrestrial communities only along realized resource-access pathways, while habitat filtering and contamination can interrupt or redirect that propagation.**

十分な観測があるなら、原Outfallごとに：
1. SOURCE：同年の真の処理場排水口と年次実測N/P/流量/処理方式、上流対照・下流環境の真のreach流下方向、別の未曝露参照河川を固定。下流TPは農業源と混合し得るので処理場由来割合を区別。
2. 水域：同じ河道・同じ季節の栄養N/P、水温、藻類、EPT形質と計数・実測羽化（可能なら生物量）。
3. 陸域：水文学的に連結する同じ河畔の**定量的植生調査区**または羽化昆虫を捕食するクモの栄養源/汚染測定。蝶を採るなら寄主植物/幼虫利用と成虫吸水の**別経路**を計測。出現記録のHNSだけでは真の摂食や送粉サービスを示せない。
4. 分析：各処理場内の真の upstream–downstream 対比を、その地域の下水を挟まない同方向の参照河道差と比較する階層BACI/傾向差設計。見たいのは「source nutrient subsidy × functional trait or access-to-riparian-channel」の交互作用、真の観測努力/年/距離/流路の調整。実際に複数年度ある場合のみ時間遅れを検討。
5. 生態学的効果推定をする前に、独立WWTPの数・国・対照河川・実測時期・river reach・河畔調査区/バイオマスのsource gateを凍結する。少数の原16 outletにGBIF 14年間の広域録を貼り付けて、因果的な"下水の特徴的な生態系"と見せない。

**現状の結論はGO（共通の真河道地点設計が可能）/STOP（同期間の局所GBIF植生・真の年別処理場N/P曝露との同時リンクが未成立）**。この結論は「下水の生態作用はゼロ」という帰無支持ではない。選択肢は、河畔植生・BMSの標準化現場データを本来の局所スケールで取得するか、元の水域・陸域の確実な結果を別々の制限された観測論文にするか。結果が出るまで水陸の媒介SEMを強制しない。

## 5. 原本再現と保存範囲

- EIV旧48サイトのscoreability source Actions #43845610859
- 真の放流口±1.5km物理河道再配置 Actions #43846527119
- 250m固定year site records source Actions #43847165094
- 500m/1km radius sensitivity Actions #43847809172
- 原本文書・source script：analysis/20261116_outfall_riparian_local_source（原正真outfall selected16およびsource 4011 EIV sourceをartifactから復元）。第三者原HydroRIVERSやGBIF観察者の元個人ID、巨大CSVをリポジトリへ再配布していない。source結果は全て探査的、レビュー前。既存 main/chocho 変更なし、Draft PR #115。
