# 2026-11-15 固定仕様：欧州H7における植物NUTRIの特異性を、同じEIVの非N形質5項目で反証

**状態：Draft #102の肯定結果と否定的時間差分を見た後の事後対抗検査。前向き事前登録でも下水因果実験でもない。研究結果を見た後の有利な変数選択は禁止。**

## 中心疑問
同じH7流域で、植物種均等NUTRIを加えると蝶の種均等寄主HNSのH6留保予測が増えた（元451 H7period、90 H6）。だが湿度・温度・光・反応度・塩分などの植生環境の総合的違いでも起こりうる。**NUTRIが同じ指標セットの他の5環境適性指標より追加情報を持つか**を調べる。

## 元ソース
- 原451 eligible H7×period / 232 H7 / 90 genuine H6 / 11 majority country groups。元パネルSHA256 `1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca`。
- 真正 original 2,808 facet JSON SHA256 `251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4`。対象外14 H7periodを除き、元454の植物speciesKeyだけを使用。
- Tichý et al. 2023, JVS e13168 / 公式EIV European元 `EIV_Europe.csv` SHA256 `ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba`。元の4,011植物NUTRI exact matchesを再検証。GBIF hybrid source-v2監査で175 hybridは採点対象外、残る6は未解決。NUTRIを補完しない。
- 同じEIVの**6つをすべて先に固定**：light, temperature, moisture, soil reaction, nutrient preference (NUTRI), salinity。生態学的予測を見て都合のよい偽対照形質だけ選ばない。CSV生列名は取得ヘッダを原本SHAと共に検査して固定する。

## 原本での候補種・サイト同一化
- 原speciesKey→現行accepted vascular `SPECIES` の厳密二名法一致のみ。曖昧・未受理・雑種・数値形質なしを欠測とする。fuzzy/属・近縁種で埋めない。
- **主要共通形質母集団**：元11,511種の中で6指標すべて数値がある受理種。6形質とも**同じ種集合を用いて**H7×期間の種均等CWMを計算する。植物観測努力の違いではなく形質自身の追加予測を比較するため、指標ごとに異なる母集団を使わない。
- その共通形質母集団で各H7×期間10種以上かつ元全Plantae speciesKeyに対する対応率40%以上を満たす条件を先に決める。
- 主たる予測検定は≥100 H7×期間、≥40 H6、≥5 majority source-country groupsを確保できたときのみ。条件不達なら原本支持率の監査結果だけを保存し、共通母集団の比較を「STOP」とする。
- 元の451 NUTRI supported行だけを出発点とし、追加の共通6形質サポートにより減った行と偏りを国/H6/時間/記録努力で記述。共通母集団の選択バイアスを明記。

## 主解析（主蝶応答は元のequal-species incidence HNS）
1. ベースライン：元`B_GEOG` （W、都市・作物土地利用・H7地理・EEA/GBIF記録努力・蝶HNS被覆率）+ 6種共通coverageと得点種数 + country×method + period。元B_GEOG内のNUTRI coverage/dict scored species数は共通6種のcoverage/scored species数に**置き換える**。同じ母集団でモデルを比較。
2. 6単独追加：各EIV CWM（LIGHT/TEMP/MOIST/REACT/NUTRI/SALINITY）を一つずつ同一ベースラインに追加、H6クラスタWLSでβ/95とH6外MSE改善/固定誤差95を測定。各形質は同じサイト標本と同じデザイン、H6 5fold、Ridge alpha10、2000回固定損失のH6クラスタ再標本化。6つの名目pはBH6補正、選択した最大のものだけを報告しない。
3. **主反証**：NUTRI以外の5形質すべてを含むモデルにNUTRIを追加して、H6-out MSE改善とWLS β/95を測定。追加NUTRI 1SDの係数、他の5形質との多重共線性条件数を計測。VIF>10、条件数>100、または列ランク不足の場合、個別のNUTRI係数を選択固有の形質効果として解釈しない。Ridgeの追加予測が出ても因果や栄養独立を意味しない。
4. 対照の相対的強さ：同じ共通6種・同じサンプルで、各1形質を他の5つに追加したモデル（6通り）のH6-out改善をすべて報告。NUTRIのみが固有と仮定しない。
5. 元H6を訓練から完全除去するleave-one-country-out（6国程度）を探索的に追加、0を含む記述区間なら「国外転送未確認」。

## 補助解析（劣後）
- 蝶記録重みHNS、植物記録重みCWMは別解析で報告し、主要の負結果を上書きしない。
- 同じH7での時間差分はすでにΔPlantN→ΔButterfly HNS β+0.0039 [−0.0227,+0.0305]、H6留保MSE改善−0.70% [−4.16,+3.10]で正の同期なし。ここで同じデータを新しい正の時間媒介に再解釈しない。

## 成否と生態学的結論
- NUTRI追加がH6外で正／他5形質制御後も再現→「一般植生環境指標だけでは説明できないNUTRI特有の**記録組成上の予測情報**」。
- 非N形質の追加予測と同程度／NUTRI追加が不確実→「一般的な環境選別（共同生物地理や記録努力）との区別は不十分」。元のW→PlantN、PlantN↔HNSの正の関連は保存。
- いずれも下水現地N、葉N、蝶幼虫利用、真の個体数、花粉媒介、河川サービスやW→plant→butterfly因果を実証しない。古い源でEQR↔蝶HNS残差相関もゼロを含み、水→陸の連鎖は未実証。

**操作ルール：** 原sourceを上書きせず、新しいDraft内にsource-only QA→共通6指標対象→数値検査を積む。肯定・否定の双方を保存。main, chocho, PR #102は変更しない。
