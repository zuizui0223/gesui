# 2026-11-20｜水道会社最終処理水TP測定が元ChemPop底生調査に先行したか：原採水日SOURCE-only

**事後探索・SOURCE-only。** 昨日の元真の2017/18 EA water company final-treated SA TP観測と英国ChemPop v2のSOURCE適格BIO調査のネットワーク対応は、50km元NEXT_DOWN上流reach候補**148異なるBIO地点／237地点×年**（[source run #44681039804](https://github.com/zuizui0223/gesui/actions/runs/44681039804)、正式API250件/ページに修正済み）。原施設の物理放流口、処理水の河川への到達、P質量kg/year、処理場の栄養から底生生物BMWPへの因果効果を直接証明しない。

## Sourceデータ・元母集団と時間順序

- **BIO側：** 元2010–18適格7,836 EIDC ChemPop v2 sampleの真のBIO_SITE_ID×BIO_SAMPLE_DATE（source-only年月日、BMWPの応答値を開かない）。2017年389source-qualified BIO site、2018年366、両年BIO union463。BIOの同一site×yearで複数元sampleの採集日があれば、日付を重複なく全て保持。
- **処理水TP側：** 元EA新WQEの処理済み下水Water Company SA type、本物の測定0348 Phosphorus Total as P、元FINAL SEWAGE EFFLUENT、単位mg/L、非負数値のphenomenonTime元採水日。2017年281原SA、2018年293、両年重複275点、8,753原測定行。現在登録されているだけのSAはサンプル日と同値ではない。
- **空間・河川：** 元HydroRIVERS v1 Europe公式ZIP SHA500da7d36ceee0aa4c82dd625f5c53d80ba2f7ac14ad022cd3a6601bf55068e5。全国5,722 SAから原生物BIOと50km以内の方向・逆方向・same reach候補を先に網羅、その後真の2017/18 TP元採水歴で選択済み。元通常500m snap、同reach unknownを上流へ含めない。
- **本検定のSOURCE対象は元sourceで上流候補だった118 BIO site（2017年）/119（2018年）**。各BIO元採水日の前に、同年の上流候補SAに実際にTP分析があったか確認する。同年上流候補がないBIOは「処理場源0」にはならず、時間に関する今回の条件付き分母から除外する。

## 解析前に固定した3つの前方窓と1つの比較

1. BIO採集日以前**30日以内**の上流SA処理水TP元測定が一つ以上。
2. BIO採集日以前**90日以内**の同条件。
3. BIO採集日以前**365日以内**の同条件（2017年は2016年のEA源測定が未取得なので左側打切り。後から放流0扱いしない）。
4. BIO採集日**以降90日以内**のTP測定は、同年の測定が生物に「先行」していない可能性を可視化するために別表で数える。これは厳密な無処理陰性対照ではない（施設や調査の季節性・source打切りがある）。

**同一暦日のTP分析は時刻順が不明なので「先行」に含めない。** 同じ日付の観測ありを独立のsame_day判定にする。生物試料と化学採水の元時刻が両方で独立に確認できなければ、厳密なTP測定先行を0日差に拡張しない。

元BIO site×yearに複数の採集日がある場合、(a)一日でも条件を満たすsource-year、(b)全採集日が条件を満たすsource-year、(c)元distinct BIO site×dateで条件を満たす行数、の3つを別々に出す。2017/18同じBIO siteの年を独立河川に足し算せず、両年union distinct BIO siteも出す。WQEのtrue samplingPoint SA source採水日、原地理、BIO採集日は元idのみで対応させ、未来の測定を過去のsourceとして補完しない。

## 品質条件と解釈

元主source当初の148独立BIO・237site×year、2017/18 281/293実TP-SA、2026全国SA5,722、HydroRIVERS原SHAを再現できなければSOURCE STOP。WQE公式最大250観測/ページ、処理水採水点のnative IDをコンマで連結（最大100）し、5000件/回など大きな未検証ページ要求・繰り返しパラメータ形式は使用しない。原API HTTP502等は一定回数再試行し、全原pageが揃わなければSOURCE STOPとし、partial sourceを全国推定へ格上げしない。

この分析は**TP measurement coverageの先行性**であって、「実際にその日に河川がPを受けた」という確認ではない。TP分析値の濃度、下水量、施設permit-outfallの地理的同定、BIOのBMWP・豊富さ・生態系サービス応答値は扱わない。事後source分析の有意性を統計的確証と呼ばない。業者の原処理場施設数をnative samplingPoint数に代用しない。

元系列のsource-only調書とnegative/STOPをすべて保持し、mainやPR #135のBMI生物係数・欧州陸域結果へ変更を加えない。
