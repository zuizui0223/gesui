# 2026-11-16 実測報告下水N/P放流源―河畔衛星NDVI―水生BQEの厳密な同一地点接続・全数監査

状態：**2022年の実測報告N/Pを持つ原出水口群は取得成功、河道源地理ゲート成功、光学画素ゲートには全48で一定の支持があるが、同じ下水系の水域化学＋水生底生動物＋陸域衛星を結んだ介入生態解析にはGOできない**。これはソース欠測・結合の限界であり、栄養塩曝露が生態効果を持たないことを示す帰無結果ではない。指標の環境・因果区別を守る。

## A. 年間報告N/Pが実測値の出水口源を中心とする設計

- genuine EEA UWWTD official 2022 source direct one-freshwater-outlet reported positive yearly measured N & P：**1,156** (ドイツ630、オーストリア526)。本来の古いHydroWASTE16出水口に2022出水口sourceを上書きできない（旧DE/AT 8処理場では最寄りの年報告正のEEA出水口が1km以内0）。
- 原公開HydroRIVERS v1 EU第三者変換の元SHA固定、元実際のNEXT_DOWN・河道Polyline/支流・流下方向・放流口からの局所河道±1,500mと250m半径のbank候補を検証：**282/1,156**源出水口が技術地理PASS、AT145/DE137；*not* true same-site effluent delivered mass in river.
- 衛星結果閲覧前に地理的分散/別施設の距離≥5kmで**48** original yearly reported NP outfalls（AT24/DE24）を選定。[真正地理 Actions #43850677284 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43850677284)。後に元282施設すべてのtrue-source-geometry派生CSVを、source値も生物値も選別せず同じ元frameから再発行した：[Source full-denominator Actions #43856642187 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43856642187)。初期48組は不変。

## B. 同じ季節と元年を使用するSentinel-2観測の品質

- Element84 Earth Search C1 L2Aからは2018候補全48だが2022候補0。公開C1の年次元登録欠落であり、植物の無観測/自然変化なしとは解釈しない。
- original same **legacy Sentinel-2 L2A** 2018と2022の両年を使う代替sourceは同じ48施設で全施設3日以上の季節・シーン存在証明、赤/NIR/SCL元COG（元の2018 band scale0.0001/offset0, 元の2022 band scale0.0001/offset−0.1）。[STAC original source Actions #43851791936 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43851791936)。ESA PB04.00以降 L2A反射率は (DN-1000)/10000 に修正しなければ、偽の2018–22年NDVI差を発生させる。
- **事前固定した12施設（AT6/DE6）** の実COG raster 10m red/NIR,20m SCL最近傍再格子、SCL=4/5のclear landのみ・両岸各≥100実画素/月の主3か月（6–8月）において、両年共通2か月を満たす施設は**6/12**、AT3/DE3。固定閾値8/12失敗。元2018 monthly gates 12/11/11 pass；元2022 monthly gates 9/6/6 pass。55 original scene source画素QA pass、37 clouds/land quality fail。よって原12-smoke判定**STOP**で変更しない。[Original sensor QA #43853111036 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43853111036)。
- これは別の**事後source機会検証**として、先に選んでいたfull48へ同じ original SCL/PB/各月最大2scene/各月両岸≥100px/両年共通2か月の品質条件を適用。[Signed original full48 raster #43854317398 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43854317398)：**26/48**で元2018・2022の6–8月の2か月以上の同暦月実画素比較が可能（**AT14/DE12**）。2018 monthly pass=47/45/44 out of48, 2022 monthly pass=28/26/28 out of48。これは source availability/tech positiveであり**生態的NDVI応答や確認的効果推定ではない**。12-smoke元不合格は独立確認に変更しない。
- 2022 EEA 年間N/P質量は2018年の同施設N/P質量を提供せず、NDVIの2018→22年差が下水源曝露の2018→22年差を表すわけではない。SCL class 4/5だけでは本来の安定植生被覆・土地利用の識別をできず、上/下流250m円の土地利用は異なり得る。ESA WorldCover21等の厳密なland coverリファレンスと上流・同河道のoutfall-free対照、受水TN/TPへの流入が引き続き必要。

## C. SAME出水口に本当の水域EQRと実測TPの観測点があるか

元水生生物原本：EEA2022 Biology本物底生EQRステーション**2,475個**、2018/22では発表源・源サポートの変動があり、2022対象内で**AT68、DE0**（ドイツ自然界の底生生物ゼロではなく**このEEA2022 BQE公開source条件下でDE0**）。
元河川化学原本：ドイツ/オーストリアのWISE6 2022 genuine総リンTP測定日が十分あった**342ステーション**、公表WGS84座標の源対応**301**（DE214/AT87）、元国・完全site scheme識別。[source original PR #30](https://github.com/zuizui0223/gesui/pull/30)。**いずれもBQEの実数EQR・TPの濃度を開く前の純粋測定機会照合**。距離は水文流下方向を証明しない。

### 48 source地理先行パイロット

[Same 48 outfalls original aquatic source-only #43856384574 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43856384574)。

| original出水口から本来WGS84地点への距離 | 河川TP地点が少なくとも1つ近いsource数 | benthic EQR地点が少なくとも1つ近いsource数 | 両方が近いsource数 |
|---|---:|---:|---:|
| ≤500m | 1/48 | 0/48 | 0/48 |
| ≤1km | 4/48 | 1/48 | 1/48 |
| ≤2km | 5/48 | 1/48 | **1/48** |
| ≤5km | 5/48 | 2/48 | 1/48 |

2km以内にある唯一の水生BQE地点も、元放流口HYRIV_IDそのものまたは直下流NEXT_DOWNの同一reachには該当せず（0/48）。同一上流/下流処理場の水と陸が両方観測されたことを意味しない。

### 48例だけがたまたま少なかったのか：282源全体の独立したチェック

原282全施設にわたるsource-only地理再照合：[Actions #43857491196 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43857491196)。元 source full282は **AT145/DE137** の完全GISゲート通過出水口、決して衛星の高NDVI・低NDVIによって選ばれたものではない。282 outfallsには同じ流域を共有する施設もあるため、282を独立生態流域数と呼ばない。

| 距離 | 実際に河川TP地点が近いoriginal outfall数 | 底生BQE地点が近いoriginal outfall数 | 両方 |
|---|---:|---:|---:|
| ≤500m | 3/282 | 0/282 | 0/282 |
| ≤1km | 7/282 | 1/282 | 1/282 |
| ≤2km | **9/282** | **1/282** | **1/282** |
| ≤5km | 31/282 | 15/282 | 14/282 |

2km以内で本来HYRIV_ID出水口本流または直下流reachに該当するBQE地点も**0/282**。5km以内まで許すと同じ/直下流reach潜在BQEは4/282であるが、離れた下流の間の支流/別排水/土地利用を無視できず、実際の受水日付や川の河道方向にまで照合したのではない。

**判定：現行EEA/WISE6源だけでは処理場の個別2022N/P＋同year物理河道のTP＋水域底生生物＋河畔NDVIの共同source設計はNO-GO。** 上下流距離を調整するだけで真の経路に変えることはできず、原France 73H7の強いTP–I2M2 EQRはAT/DEの出水口と別の生物観測資料。

## D. 明確に検証した次の元公開source

既存EEA Biologyの国別年sourceカバレッジとは独立に、**オーストリアWISA/H2Oデータ照会**は政府自身が同一システムで①水質TN/TP等H2O化学、②ベントス(MZB/付着藻類)/魚/大型水草、③EmReg-OW水生測定局に集約された排出負荷、を横断照合できると案内している：
- https://www.bmluk.gv.at/themen/wasser/wisa/daten/daten.html
- https://wasser.umweltbundesamt.at/h2odb/

このデータベースは、下水を1施設単位で源地理に一致させた後に、水域環境と生物の同時・同一観測を取り出す道として**今のEEA2022 AT68/DE0より見込みがある**。ただし**実際に2022年のMZB+TP+排水質量が同一物理stationで何行取得できるかはまだ未照合**。統合仕様に含まれることと科学的なsource support gate通過は別で、既存の282 outfallsがすぐに実証可能とは言わない。
- ドイツのHessian HLNUG公式サイトでは**2004–2022年の河川MZBモニタリング結果XLSX**を明示的に配布：https://www.hlnug.de/themen/wasser/fliessgewaesser/fliessgewaesser-biologie/ueberwachungsergebnisse/fischnaehrtiere 。公開表に2022年調査地点識別子・河道座標があり、元EEA2022報告N/Pの同出水口真正下流に結びつくかは未確認。ドイツ全国の統一公開2022 BQE観測を無検証であると断言しない。

## E. 一本の論文への生態学的意義：中止基準と次の一次仮説

**原因主張にGOする条件は**本来2022出水口負荷・同年の正しい水質TP/TN・正しい同水域水生群集・同一水文学的に到達可能な河畔で測った地上の植生or衛星NDVIが、ある程度の独立施設で、処理場を挟まない参照河道とともにそろうこと。加えて地点/季節の違い・土地利用・河川流量/他の流入源を調整する。

一本化の真正中心問いはなお：
**When does wastewater-associated resource and stress redistribution propagate from rivers to terrestrial assemblages, and when is it interrupted by habitat connectivity and resource-access traits?**

本来の水域の本物TN/TP・底生EQRと河畔NDVIは、独立同一出水口での共通曝露が未照合なら、各生態区画の別々の研究であり、生物・環境効果量を一つのmeta効果としてプールしない。下水関連地理Wや衛星のNDVI source数だけで「下水が植物・蝶・生物多様性を増加させた」または「水陸伝播がない」を導かない。

元コードとアーティファクト：analysis/20261116_uwwtd_annual_riparian_sentinel、GitHub Draft PR #116、2022 original source site and S2 original 48 signed artifacts。第三者の全衛星画像/元EEA WQ全文を転載しない。main/chocho変更せず。
