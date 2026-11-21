# 2026-11-21｜英国2018年の正式な淡水放流地点を復元：EA最終処理水採水点との対応

**SOURCE-onlyの原本照合。源別P/N年間質量やBMWPへの効果は未推定。**

## 1. 原2018年英国UWWTDに放流点は存在する

旧strict「2018年のN/P両方を実測報告していて淡水放流口が一意な処理場」集合の英国0件は、英国の全原放流地点がないという意味ではなかった。独立に archived EEA UWWTD 2017–18 v07r01 を取得し、原3科学表の独立SHAを既存2018 sourceと完全照合。

- 元英国Reporter UK5（原対象2018年）のUWWTPS 1,799施設行、DischargePoints 1,809行。
- 淡水（元FWコード）放流地点 **1,498**、全部に妥当な座標が存在、重複native outlet IDなし。
- 元処理場1,787施設が放流点1つ、11施設が複数、1施設は放流点なし。
- ただし英国元UWWTPSで正の実測uwwNDischargeMeasuredは0、正の実測uwwPDischargeMeasuredも0。これは報告された実測値の欠測であり、排出自体がゼロと解釈しない。

**原本確認：** [EEA英国UWWTD原植物・放流地点監査 Actions #44848368611](https://github.com/zuizui0223/gesui/actions/runs/44848368611)。元UWWTPS/DischargePoints/ReportPeriodのuncompressed SHAは独立に固定済み。公式 archive 元source https://sdi.eea.europa.eu/datashare/s/mGmXg3zcEqKQLTG/download

## 2. 原FW放流地点とEAの処理済み水採水地点の地理対応

[放流地点対EA採水点の原座標照合 Actions #44848570506](https://github.com/zuizui0223/gesui/actions/runs/44848570506)、元集約JSON artifact #14684576041。元FW放流地点WGS84→BNGへ変換し、EA WQEの現在登録Water Company final treated SA 5,722監視点と比較。

| 直線距離 | 1+SA候補のFW放流点 | ちょうど1つのSA候補 | さらに元施設側も1放流点だけ |
|---|---:|---:|---:|
| 100 m | 747 | 688 | 683 |
| **250 m** | **1,032** | **938** | **932** |
| 500 m | 1,170 | 1,035 | 1,030 |
| 1,000 m | 1,241 | 1,063 | 1,057 |

ここで「ちょうど1 SA」は**放流地点側から**の条件であり、同じSAへ複数の原放流地点が近い可能性を消していない。また2026年現在登録SAは2018年に実際に測定されていたことを示さない。主候補は250m以内で**UWWTD放流地点→SAが唯一、SA→UWWTD放流地点も唯一、元処理場1放流点のみ**の相互一意sourceを事前固定。

## 3. 次に必要な実測年と流向

[source-protocol](PROTOCOL_20261121_UWWTD2018_OUTLET_TO_EA_TP_TO_BIO_DIRECTION.md)に従い、上の相互一意2018放流地点候補に対応するSAから**EA実測2018年TP code0348 mg/L**が存在した場合のみ同年候補とする。そして正式2018 UWWTD元放流点の位置（SA採水点の位置で代用しない）を公式HydroRIVERS v1の500m snap/50km NEXT_DOWNへ入れ、2018 ChemPopのoriginal source qualified 366 BIO siteへの上流／逆方向placebo／同reach不明のsource関係を集計する。BIOの元採水日前90日にEA TPの採水日があったかも別sourceとして検査する。

**厳密なsource限界：** このsource作業で得る相互一意の地理候補も、行政permit/native施設IDで証明された真の同一施設ではない。2018元UWWTDの物理放流地点の位置があるだけでは処理場がその時点でN/Pを何kg/y放出したかはわからず、EA原TP mg/Lを年量へ換算できない。河川源別TP到達やBMWP因果効果を今回の成果から直接主張しない。

なおUK PRTR 2017/18元都市下水施設のwater annual totalN/totalPはそれぞれ2017年124/0施設、2018年118/8施設で、両物質同一施設年の正値は0。別制度の施設年massをUWWTDの原outletへNative IDなしで割り振ることも禁止。結果・原監査はDraft PR #137に残し、PR #136/#135/mainの生物解析結果に変更なし。
