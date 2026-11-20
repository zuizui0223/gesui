# 2026-11-20 英国の実測最終処理水TPと底生ChemPopの有向HydroRIVERS source品質監査

**SOURCE ONLY／事後探索。正式な物理放流口、下水由来N/Pの下流到達量、群集への因果効果は実証していない。**

## 固定原資料と監査

- [EA公式WQEの2017・2018年の実測総リンSA源の検証 Actions #44678735822](https://github.com/zuizui0223/gesui/actions/runs/44678735822)：現在のWater Company最終処理水採水点5,722点のうち、2017・2018年のSOURCE適格ChemPop BIO site和集合463から2kmの元BNG直線距離以内にある367現在登録SA採水点を歴史原資料で問合せ。2017年は67、2018年は71元SA採水点に真正0348総リン (TP as P) mg/L原観測あり、67採水点が両年に記録。元生物SITE ID/日付は2017年389、2018年366元SOURCE適格site。事前の直線距離2kmのみでは2017年50、2018年69BIO siteが候補だったが、源から下流とは断定しない。
- [公式HydroSHEDS HydroRIVERS v1欧州原ソース Actions #44679393760](https://github.com/zuizui0223/gesui/actions/runs/44679393760)：原ZIP 67,648,957 bytes、**938,544河道**、SHA256 500da7d36ceee0aa4c82dd625f5c53d80ba2f7ac14ad022cd3a6601bf55068e5。HYRIV_ID、NEXT_DOWN、LENGTH_KM、MAIN_RIV、HYBAS_L12フィールドを確認。英国を含むWGS84図郭lon−9〜4、lat49〜58で**28,830河道**を抽出。元EA SAと元ChemPop BIOのBNG座標（EPSG27700）で500m以下の最近接到達を確認。最近接異なるreachの同距離±1mはAMBIGUOUSとして除外、同じ河道reachesで前後位置不明の組合せも上下流比較へ入れない。
- [今回の有向NEXT_DOWNによる原SOURCE-only照合 Actions #44679903035 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/44679903035)、成果artifact ID **14604061338**。正式河道に沿いSOURCE reachからBIO reachまで最大50kmの方向性を数える。同一河道の不明な前後は除き、BIO→SAの逆向きのPLACEBOも独立に計上。生物BMWPの元数値は一切モデルに入れない。

## 位置のソースカバレッジ

2017/18の元ChemPop source適格BIO和集合463のうち **442 site**が公式河道へ500m以内で一意にスナップ、21 siteは非snap（小河川・都市水路など元HydroRIVERS未収録もありうる）。真正TP実測SA採水点和集合71のうち **66 point**が一意にriver snap、5 pointは非snap。最近接分岐tie未決定は0。スナップ距離中央値はBIO約109m、SA約117m。非snap地点を「未汚染／STW排出ゼロ」扱いしない。

| 同一年度の元SOURCE・river network比較 | 2017 | 2018 |
|---|---:|---:|
| ChemPop SOURCE適格BIO site | 389 | 366 |
| 元実測TPのあるWater Company SA point | 67 | 71 |
| BIO site一意river snap成功 | 371 | 354 |
| STW SA point一意river snap成功 | 62 | 66 |
| **SA河道→50km下流のBIO河道に到達した異なるBIO site** | **36** | **43** |
| 上流候補になった異なるSA採水点 | 31 | 35 |
| **逆向きBIO→下流SAの対照があるBIO site** | **55** | **60** |
| 同一河道で位置順序不明のBIO site | 32 | 45 |
| SA→BIOの異なるreach有向ペア | 50 | 58 |
| BIO→SAの異なるreach逆向きペア | 57 | 65 |
| 上流候補の河道距離中央値 | 13.9km | 16.05km |
| 上流候補BIOが含まれる原HYBAS_L12の数 | 30 | 33 |
| 1 BIO siteに2以上の上流SA原採水点が対応 | 7 | 8 |

**2年間のユニークな正方向候補は51 BIO site／延べ79 BIO site×year。** 元の SOURCE支持条件「50異なるBIO site」に届くが、「100 site×year」には届かない。HYBAS_L12は大きな同一本流水系内で複数になるため、5独立水文blockの証明としてそのまま用いない。

**逆方向PLACEBOも豊富にある。** 2017年55、2018年60異なるBIO siteでSA点がBIOより下流側にあり、近接した採水点を上流STWと呼ぶだけでは偽の曝露を増やす。正向き/逆向きBIO siteが複数SAに重複対応するため、55 vs 36の単純な比率を正式ランダム化pに変換しない。

## SOURCE STOPと結論上限

1. 元source frameは「現在登録SAから2017/18のBIO siteに2km以内のsource候補367点」だけであり、**本当の全国全WWTP/全上流施設**を網羅していない。50km上流のすべての施設を検索したわけではない。
2. 元SOURCEプロトコル固定閾値は異なる有向連結BIO≥50、同年の延べBIO site×year≥100、独立5真水文ブロック以上。今回は51で第一条件だけ通過、**79/100でSTOP**。結果を見てsnap500m→2kmへ変更、Same reachの順序不明を上流に追加、2010–16BIOに2017/18処理水を転用等の救済を行わない。
3. **EAのSAは最終処理水を採水する監視地点であって、正式物理outfall座標だと検証していない。** サンプリングポイントが大河川から500m以内にあることも、下水配管・放流口・小河川を経由して本当にそのreachに合流したかとは別。したがって今回の有向ネットワークは「潜在的なsource-side upstream screen」に留まる。
4. EA源は処理後水の総リン**mg/L濃度**であり、流量、同時河川水質、実際の下水由来TPの下流受水負荷kg/yearはまだない。Wessex原CSV/XLSXはHTTP403、正式WRC日流量は公開更新中。厳格EEA2018元適格実測N+P両方のsourceには英国のreporter codeがなく、**不足をゼロ源としない**。
5. HydroRIVERS v1は集水域10km²/平均流量0.1cms以下の小流路を省略するため、missing routeは測定機会不足も示す。「下水の生態影響がない」証拠ではない。

**生物学的に最も強い既存結果はPR #135にある**：同じ英国ChemPop 591site/2,624site×yearの底生科均等BMWPに対して、河川orthophosphate+nitrate+ammoniacalNの原期間集計を加えるとモデルEDFの独立した地域外予測は−0.519%、その逆の水質3種after EDFは+10.56%、LOQ1/2/3と1site1行感度でも維持。今回のSOURCE-only river routeから「本物の下水TPがBMWPを変えた」ことへの昇格は不正確。さらに正式permit/outfall identityと放流流量/監視年月の照合ができるまでは、STW→BIOの新しい生物回帰を停止する。

方法：[真正source契約](PROTOCOL_20261120_HYDRORIVERS_DIRECTED_STW_TP_SOURCE_GATE.md)、原コード analysis/20261120_wessex_wwtp/ea_SA_TP_ChemPop_HydroRIVERS_directed_route_gate.py。元のHydroRIVERS、EIDC site、EA SA native ID/座標、実TP数値などは公Gitに再配布せず、原SHAと集約件数のみ保存。Main/他PRは変更しない。
