# 2026-11-20：前回の原Bath/全国SA「2km予備抽出」後の有向河川ソース判定（SOURCE-only）

**正式実行:** [GitHub Actions #44679903035](https://github.com/zuizui0223/gesui/actions/runs/44679903035) SUCCESS、集約原JSON artifact ID **14604061338**。これは2017/18年当時に本物の総リンTP濃度が採水されたWater Company FINAL/TREATED sewage effluent EA-WQE SA公式サンプル点と、EIDC ChemPop v2の品質条件を満たす同年の底生生物調査点との、HydroRIVERSの流向に関する*source-only*近似判定。水域生物BMWPの応答・係数、河川の源別TP化学実測を開いたわけではない。

## 元ソースと重要な選択上の境界

- EA公式2026 SA登録5,722箇所、原座標BNG100%有効。ただしこの段階では、2017/18年に本当にTPが測定された元SAの抽出を**ChemPop品質適格2017/18地点から直線2km以内にある現行SA地点367点**に制限していた。その367点の原採水TPコード0348の測定は2017年67点、2018年71点（67点は両年実測あり）。これは英国の「全2017/18元STW」の母集団ではない。
- 同年原底生生物情報（結果値は閉じたまま）：2017年389・2018年366・どちらかに2017/18観測463異なるsite。
- [公式HydroSHEDS HydroRIVERS EU v1 source ZIP](https://data.hydrosheds.org/file/HydroRIVERS/HydroRIVERS_v10_eu_shp.zip)、SHA256 **500da7d36ceee0aa4c82dd625f5c53d80ba2f7ac14ad022cd3a6601bf55068e5**、938,544原reach（英国含む西欧bbox抽出28,830 reach）、本来のNEXT_DOWN、BNG座標へ変換して500m以内最寄り（同距離複数reachは曖昧）。
- 真の化学解析への統合では現地SAサンプル点と原permit/outfallの位置が一致することを証明する必要あり。ここは水文ルートを示すproxyで、実際の施設から河川への流入の原観測ではない。

## 流向に基づくSOURCE-only幾何ルート判定

| 指標 | 2017 | 2018 |
|---|---:|---:|
| 元SOURCE適格BIO site | 389 | 366 |
| EA-SAとして本当にTP採水した元採水点 | 67 | 71 |
| BIO siteが500m以内に河道へスナップ | 371 | 354 |
| SA TP採水点が500m以内に河道へスナップ | 62 | 66 |
| 元SA reach→BIO reach がNEXT_DOWNで50km以内上流側のBIOサイト | **36** | **43** |
| 逆（BIO上流、元SAは下流側）placeboで候補となるBIOサイト | **55** | **60** |
| STW SA採水点とBIOが同一reachにあり順序不明のBIOサイト | 32 | 45 |
| 真の方向が判別された元SA→BIO reachペア | 50 | 58 |
| 逆方向(placebo)元SA/BIO reachペア | 57 | 65 |

実際に「元SA reachがBIO reachの上流にある」候補の和集合は**51異なるBIO site／79 site×year**。50独立BIO siteという旧地域SOURCE下限は数の上では通るが、固定の100 site×year下限を通らない。hydrologic blocksとして出力した元HYBAS_L12の2017は30/2018は33、MAIN_RIVは12/13だが、これを**独立した5以上のHydrological blockの推定量や下水処理場介入の独立反復**とは呼ばない。

逆方向placebo元SA-BIOの地点数や同reach順序不明の数が同規模以上であることは、「上流側」と「近くにある」という選択が違うと明示するために必ず残す。逆のsourceは生態学的未影響の陰性対照と断言できない（他の水源・農地・人口・別施設が流れる可能性）。

**元判定：この2km前絞りでは100 site×year源下限に届かず、物理許可outfall位置・流量・STW別N/P年間kg質量未検証。よって生物回帰SOURCE STOP。** 500m内でsnapしなかった21 BIO site／5元TP-SAは生物種組成陰性ではなく、元HyRivが漏らす小流路/未標定として扱う。源空間距離や自己反復の条件を結果見ながら緩めない。

## 選択バイアスの独立した必要修正：all SA network-first

「生物調査点から直線2km以内の処理水採水点だけ」に源を狭めてから上流50kmを探索すると、**河道に沿って2〜50km上流の施設を研究対象から除外する**。従ってこの79 site×yearは**英国全域・真の50km上流全施設に対する上限でない**。これは元モデルの単なる統計不確実性ではなく曝露母集団の仕様制約。研究効果/標本の追加選択とは独立して、元公式全国の5,722登録SA全地点をHydroRIVERSで先に幾何判定し、下流到達・逆方向placebo・同reach順序不明候補を全て候補に含めてから、2017/18元TPの実測有無を確認する別source感度を、[事前固定プロトコル](PROTOCOL_20261120_ALL_SA_NETWORK_FIRST_HISTORICAL_TP_GATE.md)に従って実行する。

**いかなる数になっても、物理 outfall原位置の未検証や下水年度質量の欠測は残る。** ChemPopの元EDF→BMWPの既知関連、同標本に3栄養水質summaryを入れるとEDF独自予測が消えるという原結論をsource発見だけで因果確認へ昇格させない。
