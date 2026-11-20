# 2026-11-20｜EA全国最終処理水SA採水台帳とChemPop原生物地点：全国的地理機会のSOURCE-only評価

**研究対象は実際のN/P放流が届いた生物影響ではなく、全国の原採水点台帳と2017/18真正生物観測地点の地理的SOURCE機会。** 現行登録SA=2017/18年のSTW稼働・実測とは同一でない。

## 原API契約と地理

[公式API SA/SBの厳密フィルタ試験 #44677052429](https://github.com/zuizui0223/gesui/actions/runs/44677052429)：EA新Water Quality Explorer /data/sampling-point POSTの samplingPointType=SA は全100 sample membersが本当にSA、metadata totalItems=5722。対照 samplingPointType=SB ではstorm overflowだけ100件で totalItems=1356。無指定の全サンプルは総66,363、SA/SB以外も混在。samplingPointType.notation=SA はHTTP400であり使わない。キーワードSewageではSA以外も多数返り、正式フィルタの代用不可。

[原native geometry schema #44677344654](https://github.com/zuizui0223/gesui/actions/runs/44677344654)：当該SA元地点memberには id、notation、prefLabel、samplingPointStatus、samplingPointType、英国座標元geometry.asWKT等のフィールドがある。元source WKTのCRS注記・座標域をその場で検査し、WGS84の経緯度をメートルと混同しない。EA新WQE sourceに英国BNG native coordinatesがない場合は空間SCREENをSTOPして、勝手な投影に変換しない。

## 事前固定する全国SOURCE品質チェック

1. **全件ページ完走**：250件ずつskipを増やし、公開totalItems=5722と全出現member件数・一意native sampling-point notation数が一致しない場合はSOURCE STOP。最大40ページ以内に限り、API 429・source失効等は正式stopとする。
2. **座標と同名重複**：original source geometry.asWKTが英国BNG POINT(E N)であることを検証、地点位置 0<E<800,000、0<N<1,300,000 に限定。元情報にEPSG:4326と明記されている場合のみ座標変換し、CRS不明の違う単位を採用しない。重複notationで異なる位置があればSTOP。元SA点の座標確保率は**95%以上**をハードゲートとし、不達は欠測として完全全国ジオモデルをSTOP。
3. **ChemPop元1,519全siteと元2010–18 score適格707site**は prior source-only SHA、original WRC siteMatchedLocation元ZIP、BMWP reference SHAに一致させる。2017年本物の適格389site、2018年366site（合計独立サイト463）を別々に保持。全国sourceとの距離はsource site×bio siteの原英BNGメートル最短距離、0.5/1/2/5/10 kmの五固定閾値。地理的近さと河川流下を同一視しない。
4. **年代の非対称性**：全国5722 source pointsは2026年の現行台帳であり、2017/18に全てが稼働、最終TP原液を測っていたとは主張しない。全国地理マップはあくまで潜在的なSTW source sampling point opportunity上限。真の過去の曝露を得るには原pointNotationに対応する2017/18のWQE実測TP観測を再照会する追加SOURCE段階が必要。
5. **source状態と流向**：source status active/obsolete、不明を集約して保存。正式permit/outfall point identity（採水施設内位置と放流出口は異なり得る）およびNEXT_DOWNなど有向河川トポロジーがない限り「下流」ではない。nearest 500m /10kmと独立河川網のsource-to-receptor directed距離を混同しない。
6. 原WQE source名・notation・住所・観測の個票は公Gitに再配布せず、元source件数・全メンバー内容のSHA・source地点type/status/地域別件数と全国距離screenだけを成果物に保存する。EIDC/BMWP元個票も同様。

## モデルへの移行要件は以前の規約から変更しない

現地で2017/18に実測TPがあるSA点≥30、独立水文流域5以上、真正同年ChemPop生物siteで元有向河川ネットワークを経由し下流にある≥50site/≥100site×yearという基準を満たして初めてSTW source生物応答の探索的モデルへ進む。全国現行台帳の近接screenが多数あってもSOURCE GATE_PASSとはしない。測定TP原液濃度(mg/L)から年間kg/yrを勝手に計算しない。

本解析はBath50kmで実測SA点46/49、同年ChemPopの10km直線近接5/4でSTOPしたsource pilotを全国の**地理的機会**に拡張するものであり、記録の有無で試験地点・閾値を選び直さない。多国・多分類群での下水由来N/Pの普遍法則、栄養介在の因果パス、局所の個体群の変化は主張しない。
