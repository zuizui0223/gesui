# 2026-11-19｜ザールラント州Klärschlamm WFSの属性定義に基づくSOURCE-STOP

## 検証結果

[公式WFSのDescribeFeatureTypeを取得したActions #44523503343](https://github.com/zuizui0223/gesui/actions/runs/44523503343) はSUCCESS。GetCapabilitiesに正式に掲載された `Abfallbeseitigung:Klaerschlamm_Verfahren` を参照し、フィールドの構造のみを確認（GetFeature未取得、植物・蝶・水域EQR値を一切使用していない）。

- GetCapabilities：HTTP200、15,006 bytes、SHA256 `f6b7ecbabbb2b038c148c45c4ae7ada46b666c2dcf52eb26c90e4d1dc0cf5ea5`。
- DescribeFeatureType：HTTP200、2,310 bytes、SHA256 `23b6c453c44e6624cca8ba07e8fbdb427fb0480fe2470a0ee45e3a827437d3cc`。
- WFS XMLで実際に確認した全属性9個：`Klaerschlamm_Verfahren`、`SCHLAGNR`、`SCHLAGNAM`、`FLAECHE_HA`、`OIDEXT`、`OBJECTID`、`Shape`、`Shape_Length`、`Shape_Area`。
- 地理形状と区画面積を示す構造はあるが、**年・月・日・実施済みか計画か・処理場識別子・乾物量・N/P含有量の属性がない**。区画番号の外部テーブルとの正式な結合仕様がない限り、2017–2022年の実施用曝露は復元できない。
- 元WFSのフィールドに年がないことは**構造的な否定**であり、実際の施用そのものがないことを意味しない。公開GISの年次を推測したりポリゴンを時点間でコピーしたり、別年の全国汚泥処理経路統計を曝露地へ割り付けたりしない。

## SOURCE判定

**GEO_LAYER_SCHEMA_GO / YEAR_AND_REALIZED_APPLICATION_STOP / WWTP_ATTRIBUTION_STOP / NUTRIENT_MASS_STOP / BIOLOGICAL_MODEL_STOP**。

元のドイツ2016–2022年Destatis/ザールラント州資料は**汚泥を出した処理場側の利用経路集計**として独立に利用できるが、実際に汚泥を受けた個々の畑ではない。元フランスFeurs PR #129と同じ推定単位（実施用農地×STEU×年）へGIS単独で昇格させない。

## 次に広域へ進める方向

個々の処理場が発行した**bilan d'épandage réalisé（実績）**の原票があり、同一施設の処理水放流年・受水地点・流量/N/Pが確認できる場所を広く探索する。GISだけで植物や蝶へ対応づけず、SOURCEがある複数施設を独立に数える。その後も栄養同化・生物間媒介は別の検証とする。

結果ファイル：[Actions artifact #14529097444](https://github.com/zuizui0223/gesui/actions/runs/44523503343/artifacts/14529097444)。元XMLは再配布せずハッシュと属性名のみ記録。
