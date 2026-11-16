# 同じH6親流域内で植物Nと蝶HNSの空間対応を検証できるか：原本機会監査

**2026-11-15 / Source-only PASS (取得・整合) / ecological model STOP (独立H6数不足) / Draft PR #103**

前の欧州H7×期間451行・90 H6におけるPlantN↔ButterflyHNS空間的正関連（モデルW+元地理/努力調整）は、多くの大地域H6親流域の間で識別される可能性がある。厳しい検定として、同じ真正H6の中の**互いに異なるH7**同士で類似した栄養嗜好の偏りが起こるかを考えた。

## 先に固定したデータ機会条件

- 元H7×期間454、original NUTRI scored H7×期間451、真の元H6=90。元source SHAは `1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca`、EIVは `ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba`。
- 各H7で2017–19/2020–22の**両方**の品質条件を満たすものだけ採る。**219 truly distinct H7、86 truly distinct H6**。同じH7の2期を2個の地理反復として増やさない。
- 親H6内に異なるH7が**最低2つ**存在し、さらに≥100 H7、≥40独立H6、≥5国群を満たした場合のみ H6 fixed-effects 空間モデルに進む。実測前のsource-only閾値を事後変更しない。

## 元のH7/H6源資料での結果

- **86親H6のうち56親H6には、両期間の観測条件を満たしたH7が1つしかない**。
- **H7が2つ以上あるH6は30**（固定条件40に未達）。その30 H6内の互いに異なるH7は163、country_holdout groupは8。
- H6内の異なるH7種数分布：1 H7×56 H6、2×7、3×3、5×3、6×5、7×4、8×5、9×3。
- **MODEL_GATE=STOP_INDEPENDENT_MULTIH7_H6_SUPPORT**。H6固定効果の蝶HNS回帰β、permutation p、地域内予測増分はいずれも**計算していない**。

[元原本source-only GitHub Actions #43721770883 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43721770883)、コード `analysis/20261115_eiv_specificity/h6_multiH7_source_gate.py`。数値原本の再保存 `results/20261115_eiv_specificity/SOURCE_ONLY_SIX_PAIRWISE_WITHINH6_GATES.json`。

## 解釈と次の実施単位

**言えること：** 欧州のNUTRIの固有追加予測はH6を全体留保する広域分散外挿で支持されたが、**元の232 H7における同じH6親内部の空間的反復が十分ではなく、局所河川でのnutrient/hostフィルタリングまで一般化できない**。

**追加データが必要な場合：** まず返り値を見ずに、source registry/HydroBASINSの元H7から現在H7が1つしかない真正H6を選んで**同親H6内2つ目以降のH7を増やす**。元の指定2017–19/2020–22、GBIF植物種ID/蝶speciesKey/元EEA水質・taxonomy source・trait coverageを**source-only**で照合する。先に定めた≥40 H6/100独立H7条件が通らなければ生物応答をFitしない。主結果が正であっても、Wは河岸N処理水曝露ではないので因果的パスではない。

**安全な結論：** 「欧州同一H7での植物NUTRIと蝶HNSの広域空間対応と、同じH7の期間変化の非同期」。H6内の局所生態因果・実測植物被度／幼虫利用／蝶個体数／生態系サービスを主張しない。
