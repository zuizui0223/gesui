# EPA NRSA2023–24 実測河川リン × 河畔樹冠 × 底生動物：流域外予測の検証（2026-11-22）

**判定：樹冠との交互作用の広域一般化は支持されない。** 下水由来リンの結果ではない。探索的／観測研究。下水特異的な試験を原本曝露の不足でSTOPした[PR #139](https://github.com/zuizui0223/gesui/pull/139)とは、曝露概念も母集団も異なる。

## 再現可能な原本

- [SOURCE ONLY Actions #45014068137](https://github.com/zuizui0223/gesui/actions/runs/45014068137): EPA公式2023–24 SiteInfo・WaterChem・PhysicalHabitat・BenthicCount・BenthicMMIの原本SHAすべて一致。**2,139共通UID**、調査方式が一致した渡渉BERW/WADEABLE **1,307 UID**、舟型BETB/BOATABLE **825 UID**。ほか7不整合。結果値を見る前にソース成立を確認。
- [Source protocol](https://github.com/zuizui0223/gesui/blob/analysis/nrsa2324-ambienttp-canopy-benthos-20261122/docs/PROTOCOL_20261122_NRSA2324_AMBIENT_TP_CANOPY_BENTHOS_POSTHOC.md) を数値解析前にコミット。研究そのものは2023/2018既報結果を見た後の**事後的**追加解析。
- [NUMERIC Actions #45014519276](https://github.com/zuizui0223/gesui/actions/runs/45014519276): 実際のEPA原本5表・本来のUID/visit/date/site一致を独立再照合（0不整合）。5fold**全HUC8群を保持**した地域外Ridge10。検証結果のoriginal artifact **14751213267**、ZIP original digest `009b05f2de6597cfd0750a6bf523f0de7c1e9d51555e139b3a7fd6261b46034a`。公開結果は集約値のみ。
- 実測河川TP `PTL_RESULT`・NARS QA旗空欄・単位原本統一（µg/L→mg/L）、採取当日の水質。これを**処理場から届いたPと置き換えない**。TN, conductivity, 河幅, 州、生態区、緯度経度、河川次数、年と季節を調整。
- BERWでは中央樹冠密度`XCDENMID`、BETBでは河岸樹冠密度`XCDENBK`。測定定義が異なるので**別モデル**。
- 主応答EPA `MMI_BENT`、副応答`TOTAL300`のEPT（Ephemeroptera, Plecoptera, Trichoptera）個体数割合。後者は羽化バイオマスではない。

## 正真正銘の原本から推定された4組の結果

| 独立方式・応答 | UID/物理SITE_ID/HUC8 | standardized `logTP×canopy` β [HUC8クラスター記述95%] | p（補正前） | HUC8均等の地域外MSE利得 [HUC8 bootstrap 95%] |
|---|---:|---:|---:|---:|
| BERW：MMI | 1165/1046/693 | −0.03187 [−0.07737,+0.01364] | .16990 | **−0.0666%** [−0.6421,+0.5148] |
| BERW：EPT個体割合 | 1165/1046/693 | −0.07612 [−0.12654,−0.02569] | .00309 | **+0.2266%** [−1.0536,+1.4105] |
| BETB：MMI | 725/666/476 | −0.08356 [−0.15767,−0.00944] | .02713 | **+0.4220%** [−1.1797,+2.0060] |
| BETB：EPT個体割合 | 725/666/476 | −0.03241 [−0.10961,+0.04479] | .41057 | **−0.1237%** [−0.7989,+0.5618] |

固定した4件にHolm補正（探索的感度）を適用すると対応pはBERW MMI .3398、BERW EPT **.0124**、BETB MMI .0814、BETB EPT .4106。**この1件の同時点記述係数の信号は、流域外予測の改善としては安定しない。** 4件いずれも HUC8 均等 fixed OOF bootstrap95がゼロを含む。観測地点重みのMSE利得も −0.052%,+0.512%,+1.238%,+0.093%と小さい。

## 生態学的解釈と境界

- **局所的な樹冠×栄養塩の関連**は一部の応答に存在する。しかし広域普遍的な**緩和**、とりわけ「河畔林が下水の悪影響を弱める」とは支持されない。交互作用が負だからといって、単独のリン勾配の向きを評価せず「林が悪化させる」とも言えない。
- *以前の独立な2018–19 EPA解析* [PR #78](https://github.com/zuizui0223/gesui/pull/78) のChla×樹冠密度とEPT割合の交互作用も支持されなかった。ただしChla×canopy と現在の ambientTP×canopy は異なる検定である。
- 同一HUC8内の繰り返し訪問、空間共変量や未観測都市化・河床・遮光・酸素・流量などの交絡が残る。元の2018/13/08別調査期の地域外輸送検定結果とも同じ母集団ではない。
- **今後**は交互作用を主結果へ押し込むのではなく、同一方法・同一観測単位での**物理生息環境の独立予測情報 vs 栄養塩の独立予測情報**を全4条件で比較し、その後、下水という発生源を別の測定設計で検定すること。

**全モデルは posthoc 探索的、因果推定ではない。**
