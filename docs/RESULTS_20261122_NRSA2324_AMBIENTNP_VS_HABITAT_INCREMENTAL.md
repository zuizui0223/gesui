# NRSA 2023–24：河川内実測N/P と物理生息環境の独立した予測情報（2026-11-22）

**主な肯定結果：河川の樹冠＋幅は、実測河川TN/TPと地域共変量を調整しても底生動物MMIを強く予測する。** ただし水質のリンを下水由来と判定した研究ではない。交互作用が予測をほぼ改善しなかったのを受けた**事後的な探索的検定**。

## 調査方式と同じ母集団

[PR #140 原本交互作用と方法](https://github.com/zuizui0223/gesui/pull/140)、前結果[実測TP×樹冠](https://github.com/zuizui0223/gesui/blob/analysis/nrsa2324-ambienttp-canopy-benthos-20261122/docs/RESULTS_20261122_NRSA2324_AMBIENT_TP_CANOPY_BENTHOS_INTERACTION.md)、[追加検定の仕様](https://github.com/zuizui0223/gesui/blob/analysis/nrsa2324-ambienttp-canopy-benthos-20261122/docs/PROTOCOL_20261122_NRSA2324_BENTHOS_AMBIENTNP_VS_PHAB_INCREMENTAL.md)。

- **[Actions #45014778716 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45014778716)**、署名された集約JSON artifact ID **14750894606**。EPA5原本SHA・元2,139共通UIDと exact visit/date/site 0不整合を継承。
- 以前の交互作用モデルの**main effects only** HUC8均等MSEを全4群で正確に再現（BERW MMI 257.589433、EPT 0.04052410、BETB MMI 290.342764、EPT 0.03666350）。
- 比較する**同じ4つの固定モデル**：地理/河川サイズ/季節/電気伝導度B、B+河川内実測TN/TP（N）、B+樹冠密度/河幅（H）、B+N+H。Ridge α10、HUC8-disjoint fixed fivefold、各fold train-only standardization/region encoding。すべて**同じ完全観測サンプル**で比較。HUC8均等のheldout MSE改善と、fixed OOFの2,000回 whole-HUC8 bootstrap記述95%を報告。係数の因果解釈は禁止。

## 本当の追加予測情報（同じ流域を保持した外部予測）

| 河川と応答 | UID / 物理SITE_ID / HUC8 | **TN+TP追加 after 樹冠・幅** [記述95%] | **樹冠・幅追加 after TN+TP** [記述95%] |
|---|---:|---:|---:|
| BERW 渡渉：MMI | 1165 / 1046 / 693 | **+5.259%** [+2.199,+8.222] | **+5.994%** [+2.033,+9.628] |
| BERW 渡渉：EPT割合 | 1165 / 1046 / 693 | **+2.480%** [+0.353,+4.612] | +1.397% [−0.885,+3.638] |
| BETB 舟型：MMI | 725 / 666 / 476 | +0.168% [−2.159,+2.415] | **+8.475%** [+4.111,+12.636] |
| BETB 舟型：EPT割合 | 725 / 666 / 476 | −0.189% [−0.960,+0.533] | **+4.232%** [+1.043,+7.459] |

物理的環境 after B もBERW MMI +7.750% [3.633,11.500]、BETB MMI +7.595% [3.607,11.430]、BETB EPT +4.207% [0.973,7.441]、BERW EPT +2.034% [−0.370,4.428]。実測N/P after B はBERW MMI +7.029% [3.538,10.527]、BERW EPT +3.109% [0.662,5.601]、BETB MMI −0.792% [−2.200,0.589]、BETB EPT −0.215% [−1.003,0.589]。

## 科学的結論

1. **MMIでは**物理的環境の追加予測情報が渡渉・舟型の双方で正であり、特に舟型で+8.5%。この結論は栄養塩の追加情報とは独立。
2. **実測N/Pの追加情報は河川型依存**：渡渉のMMIとEPTでは正、舟型ではどちらも原本群を外した予測改善は不確か。
3. **交互作用は必要ではない**：原本実測TP×樹冠の4応答におけるheldout MSE利得は−0.12%〜+0.42%、全てCIに0。物理生息環境の**main effect information**の方が再現的な信号。

限界：河川サイズ／プロトコルによって種組成・MMI設計・生息環境勾配が異なるので、舟型で栄養塩が無影響とは言えない。今回「物理環境」は **樹冠・河幅の2指標**、河床や光の実測ではない。日本や欧州の再現ではない。元3時期（2008/13/18）NRSAのEPT群集と物理環境の予測利得に方向性は整合するが、同一の測定定義・アウトカム・母集団ではない。加重・測定誤差・生態地域間の転移も残る。

**現時点の正確な主張**：`In-stream nutrient and physical-habitat predictors add complementary but river-type-dependent information about benthic community condition across U.S. river basins.` 汚水由来の因果、下水の河畔樹冠緩和、羽化昆虫輸送、生態系サービス連鎖までは示していない。
