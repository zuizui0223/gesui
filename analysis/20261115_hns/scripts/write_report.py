"""Generate the result ledger from verified outputs, not hand-copied numbers."""
from pathlib import Path
import json,pandas as pd,platform
ROOT=Path(__file__).resolve().parents[1];O=ROOT/'output'
s=pd.read_csv(O/'model_summary.csv');g=pd.read_csv(O/'diagnostic_summary.csv');q=json.loads((O/'SOURCE_QA.json').read_text());qa=json.loads((O/'INDEPENDENT_QA.json').read_text());di=json.loads((O/'diagnostic_results.json').read_text())

def table(df):
    lines=['|解析|流域×期間|H7|H6|β / Wの1SD [95%区間]|H6予測MSE改善% [95%区間]|国均等予測改善%|','|---|---:|---:|---:|---|---|---:|']
    for r in df.itertuples():
        country='未実施' if pd.isna(r.country_gain) else f'{r.country_gain:+.3f}'
        lines.append(f'|{r.analysis}|{r.n}|{r.H7}|{r.H6}|{r.W_beta:+.5f} [{r.beta_low:+.5f}, {r.beta_high:+.5f}]|{r.H6_gain_pct:+.3f} [{r.H6_low:+.3f}, {r.H6_high:+.3f}]|{country}|')
    return '\n'.join(lines)
text='''# 寄主栄養選好と下水インフラ地理：同一流域スケールの初回HNS解析

2026-11-15。**探索的なマクロ再解析。メタ解析・現地上下流実験・因果効果の推定ではない。**

## 結論

元の454流域×期間を保った解析で、下水インフラの静的地理指標Wと、記録数で重みづけした群集HNSの正の条件付き関連が得られた。科別記録比率・土地利用・記録努力・スコア対応率を調整したβは+0.1301 [0.0626,0.1975]／Wの1標準偏差。連続的な地理5項を加えても+0.1080 [0.0389,0.1772]。

一方、全454行の未知H6への予測改善は地理補正後+5.374% [-5.901,14.780]で不確実。種を等重みとした群集HNSは地理補正後β+0.0342 [-0.0017,0.0700]。**より強い信号は種リストそのものの全面的な交代ではなく、観察記録の比率を含む組成の偏りにある。実個体数の増加や生息密度の変化は測定していない。**

## 元データと母集団

234 H7×2期間×6群の2,808取得済みfacetから、従来のEEA・土地利用との共通標本454行／232 H7／90 H6を再構築した。期間は2017–2019と2020–2022。元の303蝶キーをすべて種一覧・観察記録表に保持した。受理済み分類かつ数値スコアがある186キーを主スコア計算に使い、117の未対応・分類保留キーをゼロや無影響種に変換しなかった。

5科の原記録総数4,400,147、種キーの付いた記録4,318,269。スコア付き186キーは全原記録の91.064%を占める。流域×期間ごとの対応率は54.321–99.678%、中央値87.551%。資料の整った種のみ（既知寄主N対応率>=80%）では63キー、全原記録の13.676%まで対象が変わる。そのため主解析は全454行を保持し、この63キー結果を独立な確認と呼ばない。

元の歴史的取得時点に比べ、新しいspeciesKey再取得では454行中51行で総記録数がわずかに変わっていた。最大差6記録。今回は取得済みの新しいspeciesKey内訳を固定し、外部GBIFへ再照会していない。記録は個体数ではない。

各行のEEA報告国ラベルは15種。H7が国境をまたぐ場合等を含むため、国除外交差検証はH7ごとの最多報告国（同数なら辞書順）に固定した11群で実施し、主たる国均等評価は>=10行・>=3H6の6国。6つの独立研究という意味ではない。元の空欄の生物評価法44行はUNKNOWNのまま保持した。

## 群集HNSの定義

既存スコアは蝶の既知寄主に対応するEllenberg型Nの平均。Nは寄主生育地の栄養条件／生産性への選好であり、葉内N、蝶が摂取したN、吸塩行動、下水曝露濃度ではない。

記録重み指標：C_record(j) = Σ_s n_js HNS_s / Σ_s n_js（数値HNSのある種のみ）。種均等指標：C_incidence(j) = Σ_{s:n_js>0} HNS_s / 数値HNSのある記録種数。全種の記録を分母にした対応率を別列に残す。これは従来の蝶生息地最適値から作られたCNIとは異なる。

## 主解析・固定した感度解析

次表のβは固定した種スコアと共変量を条件とするH6クラスター95%区間。予測改善区間は3,000回の固定予測誤差の群再標本化による記述区間であり、同じ種類の不確実性ではない。MSE改善率の符号は、生物量の増減方向ではない。

'''+table(s)+'''

## 地理補正と時期間差：追加診断

位置の補正不足が新しい正の関連を説明し得るため、結果確認後の追加診断として、元のGBIFリクエスト内の照合済みpolygonからEPSG:6933の面積重み重心を計算した。緯度・経度・各二乗・積の5項を固定して追加。気候・標高・河畔微環境を実測していない点は残る。事後追加であることをPOST_RESULT_DIAGNOSTIC_ADDENDUM.mdに記録した。

'''+table(g)+'''

地理補正後の全454行の記録重み指標は、H6予測+5.374%で不確実だが、国を除外し共有H6も学習から除いた評価では国均等+6.415% [1.233,14.870]だった。元の分類同定率>=80%の390行ではH6予測+10.996% [3.171,17.102]だが、同じモデルの国均等は+1.332% [-11.798,7.833]と不確実。対象選択や検証範囲に依存するため「全欧州で頑健」とはしない。

フランス報告行のみの事前固定モデルはβ+0.0175 [-0.0367,0.0718]、H6予測−7.286%。フランス以外はβ+0.2637 [0.0997,0.4277]、H6予測+10.183% [-11.399,26.834]。地理補正を追加する前の診断であり、地理補正結果と同一仕様として混ぜない。

## 同じ流域で高HNS側へ変わったか

同じ222 H7／89 H6で期間差を取った。記録重みHNSは上昇112・低下110、平均差−0.01725、中央値+0.00355で、全体が一方向に高HNS化したとはいえない。静的Wを加えるとΔHNSのH6予測は+5.412% [0.963,8.858]改善したが、条件付き方向はβ+0.0353 [-0.0162,0.0869]、国均等予測+5.632% [-5.116,7.865]で不確実。種均等ΔHNSの予測は改善しなかった。**Wは両期間で同じであり、下水の増減が蝶を変化させた検定ではない。**

## 科内の対応と寄与種

各科の比率を調整した後も主記録重みβは正。種に割り当てられたHNSを科内で999回置換した記述的対照で、観測β0.1301は置換95%範囲[-0.0698,0.0680]を上回った（片側の有限反復比率0.001）。種均等では観測0.0528、置換95%[-0.0377,0.0372]。系統・生活史や分布の交換可能性を保証しないので確証的p値ではない。

species_slope_decomposition.csvは、群集傾きへの種別の加算寄与を保存した。高HNS種の正の記録比率係数だけでなく、低HNS種の負の係数も群集HNSの正の傾きに寄与する。個々の種の因果効果や有意差を検定した表ではない。Aglais ioを除く事前固定感度でも主傾き+0.1345 [0.0644,0.2046]が残り、一種の目録欠落だけで全関連は消えなかった。

## 欠測への重要な制約

「欠測寄主のNをすべて1」「すべて9」の2シナリオで同じ調整済みβが得られたが、これは独立な頑健性証拠ではない。2指標の差は厳密に8×(1−寄主N対応率)であり、その対応率がモデルに入っているため、線形調整後のW係数は構造上同じになる。

そこで、各スコア付き種の記録済み・未数値寄主が持ち得る既存の下限／上限を、W係数に不利な方向へ種別に割り当てた固定デザイン感度範囲を求めた。記録重みβは **−0.14956～+0.32447**、種均等βは−0.13299～+0.19106で符号が固定されない。この範囲は95%信頼区間でも起こりやすさの分布でもない。欧州EIVのx（無差別）に数値を仮定する極端な感度であり、未記載の主要寄主や全117未スコアキーを救済するものでもない。地域で実際に利用する主要寄主による補強が残る。

## 推定対象と留保

WはH7出口へ集約した静的モデル下水水量／集水面積の地理指標であり、放流口の上下流比較でも土壌・葉・蝶で測った下水N負荷でもない。Wと都市土地利用はρ=0.7806。都市性、気候、標高、記録者や種別の報告選好の交絡を除去し切れていない。河畔局所の吸塩行動、幼虫採食、送粉、個体数やバイオマスの変化は示していない。

**支持される表現：下水インフラ地理と、既知寄主の栄養選好で要約した蝶の記録組成との正の関連が、科間差だけでは説明されずに残る。ただし記録重みに強く、母集団選択・地域・欠測に依存する。**

## 数値検証・実行範囲

新規14テスト、既存スコア13テスト、既存ネットワーク7テストが通過。元の37,785正記録行をgroupbyして全454行の群集指標を再計算し、最大絶対差1.78e−15。20モデル仕様の全H6・国除外予測を、独立に構築した学習専用標準化・カテゴリ行列とsklearn SVD Ridgeで再計算（最大差5.40e−13未満）。重み付き係数・クラスター区間はstatsmodels WLSで別計算し一致。独立な生態調査を追加したわけではない。

全作業は取得済み原本上のローカル解析。GitHubへの変更・コミット・マージ、新しいActions実行はない。ソースと数値を再現するコマンドはREADME.md、入力SHAはinputs/SOURCE_MANIFEST.json、主要監査はoutput/INDEPENDENT_QA.json。HNSスコア自体の生成原本は前回のスコア作成ZIPにあり、このパッケージは固定したスコアCSV以後の新規解析を再現する。

## 出典

- Tichý et al. (2023), Ellenberg-type indicator values for European vascular plant species. Journal of Vegetation Science 34:e13168. DOI 10.1111/jvs.13168.
- 固定版原本：Zenodo DOI 10.5281/zenodo.7427088、CC BY 4.0。
- Ferrer-Paris (2014), butterfly larval-host records. DOI 10.6084/m9.figshare.1168861.v1.
- WallisDeVries & van Swaay (2017), A nitrogen index to track changes in butterfly species assemblages under nitrogen deposition. DOI 10.1016/j.biocon.2016.11.029。今回のHNSと同一指標ではない。
- 元のgesui data acquisition: GitHub Actions run43593227162、source ZIP SHA256 4b37c741d1bd28c178b0cef61c9457af0dd30797c173f9384ebe1d1d0036d07c。
'''
(ROOT/'RESULTS_20261115_HNS_MACRO.md').write_text(text)
(ROOT/'README.md').write_text('''# Reproducing the new HNS composition analysis

All sources are frozen user-provided files, not synthetic ecological observations. No new downloads are needed. Keep the source score CSV and occurrence archives unchanged. `source_speciesKey` counts are observation records, not individuals. The model does not establish wastewater-source dose or causality.

Run from this directory with Python and NumPy/pandas/SciPy/sklearn/statsmodels/shapely/pyproj installed:

```bash
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -v
OPENBLAS_NUM_THREADS=1 python scripts/build_panel.py
OPENBLAS_NUM_THREADS=1 python scripts/analyse.py
OPENBLAS_NUM_THREADS=1 python scripts/diagnostics.py
OPENBLAS_NUM_THREADS=1 python scripts/validate.py
python scripts/write_report.py
```

`prior/` contains source reference code/test subsets and the retained network matrix for count reconciliation. It is not the complete previous project; only the commands above are promised runnable from this archive. HNS trait-score generation itself is in the previously delivered host-nutrient-score package. This archive starts from that fixed CSV.

`docs/PROTOCOL.md` was written before fitting the new effects. `docs/POST_RESULT_DIAGNOSTIC_ADDENDUM.md` records explicitly post-result geography and adversarial-missingness checks. All nulls and prediction failures remain. The new data source is shared with earlier analyses; numerical verification is not ecological replication.
''')
import numpy,scipy,sklearn,statsmodels,shapely,pyproj
(ROOT/'environment.json').write_text(json.dumps({'python':platform.python_version(),**{m.__name__:m.__version__ for m in [numpy,pd,scipy,sklearn,statsmodels,shapely,pyproj]}},indent=2))
p=ROOT/'docs/PLAN.md';p.write_text(p.read_text().replace('- [ ]','- [x]'))
print('Report created:',ROOT/'RESULTS_20261115_HNS_MACRO.md')
