# ChemPop v2：英国下水曝露proxy EDFと、実測水質N/P指標を分けるsource-onlyゲート

2026-11-20。英国EDF→底生科均等BMWPの負の関係を事後的に確認後の**探索的な反証のデータ取得設計**。まだ水質調整の回帰係数、形質の環境対応、生態機構は計算しない。先にデータ機会のみを固定・監査する。

## 問いと停止条件

元のChemPop EDF_MEAN は河川中処理水のモデル推定混入比率（静的な地点特性）。栄養成分を実測した原分析ではなく、農地流出・都市化・化学汚染・河道変化を伴い得る。したがって今後、下水の「栄養選別」を訴える前に、元同一の生物試料採集日へ時点を合わせた実測水質N/P表と、EDF地域変動の両方があるかを確認する。

- 元v2 ChemPop DOI 10.5285/23f194b0-1e21-4f14-be9c-6eabdcc1feb7、326ファイルのうち287の水質項目×地域ごとの要約表（時刻・LOQ計測数など）を使う。
- 生物データは [主解析 #44548177752](https://github.com/zuizui0223/gesui/actions/runs/44548177752) の元2010–18品質適格**7,836 source sample・707 site**を固定。観測した生物のBMWP・回帰係数は見ず、元分類群のsite/date、BMWP得点可能性、地域、元化学資料のsource-only列だけでgate。
- original sourceから水質候補は**ファイル名にphos/phosph/nitrate/nitrite/ammon/ nitrogen**を含むすべて（正負の効果による選択なし）。各region×analyteで WQ_COUNT>0、WQ_MEAN_LOQ2が有限非負、BIO_SITE_ID＋BIO_SAMPLE_DATE が同一原生物SOURCE標本日と一致する組合せを数える。原水質の「window summary」を独立の同日実測地点濃度とは呼ばない。
- 分析前の最低source基準：ひとつの原analyteに対して真の**500 site×date 組合せ、150 site、4地域**、さらに対象各地域で10site以上が同じ水質と生物のsourceを持つ。WQ_MEAN_LOQ1/LOQ2/LOQ3 と実測数、LOQ以下の割合を確認し、水質原単位・分析法・検出限界に関する正式metadataが欠ければ機構推論はSOURCE STOP。
- 元サイト×日に複数BIO_ANALYSIS_IDがある場合、1つのwq summaryを人工的に独立化しない。水質表を閲覧して応答の符号が望ましい項目だけ選ぶことは禁止。
- N/P候補が品質ゲートを通過した場合だけ、正式な化学種（溶存リン酸 vs TP、アンモニア/硝酸 vs TN）の違いと単位を先に固定。主比較は同一SOURCE可能標本で B＝地理・都市・農地・水文・調査法・季節、B+NP、B+EDF、B+NP+EDF。独立地域留保でEDF独自の増分と水質の独自の増分を別々に測り、事後固定誤差区間を報告する。

**この水質接続が有意でも、下水由来のN/Pを生物が吸収して選別されたという因果媒介にはならない。** 真の処理場源別年間負荷と受水域到達量がない限り、EDFの地理proxyと栄養塩を区別した条件付き関連に限る。

データ著作権・監査：公開Gitへ元水質表やBMWP scored species個票は再配布せず、元SHA・品質件数・正式原文書だけを保持。元論文・gesui/main・他PRには変更を入れない。
