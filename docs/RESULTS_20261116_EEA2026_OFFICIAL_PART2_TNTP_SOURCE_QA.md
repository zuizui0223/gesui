# EEA Waterbase Water Quality ICM 2026：公式出版原本の取得成功と同じ生物観測局の水質照合

**2026-11-16 / PR #107 / source-only 原本品質結果／河川生物・植物・蝶の新しい効果推定なし。**

## 旧Discodata と公式2026年版の区別

[EEA公式 Waterbase Water Quality ICM 2026](https://www.eea.europa.eu/en/datahub/datahubitem-view/fbf3717c-cd7b-4785-933a-d0cf510542e1) は2026-07-07公開、対象1900–2025。別の [WISE SoE公開窓口](https://water.europa.eu/freshwater/resources/wise-soe-data-collection) は既存Discodata WISE6表が最新の公開版と限らないと明記している。元のDiscodataは yearだけ・TNだけ・TPだけではHTTP200だが、複数AND条件はHTTP403。**HTTP403を水質データ不存在と解釈しない。** [旧API source test #43797159964](https://github.com/zuizui0223/gesui/actions/runs/43797159964)。

公式ディストリビューションはNested Nextcloud `https://sdi.eea.europa.eu/datashare/` に置かれている。元の `/data/82e194cc-2506-411c-b341-9eb50b6ad5e0` 画面で公開共有 tokenを特定し、`/datashare/public.php/dav/files/ZcJ37WK6oaWonFJ/` を **PROPFIND Depth1** で照合した結果HTTP207。2026版の実ファイルの存在と正確なバイト数を得た。[公式source file list #43797712726 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43797712726)。

| EEA2026 公式公開ファイル | 原本ZIPバイト |
|---|---:|
| `WISE6_AggregatedData-csv.zip` (Part2 annual station) | **170,422,309** |
| `WISE6_DisaggregatedData-csv.zip` (Part1 individual sample) | 1,639,231,716 |
| `WISE6_DisaggregatedData_sqlite.zip` | 2,515,599,514 |
| `WISE6.sqlite.zip` | 242,604,679 |

## Part2 年別観測局の公式170MB原本を真正ダウンロードして全走査

[公式Part2 run #43798318436 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43798318436) は`WISE6_AggregatedData-csv.zip`を**170,422,309バイト**そのまま取得、SHA256 `315396800d76e46db415164851fb4f99765b7ec2654828fc42782d902e679a7a`で保存。CSV **6,450,192行**を全行検査。全外部データZIP自体をGitリポジトリやGitHub Artifactsへ再配布せず、元EEA EQR 3525 site-yearの位置・source最小数値表と品質監査だけをActions artifactで保持。

- 元の真正EEA EQR物理station×year: **3,525行／2,485固有station**。2018年1,803局×年、2022年1,722局×年。元H7×期間454・232 H7/90真正H6の場所だけを検査し、EEAオリジナルGit blob 2018/2022の現実に上流年別報告N/P正を持つstation×year **132**を独立に固定。報告なしを未曝露対照とはみなさない。
- Part2全体の2018観測年に254,728行、2022観測年に580,704行が存在する（この数は全欧州、すべての化学項目を含む）。
- 2018/22×元の正確なEEA生物局ID+schemeでの原本一致を満たすPart2の全項目行1,387（同じ観測局×年でも複数の測定項目があるため重複を含む）。河川RW水行W、TN determinand `EEA_31615-01-7`、TP `CAS_7723-14-0`、mg{N}/L / mg{P}/L、LOQフラグと最低3試料、年・元国が完全に合致すると、TNは**20 site×year**、TPも**20**、同じ局と年で両方のTN/TPが揃うのは**20 site×year**。実際に年別正の上流処理場N/P報告がある元サイトは**1 site×year**（1 H7×期間・1 H7・1 H6・1国）。
- 20 site×yearは**ベルギー**の2018年8、2022年12であり、H7×期間は合計7（2018年4、2022年3）。これは欧州でTN/TPが測られていないことを意味せず、**Part2 sourceに現れる元EQRサイト×年に集計化学が少ない**ことを示す。
- source gate: 元PR #91より固定していた、上流報告正かつ同じ観測局2018/22年の実測TN+TPがあるH7×期間≥30・H7≥20・H6≥10・国≥3には**1 H7×期間/1 H7/1 H6/1国**で未達。**STOP**。生物指標の係数や生態サービス、下水媒介モデルを作らない。

[追加の欠測構造の再検証 #43798670501](https://github.com/zuizui0223/gesui/actions/runs/43798670501)：国別・年別の20元観測局と、Part2での厳密対応ステージおよび除外理由をすべて数えた。2026公式Part2の受理項目の存在を確認したが、全欧州での**個別試料表 Part1** はまだこの段階で検証していなかった。

## 次段：公開個別試料データの原本を別途検査

元の同じ2018/22 EEA生物調査局のTN/TP個別試料記録を、公式 `WISE6_DisaggregatedData-csv.zip` で確認する。個別試料と年別集計表はEEAの仕様上別々の源であり、Part2が空でもPart1に測定がある可能性がある。採用するのは真正年・観測局ID+scheme、河川RW・水W、原本TN/TPコードと単位、定量下限、有限な測定値、**異なる3試料以上**を両方の栄養塩で満たす局のみ。原本ZIPは全長とSHAを確認し、元3,525 site-yearだけをstream filterする。

Part1が条件達成しても、TN/TPは**受容河川の水質**であり下水由来の実際の到達流量・根圏のNや葉窒素・幼虫食草・蝶個体数を直接測ったものではない。元PR #106でmainに統合した肯定結果は維持しつつ、下水特有の因果媒介の新しい推定はsource-only gateの後に別仕様で検討する。

## Provenance

- `analysis/20261116_reported_np_source/probe_eea_datashare_correct_dav_v4.py` / Actions #43797712726。
- `analysis/20261116_reported_np_source/acquire_eea_2026_part2_exact_site.py` / Actions #43798318436。
- `analysis/20261116_reported_np_source/diagnose_official_2026_part2_matched_qa.py` / Actions #43798670501。
- `results/20261116_reported_np_source/EEA2026_OFFICIAL_PART2_3525_SITE_YEAR_BRIEF_QA.json` はActionsの原本ログからSHA・達成状況・元CSVヘッダだけGitへ複製、Part2全原本の試料は複製しない。
- `analysis/20261116_reported_np_source/acquire_eea_2026_part1_exact_site.py` は別個のsource pipeline。元統合PR #106/mains/chochoは編集しない、PR #107はDraft。
