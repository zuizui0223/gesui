# 2026-11-19｜経路AとCが同じ下水処理場から分岐する実施用記録：Feurs 2018 の独立SOURCE監査

**結論：実施用を記載した行政報告書で確認できる自治体が、ROSEAU公開6年候補にない。公開候補の不在は「汚泥未施用」の証拠ではない。生物モデルは引き続きSTOP。**

## 公式記録（ROSEAUの候補抽出とは独立）

フランスDREAL Auvergne-Rhône-Alpes公開の[Feurs下水処理場の2018年施用実績報告書](https://www.auvergne-rhone-alpes.developpement-durable.gouv.fr/IMG/pdf/10_step_feurs_cas_par-cas_annexe9_epandage_2018bilan.pdf)（Suivi agronomique des épandages des boues – bilan agronomique 2018）では、処理場SANDRE **0442094S0002** の処理水が**Loire川へ放流**される一方、汚泥は農地へ運ばれると記載（原報告p.4）。対象排水は都市・産業混合であり、家庭下水のみの影響とみなさない。

2018年の**実施済み農業施用**は原報告p.5,15–16に**湿重量1,918 t、乾燥重量481.31 t（石灰を含む）、施用面積161.53 ha、12農場**と記録。p.16–17の「Tableau 7 – Cahier d'épandage ... Campagne 2018」は実施用の日付・施用畑コード・面積・湿重量・前後作物を記載し、**8つの異なる施用先自治体**が確認できる。これは申請だけの plan d'épandage ではなく bilan（実施用記録）に基づく。個別農家名・圃場IDは本監査に保存しない。

[INSEEの公式Loire県自治体台帳](https://www.insee.fr/fr/metadonnees/geographie/departement/42-loire)で、原報告の表記を8自治体のCOGコードに対応付けた：**Cleppé (42066), Poncins (42174), Saint-Laurent-la-Conche (42251), Nervieux (42155), Valeille (42319), Feurs (42094), Mizérieux (42143), Marclopt (42135)**。現行台帳のコード・名称照合であり、2018年の歴史COGとの別照合および原PDFバイナリSHAの検証は今後のSOURCE監査課題。

※元PDFのテキスト索引から本文・Tableau 7を確認したが、この時点のページ画像再取得はHTTP 403で失敗した。PDFバイナリの現物SHAは未確保であり、今後ソース再取得時に原ページ・数字を再監査する。推定値の検証ではなく独立のSOURCE-only事例。

## 先行SOURCE候補集合との真の照合

- 既存[GitHub Actions #44465289468](https://github.com/zuizui0223/gesui/actions/runs/44465289468)からFR_ROSEAU_2017TO22_PUBLIC_CANDIDATE_COMMUNE_YEAR_PRESENCE.csvを取得（artifact 14501751936）。原CSV SHA256 **971e20215a0788a229c5f0a18663e1eee9ddf70b6161fc4f70e78c8da9621333**。歴史COG厳密一致の候補和集合は130自治体、2018年は72自治体。
- 外部行政報告書で**実施用**を確認した **8自治体すべて**は、2018年72候補にも2017–22年130候補の和集合にも**存在しない**。同じLoire県の別自治体候補（2018年: 42192, 42229, 42235, 42262）があるため、単に「42県が一切なかった」のではない。ただし報告書の対象処理場が公開ROSEAUの収録対象であった保証はない。候補の全国的網羅性・収録対象の違いも未同定。
- **外部で確認した8自治体×処理場×年の陽性事例は、1つの処理場×2018年に由来する同一SOURCEクラスター**。8つの独立流域や独立処理場・年間反復として数えない。欠落率、候補データの感度、全国未施用率、正式な対照群の数は推定**しない**。本比較は公開抽出の完全性を証明するものではなく、**非収録＝無施用と扱う設計を反証**する一例。
- 再計算コード: **analysis/20261119_fr_sludge_routeC/fr_external_performed_spreading_coverage.py**。手動で原文献から設定した公的事例は **FR_EXTERNALLY_VERIFIED_ACTUAL_APPLICATIONS_2018.csv**。--candidates に上記artifactの元候補CSV、--external に事例CSVを渡す。再計算結果は結果JSONとGitHub Actions成果物に保存。単体テスト6件はローカルPASS、GitHub独立Actionsの最終成功確認は別途必要。

## 水と陸を「同じ下水起源」で結ぶ次の解析単位

**このSOURCE事例の新規性**：同じFeurs処理場が2018年に水出口（Loire川へ処理水放流）と土壌出口（8自治体161.53haへ実際に汚泥施用）の両出口を持つことを確認した。ここでの共通単位は処理場×年であり、水と陸の**生態学的応答の連鎖**が成立したという主張ではない。


**施設×年をSOURCEの共通起点**として「処理水→河川（経路A）」と「汚泥→圃場（経路C）」の二つの出口を別々に列化する。水側は実測の受水口・排出TN/TP質量と受水流量、陸側はbilanによる実際の圃場・年・乾燥汚泥量・製品N/P量。最初の共通単位はFeurs 2018のような**双方の出口が同一処理場に確実に帰属する**事例。ただし現時点でFeursの同年河川TN/TP質量、作物・河岸寄主の実際の曝露、蝶幼虫の摂食は未照合。

同年同施設の二経路マスフローと、別の応答（河川EQR／植物NUTRI／蝶HNS）を混同しない。比較可能な実施用圃場が**複数独立H6**と時点に拡大し、非施用の**確認済み**対照、別の窒素投入がそろうまでは生物への因果効果・媒介モデルを走らせない。独立の水陸経路を立証したのは**処理場出口の分岐のみ**であり、生物間の窒素移送は未証明。

## 近隣の公開GIS候補（非施用証拠、モデルに使用禁止）

[data.gouv.fr Loir-et-Cher (41) の2018年汚泥施用予定区域の圃場GIS](https://www.data.gouv.fr/datasets/agriculture-parcelle-de-perimetre-depandage-de-boues-urbaines-2018-en-loir-et-cher)と[Corrèze (19) の適地GIS](https://www.data.gouv.fr/datasets/parcelles-de-perimetre-depandage-de-matiere-organique-boues-de-stations-depuration-en-correze)は公開されている。しかし**計画・適地のジオメトリであって、実際の施用量・年月ではない**。実施用と同一視しない。1998年施用省令の[2020年改正第20条](https://www.legifrance.gouv.fr/loda/article_lc/LEGIARTI000042419149)に沿い、SILLAGE/VERSEAUのbilan d'épandageを追加照会する方針は維持する。

**既存の14 H7 / 6 H6の地理候補とFeursが同一であるとは主張しない。** 元のフランス151 H7植物SOURCE母集団、139行・72 H7の水陸共通標本、旧NUTS2/JRC SOURCE-STOPを維持。NUTRI/HNS/EQRの結果を一切開かず、PR #129をDraftに保持し、main/論文草稿PR #124は変更しない。
