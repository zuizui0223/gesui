# 2026-11-17｜真の下水→水生昆虫→河岸捕食者：公開原本の観測粒度・比較可能性と次のマクロ検証

**2026-11-17 source-first progress**. No ecological experiment was created and no field sampling was claimed. This audit distinguishes genuinely observed source/consumer material from remotely sensed terrestrial greenness and geographic butterfly host-affinity proxies. Source-only documents verified below do not themselves constitute a new causal meta-analysis.

## Original Bow River municipal effluent study: genuine specimen schema, author exposure map

2025 published study (Diesbourg, Kidd & Perrotta; doi:10.1016/j.envpol.2025.125948) links [public author original code & CSVs](https://github.com/Karen-Kidd-Lab-Ecotoxicology/Bow-River-Microbiome-2022), exact source Git commit `32b36abfcf31117ac1bd9d7f3fdbc17cdcedd795`. All raw source column and specimen/site counts came from signed [GitHub Actions #44078675291 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/44078675291), not from synthetic samples or an inferred study-size.

- Original `stable_isotope_analysis/StableIsotope2022.csv`: **184 source specimens × five actual Bow River sites**, SHA256 `e9ca82f50bc6f42a43213df8a0daa6d118b01cbaf16953e9e084d603e3097076`. δ¹³C / δ¹⁵N / %C / %N / C:N have real original instrument result fields.
- Exact original source groups: **53 aquatic larvae, 55 aquatic emergent adults, 30 riparian spiders, 30 basal resources (15 biofilms + 15 leaves), 16 terrestrial potential prey**. Each genuine river site has **6 spider specimens and 3 biofilm specimens**, and larvae/adult specimens.
- `microbiome_analysis/Raw data/metadata_noKK2078_samples_ordered.csv`: 380 rows, SHA256 `f0a77f4efb51e0f0f91b87fbabb8ffcc6ef43685bbceaeec29da5879023ba5e8`. This file contains **eight** named locations total but **five actual Bow River river sites (309 metadata records) versus three experimental ACWA stream channels (61 metadata records)**; additional 10 metadata blanks/nonassigned, not eight independent natural rivers.
- Original source `Exposure` label on the real five Bow River sites: **Cochrane, Cushing Bridge, Sunalta=Unexposed**, **Graves Bridge, Policeman Flats=Exposed**. Author site source classifications, not randomized experimental WWTP treatments. Thus **three upstream stations vs two downstream stations along one Bow River**. A 2025 paper studies three plants on the same river, not 3 independent hydrographic catchment replicates, and 184 organism samples are nested beneath five site observations.
- Original field `Month` for five Bow River natural sites = **July** in the original 2022 experiment metadata. Original stable-isotope CSV's `Date` contains February–March 2023 **analysis** date entries: do NOT treat as in-situ 2023 sampling or 2018–22 temporal replication.
- Complete source-only five-site group census [CSV](evidence/bow2022_5site_author_exposure_insect_spider_source_census.csv); code for actual isotope values `analysis/20261117_emergence_macro_source/bow_2022_isotope_site_level_3up_2down_exact_perm.py`. The *original-source numerical fit* requires the source SHA, group field QC and careful site-level normalization first. No unrestricted original specimen-level pseudo-treatment p values.

### One-site statistical/biological boundaries

A physically meaningful within-source comparison is: for each true **site**, compute original measured δ¹⁵N/δ¹³C medians separately for aquatic larvae, emergent adults, spiders, biofilm, terrestrial potential prey, and spider taxa Araneidae/Tetragnathidae. Report **downstream mean of site-level medians minus upstream mean of site-level medians**, plus **spider-minus-site-specific-biofilm** isotope offsets to distinguish trophic-level enrichment from local background shifts. No measured effluent endmember in this CSV → no precise mixing fraction. Label permutation only across the **five site-level source labels**: only C(5,2)=10 assignments, hence minimum exact one-sided p=0.1 even before spatial trend confounding. These are **source-specific exploratory checks of published local observations**, not independent multi-WWTP causal effects, not regional macro transfer rates, and not a new publication discovery if original authors have already tested the effect.

## Grand River 2022–23, eight municipal WWTPs: *two linked original public replication datasets*

The 2026 Wardlaw et al. paper ([DOI](https://doi.org/10.1016/j.envpol.2026.129116)) sampled upstream/downstream of **8 real distinct municipal WWTPs** across one Grand River watershed in 2022/23 and measured microplastics in actual aquatic insect larvae, emergent adults, and riparian Tetragnathidae spiders. Original paper reports approximate downstream relative concentrations: **Heptageniidae larvae +2.4x, chironomid adults +3.7x, riparian Tetragnathidae +2.0x**, varying by site and taxon. Most larvae→spider trophic magnification factors are **not consistently >1**, and the actual constituent is **microplastic particle**, not nitrogen, phosphorus, natural host-food abundance or carbon energy.

Public **original Dataverse Borealis replication DOIs** located by actually querying the Dataverse REST API (not invented):
- [`10.5683/SP4/U0MN9S`](https://doi.org/10.5683/SP4/U0MN9S): 2026 **aquatic insect–adult–riparian spider** original source and stage data, published 2026-09-11. Public-source file metadata, original columns/year/site keys need a separate verification before quantitative model.
- [`10.5683/SP3/TLDY6W`](https://doi.org/10.5683/SP3/TLDY6W): 2026 **WWTP effluent, sediment, biofilm and aquatic insect** original replication, published 2026-05-07. The Borealis REST version1 exposes **10 original file entries** with file IDs/checksums and published fields including `Effluent_Summary.tab`, `Biofilm_Summary.tab`, `Sediment_Summary.tab`, `AquaticInsect_Summary.tab` and full characterization tables. This means source `WWTP→river habitat→larval insect` contains genuinely available files, **not** only previously published summary ratios.
- The original publications are 2026 release within one basin; the two original source tables may differ by field year/WWTP identifiers. Do **not** assume a valid exact WWTP×year join, nor count 16 independent WWTP studies. Need original `WWTP_site_ID×year×upstream/downstream×taxon/life_stage` crosswalk from both DOIs.
- User repository contains an original source-only Borealis REST probe script `analysis/20261117_emergence_macro_source/borealis_grand8_wwtp_microplastics_source_search.py`, and [Borealis source search Actions](https://github.com/zuizui0223/gesui/actions/runs/44079825536) reached the source-index stage; a final display-only JSON TypeError occurred **after original JSON was saved**. Fixed probe adds direct second DOI and `.tab` original field schemas; no ecological effect measured until original source grade passes.

## Strict next data gate for cross-boundary macro pathway inference

1. **Original exposure:** `WWTP_ID`, year, river network upstream vs downstream must be observed in **both** published Borealis source files. Real effluent particles/L is only within-source exposure. All 8 distinct WWTP sites from one river network, not eight independent countries/basins.
2. **Real source→donor→receiver:** original measured site/year material concentrations per unit in effluent, sediment/biofilm, aquatic larvae, adult emergence, riparian spider. Check data file original units separately: effluent particles/L; sediment particles/kg; biofilm particles/g; insect/spider particles per mass or count. **Do not multiply across noncomparable units.** Compare *within-stage* paired downstream−upstream log ratios, with 8 facility units and source-specific taxon/season support; a site's absence for one stage is missing data, never a zero.
3. **Before opening any effect signs**, freeze a per-plant observed stage support ledger and a core source eligibility of at least **6/8 original plant×same-year pairs** with larval, adult and spider measured exposure contrasts. If fewer, STOP whole-chain model, retain site-/taxon-specific source facts as descriptive. If all missing source coordinate IDs and same-year/matrices, STOP and do not substitute BQE EQR/GBIF/NDVI.
4. A direct conceptual nutrient transfer proxy is **Bow stable-isotope enrichment** (food-web assimilation), whereas Grand River Borealis records directly measured **microplastic movement**. They can triangulate biological plausibility **but must not be pooled as a single numeric substance effect**. Neither addresses France122 where original 2018 positive upstream reported WWTP P exists at **0/122** source sites, nor proves butterfly host use.

The paper's cautious new *macro* claim would be: emission channel and taxa/guild identity condition which aquatic-to-riparian endpoints respond. This requires source-matched independent WWTP data, holding known original nulls and spatial gradients. A positive paper-level effect is NOT presumed.


## 6. Bow River 2022 原本同位体から実際に再計算した地点単位の数値（原著の追再現、独立発見ではない）

原著コミット`32b36abfcf31117ac1bd9d7f3fdbc17cdcedd795`をcheckoutしてSHA固定の`StableIsotope2022.csv`全184試料をsource QCした後、2022年の原著微生物叢メタデータから上流`Unexposed`（3地点）・下流`Exposed`（2地点）を照合し、**地点ごとの**δ¹³C・δ¹⁵N中央値を集計。原著の試料は独立したWWTP処理場の反復数ではない。

[実行済GitHub Actions #44080153960 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/44080153960)、元著者CSVからのderived結果artifact ID `14326522225`、完全コード `analysis/20261117_emergence_macro_source/bow_2022_isotope_site_level_3up_2down_exact_perm.py`。

| 本当に観測した元sample_group | 上流3地点平均の地点中央値 δ¹⁵N | 下流2地点平均の地点中央値 δ¹⁵N | 下流−上流 δ¹⁵N | 下流−上流 δ¹³C | site-label exact両側p（δ¹⁵N） |
|---|---:|---:|---:|---:|---:|
| Aquatic Larvae / 水生幼虫 | source artifact | source artifact | **+1.8167‰** | **+2.8500‰** | 0.1 |
| Aquatic Adult / 羽化成虫 | source artifact | source artifact | **+0.4083‰** | **+3.7083‰** | 0.8 |
| Spider / 河岸クモ | 8.5917‰ | 10.4750‰ | **+1.8833‰** | **+2.7542‰** | 0.4 |
| Spider δ¹⁵N minus matched source-site biofilm δ¹⁵N | — | — | **+2.5000‰** | — | 0.3 |

Spider の真の原著5地点中央値 δ¹⁵N（‰）：`Cochrane 8.650, Cushing Bridge 8.650, Sunalta 8.475`（上流3）、`Graves Bridge 12.750, Policeman Flats 8.200`（下流2）。**下流2地点で一貫した強い同位体濃縮ではなく、一方の下流サイトが大きい**。Original biofilm median δ¹⁵N（‰）も Cushing5.6、Sunalta5.0、Cochrane2.5、Graves8.3、Policeman−0.8 と地点差が大きい。局所biofilm基準で差を取っても、source-only `down−up` δ¹⁵N spider–biofilm差 +2.50‰は5サイトexact p=.3で、栄養移行量の定量・下水以外の交絡の除去にはならない。

すべて同一Bow River沿いの地点ラベルであり、元3up/2downの配置を入れ替える10通りの**数学的帰無の記述値**。最低1側p=.1、二側pも細かい推定は不可能。下水はランダム割付されていない。したがって**元論文の局所的な同位体的整合性を再現した**という意味が最も正確で、元2026 FR122の下水負荷陽性0/122というサンプル選択の欠落を補う「別物の全球8処理場」や独立の栄養移行メタ解析にはしない。

**次に必要**なのはBorealisの2026年8WWTPについて、著者が公開したsource DOI `SP3/TLDY6W`と受け手 DOI `SP4/U0MN9S`の真正original `WWTP_site×year×upstream/downstream×taxon/life-stage`各表の欠測/単位を照合すること。主対象をmicroplastic物質移行と明記し、N/P栄養量へ換算しない。
