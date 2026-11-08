# 2021 historical butterfly photo-composition check of a 2022 source P prediction (posthoc)

**Created 2026-11-08 BEFORE querying any 2021 butterfly-family observations for this test. However this is emphatically POSTHOC** because the original 2022 Lycaenidae photo reporting signal and the ESA WorldCover2021 adjusted WWTP P predictive gains have already been seen. Do NOT present historical records as independent treatment assignment or randomized placebo.

## Core unresolved question

Does the apparent 2022 EEA WWTP reported P source × Lycaenidae photo reporting share relation reflect an enduring prior-year pattern of butterfly photos, their habitats, tourist/urban access and observer effort, rather than a year-specific sewage source nutrient response? Public **2021** photographs can be compared with the **same pre-existing 2022** source×country×1° cells, but **a 2022 WWTP reported source may already have operated in 2021**, so this is a temporal recording-pattern diagnostic, not a wastewater-no-exposure negative control.

## Fixed input sample, before reading 2021 butterfly photo outcomes

- Take the **exact original 43 cells** in DE28, AT15 that were 2022 3-family-photo-supported **AND** independent original 2021 ESA WorldCover country-land source quality supported on [Draft PR #25](https://github.com/zuizui0223/gesui/pull/25), with 22 country×2° geographical regions. Do NOT replace outcome-poor cells with other regions after reading 2021 results. These 43 source cells previously satisfied original 2021 nonfocal Plantae reporting effort threshold; their 2022 Lycaenidae/Pieridae/Papilionidae counts are already source archived.
- Query exactly **43 × 3 = 129** iNaturalist v1 official aggregated `total_results` for **2021-04-01 through 2021-10-31** in identical country×1° bounding rectangle and official Germany place 7207 / Austria 8057, family descendant IDs Lycaenidae=47923, Pieridae=48508, Papilionidae=47223. Retain `photos=true`, `geo=true`, `geoprivacy=open`, `taxon_geoprivacy=open`, `quality_grade=research`, `captive=false`, `acc_below=500`, `per_page=1`.
- No photo pixels, individual species/observation IDs, researcher/observer names, nor sensitive transect positions retrieved. Original aggregated public total count only. Preserve API source URL and response SHA256 for all 129 calls, including zero and missing status; absent public observations are not true ecological absence.
- Before any new historical-result model, require **≥15** 2021 three-family photos in **≥30** of the **43** already frozen original source cells, including **≥18 DE** and **≥9 AT**, and **≥12** original country×2° geographic blocks, with **all three family queries successful**. These are new previously frozen 2021 photo-support thresholds not to be changed after observing 2021 values.

## Fixed response and blocked validations (only if new prior-year support gate passes)

- Source classes remain **2022** EEA original uniquely both positive measured-reported N/P mass with one FW outlet per plant, country×1° public source intensity, not direct nutrient doses. Original 2021 ESA WorldCover v200 built/crop/land area country clipping and 2021 Plantae photo counts retained exactly.
- On **identical 2021-supported country×1° cells**, create 2021 response and 2022 response as `log((Lycaenidae_photos+0.5)/(Pieridae_photos + Papilionidae_photos+0.5))`, and paired change `Y_2022−Y_2021`. These are *changes in reported photographic composition*; **not** changes in population size or butterfly use of effluent.
- Predeclared (for this POSTHOC test) descriptive 2021/2022 photo-log-odds Spearman and Pearson correlations, country-specific paired shifts and unweighted grid means. Photo opportunity in 2021 may be less than 2022 due observer sampling growth; never interpret recording count difference as population change.
- Predict 2022 photo-log-odds with fixed **leave-one country×2° block out** training-only feature scaling. Baseline **A0 = country, latitude/longitude, log1p(2021 Plantae photographs), 2021 WorldCover built/crop country-land fractions and log country terrestrial land area**; A1 = A0 + within-country standardized 2022 reported WWTP P source log1p. Replicate exact 43-model only on the NEW 2021 photo-support subset, so the sample remains the same for nested comparisons.
- Historical-adjusted baseline **B0 = A0 + Y_2021 + log1p(2021 three-family photo count)**, and B1 = B0 + same original 2022 WWTP reported P source. Compare MSE B0 vs B1 on the *same original 2021-supported sites* and geographic 2° folds. If `B1` does not add predictive information, the previously observed P-source geographic signal is compatible with **historical stable photographic composition and/or targeted observer effort**. If B1 remains positive, this still only predicts 2022 photo reporting and is not a 2022 effluent mechanism.
- Separate static-history test: predict Y_2021 with the original A0 baseline and add **2022 registered P source**. If future source metadata predicts the **prior-year** 2021 photo composition, report temporal-spatial continuity explicitly. Do not call this "no sewage in 2021" because WWTP plants likely existed before 2022.
- Exact same predictors/training-only preprocessing, unweighted 1° cells, rank testing and no new feature selection, no tuning. Also report 2021 three-family total photo-record count and cell P source rank association by DE and AT; P/N collinearity and urban/built/host/rural confounding remain.
- Do NOT compute fitted models if 2021 photo-support gate fails. Only source+data availability results may be reported.

## Interpretation firewall

- The original `chocho` 22/25 Lycaenidae added *potential* introduced larval host resource geography remains separate and unchanged.
- The original `gesui` 2018/2022 river-benthic EQR predictive gains, original 2022 local 1km 1/24 and 5km 4/17 butterfly photo support failures, and public national Norway/Sweden transect near-WWTP source failures remain separate.
- No direct juvenile/ant survival, adult protein/sodium puddling, 2021/2022 field density, actual 2021 source P discharge history or exposure at butterfly sites can be inferred.
- All outcomes will be explicitly labelled **posthoc observational photo-reporting temporal sensitivity**, regardless of sign or fit.
