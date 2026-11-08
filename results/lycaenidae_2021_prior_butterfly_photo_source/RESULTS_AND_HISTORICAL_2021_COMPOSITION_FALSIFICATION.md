# Does the previous year's butterfly photographic composition explain the 2022 wastewater-source P association?

**Date:** 2026-11-08. **Status:** complete source-backed *posthoc* public photographic reporting sensitivity; this is NOT standardized butterfly abundance, adult puddling, a sewage field experiment, or a contemporaneous phosphorus concentration/chemical mechanism.

## 1. Existing 2022 source-landscape positive, retained as a prior finding

Previous [Draft PR #25](https://github.com/zuizui0223/gesui/pull/25) original 2022 Germany/Austria iNaturalist **photographic reporting composition**, EEA official measured-reported annual N/P source mass and independently original 2021 ESA WorldCover categorical built/crop/country land-area adjustment found **P-source+13.67% 2° region-heldout MSE improvement** on **43 WorldCover-quality source-positive 1° geographical cells** (DE28, AT15). However earlier within-2° source label shuffle null had 35% equal/better improvements and 22-block bootstrap P MSE gain range was wide [−14.18%,+33.04%]. This is at most exploratory spatial prediction of photo reporting composition, not a causal wastewater–Lycaenidae effect. More precise local 1–5km river matches previously failed sampling coverage gates and cannot be rescued by broad 1° regions.

**New hypothesis (formulated AFTER the original positives were observed):** prior 2021 Lycaenidae *photographic recording composition* and 2021 three-family observer opportunity may account for apparent 2022 P-source geographic predictiveness.

## 2. Fixed historical cohort and new public source records

On the **exact original 43 source/WorldCover-supported 2022 photo locations**, source-key and country bbox fixed before querying prior year, obtain 2021-04-01 to 2021-10-31 iNaturalist v1 public research-grade georeferenced photographed observations (open geo, ≤500m, uncaptive) for descendant family taxa **Lycaenidae 47923**, **Pieridae 48508**, **Papilionidae 47223**; official country places DE7207, AT8057. No photos, observers or individual observation IDs downloaded; **129/129 public source API counts succeeded**. All 43 source IDs and original 2021 ESA land-cover validity are preserved.

- Original **43** fixed 2022 sites have a pooled **2021 photo reporting total 4,358**: Lycaenidae **2,237**, Pieridae **1,768**, Papilionidae **353**. These are source public photo *records* and not individual butterflies.
- Previously locked 2021 support criterion ≥15 three-family photo records at each of ≥30 sites, with DE≥18, AT≥9 and ≥12 source 2° regional blocks: **42/43** cells qualify, DE**27**/AT**15**, **18** 2° blocks by source request's count; note the final modeled country×2° block IDs are 22 because source ID grouping was counted with country distinct. Original 2022 sample remains unchanged; new historical supported **42-country×1°** cells are clearly a distinct matched subset and compared only to themselves.
- The source data are the original same-year 2021 historical public photographs, **not** independent randomized butterfly absence, nor a year-specific source of discharged P. Original 2022 WWTP may have been in service in 2021 or decades earlier, so 2021 is not an exposure-free control year.

[All original 2021 public source receipt URLs and checksums, with stage-A sample support](2021_HISTORICAL_3FAMILY_PHOTO_COUNTS_ORIGINAL43_SOURCE_ONLY.json) | [2021 original photo source-only workflow SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42465101895).

## 3. Historical stability of the observed Lycaenidae family photo fraction

In the 42 matched same original source cells (DE27/AT15) compute per-cell family photographic log-odds `log((Lycaenidae+0.5)/(Pieridae+Papilionidae+0.5))`. Original 2021 vs 2022 country×1° photo-logit values co-vary:

- **Pearson r +0.5143**, nominal exploratory p ≈0.000494.
- **Spearman ρ +0.5107**, nominal exploratory p ≈0.000548.
- Unweighted mean original 2022-minus-2021 logit change only **+0.0212** across matched regional cells, with country-specific means **Germany −0.0246**, **Austria +0.1038**. This is observed photo-reporting composition, not population stasis or ecological no-change.
- 2021 supported **three-family photo records** in the same original fixed geographic cohort: DE2,801, AT1,555; 2022 same selected sites DE3,321, AT2,393. This highlights uneven reporting across seasons/years and cannot be used as abundance trends.

## 4. Strict same-42-cell original 2° geographic heldout prediction

No extra source/photo queries in stage-B. Full original 2021 public land-cover and 2021 botanical record reporting, lat,lon,country used as baseline; source is original annual 2022 EEA *reported measured* unique-FW WWTP P mass, within-country training-standardized; exact unweighted same-42-cell 2° block-left-out folds, no hyperparameter tuning.

| Response and predictor baseline | Baseline heldout MSE | +EEA 2022 P source MSE | Added-source MSE improvement |
| --- | ---: | ---: | ---: |
| **2022** Lycaenidae photo family logit ~2021 WorldCover urban/crop/land + 2021 Plantae effort + lat/lon + country | 0.158674 | 0.142221 | **+10.37%** |
| **2022** Lycaenidae photo logit, same controls **plus 2021 Lycaenidae photographic composition and previous-year three-family photo reporting total** | 0.155528 | 0.149778 | **+3.70%** |
| **2021** prior-year Lycaenidae photo logit predicted from 2022 EEA reported source mass and same original environmental controls | 0.258837 | 0.280651 | **−8.43% (worse)** |

Country breakdown for **2022 after 2021 historical butterfly composition/target-group effort**: DE27 cells **+1.37%**, AT15 **+8.75%**. The photo-historical source P information **shrinks substantially** after accounting for preexisting local butterfly reporting composition and butterfly-targeted reporting effort. The source does not reliably reconstruct the prior-year 2021 photo composition using a *future-year* reported mass value, but **this does not refute the existence of prior wastewater treatment**.

**Scientific interpretation:** Much of the 2022 registered-P source geographical photo-predictive signal is compatible with **persistent year-to-year habitat/observer/taxon composition**, not a newly identified 2022 sewage-specific phosphorus mechanism. The remaining +3.70% is a small posthoc descriptive conditional photo-reporting prediction. No direct WWTP source-to-butterfly-exposure link, actual biologic butterfly abundance, adult puddling protein/Na, native/introduced larval host quality or ant symbiosis was measured.

[Source-frozen complete 2021→2022 historical photo CV GitHub Actions SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42465437203) | [full same-42 region predictions and historical controls](2021_PRIORYEAR_LYCAENIDAE_COMPOSITION_VS_2022_P_SOURCE_BLOCK_CV.json).



## 4b. Explicitly posthoc geographical prediction-error uncertainty

A separate follow-up was frozen **after** observing the +3.70% historical-adjusted improvement, and resampled **entire original 22 country×2° geographical heldout test groups** 1,999 times from the exact saved original 42-cell predictions (seed 20261108). The original predictions, photo/source data and original 2021 support selection were never modified. These are **descriptive spatial resamples of correlated heldout prediction errors**; they are NOT intervals for a causal P effect and do not replicate independent monitoring datasets.

| Model: source-P incremental 2° heldout MSE improvement | Observed gain | 1,999-block-resample descriptive 95% interval | Fraction source MSE gain ≤0 |
| --- | ---: | ---: | ---: |
| 2022 with original WorldCover/2021-Plantae/geography baseline | +10.37% | **[−27.85%, +31.60%]** | 23.76% |
| **2022 after 2021 butterfly photo composition + target-group photo count** | **+3.70%** | **[−31.10%, +28.27%]** | **38.12%** |
| 2021 historical photo composition predicted by the original reported 2022 source-P | −8.43% | [−20.99%, +1.00%] | 95.65% |

The key new uncertainty result is that **+3.70% historical-adjusted gain is NOT regionally precise or clearly different from zero in this descriptive geographic error resample**. The result should be described only as a small, geographically uncertain posthoc apparent increase in prediction for public butterfly *photo reporting composition*. A 2021–2022 causal source-phosphorus effect is not established.

[Regional bootstrap reproducible GitHub Actions SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42465827581) | [full frozen original-block error resamples and three-model uncertainty](2022_P_SOURCE_AFTER_2021_LYCAENIDAE_HISTORY_2DEG_ERROR_BOOTSTRAP.json).

## 5. What could actually distinguish mechanisms?

Original 2022 reported effluent mass is neither actual organic P/N input to the 1° terrestrial area nor a butterfly mud-puddling substrate resource. More observational controls alone cannot replace:
- **Coordinator-approved standardized butterfly per-visit sampling and observer effort**, matched to original 2022 source-positive and source-free reporting geography; privacy-safe eBMS/UFZ/Viel-Falter source-only feasibility requests remain UNSENT in [Draft PR #23](https://github.com/zuizui0223/gesui/pull/23).
- **Actual location/year-specific nutrient and water chemistry**, treatment timing, outfall-to-habitat connectivity, water quantity, dilution, landscape built/crop and larval host/ant ecology.
- For life-stage divergence, actual adult substrate drinking behavior and larval growth/survival at matched sites must be observed. Reconstructed `chocho` **22/25 Lycaenidae introduced host *potential resource geography*** is neither direct local host feeding nor sewer exposure and stays a separate positive.

The previous explicit 2018/2022 aquatic EQR out-of-block P/N predictive improvements are separate water-community findings. **Retain both positive photographic spatial-information results and evidence weakening the sewage-specific interpretation; do not blend them into a single causal ecology claim.**

## Provenance
- [Original 2022 EEA source 1156 Germany/Austria public FW outfalls](https://github.com/zuizui0223/gesui/pull/23).
- [Original 2021 WorldCover country-clipped built/crop and 2022 photo source P falsifier](https://github.com/zuizui0223/gesui/pull/25).
- [Historical 2021 photo source protocol, written before opening 2021 observations](../../docs/LYCAENIDAE_2021_VS_2022_HISTORICAL_PHOTO_REPORTING_POSTHOC_LOCK.md).
- [All original source photos counts/receipts](2021_HISTORICAL_3FAMILY_PHOTO_COUNTS_ORIGINAL43_SOURCE_ONLY.json).
- [2021 history-vs-2022 region-heldout source models](2021_PRIORYEAR_LYCAENIDAE_COMPOSITION_VS_2022_P_SOURCE_BLOCK_CV.json).
