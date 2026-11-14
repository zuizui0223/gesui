# Three EPA NRSA survey cycles: independent-period, zero-shared-HUC8 transport of ambient nutrients versus physical habitats for benthic EPT

**2026-11-14 | Original US EPA source data rebuilt three times, zero HUC8 training–test overlap, source-locked before new results, independent NumPy matrix audit SUCCESS. DRAFT observational macroecology, NOT causal wastewater or ecosystem-service evidence.**

## Actual original data and source identity

Reconstructed directly from previously publisher-SHA-verified and source-frozen original raw EPA CSVs in one GitHub Actions runner rather than copying coefficients from earlier READMEs:

| Authentic original NRSA survey cycle | Original source-valid native survey UID visits | Distinct original HUC8 | States | Benthic EPT operational response |
|---|---:|---:|---:|---|
| 2008–09 | **586** | **437** | 48 | wadeable **BERWW** original `TOTAL300` counted EPT share, directly reconstructed using publisher `TAXA_ID→ORDER`; >=95% of counted individual original ORDER linked per UID; conservative known-EPT/ALL-original-count lower bound |
| 2013–14 | **1,116** | **706** | 48 | original **BERW** wadeable `TOTAL300` EPT fraction, original source sites also have measured periphyton Chla |
| 2018–19 | **946** | **603** | 47 | original **BERW** wadeable `TOTAL300` EPT fraction, original source sites also have measured periphyton Chla |
| **Combined real source native UID visits** | **2,648** | **1,030 distinct union HUC8**, not sum 1,746 | — | Original visit fraction, **not species richness, benthic biomass, adult emergence or ecosystem service** |

Native **HUC8** is not a decimal numeric measurement. Original EPA 2008 raw HUC8 codes sometimes omit the leading geographic zero (7 digits), and 2013 raw codes are sometimes prefixed by `H` (e.g. `H10270102`); [source-first QA](https://github.com/zuizui0223/gesui/actions/runs/43566511111) only normalizes literal 7/8 digit codes or `H`+8 digits to a canonical exact eight-digit code, preserving original raw source-length/prefix counts and enforcing **no collapse of distinct original HUC8 values**. HUC IDs are NEVER inferred using nearest sites, counties or geocoded buffers. Source original survey samples remain locked to prior real EPA source panels and native original UID/SITE_ID.

**Real overlap across original cycles**, after canonical identity:
- 2008 vs 2013: **260** HUC8.
- 2008 vs 2018: **203** HUC8.
- 2013 vs 2018: **419** HUC8.
- Shared by all three: **166** HUC8. Total distinct union **1,030**.

All three original source eco9 fields resolve to the SAME nine original categories: CPL, NAP, NPL, SAP, SPL, TPL, UMW, WMT, XER. These are comparable labels, but source method BERWW vs BERW and sample year differences remain, and 2013/2018 samples were restricted by *periphyton Chla source eligibility*, unlike 2008. This is an **unweighted observational quality-conditioned survey subcohort**, not nationally weighted river-mile condition or a single identical historical panel.

### Fixed identical source predictor contract, no leakage

[Protocol fixed BEFORE reading any new external transfer MSE](https://github.com/zuizui0223/gesui/blob/analysis/nrsa-threecycle-ept-crossperiod-huc8-20261114/docs/NRSA_THREE_CYCLE_EPT_AMBIENT_NP_PHAB_HUC8_TRANSPORT_PROTOCOL_20261114.md):

- Environmental baseline B: original latitude/longitude, **native documented Strahler order** (2008 source '1st' to 1 through '7th' to 7; ambiguous 8th+ excluded), field sampling date season sin/cos, log10 real field conductivity, original state and original nine-level eco9 (one-hot **trained only on fitting cycles**).
- **True physical field habitat PHab**: original streambed embeddedness XEMBED, mid-channel canopy density XCDENMID, riparian canopy existence XPCAN, overhanging vegetation XFC_OHV, original human disturbance W1_HALL, log10 original wetted stream width XWIDTH. Source 2013 XGV and 2018 XGVEG were deliberately NOT substituted for one another.
- **True ambient stream nutrient N+P**: measured TN and TP source-aware converted to log10 mg/L (2008 CHEMW NTL & PTL originally µg/L; 2013/18 NTL mg/L and PTL µg/L) with each original survey's disallowed CHEM flags excluded; not identified as municipal WWTP effluent loads.
- Four fixed nested Ridge α=10 models: B, B+NP, B+PHab, B+PHab+NP. No year as predictor of an *unseen* era. No posthoc feature selection or tuning after result signs. Scaling and eco/state vocabulary fit **only on train**.
- Each training HUC8 receives **equal total regression sample weight** regardless of multiple native survey visits; primary validation outcome **equal-weight original HUC8 MSE**, not a summed-visit precision claim. Equal visit training and test error are separately frozen sensitivity analyses.
- For each heldout ORIGINAL cycle, REMOVE from all other original cycles **every training UID whose native canonical HUC8 appears anywhere in the heldout target cycle**. Test original UIDs retained, with no watershed leaking between training and test. Target test sample HUC8 bootstrap 2,000 block resamples of *fixed heldout errors* yields descriptive95, not standard errors from independent retraining.

## Main real external survey-period result: all three test periods

[Original three-cycle direct source/EPT transfer SUCCESS #43566511111](https://github.com/zuizui0223/gesui/actions/runs/43566511111). **Interpret the two conditional model gains separately**: physical added AFTER actual nutrients, and actual nutrients added AFTER physical.

| Fully untouched survey-cycle TEST | Original target UID /HUC8 | Train original UIDs /HUC8 AFTER removing target HUC8 | PHab added AFTER NP, HUC8-equal MSE gain [2,000 block bootstrap95] | NP added AFTER PHab, HUC8-equal MSE gain [bootstrap95] |
|---|---:|---:|---|---|
| **2008–09** | **586 /437** | **1,199 /593** | **+7.4519% [ +1.585%,+13.270% ]** | **+2.5863% [ +1.058%,+4.024% ]** |
| **2013–14** | **1,116 /706** | **444 /324** | **+7.9790% [ +3.616%,+12.495% ]** | **+0.1992% [ −2.152%,+2.640% ]** |
| **2018–19** | **946 /603** | **675 /427** | **+11.1927% [ +8.341%,+14.129% ]** | **+0.8880% [ −0.735%,+2.518% ]** |

**Exact primary source test HUC8 equal MSE** (source EPT 0–1 share squared units):

| Heldout ERA | B | B+NP | B+PHab | B+PHab+NP | Train-only intercept reference |
|---|---:|---:|---:|---:|---:|
| 2008–09 | 0.0420345 | 0.0398598 | 0.0378689 | **0.0368894** | 0.0534264 |
| 2013–14 | 0.0462573 | 0.0460751 | 0.0424834 | **0.0423988** | 0.0548201 |
| 2018–19 | 0.0482450 | 0.0467949 | 0.0419297 | **0.0415573** | 0.0587138 |

**Descriptive within-cycle geographic fivefold CV** on the SAME frozen original native sites, NOT the main external-period test:
- 2008: PHab after NP +6.680% [2.209%,10.759%]; NP after PHab +1.615% [−1.098%,4.223%].
- 2013: PHab after NP +9.333% [5.303%,13.140%]; NP after PHab **−0.037%** [−1.075%,0.883%].
- 2018: PHab after NP +11.160% [6.740%,15.627%]; NP after PHab +1.169% [−0.392%,2.886%].

**Important supports and weaknesses:** heldout 2013 cycle has only **444 train visits across 324 external HUC8** after strict source-HUC exclusion, making the 2013 domain-transfer model small and potentially fragile. Source state category CT/RI absent from its train and mapped to OneHotEncoder unknown-category all-zero, without test leakage; heldout 2018 train has CT/MA/RI missing. Across-cycle sampling years and BERWW/BERW source protocols differ; no repeated measures or experimental effects identified.

## Even more stringent: target HUC8 never present in ANY other cycle

This is a **subset** of each already strict zero-leak target, defined entirely by original source geography, NOT selected by favorable effects. Target-only HUC8 clusters with all other survey cycles absent:

| Original target-only HUC8 (no other survey-cycle HUC8) | Actual target UIDs / unique HUC8 | NP after PHab, HUC8 equal MSE gain [bootstrap95] | PHab after NP, HUC8 equal MSE gain [bootstrap95] |
|---|---:|---|---|
| 2008–09 | 172 /140 | **+4.696% [ +1.692%,+7.404% ]** | +3.577% [−8.796%,+15.437%] uncertain |
| 2013–14 | 248 /193 | −2.283% [−7.595%,+2.824%] uncertain | +7.802% [−0.612%,+16.101%] uncertain |
| 2018–19 | 184 /147 | +0.998% [−2.190%,+4.009%] uncertain | **+15.669% [ +10.017%,+21.086% ]** |

[Separate closed-form, manual weighted NumPy ridge normal equations, independent train-only categories, all primary model MSEs and source-HUC8 group bootstrap reproduced to numerical tolerances SUCCESS #43566835585](https://github.com/zuizui0223/gesui/actions/runs/43566835585). It also independently computed novel-HUC8 fixed-OOF bootstrap bounds (not independent ecological field samples). Original upper-bound source EPT if ALL unknown-order counted organisms were sensitive taxa and visit-equal training also leave the main PHab-vs-NP conclusion qualitatively unchanged. Exact original raw original field files were re-read from EPA original sources, their source SHA verified, but the code implementations share the same actual EPA observations, not independent replication.

## Ecological result, legitimate interpretation and competing explanations

**Evidence-backed macroscale association/predictive regularity:** after holding an entire national original survey cycle AND every matching HUC8 out of training, *measured river physical habitat* contributes useful generalizable information about the relative share of EPT counted organisms in all three survey periods. The additional **in-stream N/P** contributes independently detectable geographic predictive information on 2008 test rivers only under the fixed source-locked model, whereas the 2013 and 2018 increments have block uncertainty intervals containing zero. This is more stable across historical sampling cycles than the previously posthoc 2018 EPT×periphyton Chla×embeddedness interaction, which failed on never-overlapped independent 2013 catchments in [PR #79](https://github.com/zuizui0223/gesui/pull/79).

Potential ecological hypotheses, **NOT proven mechanisms**: substrate embedding influences interstitial habitat and assemblage detectability; physical refuge, canopy/light context and channel size condition the benthic community; chemical nutrient sources are partly confounded with runoff and habitat. But the joint physical predictor group includes **W1_HALL human-disturbance index**, width and vegetation, which are not necessarily independent riverbed mechanism variables or measured daylight. **SHAP/predictive importance is not ecological causality**; the model does not isolate a riverbed's direct effect from any other physical predictor, nor prove wastewater-origin nutrients.

### Explicit NON-claims
- Not independent experimental manipulation of sewage, stream nutrient delivery or measured municipal POTW discharge: original [PR #74](https://github.com/zuizui0223/gesui/pull/74) actual historical WWTP exact-source dual N/P report availability failed the prespecified site gates; all here is **ambient real measured river N/P**.
- No 2008 periphyton AFDM/Chla in these original publisher sources. Therefore external prediction of the relative EPT **counted fraction** is not a measured multi-trophic bottom-up nitrogen→algae→insects chain, no trophic biomass/effective insect export rate, no water purification process rate, no butterfly pollination.
- Relative EPT fraction ≠ absolute aquatic insect biomass or adult emergence; original 2008 unknown taxon order bounded, later BERW and 2008 BERWW protocols and source flags need cross-method caution.
- All HUC8-group resample intervals are fixed OOS error quantiles, **not confirmatory p-values, causal-effect confidence or independent refitting**. Distinct original 1,030 HUC8 rather than 2,648 biological independent replicates; groups may have repeat visits, sampling is unweighted by national spatial survey weights.
- No data from eBird, no TinyFish. No main merge. Keep all positive and negative old PRs #73–#80 and new PR Draft status.
