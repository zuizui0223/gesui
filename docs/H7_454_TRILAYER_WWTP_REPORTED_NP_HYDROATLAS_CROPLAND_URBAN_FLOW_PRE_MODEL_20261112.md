# Same-basin macro test: do original WWTP spatial and annual N signals survive cropland/urban/natural-flow controls?

**2026-11-12 | New ORIGINAL landuse data, numeric coefficients not yet seen; posthoc design.**

Previously observed: 232 European H7 basins /454 basin-periods, original static HydroWASTE outlet-pressure extra EQR heldout MSE gain **+7.328%** but H7 bootstrap wide and country-equal-weight EQR coefficient interval spans zero, with plant/5-family butterfly source index not comparably predictive. Original 2018/22 reported N/P load corresponds to EQR extra **+0.928%** beyond static pressure but only **56 positive matched reported H7-period sources** (38 basins); the report-source EQR gain reverses in French-only source coverage (4 source-positive of 299 France rows). All of these are posthoc exploratory, not identified pollution effects.

**Now genuinely obtained independent source:** official Figshare HydroATLAS BasinATLAS v1.0 Level7 **raw DBF** using verified Range read from original 4.276GB publisher ZIP without downloading it all; original publisher global ZIP Figshare file ID20087237, MD5 metadata `ce5013294bdcb6f884ff6cdf1a1acae0` (full compressed ZIP not downloaded/MD5 reverified), native only Level7 member CRC32 `ef18ec36` and SHA256 `1c73db5c07cb02ee6d10b963a0fc9fa4321fab55ff43721e8d6cd5c1054fb113`. Original 57,646 Level7 records and exact **1,219/1,219 study-basin unique HYBAS_ID matches**. Source-only covariate manifest SHA256 **`ce6780e01d2b98ed8e852574adb827cb3125be792a164e8be9e2ebb917b3a686`**, [Actions original true DBF SUCCESS #43225986249](https://github.com/zuizui0223/gesui/actions/runs/43225986249). Canonical native field codes:
- `crp_pc_use`: cropland percentage upstream cumulative catchment, 0–96 original integer %
- `urb_pc_use`: urban extent percentage upstream, 0–100 %
- `dis_m3_pyr`: natural modeled discharge **at H7 basin outlet** m³/s, original float 0–21,894.457
- `crp_pc_sse` and `urb_pc_sse`: cropland/urban in the **local** H7 polygon, and `glc_pc_s22` artificial landcover.
- All are STATIC original historical products, not true time-varying land use for 2018/2022. Source-values 0 are legitimate recorded zero areas but no source/measurement dates are imputed.

## Exact fixed evaluation cohort and predictors

No change to **454** original source-eligible native `HYBAS_ID_L07×PERIOD` rows, **232 distinct H7 basins** (227 P1/227 P2), 15 original study countries, shared exact outcome measurements:
1. original original EEA 2018/2022 benthic EQR station-equal mean;
2. GBIF original native-basin Plantae speciesKey index 50-draw rarefied expected richness `S50`;
3. original GBIF Lycaenidae log odds out of five accepted butterfly families.

Original EEA/GBIF country and assessment procedure, original river-site modeled `DIS_AV_CMS`, H7 local area, native original EEA site sample size, native Animalia/Plantae/five-family indexed recording opportunity, time bin. All source-relative covariates are same rows; no change in output type or fauna and plants.

**Source exposures fixed**:
- `W_static`: original modeled HydroWASTE cumulative discharged volume to the **basin outlet** per upstream area, NOT year-specific actual load or delivered dose at original internal monitoring station;
- `N_reported`: original 2018/2022 reported upstream N source tonnes/year diluted by **modeled river DIS_AV_CMS at original site**, converted to mg/L-equivalent hypothetical full source delivery and log1p, with source-matched zero meaning **no positive reported match**;
- `P_reported`: same source P equivalent, converted to µg/L-equivalent and log1p. N/P source logs correlation 0.936. **Do not claim independent nitrate vs phosphorus mechanism**.

## A priori nested models on same rows / outcomes BEFORE new land-use adjusted numeric coefficients

B0: original full country + EEA method + period + site modeled river flow + source EEA station count + basin area + GBIF three source recording opportunity counts; no wastewater source predictor.
B1: B0 + W_static (must reproduce original PR #65 numerical EQR/Plant/Lycaenidae MSE exactly BEFORE interpreting new values).
B2: B0 + **upstream cropland `crp_pc_use` + upstream urban `urb_pc_use` + `log1p(natural H7 outlet DIS_m3_pyr)`**.
B3: B2 + W_static (tests whether modeled wastewater *geography* adds predictive information after agricultural/urban and flow).
B4: B3 + N_reported (tests annual reported N source incremental value after land use and static W).
B5: B3 + P_reported (same for P).
B6: B3 + both reported N and P (source collinearity diagnostic, not causal attribution).

**Source-supported secondary sensitivity, outcome-blind**: replace upstream cropland and urban percentages in B2 with local H7 `crp_pc_sse` and `urb_pc_sse`, retaining the same original natural modeled outlet discharge; compare B2local, B3local and B4local to the corresponding upstream versions. No selecting a more significant spatial scale after seeing results.

All biological response families use the exact same original **454** H7 basin×period rows and identical original outcome missingness, with 5 fixed hash-of-H7 basin-disjoint folds (both periods of same basin in same fold), Ridge α10 (intercept unpenalized) and predictor standardization fitted solely in training. Report ALL signed out-of-basin MSE gains for B0→B1, B0→B2, B2→B3 and B3→B4/B5/B6; 2000 original H7-basin resample descriptive 95 intervals; plus H7 clustered descriptive OLS coefficient for W_static and N_reported when included, country-equal-total-weight sensitivity, warn on source collinearity.

**Stop criteria before fitting effects:** original 1219 HydroATLAS H7 source SHA matches, all 454 basin IDs mapped, original land-use percentages in [0,100] and natural discharge >=0 and finite for every row, source-pinned annual NP manifestation SHA matches, original B0/B1 baseline CV reproduced. If fails, mark source-quality stop, no guessed replacements. Original `hft_ix_s09` raw human footprint 2–388 is NOT used because scaling and overlap with urban source need independent QA.

## Scientific interpretation

BasinATLAS land cover is based on historically static baseline and may not represent 2017/22 actual land use. Cropland, urban cover and natural basin outlet discharge do not perfectly capture agricultural nutrient sources or actual receiver water chemistry; source and response reporting remain highly selective. Regression beta, even negative, is not a sewage-specific counterfactual. Unequal GBIF record observation process and France-dominated EQR assessment methods remain unresolved.

Any loss of positive extra N or W prediction after land-use adjustment is a substantive corrected result, NOT an excuse to change cohort, model scale, predictor family or response. Do not infer measured pollination, insect prey emergence or any service benefit, and do not merge any working branch/main/chocho. No eBird.
