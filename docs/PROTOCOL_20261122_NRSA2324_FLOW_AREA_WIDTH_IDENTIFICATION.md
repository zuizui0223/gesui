# NRSA2023–24 modeled flow and upstream area stress test — frozen before new numeric effects (2026-11-22)

## Why this is necessary
Original [EPA 2023–24 survey design](https://www.epa.gov/system/files/documents/2025-08/nrsa-2023-24-survey-design-documentation-20220209_0.pdf) explicitly states `MEAN_FLOW` is **modeled mean annual stream flow (cfs)** from NHDPlusHR, and calls Strahler order a potentially poor proxy for stream size; river-size sampling strata were based on modeled flow. Original SiteInfo also includes `WSAREASQKM` (source upstream area, km²), `FLOW_SIZE` (design size class) and `URBN_NRS23` (2021-based urban stream label). The previous positive wet-channel-width result controlled stream order but did not directly control modeled mean annual discharge.

**Question**: Does actual field-measured channel width `XWIDTH` add independent heldout information after NHDPlusHR modeled flow and drainage area are already accounted for? A positive width signal after these variables is evidence of **additional predictive information**, not of an intervention or mechanical cause. The fully observational biological outcome and the concern were chosen after prior results, so this remains **POSTHOC**.

## Fixed source gate, before opening new effect directions
Use exact SHA-pinned official NRSA2023–24 SiteInfo and physical habitat metrics, original UID/SITE_ID/VISIT_NO/DATE_COL match. Inspect only availability and source covariate distributions, NOT MMI/EPT numeric outcomes. Report missing/sentinel/positive `MEAN_FLOW`, `WSAREASQKM`, `XWIDTH`, and `URBN_NRS23` by BERW/WADEABLE and BETB/BOATABLE separately; report log10(width) versus log10(modeled flow) Pearson correlation. Source GO per method if the original outcome-complete cohort retains >=250 visits and >=100 HUC8 for primary tests. Do not pool wadeable and boatable. Do not treat missing flow as zero.

## Numerical model once source gate passes
Reuse exactly prior five original EPA source files and numeric original source checks, with original four response strata: BERW MMI/EPT and BETB MMI/EPT. Compare **same complete-case observations within each condition** and train/heldout HUC8 groups. Controls always use:
- geography (latitude, longitude), original state/eco9, original Strahler order, year and season, original measured conductivity;
- actual in-stream measured TN and TP, positive source QA flags and original unit conversions;
- measured canopy (BERW `XCDENMID`, BETB `XCDENBK`);
- source type indicator `URBN_NRS23` where its valid categories are present.

Test paired nested heldout Ridge(alpha=10) predictors:
A. B + log10(modeled mean annual FLOW) + urban + N/P + canopy versus A + log10(measured channel width).
B. B + log10(upstream watershed AREA) + urban + N/P + canopy versus B + log10(measured width).
C. B + both log10(flow) and log10(area) + urban + N/P + canopy versus C + log10(width).
D. Same C, but additionally include flow-size design class `FLOW_SIZE` categorical (when present) and original in-stream measured sediment embeddedness `XEMBED` where valid; compare only on the resulting identical complete-case cohort. D is exploratory sensitivity, not a gate-rescue or confirmation.

E. On the **same complete-case cohort as C**, additionally compare `B+modeled flow+watershed area+urban+canopy+width` versus the same plus **actual measured TN and TP**. This tests the previously observed wadeable-versus-boatable nutrient contrast after controlling river hydrological scale, rather than assuming the earlier contrast is intrinsic to sample method.

Primary is A and C; B is a second size axis, D explicit physical-sediment sensitivity, E is the nutrient-after-hydrologic-size falsifier. All additions are fixed here before reading the new numerical results. Preserve fixed five whole-HUC8 disjoint folds based on original hashed HUC8; train-only standardization and category encoding; HUC8-equal heldout MSE gain and 2,000 fixed out-of-fold whole-HUC8 bootstrap descriptive 95%. Also show heldout MSE of base+size without width, diagnostic Pearson correlation width~flow, numeric cohort n UID / physical SITE_ID / HUC8, and descriptive sign/uncertainty for all four strata. **All failed and null groups must appear**.

## No overclaim
`MEAN_FLOW` is modeled long-term flow and not measured low-flow water chemistry discharge on the sampling day; `XWIDTH` is a physical width proxy rather than a manipulable restoration treatment. `WSAREASQKM` source metadata and missingness must be audited before interpreting it. The wadeable/boatable sample-method distinction is inseparable from real hydrological and ecological size differences in a cross-sectional survey. No wastewater effluent exposure, sewage downstream reach dilution, river-to-land nutrient movement or causal restoration effect is being tested in this sensitivity analysis.
