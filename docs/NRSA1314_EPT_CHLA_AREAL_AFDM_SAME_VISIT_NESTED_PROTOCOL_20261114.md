# Same-visit 2013–14 source-first AFDM–periphyton chlorophyll–EPT standing-stock comparison

Date 2026-11-14. Source and first cross-cycle n1116 Chla × XEMBED effect **already inspected**. This next same-sample comparison is explicitly exploratory/posthoc. Do not reinterpret as independent pre-registration.

**Fixed original input:** n1,094 original probability-sampled NRSA 2013/14 wadeable `BERW` station UIDs, real `TOTAL300` standardized EPT fraction, exact same UID water total N/total P and original site measured `PCHL CHLA_RESULT` in µg/cm², `PBIO AFDM_RESULT` in mg/cm², original NARS lab quality flags no disallowed code, original physical habitat XEMBED, XCDENMID, XPCAN, XFC_OHV, W1_HALL and XWIDTH with source native UID and visit/date checks. These 1,094 sites were selected solely by measured-source quality and unit availability, not EPT effect sign. Different 2013 AFDM `MG` (not `MG/CM2`) source records expressly excluded. The current data are not sewage-derived nutrient dose.

Fixed cross-biota tests on all same 1,094 EPA UID:
- EPT counted *fraction* ~ geography/state/eco9/year/season/stream order/field conductivity/observed TN/TP plus real PHab; add original measured Chla alone, areal AFDM alone and both. Whole HUC8-disjoint 5-fold Ridge10, same exact folds and samples for every model.
- HUC8-cluster descriptive OLS original standardized +1SD log10 AFDM and +1SD log10 Chla conditional on each other and all geography/river/physical/nutrients. Report sign, CI and effect magnitudes separately.
- Nested predictive gains `base→+AFDM`, `base→+Chla`, `+Chla→+both`, `+AFDM→+both` all with fixed heldout HUC8 2000-resample descriptive95. Report negatives/nulls.
- Secondary measured ambient nutrient geography: log10 AFDM or log10 Chla ~ geography + PHab baseline then +logTN+logTP, same native EPA sites; HUC8 OOS gain and partial TN/TP coefficients with cluster95.
- No mechanistic mediation estimated, because sampled same-visit observations are nonrandomized and do not identify nutrient origin or intervention, and adult insect emergence/biomass or trophic services absent. EPT proportion is not total insect biomass, AFDM standing stock ≠ net primary production, Chla pigment ≠ autotroph mass flux. In this posthoc derived context, Chla/AFDM ratio may be described but is not independent field replication.
- Spatial sampling is unweighted probability-site *subset*, not nationally weighted population estimate. Report raw sample n and geography; no source sign-selected data or mixed BERW/BETB.
