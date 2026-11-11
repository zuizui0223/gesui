# Same-river riparian vegetation source-first coupling test (NRSA 2018–2019; posthoc after nutrient–EPT)

**Source/metric definitions fixed before examining newly joined riparian measurements.** EPA's original 2018–2019 NRSA physical-habitat full-metrics CSV [official EPA catalogue](https://www.epa.gov/national-aquatic-resource-surveys/data-national-aquatic-resource-surveys) and [original metadata](https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_physical_habitat_larger_set_of_metrics_-_metadata.txt). Source original URL:
`https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_physical_habitat_larger_set_of_metrics_-_data.csv`

Same EPA native unique site–visit `UID` must join to frozen original NRSA 2018–2019 benthic BERW + TN/TP source manifest, SHA256 `fbf53b8379d4e6564dbb51b3705b02fdd03e4a2caacad06a2314d643c8412e91` (source frozen before seeing nutrient and EPT effects). Keep sample grain `UID` and *physical site* `SITE_ID` separate, never invent same-river field values.

### New source-only endpoints locked before new numeric habitat values

1. `XGveg` — field-estimated riparian **ground-layer vegetated areal cover** (not botanical richness, butterfly larval hosts or flowering resources);
2. `XPCAN` — reach proportion with **canopy present** (not plant species richness);
3. `XFC_OHV` — stream channel areal fraction with **overhanging vegetation under 1 m above water** (potential terrestrial-to-aquatic habitat/subsidy proxy; not measured litter, insect emergence or pollination).
4. `W1H_PIPE` — proximity-weighted index of observed pipes entering/leaving water; **human disturbance only**, *NOT* confirmed municipal WWTP effluent;
5. `W1_HALL` — distance-weighted overall human disturbance indicator, a potential confounder.

**Unit QA:** verify all metadata parameter names exist exactly in original CSV. Values must not be inferred from labels or visually sampled; treat 0 as an observed zero when valid, not missing. For XGveg/XPCAN/XFC_OHV require 0–1, for W1H_* nonnegative; if measured metric units differ, stop. Duplicated UID rows are ambiguous until original sample-type metadata are reconciled without inspecting outcomes.

### Phase 0 — before new physical habitat numbers

Re-download exact public original CSV, record SHA256, headers, sample types, per-UID duplicates, source-matched native UID support across **1,221 original UID visits**, state and ecoregion counts, share nonblank for XGveg/XPCAN/XFC_OHV. No raw vegetation numeric values are printed or modelled.

**Fixed source support gate:** >=650 distinct eligible source-matched UIDs with all three primary riparian metrics nonblank; >=25 states, >=6 `AG_ECO9`; all selected join UIDs unique; source-matched original `SITE_ID` and visit agree. If fail STOP; do not rescue by lowering support, changing field metrics, or selecting river types after effects.

### Phase 1 — if Phase 0 passes

Freeze original UID member list, source SHA and declared metric meanings before reading riparian numbers. Then on same original river visits, compare EPT fraction and predeclared fractionally allocated FFGs against TN/TP and riparian metrics, adjusting for state, eco9, latitude/longitude, stream order, year, conductivity and original NARS flags. Include a *descriptive* `log10(TP) × XGveg` interaction to test whether EPT nutrient association varies with riparian ground vegetation; not causal buffering nor measured service. Use original SHA256-HUC8 disjoint 5-fold and holdout reference model, plus original source flags and site-unique sensitivity; do not claim out-of-sample generality where folds fail.

Do not connect to unknown plant species, flower visits, adult insect emergence, wastewater source loads, or ecosystem service magnitude. This is a **measured habitat–aquatic functional association at one place/time** and an important intermediate empirical bridge only. Any result is **exploratory** because source test was designed after original river nutrient–EPT signs were known. No use of eBird.
