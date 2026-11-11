# Same-river water chemistry × benthic EPT × riparian vegetation: outcome-locked exploratory bridge (2026-11-12)

**Purpose.** Test a measured, same-place **aquatic–riparian ecological coupling** rather than presenting global sewage × GBIF photos, US rivers and Europe EEA as one observed chain. No new river sample choice based on habitat/EPT outcomes. POSTHOC relative to previously seen NRSA TN/TP–EPT results.

**Already frozen sources, before this study's riparian numerical outcomes**
- Original 2018–19 EPA NRSA benthic BERW + TN and TP + site source manifest: SHA256 `fbf53b8379d4e6564dbb51b3705b02fdd03e4a2caacad06a2314d643c8412e91`, exact 1,221 UID original visits.
- Original NRSA `XGVEG` ground vegetation cover; `XPCAN` canopy presence; `XFC_OHV` overhanging vegetation; `W1_HALL` field-observed human disturbance; `W1H_PIPE` pipe influence/effluent proximity **NOT a WWTP identifier**. Official 2018–19 EPA PHab source SHA256 `69b3b56962c04492423bd7eaf4d074c6ae8d253924d683ba12eb25fd857ce063`.
- Phys-hab native source-only original UID pass: 1,220 matching 3-cover-field visits / 48 states / 9 eco9; frozen before new numeric PHab values SHA256 `23cf9546b751d1eaedb9aeefcdb63daed1edf14d41bae9f551240252c90f2012`, source run [43048691795](https://github.com/zuizui0223/gesui/actions/runs/43048691795).
- 2018–19 original EPT result `TP partial β −0.04775`, n1,220, HUC8 heldout joint TN+TP +2.35% vs geography alone. This **known result** is why bridge was designed; never claim a preregistered new prospective EPT hypothesis.
- GBIF Lycaenidae corrected original key 5473, original 360 fixed river site 5km 2017–2021 source support **only 9 joint sites with plant≥50 and Lycaenidae≥20**, so **NO butterfly mediation or pollination model** under the existing thresholds. [SOURCE correction](https://github.com/zuizui0223/gesui/actions/runs/43048245515).

**Authoritative definitions:** EPA [PHab metadata](https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_physical_habitat_larger_set_of_metrics_-_metadata.txt). Ground, canopy and overhanging vegetation are structural cover, **not plant species richness, flowers or butterfly larval food plants**.

## Frozen same-UID model family, before inspecting numeric PHab coefficients
Outcome: original BERW `EPT = (Ephemeroptera+Plecoptera+Trichoptera) counts / TOTAL300 subsample counts`. Neither adult emerging insects, aquatic insect biomass, or a quantified ecosystem service.

Predictors / filters:
- Original UID in both independently frozen source manifests, strict original SHA files, same native site identity and 2018–19 time. Exclude invalid/absent values, do not replace missing with 0. Validate XGVEG, XPCAN, XFC_OHV in [0,1], human disturbance and W1H_PIPE nonnegative. TN `mg N/L`, TP original `ug/L` then /1000 to `mg P/L` before taking log10.
- Joint complete case: year, lat/lon, stream order, conductivity, TN/TP, ground, canopy, overhang and W1_HALL; must retain >=650 visits, >=25 states, >=6 eco9. No post-result dropping countries/folds or adjusting radii.
- Original five spatial HUC8 folds by SHA256(huc8) modulo 5, Ridge α=10 training-only standardization of features and no intercept penalty, all models on exact same cases. HUC8 group bootstrap 2,000 descriptive resamples for **out-of-fold** site mean MSE gain; original HUC8 are correlated by river basin and not random treatment units.
- Baseline **B**: latitude + longitude + year + stream order + log(conductivity) + state and eco9 + W1_HALL. **N**=B+logTN+logTP. **H**=B+XGVEG+XPCAN+XFC_OHV. **NH**=B+TN/TP+3 riparian variables. **NHI**=NH + logTP × (XGVEG−0.5). NHI is **one prespecified exploratory effect modification**, not a mediation estimate.
- Evaluate ΔMSE B→N (increment nutrient); B→H (habitat); N→NH (extra habitat after nutrients); H→NH (extra nutrients after habitat); NH→NHI (vegetation modification). **Report all positive and negative** improvements as they are, plus leave-one-state-out slope sign and clean-flag / >=250 count / one-site-one-visit sensitivity when feasible.
- OLS partial slopes with HUC8-cluster descriptive 95% intervals for N, NH, NHI, and interaction. No claim of causal buffering because riparian cover, chemistry, and EPT are measured at the same time and share confounders (light, substrates, hydrology, land use, geographic clustering).
- If an effect improves predictive fit but 95% block resample crosses zero or country/region transfer fails, report **exploratory** evidence only.

## Layer-by-layer scientific boundary
1. **WWTP ↔ benthic EQR Europe**: measured source mass has small 2018 and 2022 out-of-region prediction gain, NOT tested in EPA US site nutrient source; no causal WWTP dose.
2. **Measured stream nutrients ↔ EPT US**: genuine same-visit negative EPT fraction vs TN/TP (not organism number or emergence).
3. **Riparian structural cover ↔ aquatic EPT US**: this planned new same-UID model can establish descriptive **co-organization** only. This is the only new empirical edge directly measurable now.
4. **Aquatic EPT ↔ adult insect export ↔ predators**: emergence and predator ingestion not measured; no mediation or ecosystem service.
5. **Riparian flower/host plants ↔ Lycaenidae ↔ pollination**: no same-site standardized flora, butterflies, interaction or seed-set. GBIF source gate fails; no causal chain can be fitted.
6. **Port wastewater generated ↔ GBIF Plantae/host photorecords**: separate global *recording opportunity* signal, reduced by rarefaction and Animalia/Fungi counts; never pool with US river/Europe EQR as sequential data.

Keep all original source and result files unchanged. New files only in a separate Draft branch stacked on `analysis/aquatic-terrestrial-ffg-gbif-support-audit-20261112`; no `main` or `chocho` changes.
