# 2022 reported WWTP N/P masses vs measured river chemistry — EU study

Original source-frozen real results from 2026-11-07, [source geometry run](https://github.com/zuizui0223/gesui/actions/runs/42374105184), [year-aligned measured chemistry run](https://github.com/zuizui0223/gesui/actions/runs/42375542151), [independent holdout audit](https://github.com/zuizui0223/gesui/actions/runs/42376191822), [saved-fold cross-version replay](https://github.com/zuizui0223/gesui/actions/runs/42376519412). All retained files carry exact SHA256 receipts, including original source and original returned 2022 station median datasets.

**Sources:** EEA UWWTD 2022 public treatment plants and official freshwater discharge points (https://www.eea.europa.eu/en/datahub/datahubitem-view/6244937d-1c2c-47f5-bdf1-33ca01ff1715), reported **measured annual load** uwwNDischargeMeasured/uwwPDischargeMeasured in **metric tonnes/year** (both-positive WWTP selection). Peters et al. Global Nutrient Yields Figshare v5 (https://doi.org/10.25400/lincolnuninz.30546371.v5), original measured site-date TN.2 and TP.2 harmonized concentrations (mg/L) and metadata under CC BY 4.0. This repository retains only derived source→station features, matched **2022 station median** values, frozen regional folds and source citations; original large raw datewise archives remain with the public primary publishers.

## Source crosswalk
5,492 EEA 2022 plants report both positive measured N/P; plant code joined to public unique discharge point, restricted to EEA FW direct-freshwater codes; 2,969 qualifying discharge points snap to HydroRIVERS EU-derived graph <=500m. Same-reach source excluded (unknown order). Directed reach-to-reach source loads accumulated only to stations within 50km downstream along the river network. Official HydroRIVERS original has not been independently verified against the hash-locked third-party GeoParquet conversion.

## Empirical year-aligned outcome
Original year 2022 measured station records only, >=4 positive source observations from >=3 months per site. **TP n=1,291** (53 with EEA-reported 2022 upstream measured P load; 27 spatial blocks). **TN n=159** (26 upstream-measured-N loaded; 27 spatial blocks).
Raw source load vs measured in-river concentration rho: **TP -0.06044**, **TN -0.18594**, using all original stations, descriptive, not causal.
Spatial heldout 5-fold Ridge(alpha=10), same station denominator across model families: baseline lat/lon/mean riverQ/drainage area; M1 adds true reported 2022 upstream source load, M2 negative-control downstream WWTP load; all predictors held constant per fold.
- TP baseline MSE **0.721242** vs +upstream **0.735449** (1.97% worse), +downstream **0.720317** (+0.13% improvement).
- TN baseline MSE **0.786210** vs +upstream **0.812399** (3.33% worse), +downstream **0.819693** (4.26% worse).
- TP with published crop and population covariates n=446 (46 exposed): baseline MSE **0.583402** vs +upstream **0.591134** (1.33% worse), downstream **0.588351** (0.85% worse). TN crop/pop n13/exposed2 fails source support.
No new source-specific positive prediction; this does not establish no sewage ecological impact.

## Physical dilution gate, exploratory
Convert annual EEA load tonnes/year to hypothetical uniform annual mass concentration at HydroRIVERS mean stream flow: [t/yr] x 1e9 [mg/t] / ([m3/s] x 1e3 [L/m3] x 365.25 x 86400 s/yr).
Among the **53** P-source-supported stations, median potential is **0.00228 mg P/L**, median site-specific ratio to observed TP station median **4.07%**, with **14** ratios >10% and **0** ratios >100%. Among **26** N-source-supported stations: median potential **0.04866 mg N/L**, median relative ratio **2.39%**, **3** >10%, **0** >100%.
These are **NOT true fraction contributions, not rigorous seasonal upper bounds, not effluent measurements**, and exclude unreported sources. Annual wastewater source load and long-term mean river flow are incomparable to seasonal/low-flow pulses; retention, flow, agricultural flux and monitoring selection remain major uncertainties.

## Reproducibility
2026-11-07 GitHub runner NumPy 2.4.6 / sklearn 1.9.1 GroupKFold tie-handling differs from some other versions. Exact original test-block memberships are frozen in AUDIT_2022_SOURCE_REPLAY.json and permit independent cross-version replays. Script scripts/replay_eea2022_saved_folds.py is the portable verifier. All MSE values checked to 1e-7 and file bytes hashed. Tiny Spearman changes <=1e-4 can arise when 10-digit exported loads merge near-tied ranks; original unrounded association is retained in source JSON.

**Scope:** This supports observational wastewater-load model testing and realistic source-to-receptor dilution hypotheses. A robust causal study needs matched actual receiving reach water flows and full source accounting (agriculture, land-use, other wastewaters) across multiple seasons and catchments. Preserve previously positive chocho GEB/port/island resource findings and separate positive German six-river TP exploratory trends without forced integration.
