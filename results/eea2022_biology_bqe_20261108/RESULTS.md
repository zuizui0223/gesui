# Observed aquatic biological condition versus official reported wastewater N/P (2022)

Original real result: https://github.com/zuizui0223/gesui/actions/runs/42432384196 .
Original source-only 2022 river-EQR site registry and deterministic geographical holdout folds: https://github.com/zuizui0223/gesui/actions/runs/42432030617 .
The scientific data and all effect directions are preserved unchanged; INDEPENDENT_AUDIT.json recomputes all four fixed models for P and N using a separate implementation and compares MSE from the same original 3-degree geographic blocks.

## Actual EQR, independent source support
- River benthic macroinvertebrate normalized EQR: **2,475** sites from 13 reporting countries, 45 spatial blocks, including **147** sites with 2022 measured reported EEA source N/P on directed upstream reaches. Original 2022 source-only geography/biology rows were 3,282; 807 Polish source rows are absent from the original numeric response panel, important coverage limitation.
- Phytobenthos normalized EQR: **2,784** sites, 12 countries, 44 blocks, **110** exposed. Original source geography/biology rows were 3,663; 802 Polish and 77 Austrian source rows are absent from original numeric panel.
- Fish EQR: **1,243** source-linked numeric sites across only four reporting countries, **fails preregistered >=5-country ecology inference threshold**; no outcome model fitted.

## Source-frozen heldout ecological prediction, descriptive only
Original fixed five spatially held-out regional folds, Ridge(alpha=10), baseline country/classification method + latitude/longitude + river flow + upstream drainage area. Compare upstream 50km measured EEA P or N source mass to a negative-control physically DOWNSTREAM WWTP source mass at the same sites.
- Benthic invertebrates: original baseline MSE **0.058338**, plus upstream reported P mass **0.057695 (1.10% better)** and N mass **0.057523 (1.40% better)**. Downstream placebo is **0.058541 (0.35% worse)** for P and **0.058590 (0.43% worse)** for N. Original Spearman EQR vs source mass rho ~**−0.072**, suggesting EQR is *lower* at high observed source load on average, not that wastewater is beneficial.
- Phytobenthos: baseline MSE **0.035934**, upstream P **0.035876 (0.163% better)** and upstream N **0.035871 (0.176% better)**; both near negligible. Downstream placebo does not improve; raw rho ~+0.0026.
- **N and P source masses are nearly perfectly collinear on full panels (mostly joint source absence), and strongly correlated even among exposed sites**. N vs P driver cannot be separated by this source-only both-positive EEA sample.
- The independent audit reports positive geographic blocks/countries, 5000 fixed-block bootstrap draws, source-coverage attrition and load collinearity. It is a **posthoc prediction heterogeneity diagnostic**, not a confirmatory interval for causal sewage effect. P/N nutrient source loads, EQR and coordinates remain from original reports, not posthoc selected.
- The EEA biological source has national calibrated EQR methods; OneHot encoding of country + country-method on training folds mitigates *some* scaling differences but cannot remove source-network and land-use confounding. River morphology, agriculture, other discharges, data reporting and seasonal flow may account for a small prediction increment.

## Scientific conclusion
**Positive exploratory physical direction-specific prediction for invertebrates; no demonstrated N-vs-P specific mechanism and no sewage-caused biodiversity improvement/deterioration.** An actual causal wastewater effect needs independent agricultural/urban sources, actual treatment load mass/dilution through the receiving reach, standard biological BQE quality classification, and region-robust replication. The separate original gesui EU 2022 ambient river TP/TN concentration source prediction null, positive single-river TP gradients and chocho island/port/butterfly host positives are preserved, not erased.
