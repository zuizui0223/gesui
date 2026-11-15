# Source-only protocol: Within-H6, between-distinct-H7 comparison of plant N and butterfly HNS

**2026-11-15; started AFTER PR #102 geographic positives, temporal nulls and PR #103 five paired EIV controls. This is a new posthoc adversarial design, not preregistration.**

### Question
Does the original plant-N / butterfly-HNS *spatial* covariation persist between different genuine H7 watersheds **inside the same H6 parent**? Existing H6-out transfer mostly identifies *between-H6* geographic differences, and 2017–19→2020–22 Δ within the same H7 gave no temporal synchronization. We must separate between-region association from within-parent heterogeneity before writing "linked terrestrial ecological response".

### Strict source-only opportunity gate before outcome
- Original SHA-pinned panel `1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca` 454 H7×periods/232 H7/90 true H6.
- Verified source-v2 NUTRI exact-score eligibility: 451/454 original H7×period; 4,011 numeric matched plants, no hybrid score imputations. Use unchanged original 451 H7×periods.
- Only distinct H7 with **both** 2017–2019 and 2020–2022 NUTRI-supported periods (nominal 219 H7/86 H6) to avoid observation-quality asymmetry and avoid treating periods as independent geographic units.
- Only H6 parents with **at least two distinct eligible H7 watersheds**, not two periods of one H7. Freeze lower bound: ≥100 distinct eligible H7 watersheds, ≥40 genuine distinct H6 parents, ≥5 source majority-country groups. If not met, STOP and report the exact support failure. Never lower minimum based on outcome.

### Models if source gate passes
1. For each eligible H7, take mean over exactly two periods of its source plant species-equal NUTRI, butterfly species-equal HNS, W (must be static across periods), and original source/effort/trait coverage covariates.
2. H6-parent fixed intercepts and source W, urban/crop/hydrologic geography, survey effort, NUTRI/host trait completeness. **Omit country×method dummies if nested/collinear with H6 FE**. Estimate PlantN β standardized per source H7 within cohort with H6-clustered descriptive95.
3. Include alternative record-weighted HNS in a strictly subordinate check. Only original 2-period site means, no temporal-causal inference.
4. Separate H6 cluster block residual permutation: permute under the fixed no-PlantN model among **distinct H7 within each H6** to obtain a descriptive null, retaining the original H7 mean and source H6. Use 9,999 draws seed=20261115. This is an approximate exchangeability test, not an experiment or an H6-out prediction. 
5. To avoid comparing in-sample H6 FE to out-of-H6 prediction, label FE as *within-H6 conditional association* only. Report source H6 cluster SE, 9,999 blocked residual-label test, leverage by H6, trait scoring heterogeneity, and rank diagnostics. A null cannot prove no local effect.

### Claim rules
- If within-H6 persists, ecological sorting is not only between broad H6 regions. This remains a GBIF-recorded trait composition signal with shared EIV scoring sources; no measured discharge N.
- If within-H6 does not persist, the strong H6-out geographic prediction is potentially driven by broad regional plant/butterfly turnover, not watershed-scale local sewage-linked mechanisms.
- Neither outcome is true wastewater-to-riverbank nutrient flux, observed host ingestion, butterfly abundance or ecosystem-service production. Original temporal non-synchrony remains a separate negative result.
