# Post-interaction taxon mechanism test, SOURCE frozen before decomposition outcomes

**Date 2026-11-14. NEW exploratory posthoc after seeing [PR #78](https://github.com/zuizui0223/gesui/pull/78) embeddedness×Chla association. Not independent confirmation.**

Question: the observed EPT *counted fraction* from a standardized benthic sample is compositional. The same relative-fraction association with periphyton chlorophyll may arise because other insect groups increase, not because EPT biodiversity/actual abundance declines. Decompose counted community *on the same original 946 NRSA2018/19 UID visits* without changing selected samples or PHab model.

Three pre-fixed distinct **source-observed** outcomes:
- `chironomid_counted_fraction = sum(BERW TOTAL300 positive counts where native FAMILY = Chironomidae) / original total counted subset`; this is not absolute in-stream density or mosquito biomass.
- `log1p(EPT native taxon-ID richness)`: unique native TAXA_ID entries with positive counts in original Ephemeroptera/Plecoptera/Trichoptera orders, with original total counted-subsample number controlled.
- `log1p(nonEPT native taxon-ID richness)`: distinct original native TAXA_ID entries with positive counts in other orders, with source total counted-subsample number controlled. Original taxon IDs are not species IDs.
- Background given reference: EPT original counted fraction, whose previous full n946 embeddedness interaction was β+0.0192, Holm5 p≈0.0108 and geographically heldout +0.74% MSE gain of uncertain reliability.

Fixed observational model for EACH new response: original state/eco9/geography/year/season/Strahler order/conductivity, original ambient TN and TP, all seven PHab continuous controls, original areal logChla, original Chla², and centered 1SD logChla×centered 1SD XEMBED. For richness outcomes include log(total source BERW counted subsample size). HUC8 cluster descriptive OLS interaction and 5fold whole HUC8 heldout alpha10 Ridge; compare quadratic reference Q vs Q+embeddedness interaction on identical source cohort, 2000 HUC8 fixed OOF error bootstrap. Apply Holm3 correction across the **three** new chosen outcomes. Report all nulls and no invented insect dry mass.

Potentially useful negative/control logic: a chironomid positive interaction or EPT richness null suggests a relative-composition mechanism, not observed insect density loss. Conversely taxon-ID richness interactions would indicate a taxonomic component, still without causal attribution. The decomposition is observational: wind, nutrients, habitat, sampling and detection may confound. The original DMR historical sewage N/P gate [PR #74](https://github.com/zuizui0223/gesui/pull/74) fails; no wastewater-specific cause, adult insect export, riparian predator flux or ecosystem service is assigned.

Pre-numeric sample guard: exactly the source-verified 946 native UID /603 HUC8 from [PR #77](https://github.com/zuizui0223/gesui/pull/77) and raw official EPA benthic count file SHA256 `4e9b3a0798b25f67c7a96d4800b94f86544bbd623abaf119358e067b158324a7`; reject if any original EPT count fraction disagrees with independently reconstructed old panel. All quantities are same original EPA source rather than other literature's summarized estimates.
