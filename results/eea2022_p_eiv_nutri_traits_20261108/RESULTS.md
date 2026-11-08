# EEA2022 measured P-source matched river plant nutrient affinity — audited

**Original EIV run:** https://github.com/zuizui0223/gesui/actions/runs/42378520227

The frozen 18 pairs / 36 sites / 1,006 original GBIF Plantae speciesKeys were independently crosswalked to European Tichý et al. (2023) EIV NUTRI. Exact accepted vascular species-level binomial matching only; no genus/fuzzy substitutions. **911 accepted vascular** GBIF speciesKeys, **624 numeric exact EIV NUTRI matched**, and **25/36 sites** pass original >=10 rated taxa and >=40% matched taxa threshold. Only **9/18 whole pairs** pass, below original mandatory >=10 support gate.

**No prespecified positive or negative wastewater vegetation inference is possible.** Within the nine descriptive eligible pairs, source-positive minus source-zero species-equal mean NUTRI = **−0.2552**, 3 higher vs 6 lower (posthoc exact sign p 0.5078); a record-weighted sensitivity has mean −0.1383. Neither is causal evidence that wastewater lowers plant N preference.

**New independently verified missing-data result:** the source-positive member has LOWER EIV matching fraction in **15/18** pairs (exact two-sided sign p **0.00754**, a posthoc *data-quality* diagnostic, not the preregistered environmental effect). Source-positive has more plant GBIF record counts in 13/18 pairs (posthoc sign p 0.0963), while site-level rated-species quality is 14/18 source-positive and 11/18 controls. Uneven GBIF effort and trait-missingness can change which whole pairs are eligible. See INDEPENDENT_EIV_QUALITY_AUDIT.json for all pairs and source reproduction.

The 2022 EEA source P mass is a measured facility annual load (t/year), NOT observed soil or river-edge P at GBIF plant coordinates. Source-zero controls are only *reported-source-zero* within a 50km river graph. EIV NUTRI is an ecological species indicator, not effluent N concentration or plant tissue phosphorus. No causal sewage N/P→plants mechanism is established.

To preserve data reuse without unnecessarily republishing a third-party trait catalogue, original numerical EIV-by-species source remains at its own published URL/DOI. This repository retains source-verified GBIF taxonomy IDs, exact match flags, per-site EIV summaries, original source output hashes and the independent offline audit, **not** the source full EIV table nor the 624 species-level original numerical NUTRI scores.

This project remains independent of original chocho port/island butterfly host opportunities, which retain their positives unchanged.
