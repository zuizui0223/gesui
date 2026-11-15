# HNS taxon influence audit — exploratory, original-data only

Local/CI reproduction for 2026-11-15 posthoc HNS source-score influence. **No causal SEM, field effect, species extinction model, or wastewater dose is inferred.** The input matrices are derived aggregate GBIF speciesKey source occurrence counts and documented larval host plant habitat N indicators, **not raw GBIF occurrence IDs**.

The base original H7 panel is stored in stacked Draft PR #87 (`results/20261115_sync/hns_macro/H7_454_HNS_panel.csv`). The model implementation is the unchanged `analysis/20261115_hns/scripts/model_core.py` from PR #87. No source-only controls are recomputed for each taxon exclusion; this is sensitivity to the response definition only.

Source EIV attribution: Tichý et al. (2023), Journal of Vegetation Science, https://doi.org/10.1111/jvs.13168; versioned original https://doi.org/10.5281/zenodo.7427088 (CC BY 4.0). Host links: Ferrer-Paris et al., https://doi.org/10.6084/m9.figshare.1168861.v1. Source GBIF API and source SHA metadata in Draft PR #85–#87.

Data source hashes are checked before any fit. `inputs/species_records_with_HNS.csv.gz` is an aggregated index, not source taxon specimen records; species-specific adult diet, sex, body condition and larval feeding are unknown. External validation means H6 geographical holdout of the SAME study, not a new field campaign.
