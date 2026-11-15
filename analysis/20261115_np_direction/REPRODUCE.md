# Replicate annual reported wastewater-source upstream–downstream direction check

This directory is an exploratory source-index negative control for 454 genuine H7×period records. The station-year EEA original source population, baseline H6-fold definitions and HNS/plant/EQR outcomes remain unchanged. It is NOT a causal experiment or measured TN/TP mediator model.

Use the frozen original content, not a new recent GBIF re-fetch. Inputs are retrieved from GitHub sources:

- Original EEA exact biology station/H7 mapping: Actions #43165326244 (SHA256 `a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89`).
- EEA original 2018 and 2022 site-level reported WWTP source fields, Git blobs `0f30bfbc3dac13f6a182c8a910e9899c04778322` and `325f199e985d9862ffee5b0edb205c61a7be3833`.
- Original annual upstream-NP source-only artifact: Actions #43222864099, source JSON SHA256 `012ad79feabc8c11e0d74f803730b5b7913b96a64a6af26fb0c97263d9aa7c26`.
- Frozen 454 H7 biological cohort original Action #43215754868 and original host-HNS panel source hash `1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca`.
- Source-direction-recovered artifact: Actions #43679382449, full artifact digest `sha256:c07815e270f7f48f6b24629490b4827e31c4448ae421afc47f9ff0bf85d65e90` (note that this digest is the artifact ZIP, not the CSV itself).

For local replication see code checked in as `analysis/20261115_np_direction`, plus the prior code from `analysis/20261115_reported_np_hns/scripts/source_join.py` and `analysis/20261115_multilayer/scripts/multilayer.py` already in the parent source branch. GitHub Actions automatically fetches the two source-only artifacts and re-fits the H6-out models.

Run `python -m unittest discover -s analysis/20261115_np_direction/tests -v` for original unit tests. The source prediction code also re-fits 27 tested models via a separate sklearn implementation. All reported MSE gains and intervals are conditional descriptive diagnostics and should never be reinterpreted as biodiversity, insect biomass, larval feeding, actual adult puddling or WWTP-specific treatment effects.

No individual third-party EEA station value CSV, full host original roster, raw GBIF occurrence or third-party original nutritional indicators are redistributed in this PR. Only derived, auditable statistical source indices are committed with original source hash receipts.
