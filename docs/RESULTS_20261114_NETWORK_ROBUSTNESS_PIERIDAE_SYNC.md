# Same-H7 butterfly–host network robustness, Pieridae and transfer

Exploratory posthoc source-replay using the SAME GBIF source 454 basin×period records, 232 H7 and 90 H6 from PR86. Data are potential host opportunities computed from occurrence indices and a historical host catalogue, NOT observed herbivory.

| Exclusion of source-ranked host plant contributions | Observed link index | Spatial constrained null mean | Index excess over null |
|---|---:|---:|---:|
| None | 5.4393 | 4.6568 | +16.80% |
| Top 1 (Urtica dioica) | 4.9168 | 4.2060 | +16.90% |
| Top 5 | 4.0899 | 3.6358 | +12.49% |
| Top 20 | 2.6975 | 2.5855 | +4.33% |

Top-20 plants accounted for 50.41% of original source candidate-link index. At top-20 exclusion, observed 2.6975 falls inside spatial constrained permutation 95% interval [2.4510,2.7303]; do not cherry-pick bootstrap intervals to claim general structure confirmed.

## Five-family WWTP-geography predictor; source coverage ≥80% (n390/217 H7/83 H6)
| Potential host-support index | H6 heldout MSE gain % | H6 fixed-error bootstrap95 |
|---|---:|---|
| Lycaenidae | -2.2508 | [-9.8612,+5.0250] |
| Nymphalidae | -2.8159 | [-6.5333,+0.7269] |
| Pieridae | +12.9552 | [+3.4053,+20.1493] |
| Papilionidae | -5.1663 | [-9.4374,-1.7455] |
| Hesperiidae | +1.0391 | [-8.6973,+10.1352] |

Pieridae additional model controls for within-family recording effort and its 33 butterfly source species composition: +12.1748% [+3.8957,+17.7877] H6 MSE gain, but **six-supported-country equal-transfer +0.935% [-9.866,+9.758]**, and FR-only +0.1665% [-1.236,+1.694]. Thus no stable continent-wide wastewater effect. Full network geographic null link-excess W adds -1.124% [-3.055,+0.571] H6 MSE gain (no support).

Statistical scope: many outcomes/model choices posthoc; these are MSE predictive improvements NOT signs or percent changes in butterfly abundance, food biomass or interaction rates. W is static modelled basin-outlet HydroWASTE geography; its confounding with urban land use remains.

Original local output SHA provenance is in ../results/20261115_sync/LOCAL_ARCHIVE_SHA256.csv, linked earlier PR86 original source and local ZIP SHA256 3b37ab244cb3c7755e75a2a70498a4314ca5c885517c1c88e7017f40823e8846.

2026-11-15 new verification: 8/8 robustness source/unit tests pass; 1,709 recorded files verify with original byte+SHA manifest; no independent ecological replication.
