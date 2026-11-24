# Post-sync follow-up: early-only reference, country transfer, genus composition

Exploratory posthoc verification after the original snapshot was published at commit `305d0936943136aa886c486740cd5dd2502f2bb8` / PR #153. No new field samples. See `PROTOCOL_EARLY_REFERENCE.md` (fixed locally before the additional fits). The early-only reference uses 227 early H7-periods/89 H6, 500 nonzero host profiles, 278 family-frequency-geography strata, and the original 499 permutations; the same reference is used in both periods. The early-only null still includes early records from heldout basins and countries, so the forecast is retrospective/transductive. Country holdout applies to regression estimation, not the whole reference construction.

## Main result

For the quality cohort 173 H7/77 H6, the EQR predictor added after the previous network index and fixed geography/effort covariates predicts the plant-side change term by +13.658% H6-out MSE gain [6.000,21.835] with the early-only null, compared with +14.663% [7.123,22.235] with the original full-period null. The consumer term is -1.070% [-3.579,-0.166]. All 222 repeated H7 give +9.483% [3.513,16.336] on the plant side. These are descriptive fixed-H6 error intervals, not independent causal effect sizes. All fitted EQR coefficients for the quality plant-side response are negative.

Country-heldout regression with every shared H6 purged (no country fixed effects), early-only reference, quality panel: FR +12.043%, IT +16.503%, LT +10.894%, country-equal gain +11.706%. All 222: FR +10.775%, IT +22.469%, LT +10.789%, NO +15.523%, country-equal gain +13.336%. Only eligible countries are shown (fixed >=10 H7, >=5 H6, training >=30 H7 and >=15 H6). These are conditional descriptive transfers, not fully independent country-level studies or strict source-heldout prediction.

## Nested taxonomic decomposition

All 990 host keys have consistent exact catalogue family and genus (71 families, 437 genera). Three terms add exactly to the plant-side change: family mass, difference of genus/family mass terms, and within-genus identity. Mass refers to expected host species inclusion after record rarefaction, not abundance.

| EQR MSE gain, quality panel | Full-period null | Early-only null [fixed-loss 95%] |
|---|---:|---:|
| Family mass | +5.034% | +6.113% [1.567,10.500] |
| Within-family identity (combined) | +9.610% | +6.946% [0.166,14.506] |
| Genus mass within family | +6.439% | +4.968% [-0.207,10.512] |
| Within-genus identity | +1.312% | +0.592% [-2.289,3.961] |

No robust within-genus species-specific prediction; the combined within-family composition signal survives but its genus component is uncertain with the early-only reference. Do not turn this decomposition into measured functional traits or direct causal mediation. The official original full EIV six-trait CSV could not be downloaded in this session; grassland adaptation, agricultural use, alien status, life form and nitrogen fixation remain untested. Do not infer them from plant genus.

## Same physical French river measurements

Fixed 73 H7/24 H6 and 263 original stations: early-only reference, measured TP added to plant-side change -2.303% [-5.866,1.441]; same-station EQR +3.193% [-2.561,8.609]. The broader aquatic–terrestrial correspondence does not establish a specific wastewater or measured river P link.

40 model fits, independent Ridge predictions max absolute difference 1.832e-15; decomposition identity max error 3.331e-16. Preserve all country results, both reference variants, all nonpositive results. Exact code and derived outputs in this directory; `FOLLOWUP_MANIFEST.json` lists byte identities. To replay, after restoring original archives and the preserved previous script/panel, run `python gesui_early_reference_falsifier.py <workdir>`.

Earlier snapshot numerical replay: `LOCAL_REPLAY_VERIFICATION.log` confirms 16/16 original bytes and all eight output tables numerically identical. Original host-exclusion generator was unavailable; the earlier aggregate results and reconstructed panel were archived verbatim. No raw source rehosting.
