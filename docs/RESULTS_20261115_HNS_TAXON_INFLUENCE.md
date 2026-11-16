# HNS geographical association: taxon influence and conditional robustness

15 November 2026 · original sourced, **exploratory POSTHOC**. This sensitivity uses the **same 454 European H7 basin×period rows (232 H7, 90 H6)** and does not estimate environmental causal effects or a true species-extinction response.

## Source, metric and ranking

- Source original panel H7 SHA256 `1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca`. All original species GBIF source records and host HNS species traits were previously recovered and verified in Draft PR #87 / originating Actions #43593227162. Independent source species × period long-format index is included as a compressed **derived** record-level table, not original photographed specimens.
- HNS is the mean of the **habitat nutrient indicator of historically documented larval host plants**, weighted by **GBIF butterfly source occurrence records**. It is NOT leaf N, adult sodium/protein intake, butterfly abundance, actual larval feeding, or an experimentally observed wastewater effect. Source historical taxon score is fixed, with missing labels excluded from the denominator, never zero-imputed.
- Deletion order is ranked **after observing** the original non-geographic HNS conditional W β: absolute additive species contribution, using the SHA-pinned original `species_slope_decomposition.csv`. This is **post-result cherry-selection for influence stress-testing**; a deleted taxon is removed from the definition of the response index, not from an actual ecological community. Source counts, country, period, family source shares, effort covariates and H6-group memberships stay unchanged.
- All models add original static HydroWASTE basin-outlet geography W after recording effort, urban/crop and flow, country, method and second-order spatial location terms. Clustered original observational W coefficient is distinct from H6-group-heldout predictive added value. H6 bootstrap intervals are descriptive fixed-error resamples, not formal causal confidence bounds.

## Exact source result (454 H7×periods; source score-eligible original species only)

| Posthoc excluded high-contribution butterfly keys | Adjusted HNS association per +1 SD modelled W [H6 descriptive 95] | W addition H6 heldout MSE improvement [descriptive 95] |
|---:|---|---|
| 0 | **+0.10802** [+0.03888,+0.17717] | +5.374% [-5.901,+14.780] |
| 1 | +0.11149 [+0.04668,+0.17630] | +4.034% [-10.769,+15.578] |
| 5 | +0.09948 [+0.03930,+0.15965] | +2.962% [-13.561,+16.148] |
| 10 | +0.07118 [+0.01104,+0.13131] | **−1.876%** [-15.037,+8.256] |
| 20 | **+0.04440** [+0.00029,+0.08850] | **−2.428%** [-14.681,+6.019] |
| 30 | +0.04761 [+0.00732,+0.08790] | +11.135% [+2.812,+18.107] |

**Interpretation:** The original positive HNS association is not wholly one species, but becomes much smaller after excluding the top 10–20 posthoc influential taxa. Out-of-H6 predictive gain changes sign across choices, even rebounding for 30 removed taxa. This is instability/response-definition sensitivity and cannot be interpreted as species displacement, pollution response, or selection on larvae. Inferences from the selected removal count are not multiplicity-corrected or independent validation.

## Reproduction and numerical QA

- Six recalculated candidate HNS indices retain 454 original rows, 232 unique H7 and 90 H6; zero basin×period receives a made-up zero denominator.
- Original unablated HNS exactly reconstructed from the frozen 454×303 GBIF original speciesKey record matrix; maximum absolute numerical discrepancy < 1e−10.
- Independent `pandas.groupby` on the **original long-format species×H7×period source table**, not the optimized NumPy species-by-basin matrix, reproduces all six HNS index definitions on 454 groups; original unablated max error 1.78×10−15. This is different code on the SAME source, not an independent ecological study.
- Four synthetic contract tests pass; full source validation also checks all six fitted W betas against signed original result and no changed H7/H6 grain.

## Remaining required falsifiers

- Keep trait database completeness/missing-host uncertainty as a separate major source of sign fragility (existing adversarial plausible coefficient bound includes zero).
- Verify any response pattern against true wastewater **received concentration/dose**, not the static outlet geography W. Compare actual wet reference/stormwater sources. Use species-level adult puddling and host larval occupancy observations to distinguish adult mineral intake and host plant pathways.
- Do not concatenate unrelated studies into a pseudo-mediation SEM. Aquatic EQR, plant S50 and butterfly HNS cross-realm coresponse is distinct from measured insect biomass and consumer food transfer (see stacked PR #88).

Derived host affinity source attribution: Tichý et al. (2023), *Journal of Vegetation Science* DOI `10.1111/jvs.13168`, Zenodo `10.5281/zenodo.7427088` (CC BY 4.0). Ferrer-Paris et al. historical butterfly larval host compilation DOI `10.6084/m9.figshare.1168861.v1`. GBIF species index from original frozen source retrieval run; see PR #87 / #86 for original acquisition and validation.
