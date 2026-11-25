# Common-support attenuation audit and Clarke publisher-source retry (2026-11-24)

Exploratory, same underlying observations; no independent field replication. No Clarke rows were available for ecological fitting at this snapshot. Source-only retry 45386260304 queries the actual five XLSX keys in the official PMC OA inventory obtained in run 45385384325 (artifact 14915862892, ZIP SHA256 7af61933a1bd5775e4684239f6ba7aa0f498bfb87a3075a045778998fc578281). Attempt 2 receipt preserved unchanged; direct CSVs rejected (401/403), prior JATS XML returned HTML, original PMC OA file inventory lists the actual supplement keys. The publisher and Dryad README establish that Appendix S1 is a view of the accepted host catalogue; workbook content and evidence labels must be audited before use.

## Factorial audit of original source
All six original source ZIPs and original 454 network were replayed by the existing strict-reference code. Using only the **old catalogue edges**, compare 52 mapped Lycaenidae keys with the 40 accepted common keys in the prior non-Wikipedia comparison, on both the full original cohort and the previously fixed >=80% common inclusion-mass cohort. Use the same early training H6-only strata, ridge alpha 10, paired H6 losses and fixed-loss descriptive 95% intervals. No new source thresholds or outcomes chosen based on results. This is an exploratory selection audit, not a causal mechanism test.

| Cohort | Old 52 gain | Common 40 gain |
|---|---:|---:|
| quality 173/77 H6 | 7.328 | 9.455 |
| quality common 111/52 | 2.095 | 3.214 |
| all 222/89 | 8.228 | 9.843 |
| all common 146/62 | 1.771 | 2.357 |

Restricting taxa alone does not attenuate prediction. Restricting the sites used for both fitting and scoring does. This does not identify a biological coverage effect.

## Fixed all-eligible training / stratified holdout scoring
Fit the same full-eligible H6-disjoint models, then score their fixed predictions on retained or excluded H7 rows. H6-equal averaging is separate for each subset. In common40 quality, retained +11.546% [-1.951,20.943], excluded +7.646% [2.634,15.146]; all retained +11.897% [4.426,20.211], excluded +5.533% [-1.733,14.463]. The full-trained model carries information even in the retained subset; re-training on restricted support changes performance. Shared-H6 comparisons (15 quality /20 all) have wide intervals crossing zero. This does not distinguish sample size, geography, covariate distribution, changing training-defined reference, or model fitting effects. The retained cohort has lower mean early EQR (quality 0.640 vs 0.693, all 0.623 vs 0.681; raw H7-weighted descriptive means). Coverage involves both periods, so is not a prospective sampling rule. Do not treat retained vs excluded as independent groups or p-values from difference of significance.

Both scripts' independent ridge prediction checks agree to 1.34e-15. Code, every result and loss, coverage summary and shared H6 contrasts are archived. No independent catalogue result, no wastewater causal effect, and no new actual TP or field validation.
