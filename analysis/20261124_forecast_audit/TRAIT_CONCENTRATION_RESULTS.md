# Species concentration and host-topology audit — 2026-11-24

Posthoc diagnosis of the unchanged third-period forecasts from run 45394379127, not a new field sample or model fit. All eight source archive hashes checked against the original receipt; independent hypergeometric calculation reconstructs 297 species-level structured-null loss rows, max difference 9.36e-17. Original quality 143 H7/70 H6, common Lycaenidae 133/67. The Clarke 46-key topology table is computed directly from the fixed expanded Clarke graph (status <=3, blob 1bb8c445d9c83e523000effb0baddf940013be52) and rounded to eight decimals; its Git blob 64c8c0f2d681b4596a7f0264d6501e0e23c25ba9 matches the local table exactly. Source key lists are unchanged.

## Species-count standardization
For each family and graph, sample 10 taxa uniformly without replacement 5000 times, rank positive absolute MSE reductions vs structured null, compute top contributor share and effective contributor count (inverse HHI). Exclude no-positive draws (three Nymphalidae draws); these are descriptive conditional resamples, not a formal between-family effect test. Papilionidae has exactly 10 species so its value is constant. Each species loss averaged over H7 then equally over H6. All 205 mapped accepted original species; common panel 46 taxa each graph. Prior outcome-derived species ranking is not used for selection or model refit.

| Family / catalogue | n | positive | top share all (%) | median top share 10 (%) | draw 2.5–97.5% | effective count median 10 | baseline loss top share median 10 (%) |
|---|---:|---:|---:|---:|---|---:|---:|
| Hesperiidae | 26 | 13 | 27.84 | 54.92 | 34.13–88.53 | 2.489 | 43.25 |
| Lycaenidae | 52 | 33 | 14.70 | 42.52 | 24.85–74.16 | 3.477 | 42.83 |
| Nymphalidae | 97 | 45 | 15.44 | 61.99 | 34.31–98.75 | 2.095 | 44.74 |
| Papilionidae | 10 | 6 | 82.19 | 82.19 | constant | 1.455 | 53.96 |
| Pieridae | 20 | 10 | 38.44 | 57.38 | 37.43–89.11 | 2.335 | 36.26 |
| Common Lycaenidae, old | 46 | 29 | 16.46 | 41.58 | 23.55–71.52 | 3.627 | 38.46 |
| Common Lycaenidae, Clarke | 46 | 31 | 16.50 | 48.15 | 29.18–80.55 | 2.991 | 38.54 |

Smaller concentration in Lycaenidae persists descriptively after fixing species count, but draw distributions are broad, depend on sampling frame and baseline error. Original five-family Lycaenidae median concentration is similar to its baseline-error concentration; Papilionidae is concentrated even before host addition and more concentrated for positive host gain. No claim of statistically demonstrated distributed ecological mechanisms or interaction redundancy.

## Matched species consistency across host catalogues
For 46 common Lycaenidae, Spearman correlation of signed per-species absolute loss reductions old versus Clarke: +0.73136 vs structured null, +0.70583 vs baseline. 24 species improve under both structured-null comparisons, 10 worsen under both, 12 discordant. Top five contributor species are the same (different order): Lysandra coridon, Cyaniris semiargus, Plebejus argus, Lysandra bellargus, Elkalyce argiades (original source canonical labels). This is a descriptive shared-outcome correlation, not independent field replication or significance test. Degree rank correlation +0.58981, Jaccard +0.72253, indicating catalogue differences.

## Trait mechanism tests
Within-family standardized ranks of signed absolute loss reduction regressed on rank of prior mean and prior SD of expected inclusion, plus log-order-equivalent catalogue degree and mean pairwise host Jaccard to other modeled species of same family. Both degree and Jaccard included jointly; fixed family intercepts. 5000 genus-cluster resamples within family, fixed ranks, conditional on historical/heldout fixed predictions and errors. Six primary structured-null comparisons (two features × three panels) Bonferroni intervals; baseline sensitivities also reported locally. No explicit phylogeny, no actual host physiology.

| Panel | n / genera | Degree beta [Bonf6] | Jaccard beta [Bonf6] |
|---|---|---|---|
| Five families, old | 205 / 98 | +0.0893 [-0.0909,0.2506] | -0.0121 [-0.1848,0.1886] |
| Common Lycaenidae, old | 46 / 32 | -0.3016 [-0.7893,0.1846] | -0.0031 [-0.4055,0.5161] |
| Common Lycaenidae, Clarke | 46 / 32 | +0.0653 [-0.5896,0.8501] | -0.1334 [-0.6684,0.4714] |

No primary trait contrast supported; even nominal 95% intervals include zero. Do not assert that specialization or host sharing explains species-level response. Degree does not equal physiological breadth and genus is not phylogenetic correction. No re-tuning for a positive. Consider future independent temporal change-on-change design (fourth period or a separate geographically held-out three-period estimand), and direct EQR/wastewater test, not new claims on this same holdout. Main unchanged. Local reproducibility archive contains full code and tables; remote CI for this specific new diagnostic not executed.
