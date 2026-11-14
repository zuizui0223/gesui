# Six-family butterfly larval host breadth and chocho cross-catalogue audit — 2026-11-14

**NEW exploratory analysis of external larval-host records, NOT an observed wastewater effect, butterfly feeding event or original river-basin species-level network. PR #84 remains Draft/unmerged.**

## Original source and reproducibility

Original Ferrer-Paris et al. (2014) butterfly larval-host compilation: DOI https://doi.org/10.6084/m9.figshare.1168861.v1; verifiable mirror https://github.com/globalbioticinteractions/ferrer-paris; original CSV Git blob SHA1 ea93a5c0fdb17e361a9b6a9637a497a12e90493c.

Parsed **43,645 source records** with exact original 7-column CSV schema. After requiring a source butterfly name resembling a binomial, a plant name resembling a binomial and a non-missing plant-family label, **3,785 butterfly names in six recognized butterfly families** were available: Nymphalidae 1,447; Lycaenidae 880; Hesperiidae 748; Pieridae 325; Papilionidae 271; Riodinidae 114. Names are source-reported, **not independently resolved to accepted taxonomy**. Same butterfly × plant name duplicates collapsed; raw rows are not independent ecological measurements. Four original rows with missing butterfly family were excluded from family analyses.

For each butterfly name with N >= m distinct accepted-*looking* binomial host plant names, compute the hypergeometric expected host-plant family number among exactly m distinct recorded host taxa:

E(F_m) = SUM over plant families j [ 1 - C(N - n_j, m)/C(N,m) ].

Here the N are **source-known distinct host-name strings**, not plant individuals or in situ host biomass. This analytically standardizes recorded host-list length, but does NOT eliminate effort/synonym/phylogeny/geographic source bias. For descriptive uncertainty, 4,000 resamples of complete butterfly genera within each family. Inter-family difference intervals use **paired unsorted independent bootstrap draws** (not differences of sorted marginal quantiles).

## Results: 5 recorded host plants, then 10 recorded host plants

| Family | Species n >=5 | Mean expected host families in 5 | Species n >=10 | Mean expected host families in 10 | Genus-resampled descriptive95 of 10-host mean |
|---|---:|---:|---:|---:|---|
| Lycaenidae | 111 | **2.1626** | 41 | **3.3535** | [2.5333,4.2339] |
| Nymphalidae | 307 | 1.8262 | 100 | 2.6803 | [2.1310,3.1901] |
| Pieridae | 67 | 1.4871 | 35 | 1.8113 | [1.4429,2.4109] |
| Papilionidae | 75 | 1.5640 | 38 | 2.3954 | [1.1190,2.6347] |
| Hesperiidae | 111 | 1.4710 | 34 | 1.8648 | [1.4952,2.3611] |

Source Riodinidae had only one butterfly name eligible for 10-host comparison and is NOT placed in this five-family inferential contrast.

**Lycaenidae minus each other family, at 10 host names:**
- Pieridae +1.5422 [genus-resampled descriptive95 +0.5477, +2.4996].
- Hesperiidae +1.4887 [+0.5279, +2.4548].
- Papilionidae +0.9581 [+0.1128, +2.5076]. Only **8** genera among the 38 10-host-eligible Papilionidae, 24 named Papilio.
- Nymphalidae +0.6732 [−0.2654, +1.7213]: **not conclusively different**.

At 5 known host names, Lycaenidae minus Pieridae is +0.6755 [+0.2436, +1.0979] and minus Nymphalidae is +0.3363 [−0.0518, +0.7486]. The higher-support 5-host contrast retains the same broad pattern with a null-intersecting Nymphalidae comparison.

Remove all source rows attributed to Wikipedia as posthoc source sensitivity: expected distinct host families among 10 are **Lycaenidae 3.3797 (n38), Nymphalidae 2.6600 (n96), Pieridae 1.7910 (n33), Papilionidae 2.3997 (n38), Hesperiidae 1.8694 (n34)**. Lycaenidae vs Nymphalidae still uncertain.

Do NOT equate observed mean catalog breadth with any specific species-level local breadth, pollution tolerance, adult puddling or evolutionary specialization. Selection of species with >=10 documented hosts favors intensively studied and/or generalist taxa and is not a random global butterfly sample.

## Actual cross-catalogue comparison with frozen chocho resource geography

Connected GitHub repository zuizui0223/chocho fixed data/frozen/figure_sources/anthropogenic_species_metrics.csv contains **239 names** with reconstructed known larval host geography. Exact case-normalized source-name matching yielded **228 species** with the Ferrer-Paris strict source-host list:
Hesperiidae 51, Nymphalidae 94, Pieridae 42, Riodinidae 1, Papilionidae 16, Lycaenidae 24.

Spearman correlations in these **same 228** source-matched butterflies:

| Diagnostic | Spearman rho |
|---|---:|
| Original Ferrer host species names vs chocho resolved host species | +0.9442 |
| Original Ferrer host family count vs chocho resolved host family count | +0.6574 |
| Ferrer host family count vs chocho introduced-host geographic resource expansion (log ratio) | +0.0331 |
| chocho host family count vs same geographic expansion | +0.0020 |

**Not independent validation:** the Ferrer database itself includes the historical HOSTS DB and the chocho pipeline uses a pinned HOSTS mirror, so source overlap can drive correlation. Matching 228 source names does not identify which are locally present in sewage-exposed watersheds. Geographic resource expansion is chocho's **modeled historical host redistribution**, not measured wastewater treatment or a local population effect. Near-zero Spearman is not formal absence/equivalence.

## What is and is not ready for the original 454 H7 basin-period comparison

**Available now:** global source-derived larval host family breadth of 3,785 butterfly-name strings; exact 228-name cross-catalogue QA with existing chocho host resource work; EU same-H7-period aggregated EQR, plant 50-record S50 and five butterfly-family recording composition in 454 periods.

**Missing and cannot be imputed:** actual plant species identities per H7×period, butterfly species identities per the same H7×period, source-vetted geographic host links with accepted keys, species-specific feeding, adult puddling, ant attendance, actual flowering and insect biomass. The old source GBIF facet script computed species facets but retained **summary counts only** in the original persisted artifacts; source request URLs and hashes do not restore species identities.

Next spatial source acquisition contract: use **the same originally frozen 234 EEA H7 catchments, two fixed date windows, identical original polygons and GBIF quality settings**, but retain **speciesKey histogram identities** for both plants and five butterfly families; paginate all facets, retain exact GBIF response hashes, missingness and unresolved taxonomy, perform original-source keys and species-level host name accepted crosswalk. Never substitute chocho's historical ports or islands for the river basins; never infer realized larval-feeding edges from mere co-indexing. Country/H6 out-of-area validation and site-level effluent chemical verification remain required.

## Inspection/reproduction

- Species-level derived 3,785-butterfly host breadth source table: [CSV](../data/derived/FERRER_PARIS_2014_BUTTERFLY_SPECIES_HOST_FAMILY_RAREFACTION_SOURCE_DERIVED.csv)
- Exact 228-butterfly Ferrer x chocho source comparison: [CSV](../data/derived/FERRER_PARIS_2014_X_CHOCHO_239_BUTTERFLY_HOST_CATALOGUE_MATCHED_228.csv)
- Offline six-family source rarefaction and bootstrap script: [Python](../scripts/gesui_larval_host_rarefaction_sixfamilies_20261114.py)

The script requires the original external raw CSV and has not been claimed to run successfully in a new GitHub Actions job. Both original source datasets are separately citable. Keep all older successful/null wastewater PR results intact, no main/chocho merge, no eBird or TinyFish.
