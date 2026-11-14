# Same-grain butterfly–plant candidate host network: actual exploratory results

2026-11-14. New original-source calculations, not observed feeding, sewage chemistry, pollination, biomass or a meta-analysis. Prior positive/null results retained; main and chocho unchanged.

## Source population and actual work completed
The same original 454 H7×period rows, 232 H7 catchments, 90 H6 parents and 15 countries were rebuilt from the signed original EEA/HydroWASTE, H6 and HydroATLAS archives. Periods are 2017–2019 and 2020–2022. Every comparison below uses the same H7×period biological grain; different periods in one basin are not independent catchments.

Original species source: run 43593227162, artifact 14113695873, ZIP SHA256 4b37c741d1bd28c178b0cef61c9457af0dd30797c173f9384ebe1d1d0036d07c; 2808 successful source histograms across all 234 originally selected H7 catchments. The 454-row inference cohort is unchanged.

New source-only taxonomy and host retrieval: GitHub Actions run **43595715491 SUCCESS**, artifact **14114436781**, ZIP SHA256 **1cf6ed4ee41282a7f1be2e5cc7e1542a186b4fef5e1a46adf5861543c0ebbd94**. Original Ferrer-Paris 2014 catalogue DOI 10.6084/m9.figshare.1168861.v1: **43,645 source rows**, Git blob ea93a5c0fdb17e361a9b6a9637a497a12e90493c, SHA256 a7369d0da02e4e4a8559b5bb9f59b9f2650a4ca790a6792117e21ab53f377f64.

All **303 recorded butterfly speciesKeys** were queried through official GBIF legacy-key usage and synonym services. Exactly matched catalogue aliases, with ambiguous aliases excluded and family identities checked, supply hosts for **209 butterfly keys**. Of **1020 plant names queried**, **996** pass exact species-level Plantae accepted-key matching; **24** are unresolved/rejected and **zero** API requests are missing. Collapsing names mapped to the same accepted plant key yields **990 host plant keys** and **1565 unique catalogue butterfly–plant pairs**. **997 of those distinct pairs** have both partners recorded together in at least one original H7×period. Across the 454 rows there are 140128 candidate pair occurrences; this is NOT 140128 independent observations.

## Coverage is not the same as missing-trait exclusion
| Butterfly family | Recorded speciesKeys | Keys with resolved catalogue hosts | Share of original butterfly records represented |
|---|---:|---:|---:|
| Lycaenidae | 75 | 52 | 96.79% |
| Nymphalidae | 149 | 100 | 94.74% |
| Pieridae | 33 | 20 | 92.38% |
| Papilionidae | 15 | 10 | 99.65% |
| Hesperiidae | 31 | 27 | 99.52% |

Overall, mapped butterflies represent **95.01% of the five-family species-assigned source records**, but only 209/303 source taxon keys. Catalogue-unknown butterflies remain in the source population and the coverage denominator; they are not assigned zero hosts. No all-traits-complete species filter is used. Record representation is not ecological representativeness.

GBIF marks 13 of the 303 butterfly keys DOUBTFUL, including four keys contributing 22 mapped host edges. The primary analysis retains original source taxon-ID semantics; a strict accepted-only sensitivity removes those edges without changing the source community or catchment population. Twenty catalogue butterfly name strings had ambiguous aliases and were not assigned to multiple keys. These are source-key counts, not independently validated numbers of biological species.

## Standardized opportunity metrics
For each H7×period, retain the exact original speciesKey histograms. Compute butterfly inclusion probabilities in 20 sampled butterfly records and plant inclusion probabilities in 50 sampled plant records using the hypergeometric distribution. Counts are GBIF records, not censused organisms.

Expected candidate links = sum over catalogue pairs of Pr(butterfly included in 20 records) × Pr(plant included in 50 records).

Host-support index = the inclusion-weighted fraction of catalogue-mapped butterflies expected to have at least one known host represented in 50 plant records. A butterfly's host probability uses the summed records of all its distinct accepted host keys. This is a standardized record-index opportunity, not the fraction of larvae actually feeding. Median host support over all 454 rows is **0.32399**; median raw co-recorded link count is **314.5**.

## Positive structural result: host correspondence beyond fixed margins
Under 499 permutations of host-plant identities **within plant family × broad source occurrence-frequency bins**, each butterfly's catalogue host count, the multiset of plant degrees, observed communities and rarefaction depths are unchanged. H6-equal observed expected links are **5.4393**, versus randomized mean **4.0680**, a **33.71% excess**. The middle 95% of null draws is **[3.6187, 4.5637]**; none of 499 null means reaches the observed mean (finite Monte Carlo indicator 0.002, not a causal p-value). 912/990 host keys belong to non-singleton permutation strata.

A **posthoc stronger geographic sensitivity** additionally restricts swaps to one of six fixed-seed clusters of normalized host occurrence patterns across H6 parents. The null mean becomes **4.6568**, middle95 **[4.3259, 5.0135]**, and the excess remains **16.80%**. Here 826/990 host keys can move. This is compatible with structured resource matching beyond these chosen margins, but shared geography, habitat and observation bias are not fully eliminated. It is not proof of local feeding or wastewater mediation.

Source sensitivities retain the basic structured-null excess: accepted-only butterfly edges **1543 pairs, +33.06%**; excluding edges whose only source is Wikipedia **1443 pairs, +33.38%**. These are not independent ecological replications.

## Wastewater geography does not provide stable incremental prediction
Primary source-coverage cohort: ≥80% of original records identified to speciesKey in all six groups, **390 rows /217 H7 /83 H6 /14 countries**. Full source cohort 454/232/90 and ≥90% cohort 293 rows are sensitivities. Do not equate these quality screens with measured ecological sample coverage.

Baseline: original area, EEA and GBIF recording opportunities, flow, upstream crop/urban land use, period, country and monitoring method; additionally butterfly S20, plant S50 and catalogue coverage. Add the original static modelled HydroWASTE basin-outlet index W only after these terms. Train-only transformations, H6-equal Ridge alpha10, fixed H6-hash fivefold partitions, 3000 fixed-heldout-error H6 bootstrap draws. Country holdout omits country effects and also deletes training H6 groups shared with the heldout country. Original source method is blank for 44 rows; it remains an explicit unknown category, not a measured method.

| Endpoint, same primary 390 rows | H6-disjoint MSE gain from W | H6 descriptive95 gain interval | Supported-country equal MSE gain |
|---|---:|---|---:|
| Standardized host-support index | **−3.327%** | **[−7.155,+0.496]%** | **+9.085%** |
| log1p(expected candidate links) | **−2.653%** | **[−5.115,−0.206]%** | **+6.097%** |
| Candidate link excess over structured null | **−2.071%** | **[−3.918,−0.416]%** | **+0.763%** |

Country transfer and within-region H6 holdout therefore differ; do not erase the country-positive diagnostics or claim W has no information anywhere. They are different generalization tasks/baselines. The H6 primary test does not support a stable wastewater increment. Predictive deterioration is NOT evidence wastewater reduces feeding or host supply.

Adjusted within-sample W slope for host support is **−0.01260 per +1SD W**, H6-cluster descriptive95 **[−0.02748,+0.00227]**. Its direction is uncertain. W and upstream urban cover remain strongly correlated (**Spearman rho +0.7806**). The source is modelled basin-outlet wastewater geography, not measured effluent at a plant or butterfly.

All-454 host-support W gain **−3.112% [−6.295,−0.625]**; ≥90% coverage gain **−2.635% [−7.522,+2.349]**. In the primary390 cohort, accepted-only W host-support gain **−3.190% [−6.719,+0.308]**, no-Wikipedia gain **−2.425% [−4.655,−0.484]**. Thus none of these H6 variants supports the proposed source increment.

## Which candidate resources contribute most? Descriptive, not sewage-specific
H6-equal contributions to standardized expected links: Urtica dioica **9.606%**, Dactylis glomerata **5.292%**, Humulus lupulus **3.480%**, Alliaria petiolata **3.374%**, Prunus spinosa **3.055%**. The top20 recorded host plants account for **50.41%** of the candidate-link index. This is concentration in reconstructed opportunities, NOT measured feeding mass, nutrient transfer, pollination or a claim these plants increased because of wastewater.

## Verification and status
All **1630 original API response hashes** were verified after artifact download. Independent SciPy hypergeometric probabilities and sparse edge enumeration reproduce all454 opportunity metrics: maximum expected-link discrepancy **8.51e-10**, host-support discrepancy **1.67e-11**. Weighted Ridge predictions were checked against separately fitted scikit-learn weighted-centering/Cholesky models; numerical agreement is near floating-point precision. Seven local test cases cover exact taxonomy gates, rarefaction by exhaustive enumeration, unknown-catalogue handling and null degree margins. All spatial training/test intersections are explicitly checked empty.

The source acquisition ran successfully on GitHub Actions. **The new network fits and independent numerical audit were run locally on downloaded signed original artifacts, not claimed as a new successful GitHub numerical workflow or a second field study.** Complete inputs, scripts, tests and outputs are delivered in the conversation reproduction archive.

European 2024 host checklist DOI 10.5061/dryad.1vhhmgr12 / article 10.1002/ece3.10834: official metadata and file inventory were obtained, including evidence-status and country fields. Its data rows have **not yet been integrated** into these results. The initial Ferrer catalogue shares HOSTS history with chocho and is not an independent validation. Next source upgrade is evidence- and country-filtered regional host use, rather than treating catalogue assertions as observed local interactions.
