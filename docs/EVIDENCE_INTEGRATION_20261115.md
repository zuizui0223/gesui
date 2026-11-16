# gesui evidence integration — 2026-11-15

Status: exploratory, source-backed audit; PR draft; not peer reviewed. This document indexes all existing gesui evidence families plus the previously local 2026-11-14/15 analyses. It does not assert that every prior raw binary or 3,508 files is stored as Git blobs.

## Existing GitHub results (keep original PRs and negative evidence)
| Evidence | Original tracking | Observed claim; do not expand |
|---|---|---|
| European upstream reported wastewater load × benthic EQR (2018 and 2022) | [PR #14](https://github.com/zuizui0223/gesui/pull/14), [PR #17](https://github.com/zuizui0223/gesui/pull/17) | Regional-out EQR prediction improved ~1% when adding reported upstream N/P in two years; EQR is not taxa richness, years share 1,167 sites and exposure is confounded |
| Wastewater dose actually delivered to water chemistry | [PR #8](https://github.com/zuizui0223/gesui/pull/8), [PR #9](https://github.com/zuizui0223/gesui/pull/9), [PR #74](https://github.com/zuizui0223/gesui/pull/74) | Upstream modeled water quantity and even reported mass often did not improve nutrient concentration generalization. US matched effluent N/P source gate failed; never infer site dose |
| Plant N preference and whole flora | [PR #10](https://github.com/zuizui0223/gesui/pull/10), [PR #12](https://github.com/zuizui0223/gesui/pull/12) | Matched plant comparison underpowered; EIV trait missingness asymmetric; no confirmed wastewater-induced floral selection |
| Same-station aquatic algae/periphyton and EPT | [PR #73](https://github.com/zuizui0223/gesui/pull/73) to [PR #83](https://github.com/zuizui0223/gesui/pull/83) | Three EPA eras show habitat association with EPT counted share; rarefied EPT richness transfers in only two of three eras; EPT clinger share decreases and tolerance value increases with embeddedness in all three. None identify wastewater-specific dose or emergence biomass |
| Same-H7 butterfly family composition and turnover | [PR #84](https://github.com/zuizui0223/gesui/pull/84) | H6-local prediction gain +10.486% for GBIF family Hill-Shannon index, but country transfer and temporal coupling fail; Lycaenidae-specific effect not established |
| Source-identical species counts and candidate host-network | [PR #85](https://github.com/zuizui0223/gesui/pull/85), [PR #86](https://github.com/zuizui0223/gesui/pull/86) | 2,808 source facet opportunities recovered; host network 997 candidate butterfly–plant pairs co-indexed at least once. Geography-constrained candidate links exceed a constrained null by 16.80%: not observed larval feeding |

## This integration: previously locally computed analytical stages
1. [2026-11-14 network robustness and Pieridae](RESULTS_20261114_NETWORK_ROBUSTNESS_PIERIDAE_SYNC.md): top-20 known host removal reduces constrained-network excess from 16.80% to 4.33%; modelled WWTP geography improves Pieridae host-support H6 MSE 12.96%, but country transfer uncertain, and full network W effect not supported.
2. [2026-11-15 butterfly HNS source/trait scores](RESULTS_20261115_HOST_NUTRIENT_SCORE_SYNC.md): 303 gesui speciesKey rows (190 with at least one numeric N host), 230 chocho roster names (147 numeric HNS); different catalogues overlap and are not independent studies.
3. [2026-11-15 HNS basin macro analysis](RESULTS_20261115_HNS_MACRO_SYNC.md): 454 H7×period, 232 H7, 90 H6, score-eligible 186 accepted butterfly keys, 91.06% of 4.4M raw five-family indexed records. W-associated record-weighted mean HNS beta +0.1301 [0.0626,0.1975], geographic adjustment +0.1080 [0.0389,0.1772]. H6 predictive generalization and country validation vary; butterfly counts not individuals.
4. [Structural pathway/SEM evidence contract](PATH_ANALYSIS_IDENTIFIABILITY_20261115.md): only same-unit observed arrows can be modeled; no mediator claiming confirmed food-web energy/biomass from disconnected regional indices.

## Important negative and selection-sensitive results retained
- The all-454 H6-out prediction MSE improvement of the geographic-adjusted HNS model is +5.374% with descriptive 95 interval [-5.901,+14.780] (crosses zero); above-zero conditional W coefficient is not proof of transferable prediction or causality.
- Species-incidence-weighted geographic-adjusted HNS beta +0.0342 [-0.0017,+0.0700] spans zero; change is stronger in GBIF record-weighted composition than equal-species incidence.
- Across 222 repeated H7 periods, HNS increased in 112 and decreased in 110. Static W cannot prove increasing effluent caused temporal change.
- In adversarial species-level unscored-host assignments the W beta can range from -0.1496 to +0.3245; source trait completeness is a major limitation.
- In the network, adding W for geographic constrained null-excess links has H6 MSE gain -1.12% [-3.06,+0.57] in quality cohort: no universal network enhancement.

## Provenance and repository safety
- Current evidence is in unmerged Draft PRs. **Do not merge main or chocho based on these exploratory results.**
- The large original EIV Excel, original chocho host roster, provenance API response bundles and raw binary zip files are not duplicated into public Git history; their original publisher/Actions URLs and SHA receipts are in [archive_manifest](../results/20261115_sync/LOCAL_ARCHIVE_SHA256.csv). Public redistribution status of chocho original host roster must be reviewed before publishing raw relationships.
- This is a synchronized **result ledger, compact numerical outputs, and integrity audit**, not a full byte-for-byte GitHub upload of all local binary archives. The user retains locally supplied full reproducible ZIPs.
- HNS is mean host-plant habitat N preference (Ellenberg-type N), not leaf tissue nitrogen, adult mineral intake, delivered wastewater N/P, or an established existing CNI metric.
