# Three US EPA NRSA cycles: original-count EPT rarefaction depth and genuine autecological trait composition

**2026-11-14 | Direct EPA native site/count source; rare100 m2209 vs expanded m2403 original source run SUCCESS, native 2008/2013 EPA autecology plus 2018 strict taxonomic identity source and functional model SUCCESS. New analyses remain POSTHOC to earlier [gesui PR #82](https://github.com/zuizui0223/gesui/pull/82). NO sewage-specific causal or ecosystem-service claims.**

## I. Original source individual-count depth sensitivity — 100 compared with 200

Genuine original EPA 2008/09 wadeable BERWW, 2013/14 and 2018/19 BERW original TOTAL300 taxon-ID counts from exactly the original source cohort n2,648, original source-verified order and HUC8/chemistry/physical habitat; only count N>=100 or >=200 and >=95% counted original individual taxonomic ORDER+TAXA_ID coverage. The source m100 cohort is **n2,403** versus the prior m200 n2,209, adding **194** source visits: 2008 476→523 (+47); 2013 957→1,035 (+78); 2018 776→845 (+69). [Source-only original EPA n100 population before effect SUCCESS #43570309532](https://github.com/zuizui0223/gesui/actions/runs/43570309532).

For the **SAME original native m200 cohort n2,209** the *only* change is hypergeometric individual-based expected EPT taxon-ID richness E[S100] versus pre-established E[S200]; source original specimen counts, observed river TN/P, six physical habitat predictors, fully heldout cycle and all original target-period HUC8 excluded from other-two-period training are unchanged. Response is log1p(expected number of ORIGINAL source EPT taxon IDs from 100 counted specimens), not species density nor adult insect emergence.

[Authentic EPA original count and whole HUC8 external m100 model SUCCESS #43570830665](https://github.com/zuizui0223/gesui/actions/runs/43570830665), source period external test HUC8-equal Ridge alpha10 after original measured stream N/P:

| Whole untouched survey-cycle target | Original m200 EPT taxa: PHab external MSE gain (prior PR82) | SAME native sites at m100: PHab gain [HUC8 fixed OOF 95] | Expanded original m100 site panel: PHab gain [HUC8 fixed OOF95] |
|---|---|---|---|
| 2008–09 | +15.773% | **+15.805% [+9.196,+21.910]%** | n523 **+15.106% [+9.175,+20.489]%** |
| 2013–14 | +3.329% **UNCERTAIN** | +5.364% **[−1.259,+11.468]%** | n1035 +3.903% **[−2.520,+10.021]%** |
| 2018–19 | +14.372% | **+15.188% [+11.016,+19.264]%** | n845 **+13.353% [+9.600,+16.902]%** |

**Both precommitted universal external three-cycle richness PHab criteria FAIL**, at the SAME sites and at expanded sites: 2013 always has fixed HUC8 resample interval spanning zero. This is not an artifact of m200 rarefaction depth, and does not erase original PR81 three-cycle *relative EPT counted fraction* predictive gains.

Original measured ambient TN/TP added AFTER PHab for EPT rare100 richness: on same 2209 source n, 2008 **+8.527%** HUC CI [+4.556,+12.436], 2013 **+5.411%** [+0.937,+9.692], 2018 +4.408% [+0.031,+8.840] (lower bound near zero); on expanded 2403 source visits, 2008 **+7.306%** [+4.397,+10.202], 2013 **+4.947%** [+1.186,+8.569], 2018 +3.158% [−1.294,+7.505] UNCERTAIN. Response itself changes with rarefaction depth, so predictive *percentages* are not directly comparable physical effect sizes or nutrient coefficient directions. The expanded population adds short samples and changes sample-selection ecology; do not substitute its larger n for independent field replication.

Source original 2013 heldout model fitting with 100-depth same-site source has **n355 original other-period training UIDs**, with expanded 100-source n387: both below prior PR81 larger n400 training rule but meet new explicitly frozen m100 train >=300/200 HUC8 source support. No unreported source missingness imputation. Fixed HUC8 bootstrap is resampling fixed out-of-period errors, not independent refit uncertainty.

## II. Official EPA taxon-level autecology identity, measured-count-functional composition

[EPA NRSA original 2008 Taxa List and Autecology](https://www.epa.gov/sites/default/files/2016-06/nrsa_0809_benttaxa.csv), original [2013 benthic taxa and FFG/HABIT/PTV metadata](https://www.epa.gov/sites/default/files/2019-04/nrsa1314_benttaxa_meta_04292019.txt) and original 2018/19 EPA BERW taxon-level counts: exact source native TAXA_ID and normalized native TARGET_TAXON name AND ORDER matching (no posthoc inference or same numeric ID alone). Original [taxon ID, source, field-name and trait-code audit SUCCESS #43570200596](https://github.com/zuizui0223/gesui/actions/runs/43570200596) and [exact true native-UID site-level trait availability before effects SUCCESS #43571110513](https://github.com/zuizui0223/gesui/actions/runs/43571110513). EPT counted-individual trait assignment coverage:
- **FFG functional feeding group:** 98.34% /96.86% /96.13% in 2008/2013/2018.
- **HABIT substrate/microhabitat role:** 99.17% /99.03% /97.85%.
- **PTV source pollution tolerance attribute:** 94.74% /93.56% /93.07%.
Each EPT-positive survey UID needs >=80% of **actual EPT counted individuals** assigned the specific original trait; original no-EPT sites have an **undefined**, not zero, trait composition. All original site and strict external HUC8 support criteria passed. True source n site UIDs: scraper EPT guild 454/892/701, clinger habit 462/918/719, mean PTV 440/851/674.

Original EPA FFG and HABIT codes may list MULTIPLE modalities. A source label such as SC,SH is given fractional 0.5 scraper weight, rather than each mode given a whole count. This is a model of taxon-role assignment, not actual field ingestion/functional flux. PTV is original taxon-assigned pollution tolerance, 0–10, not an experimental sensitivity measurement.

### III. New genuine EPT functional trait environmental results

[Original source-count, exact ID official trait and 3-era heldout functional model SUCCESS #43571322998](https://github.com/zuizui0223/gesui/actions/runs/43571322998). Every comparison trains the *other TWO* original survey periods after excluding ALL HUC8 in the heldout original full survey-period site cohort, even if that site lacks an original EPT trait. Same fixed four baseline/ambient stream N/P/PHab/both models; train-only one-hot and scaling, equal HUC8 training and test MSE; 2,000 original HUC8 fixed OOF error bootstrap descriptive95.

**Clinger habit within original counted EPT organisms** (defined by EPA HABIT code CN, fractional for multiple listed modes):

| Original untouched external survey cycle | n site UIDs / HUC8 | PHab added AFTER real measured ambient N/P, external HUC8 MSE gain [descriptive95] | HUC8-clustered descriptive β per +1 SD actual XEMBED, controlling TN/TP, other five PHab, site and season [95] |
|---|---:|---|---|
| 2008–09 | 462 /360 | **+9.507% [+5.676,+13.436]%** | **−0.10416 [−0.14263,−0.06569]** |
| 2013–14 | 918 /605 | **+7.737% [+3.789,+11.596]%** | **−0.10849 [−0.13388,−0.08310]** |
| 2018–19 | 719 /481 | **+3.811% [+0.063,+7.060]%** | **−0.07676 [−0.10843,−0.04508]** |

**This original field trait composition is the clearest newly observed three-period consistent *predictive/associational* result**: original EPA EPT community with lower fraction of taxa assigned clinger habit under higher riverbed embedding, with useful six-PHab source environmental information in heldout original cycles and nonoverlapping training catchments. It is a difference in taxon composition, **not the same individuals evolving a new habit** or the response to experimentally adding sediment. Multiple habitat covariates, flow, shading, nutrients and taxa detection may confound the embeddedness slope. The 2018 habitat external CI lower bound +0.063% is marginal, not robust proof across all settings.

**EPT pollution-tolerance community-weighted score (source PTV)**:
- Conditional original XEMBED gradient positive in all three periods: 2008 β **+0.22283 [+.10249,+.34317]**, 2013 **+.14369 [+.05958,+.22779]**, 2018 **+.20107 [+.10528,+.29687]**. Higher PTV denotes taxon assemblages with higher assigned pollution tolerance, NOT evolution or realized toxicity exposure.
- PHab added after N/P external prediction is **+4.267% [+1.695,+6.640]** only in 2008, **−0.637% [−4.298,+3.034]** in 2013 (uncertain), **+0.0277% [−2.540,+2.586]** in 2018 (uncertain). Hence observed source partial gradient sign may be stable while the grouped environmental *predictive extra contribution* is not.
- TP partial within-era PTV β positive 2008 **+.15431 [.03365,.27498]**, 2013 **+.16279 [.06993,.25564]**, 2018 **+.05579 [−.06264,.17422]** crosses zero. It is not verified WWTP-derived or chemical cause of altered tolerance.

**EPT scraper functional feeding FFG role (source SC, multimodal coded fraction):** PHab source external added values 2008 −1.236% [−5.265,+3.345], 2013 −1.038% [−4.664,+1.974], 2018 +1.194% [−1.927,+3.666]. All uncertain/worsening. **Do not infer a universal increased or decreased benthic algae grazing ecosystem service** or treat scraper metric as observed periphyton consumption.

### Scientific synthesis, important remaining uncertainty

The measured original EPA data now separate **four different ecology indicators** on actual original stream sites: EPT relative counted share, count-standardized EPT original TAXA_ID richness, original role-assigned clinger/scraper composition and the source taxon-assigned weighted PTV pollution tolerance. Original physical habitat can transport across cycles for relative EPT share and clinger composition, while independent richness and scraper role are less stable. Source functional PTV vs original embedding shows a consistent conditional gradient but not consistently useful out-of-period PHab prediction. **Function and diversity respond to different dimensions of the environment**; a single taxa count or one EPT fraction is not a sufficient proxy for ecosystem service, ecological mass production, or functional transfer.

All results are observational source taxon composition in the same EPA field cohorts. 2008 wadeable BERWW differs from 2013/18 BERW protocols, the sampled 2013/18 source cohorts were originally filtered by valid periphyton Chla even when these trait outcomes do not require Chla, and 2018 traits are inherited from published 2008/2013 autecology through rigorous identity matching (not independently measured 2018 physiological roles). Multiple posthoc comparisons, HUC8-correlation and source missingness may inflate apparent universality. No field habitat manipulation or evolutionary observation.

**The original WWTP-specific historical DMR actual dual N/P source delivery gate failed** in [PR #74](https://github.com/zuizui0223/gesui/pull/74). Original 2013 new-catchment EPT×Chla slope generalization failed in [PR #79](https://github.com/zuizui0223/gesui/pull/79), and original PR82 rare200 3-era universal EPT taxon richness PHab gate failed. This document retains those nulls. No observed algal primary production rate, dry emerging insect biomass, actual terrestrial food-web flux or realized water purification/pollination service. eBird/TinyFish excluded. DRAFT unmerged; do not alter gesui/main.


## Independent source-count and numerical cross-implementation validation

[Full original EPA EPT FFG/HABIT/PTV functional-value and cross-period numerical audit SUCCESS #43571824349](https://github.com/zuizui0223/gesui/actions/runs/43571824349):

- Reconstructed **every actual source 2,209 eligible UID** EPT scraper SC, clinger CN and PTV composition weighted by exact source individual TOTAL300 counts, including correct exclusion of sites with original EPT=0 or <80% source role coverage; all independently calculated source response values matched the separately implemented original source result (raw numerical tolerance 1e-10).
- Separately reimplemented **all three eras × three functional traits × four out-of-era models = 36 HUC8-equal external Ridge10 MSE calculations**, with a manual normal-equation weighted NumPy implementation (no use of original sklearn predictor helper), original full target survey-period HUC8 excluded from train.
- Independently reconstructed each of the nine outcome-era PHab and NP gains and their 2000 HUC8 block-error fixed-holdout descriptive 95% quantiles to source tolerances; independently checked all **27 XEMBED/TN/TP adjusted source descriptor coefficients and HUC8-cluster intervals**. A ±0.3% half-CI-width numerical tolerance handles differences in finite-sample rank correction with source state+eco categorical matrices (not 0.3% biological uncertainty or an ecological signal).
- Validation is **independent code over the exact same original field specimens**, not a separate empirical confirmation, ecosystem service measurement, wastewater intervention or novel field trait measurement. Missing original taxonomy remains missing; role assignments are source autecology, not realized diet.
- The m100 original empirical results in this document have their own successful original GitHub run #43570830665, and source-only count gate #43570309532. The high-level m100 external result is therefore **source-backed and reproducible via repository workflow**, but no separately rewritten numerical matrix-implementation QA for m100 is claimed in this ledger.
