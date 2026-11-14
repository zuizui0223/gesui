# EPA three-cycle EPT **taxon-ID richness** versus EPT **counted fraction**: source-blind standardization and cross-era transfer

2026-11-14. This is deliberately **posthoc** to positive full-six-PHab external predictive value for counted EPT fraction in [PR #81](https://github.com/zuizui0223/gesui/pull/81), but the **new measured-response definition, effort standardization, original source-cohort gates and fixed outcomes are frozen before opening richness model results**.

## Research test

The previously supported generalization is specifically **relative count fraction** of Ephemeroptera, Plecoptera, Trichoptera (EPT), not a number of taxa, absolute insect numbers or insect dry biomass. Test whether the same original six physical field habitat features and actual ambient river TN/TP also generalize as **richness**, independent of differences in original source counted individuals per reach.

### True original source contract

Recompute raw taxon-by-survey visit counts from the exact source-original EPA CSVs and SHA:

- **2008–09**: EPA original wadeable `BERWW` + `SAMPLE_CAT=P`, original published `TOTAL300` positive counted subsample and `TAXA_ID`; **original national taxa `TAXA_ID→ORDER` crosswalk** source-SHA `6165d5a7ec8d4ee3e85149255fd967f6268d870d49c57de6302ca4e124673b26`. Original Benthic Count publisher CSV SHA `8251fbe2d2bfc699225873ec43231feda0b110bd112badb61d6e8de486870d31`. Restrict to the **same original 586 native UID** frozen in PR81.
- **2013–14**: EPA original `nrsa1314_bentcnts_04232019.csv` `SAMPLE_TYPE=BERW`, original source SHA from [original 2013 publisher six-table frozen manifest](https://github.com/zuizui0223/gesui/actions/runs/43560769667); positive `TOTAL300`, `TAXA_ID` and **native `ORDER`**, restrict to exact frozen **1,116 original UID** from PR81.
- **2018–19**: EPA original `nrsa_1819_benthic_macroinvertebrate_count_-_data.csv` `SAMPLE_TYPE=BERW` and same original `TOTAL300`, `TAXA_ID`, `ORDER`, SHA `4e9b3a0798b25f67c7a96d4800b94f86544bbd623abaf119358e067b158324a7`; original **946 native UID** from PR81.

Only verified source species taxonomy `TAXA_ID` and **original order**; no automatic genus-to-species assignment, no inferred order for unknown original `TAXA_ID`. Standardized original `TOTAL300` counts are individual counts *in a laboratory subsample*, NOT organisms/m². If duplicate `UID×TAXA_ID` appears, aggregate original counts within that taxon only after disclosing duplicate row number; never add duplicates as new richness. Audit integer individual counts and original `TOTAL300` sum, source fraction `known-EPT/all counted` must **exactly match frozen source EPT share** for all 2,648 real UID visits (floating tolerance 1e-9). At each site original source counted organisms with a recognized `ORDER` must be at least **95%**; 2008 original already enforces this. Source missingness/unknown taxonomy is **not biological zero**; no source data/UID or sampling rules are retroactively loosened based on signs.

### Fixed sample-effort standardization prior to inspecting the new model

Original counts can be below 300 even though the field is called `TOTAL300`. Define primary **analytical individual-based rarefaction to exactly m=200 counted individuals**, only for sites with source total N ≥200. For each true original `UID×TAXA_ID` count `n_i`, calculate expected number of taxon IDs in a random subsample without replacement of m organisms:

`E[S_m]=Σ_i (1 - C(N-n_i,m)/C(N,m))`.

For **EPT taxon-ID richness**, sum only records whose original `ORDER` is one of Ephemeroptera/Plecoptera/Trichoptera, while N includes **ALL original counted organisms**, including known other orders and unknown-order specimens. For **non-EPT original taxon-ID richness**, sum ONLY known non-EPT records, not the unknown-order fraction. **Total source taxon-ID rarefied richness** includes any original identified TAXA_ID even if order unknown, but an empty original taxon ID is not an organism-to-species assignment. Count-source flags and empty IDs gate, never impute taxa. Outcome richness is **not species richness, relative count share or insect biomass**.

Use `scipy.special.gammaln` or equivalent exact hypergeometric log combination, bound probability [0,1], not Monte Carlo stochastic rarefaction. Store source-only count N distribution per cycle and number of complete source UID/N≥200/known order coverage; aggregate response results to JSON (no per-visit private coordinates or raw individual field records committed).

**Pre-fixed source viability gate at m200:** each cycle ≥250 genuine eligible UID, ≥150 canonical HUC8, ≥20 original states; source matched EPT fraction for all site UIDs, count integer, no unsupported taxonomy. If any cycle falls below, do **not** fit the claimed 3-cycle 200-individual comparative result. A secondary **m100** source-standardization may be run only if m200 gate fails in a separate explicitly exploratory protocol, not relabeled as the original 200 model or chosen by ecological sign.

Also retain original source **EPT counted fraction** as a CONTROL outcome on **the exact same rarefied-quality eligible UID sites**. This lets us distinguish a source *sample-selection / counted-effort* effect from a genuine difference in whether physical habitat explains compositional EPT share versus standardized EPT taxon count.

### Cross-survey-year out-of-HUC8 test (fixed)

Reuse exactly [PR81's source-frozen four baseline models](https://github.com/zuizui0223/gesui/pull/81): B=lat/lon, original stream order, sine/cosine field season, log conductivity, original state/eco9; +NP=source-measured real river log TN/TP; +PHab = six strictly matching original PHab field canopy, riverbed embeddedness, overhanging cover, human disturbance and log channel width; +both=all. Ridge alpha 10 fixed; original training HUC8 **equal total fitting weight**, original target HUC8 equal MSE. Each target full cycle held out while **ALL shared canonical HUC8 are removed from the other two training periods**; train-only standardization and category vocab. No year extrapolation or test-target response used in training.

For each of three original response types (primary **rare200 EPT taxon-ID richness**, negative-control **original EPT counted share** on same rare200 cohort, secondary **rare200 non-EPT taxon-ID richness**) calculate all 4 model MSEs, NP added *after* PHab and PHab added *after* NP; 2,000 fixed HUC8 heldout error bootstrap descriptive95. **All nulls, negative and heterogeneous effects retained.** Richness is a nonnegative expected number of original taxon IDs rather than integer species number; put a supplementary `log1p` transform only if fixed before fitting (primary will use `log1p(expected rich200)` to reduce influence of rare high count taxa, with raw expected richness as source descriptive QA).

Predeclare primary novelty threshold: **all three survey eras** external PHab gain has positive 95% HUC8-block descriptive bootstrap for rare200 EPT richness *after measured ambient N/P*. If not, reject universal predictive richness claim. Predictive gains do not identify cause/biological abundance. Distinguish 2008 `BERWW` from 2013/18 `BERW`, and source 2013/18 samples are original periphyton-Chla-eligible while 2008 is source-only chemistry eligible; national river lengths are not sample-weighted.

**Do not infer:** wastewater source-delivered total N/P (failed original DMR gate in [PR #74](https://github.com/zuizui0223/gesui/pull/74)), full energy subsidy, primary production, insect emergence mass, predator feeding, pollination or water purification. Keep PR81's robust EPT counted-fraction transport and 2013 unshared-basin EPT–Chla interaction non-replication. Branch Draft /unmerged, no eBird, no TinyFish.


## Source-only census outcome and external training support gate (no new ecological model inspected)

The [raw EPA 3-cycle count source-only audit SUCCESS #43568105806](https://github.com/zuizui0223/gesui/actions/runs/43568105806) reported original native count>=200 and >=95% order+taxonID eligibility, all independently matching the already frozen original counted EPT fraction:
- 2008: **476** original UIDs, **369** source HUC8, **46** states; 2008 original TOTAL300 median 300, 2008 original raw positive count integer source n586 broader before m200.
- 2013: **957** original UIDs, **623** source HUC8, **47** states.
- 2018: **776** original UIDs, **512** source HUC8, **47** states.
- Total **2,209 original native UID** source opportunities; these source HUC numbers are sums by cycle, **not distinct union catchments** (HUC8 identity overlap is re-audited before model).
- Each survey original taxon-count native UID shares exactly match frozen PR81 EPT fraction: worst difference ≈1.1×10^-16; original positive `TOTAL300` values integer and no repeated original UID+TAXA_ID rows across selected original source cohorts.

**New rarefied-subcohort model-fit support gate, fixed before any newly calculated EPT richness or external-prediction outcomes:** for a left-out period, target after rare200 ≥250 original UIDs and ≥150 canonical HUC8, and other-two-period training AFTER removing every source original HUC8 present in the target entire wider original cycle must retain ≥300 original rare200 UIDs and ≥200 distinct canonical HUC8, ≥15 states. This is explicitly a **smaller effort-standardized subset** with its own transparent a-priori support gate; it is not a claim that the original PR81 larger-sample ≥400 UID gate was met. If a period fails, export its support QA and do NOT fit it or weaken this threshold after seeing rich EPT effects. All four baseline models and HUC8-equal weights remain identical to PR81.
