# Three original EPA NRSA surveys: EPT taxon-ID richness versus counted EPT fraction at equal 200-individual survey effort

**2026-11-14 · New source-frozen empirical macro analysis / posthoc to previously established PR #81 EPT fraction findings. Original source audit, actual EPA field count reconstitution, and independent mathematical + model QA SUCCESS. DRAFT, not sewage-specific causality or realized ecosystem-service responses.**

## Question

Does the historically transportable physical habitat information for the **proportion** of counted mayfly, stonefly and caddisfly (EPT) organisms also transport for EPT **taxon-ID richness** when original taxon count effort is standardized to the same number of sampled organisms? These are **different functions of the same original benthic count samples**; they do not comprise independently measured ecosystem layers, and original `TAXA_ID` does **not automatically equal species**.

### Genuine EPA native field source and before-effect quality gates

Directly reconstructed original EPA published raw wadeable `TOTAL300` counts and native taxonomic identity on source-locked PR #81 field visits, retaining **exact original numerical field measured river TN/TP, six physical habitat metrics, UID+HUC8 source geography**, 2008 official `TAXA_ID→ORDER` mapping as originally fidelity-validated with EPA's published EPT share, and 2013/18 original native `ORDER`:

- [Raw source **before any new richness model** SUCCESS #43568105806](https://github.com/zuizui0223/gesui/actions/runs/43568105806): original source count distributions, `TOTAL300` integer field, native 2008 BERWW versus 2013/18 BERW, 100/200/300-organism source availability, 95% known original ORDER+taxonID per site, raw publisher SHA256, no duplicated native UID×taxonID record and no surprise source joins.
- All original counted source EPT **proportions recalculated from raw specimens matched the previous PR #81 frozen original ratios** to max difference ~1.1×10^-16. Thus the new richness comparison shares exactly the prior trusted original field count cohort before m200 effort standardization.
- Original biological samples may have **fewer than 300 individuals** although field is named `TOTAL300`; median source N300, original SOURCE range N1–300. Crucially do **not** regard short original counted samples as true absent EPT or as records sampled with 300 individuals.
- Predefined **m=200** hypergeometric individual-based rarefaction only when source original `N>=200` and >=95% of all counted individuals have actual original taxonomy+order, with each EPT taxonID counted once after original UID×TAXA_ID grouping. The expected number of EPT taxon IDs in 200 original counted organisms is:

\[
E[S_{200,\mathrm{EPT}}]=\sum_{i\in\mathrm{EPT}}\left[1-\frac{\binom{N-n_i}{200}}{\binom{N}{200}}\right].
\]

Here **N includes ALL counted original organisms**, including genuinely unknown-order ones. Only original source EPT `ORDER` contributes to EPT richness; source known non-EPT IDs alone contribute to `nonEPT`, unknown source taxon IDs are NOT invented as species. Use source-only `TAXA_ID`, not third-party inferred classifications. Results are exact hypergeometric **expectations**, not stochastic resampling or true plot-level insect biomass.
- Source-only fixed site/region/effort gate >=250 native UIDs, >=150 HUC8 and >=20 states PER ERA **PASSED for all**. Since retaining only N≥200 inevitably shrinks training samples, a further externally locked pre-numeric reduced-sample train gate was **>=300 real source UIDs /200 HUC8 /15 states** with every original target-period HUC8 removed; this is a separate declared rare200 subset support gate, NOT an assertion that original PR81 n400 train gate remained satisfied.

| Genuine original EPA period | Original prior EPT source UIDs | Source m200 valid UID /distinct HUC8 /states | Mean expected EPT taxonIDs per 200 counted organisms | Spearman original EPT fraction vs rarefied EPT taxon richness |
|---|---:|---|---:|---:|
| 2008–09 | 586 | **476 /369 /46** | **9.007** | **+0.731** |
| 2013–14 | 1,116 | **957 /623 /47** | **9.041** | **+0.776** |
| 2018–19 | 946 | **776 /512 /47** | **9.772** | **+0.742** |
| **Total** | **2,648** | **2,209 source UID visits** | — | — |

The 369+623+512 period-specific HUC8 are **NOT all independent union basins**; the original three-period source HUC8 overlap and geographic canonicalization from PR #81 are preserved. The native US EPA 2008 method BERWW may differ operationally from later BERW and the EPA taxonomic resolution can vary with survey cycle. These are original survey **quality-conditioned** subcohorts, and conditioning on N≥200 likely excludes some low-abundance/difficult-to-sample sites: a potential ecological missing-not-at-random / selection problem. It does not prove a community of 200 organisms naturally existed at every source reach.

## True cross-era, zero-overlapping-original-HUC8 external test on SAME rarefied source UID sites

[Original source measured m200 and full external model SUCCESS #43568477035](https://github.com/zuizui0223/gesui/actions/runs/43568477035). Three distinct responses with **the exact same m200 eligible site/visit IDs within each period**:

1. **Primary:** `log1p(E[S_200, EPT native taxonIDs])` — counted EPT diversity under equal organism sampling effort, not true insect density or number of species.
2. **Same-site control:** original `EPT counted relative fraction` (unchanged source PR81, restricted to same m200 sites).
3. **Secondary:** `log1p(E[S_200, source known nonEPT native taxonIDs])`, with original unknown-order individuals excluded from nonEPT labels but included in the sample denominator. This is a *taxonomic contrast*, not automatically tolerant or generalist insects.

Fixed original PR81 covariates B (lat/lon, original Strahler, season, measured conductivity, state + original eco9); B+measured ambient TN+TP; B+same six original field PHab (XEMBED, XCDENMID, XPCAN, XFC_OHV, W1_HALL, log original wetted stream width); and B+both. Ridge alpha10, categorical levels and numerical scaling trained only on original other survey cycles; **HUC8-equal fitting weight and HUC8-equal heldout MSE**, 2,000 fixed heldout HUC8 block-error resamples (descriptive95). For target heldout cycle, all other-period training UIDs with ANY canonical HUC8 found in that target's **entire original PR81 586/1116/946 cohort** are excluded (more conservative than removing only rare200 eligible target HUC8). Zero train/test HUC8 overlap verified. No original test responses used in feature selection or training.

### Main comparison: physical habitat added AFTER measured TN/TP

| Heldout original survey period | Original m200 test UIDs /HUC8 | Original training UID/HUC8 after strict HUC8 exclusions | **EPT counted fraction** habitat MSE gain [HUC8 bootstrap95] | **Rare200 EPT taxon-ID richness** habitat gain [bootstrap95] |
|---|---:|---:|---|---|
| 2008–09 | 476 /369 | 1,014 /537 | **+10.125% [+3.816,+16.400]%** | **+15.773% [+9.403,+21.521]%** |
| 2013–14 | 957 /623 | **355 /272** | **+8.881% [+4.010,+13.547]%** | +3.329% **[−3.046,+9.397]%** (uncertain) |
| 2018–19 | 776 /512 | 555 /375 | **+10.895% [+7.386,+14.378]%** | **+14.372% [+10.477,+18.233]%** |

**Predeclared universal EPT rarefied taxon-richness transport gate was FAILED**: the 2013 original survey period interval includes zero. The earlier PR81 *counted EPT fraction* generalization still independently supports all three periods even on the m200 site subset, but the newly derived count-standardized **EPT taxon-ID diversity** does **NOT** demonstrate this three-cycle same-sign generalization with stable uncertainty. It is scientifically incorrect to equate EPT relative abundance, taxon richness and absolute insect biomass.

### Genuine ambient river N/P added AFTER physical habitat

| Heldout period | EPT fraction added N/P MSE improvement [bootstrap95] | Rare200 EPT taxon-ID richness added N/P [bootstrap95] | Rare200 **nonEPT** taxon-ID richness N/P [bootstrap95] |
|---|---|---|---|
| 2008–09 | +3.278% [+1.146,+5.279]% | **+9.287% [+5.066,+13.586]%** | −2.237% [−4.492,−0.046]% (predictive worsening) |
| 2013–14 | +0.709% [−2.166,+3.698]% | **+6.158% [+1.677,+10.431]%** | +0.698% [+0.084,+1.354]% (small) |
| 2018–19 | +0.952% [−1.612,+3.471]% | +4.454% [−0.216,+9.125]% (uncertain) | +0.351% [−1.249,+2.044]% (uncertain) |

**Critical biological distinction:** measured source ambient TN/TP holds more cross-period predictive information for **rarefied EPT taxon-ID richness** than for EPT counted fraction in the 2008 and 2013 target tests, even after field physical habitat is measured. This is NOT evidence of the direction or causal origin of effects: MSE reduction is prediction, not an N/P coefficient sign or experimental manipulation. NonEPT taxon-ID richness gains are weaker and sometimes negative, suggesting indicator-specific environmental information (but multiple exploratory outcomes and regional source taxa resolution remain).

**Additional physical habitat prediction for nonEPT** after measured N/P: 2008 +0.914% [−3.779,+5.402] uncertain; 2013 +2.953% [+0.555,+5.080]; 2018 +3.634% [+1.319,+5.927]. EPT and nonEPT are computed from the SAME original benthic specimens so these are not independent community datasets.

### Context: Within-cycle only, original 5fold disjoint HUC8

When each source era is both training and test but different HUC8, EPT rarefied richness PHab added after NP yields +17.654% (2008), +6.631% (2013), +11.725% (2018). Same-site EPT fractions +8.522%, +8.850%, +8.321%, respectively. These **within-era** results are descriptive comparisons and **not replacements** for externally heldout national survey-period transfer, which can fail under domain shifts.

## Truly independent numerical validation

[Independent fresh original EPA field source rebuild + SciPy `hypergeom.pmf(0)` mathematically distinct expected 200-organism taxon richness calculation, and manual weighted NumPy Ridge normal-equation reproduction of ALL 3 cycles ×3 outcomes ×4 model external MSEs and HUC8 bootstrap intervals SUCCESS #43568761030](https://github.com/zuizui0223/gesui/actions/runs/43568761030). Every original 2,209 rarefaction calculation and primary model comparison reproduced within source numeric tolerance. This is independent **code**, not independent field replication; all source data are the same real EPA survey organism counts.

**Original unknown-order taxon-ID sensitivity**: for ALL cycles, a maximally conservative hypothetical EPT **upper bound** includes every genuinely unknown-ORDER original source `TAXA_ID` as though it were EPT; this is not a claim those original taxa were sensitive EPT orders. The original rarefied EPT richness PHab predictive gain remains:
- 2008 **+13.900%** [ +6.921,+20.061]% instead of lower-source +15.773%; N/P +8.993% vs +9.287%.
- 2013 +3.853% [−2.372,+9.710]% (still uncertain) instead of +3.329%; N/P +6.576% vs +6.158%.
- 2018 **+13.894%** [+10.025,+17.703]% vs +14.372%; N/P +4.080% vs +4.454% (still uncertain).

Thus the predeclared **three-period universal rarefied EPT taxon richness hypothesis still fails**; the 2008 and 2018 positive PHab results cannot be attributed merely to unknown source order identification. No posthoc order assignment used for the primary model.

## Discussion, source-selection limits and ecological interpretation

The most defensible macro conclusion now combines [PR #81](https://github.com/zuizui0223/gesui/pull/81) with this new measured-source taxonomic contrast: **physical habitat carries historically transportable predictive information for the *relative EPT organism share*, but the transfer of count-standardized EPT *taxonomic diversity* is more context-dependent; ambient nutrients often provide stronger incremental information for taxon diversity than for organism share**. The two EPT endpoints are moderately to strongly correlated (source Spearman rho +0.731 to +0.776) but are different data functions. This indicates response-grain dependence of an environmental predictor, not a verified food-web production/consumption mechanism.

Important potential biases:
- **Taxonomic resolution drift and BERWW vs BERW protocols:** source `TAXA_ID` is not guaranteed species-level or identical rank across cycles. Even exact m200 rarefaction does not correct a different taxon-ID definition or sampling habitat protocol. 2008 source raw taxon+ORDER mapped >=95%, 2013/2018 same threshold, but 5% unknown source taxon orders yield bounded rather than exact classifications.
- **Source selection by sample individual counts:** original m200 source excludes 439 of the earlier 2,648 original source visits. Some may have been low individual counts because of actual ecological stress, so conditioning on N≥200 may preferentially exclude the most disturbed reaches. This is **not an absence of those taxa**; both EPT fraction and richness results are explicitly for this source-measured subset.
- **Small strict cross-period training in 2013 holdout:** only 355 original native UIDs /272 HUC8 across older 2008 and newer 2018 after removing all overlapping target HUC8, lower than previous PR81 train n400 threshold but passes the **new explicitly pre-declared rarefied-subset >=300** threshold. Generalization intervals still wide.
- **Predictive HUC8 bootstrap is descriptive fixed-error resampling, not repeated retraining, study-independent p-value or cause-specific effect CI**. Three outcomes across three cycles are exploratory multiple comparisons; don't selectively emphasize one positive while suppressing the failed universal gate.
- **Not wastewater:** actual historical same-source municipal effluent dual N+P monitoring at original NRSA visits failed old strict [PR #74](https://github.com/zuizui0223/gesui/pull/74) eligibility; ambient stream N/P are real measurements but cannot be assigned sewage origin. None of these data measures aquatic insects' emerging dry mass, aquatic insect resource export, genuine ecosystem services, phototrophic net primary production, butterfly pollination or water purification flux.
- Preserve earlier weak 2013–2018 riverbed×Chla effect **failure at completely new HUC8** in [PR #79](https://github.com/zuizui0223/gesui/pull/79). The 2008 source period has no original matching periphyton pigment/AFDM from official table and is not a direct measured multi-trophic cascade.

Keep DRAFT PR unmerged, do not alter main, no eBird and no TinyFish.
