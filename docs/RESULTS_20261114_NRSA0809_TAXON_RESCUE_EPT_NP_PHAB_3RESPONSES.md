# Third national survey NRSA2008–09: original rescued EPT counts, ambient nutrients and physical river habitat

**2026-11-14 | Genuine original source data and independent matrix QA COMPLETE. This is exploratory macroecology, no sewage-specific N/P source effect or ecosystem-service loss established.** Extends 2013/14 + 2018/19 true biological/biomass evidence [PR #79](https://github.com/zuizui0223/gesui/pull/79) without rewriting positives and failures.

## EPA published original data and false absence prevented

The official [EPA National Aquatic Resource Surveys source table](https://www.epa.gov/national-aquatic-resource-surveys/data-national-aquatic-resource-surveys) contains real **2008–09 site information, CHEMW, PHab, Benthic Count, Benthic MMI**, but NOT original matching 2008 periphyton Chla or AFDM stand-alone CSV. Thus **this period cannot independently repeat the 2013–14 /2018–19 EPT×Chla or AFDM×EPT result**. No water-to-land ecosystem service flux measured.

- [Original source-only retrieval and native-UID support SUCCESS #43562414140](https://github.com/zuizui0223/gesui/actions/runs/43562414140): official SiteInfo/water/PHab/MMI source n2,320 original native UID, original benthic count n91,342 rows across n2,121 UIDs, exact five-original-table overlap **2,121** UIDs. Raw publisher bytes original SHA256 pinned. *Not* 2,121 numerically valid ecology observations.
- [Original primary probability wadeable source-only method gate SUCCESS #43562599406](https://github.com/zuizui0223/gesui/actions/runs/43562599406): **1,299** eligible raw BERWW primary benthic source + CHEMW water primary + PHab + original site + published MMI native UID, **796 original HUC8 /48 states**. Original BOATABLE `BELGB` source and CHEMB chemistry explicitly excluded, original sample duplicates not merged.
- Published derived EPA 2008 MMI source field `EPT_PIND` (0..100 percent) only numeric in **286** of these original 1,299 source UIDs. Original strict published-only EPT model gate min300 and original unparsed categorical Strahler field both caused [failed strict numeric run #43562951528](https://github.com/zuizui0223/gesui/actions/runs/43562951528). This failure was **not hidden or reinterpreted as EPT=0**.

## Actual original source taxonomy recovered — source and independent publisher fidelity before effects

EPA authentic [2008 original Taxa List/Autecology](https://www.epa.gov/sites/default/files/2016-06/nrsa_0809_benttaxa.csv) contains **TAXA_ID and ORDER** for 1,323 unique original taxon IDs. [Independent source-only taxonomy join SUCCESS #43563239544](https://github.com/zuizui0223/gesui/actions/runs/43563239544):
- ZERO original TAXA_ID→ORDER conflicting identities, **98.9328% of original 333,488 positive BERWW TOTAL300 counted organisms** map to known EPA original order.
- **1,232** original primary BERWW source UIDs each have >=95% of originally counted individual samples assigned to unique source order. Source-nonidentified specimens are **unknown**, not artificially called non-EPT.
- [Exact original 2008 published `EPT_PIND` vs newly reconstructed `TOTAL300×ORDER` direct numerical concordance SUCCESS #43563366026](https://github.com/zuizui0223/gesui/actions/runs/43563366026): **279** same-site original valid EPT percentage records, median absolute discrepancy **0.00296 percentage points**, **none** with discrepancy greater than 1 point, **100%** inside the original unknown-order logical lower/upper count bounds. This is **a validated alternative ORIGINAL EPT counted-fraction construction**, not fabricated values or use of a late proxy variable.
- Source original `STRAHLERORDER` was actually coded categorically as `1st, 2nd,...,7th,8th+`, not absent, in both original site/PHab source. Source-verified monotonic labels first through seventh converted to 1..7; the single `8th+` has no exact finite order and is NA (not guessed 8). These source corrections were fixed before inspecting new ecological coefficients.

The original source prespecified separately an exploratory newly reconstructed EPT counted share primary, published EPT_NTAX secondary, and original MMI_BENT secondary. Each retains its own genuine original UID/quality support; whole HUC8-disjoint fivefold Ridge10, train-only scaling, baseline native stream geography/state/eco9/Strahler/season/conductivity plus measured physical habitat (XEMBED, midchannel canopy, canopy presence, overhang vegetation, human disturbance and true channel width), and actual unalerted source CHEMW TN/PT (both original µg/L converted to mg/L).

## Genuine 2008–09 measured ambient N/P and benthos: heldout spatial effects

The original [EPA 2008 direct native source model SUCCESS #43563594143](https://github.com/zuizui0223/gesui/actions/runs/43563594143) found strict **unalerted CHEMW chemistry + field PHab + native BERWW measured** n**586** original visit UIDs across **437 HUC8/48 states**, 560 physical SITE_ID. Original EPT_NTAX response was available on a smaller **n457 /357 HUC8/43 states**. Nutrient field N/P are **ambient measured instream**, not actual POTW wastewater loads or delivered source concentration.

| 2008 genuine source benthic response | Strict native visits /HUC8 | Adding log TN+TP AFTER actual PHab, HUC8 heldout MSE gain | Fixed OOF 2,000 HUC8 bootstrap descriptive95 | Actual PHab added AFTER TN/TP |
|---|---:|---:|---|---:|
| **Raw counted EPT individual fraction** (EPT / all original TOTAL300 counted) | **586 /437** | **+1.9252%** | [−0.3875,+4.0680]% (uncertain) | **+8.1526%** |
| **log1p(original reported EPT taxon count)** | **457 /357** | **+8.0350%** | **[+2.5376,+13.1520]%** | **+14.9493%** |
| **Original published benthic MMI score** | **586 /437** | +3.1022% | [−0.3581,+6.3079]% (uncertain) | **+17.1748%** |

Clustered adjusted observational partial slopes per +1 SD source log10 original ambient nutrient:
- 2008 **EPT counted fraction**: TN **−0.02255 [−0.04934,+0.00424]**, TP **−0.01670 [−0.04208,+0.00867]**, both interval include zero.
- 2008 **log1p(EPT taxon count)**: TN **−0.10907 [−0.18848,−0.02965]**, TP **−0.12010 [−0.20090,−0.03929]**. Both negative in this source subset after controlling measured bed/vegetation/space, with positive HUC8 heldout nutrient incremental prediction. This **does not** establish nutrients *caused* EPT taxa to vanish.
- 2008 **benthic MMI (0–100 scale)**: TN **−1.9271 [−4.0870,+0.2328]** uncertain; TP **−2.4774 [−4.4389,−0.5159]** conditional negative.
- Physical habitat predictive information is substantial across all three endpoints, so a one-factor universal nutrient causal account is not supported.

### Numerically independent source bound and QA

[Independent fresh raw EPA 2008 source rebuilt and manual NumPy Ridge/HUC8 sandwich QA SUCCESS #43563980256](https://github.com/zuizui0223/gesui/actions/runs/43563980256) independently reproduced all three true original endpoint model MSEs, nutrient adjusted coefficient values, HUC8 resample errors and original source conditions. For one smaller n457 categorical model the finite-sample HUC8 sandwich interval had a relative absolute discrepancy ~0.13% of CI half-width due to independent rank/correction convention, under documented **0.3% relative-CI tolerance**; do not falsely call interval bits identical.

Conservatively unknown-order BERWW specimens were **not assigned as non-EPT**. In the strict same original n586 cohort, *lower-bound* mean EPT share from known-EPT / all source counted individuals = **0.28921**; *upper-bound* if every original unknown-order specimen were EPT = **0.29309**. Under the upper bound, additional TN/TP OOS MSE gain +1.9006% vs lower bound +1.9252%; measured TN β−0.02167, TP β−0.01678 (both still with intervals spanning zero). The biological conclusion is not produced by assigning missing taxon order as zero.

## Three-cycle interpretation and limits

Genuine original 2008/09 evidence now adds a **third distinct field survey period** to prior independently sourced 2013/14 and 2018/19 EPA panels. Across periods, habitat indicators often provide useful prediction of relative benthic community condition and nutrition/physico-chemical co-gradients can be present, but **2008 lacks periphyton Chla/AFDM**, 2013/18 original EPT×Chla riverbed effect **fails on never-before-shared catchments** ([PR #79](https://github.com/zuizui0223/gesui/pull/79)) and 2008 EPT_NTAX has response-specific missingness. The three survey panels differ in protocols `BERWW` vs `BERW`, time, state/eco source and habitats; they are not three independent experimental replicates or national populations estimated without weights.

No measured WWTP-specific source N/P was joined, and the original historical discharge dual reporting in [PR #74](https://github.com/zuizui0223/gesui/pull/74) did not meet prespecified eligibility. No actual algal net primary production, dry emergent aquatic insect biomass, water treatment ecosystem service, terrestrial consumers or pollination. **Do not claim wastewater-induced taxa declines or a proven trophic ecosystem service cascade.** Keep positives, nulls, methods-fidelity and geographic transfer failures visible, DRAFT branch unmerged, no eBird.
