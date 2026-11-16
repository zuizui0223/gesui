# US NRSA2018–19: aquatic functional composition along real river TN/TP — source-locked exploratory macro result

**2026-11-11 | Draft PR #58 | POSTHOC independent original-data analysis.**

## Scope (what has *not* been shown)
This study **does not establish wastewater causing aquatic-insect changes or terrestrial ecosystem services**. EPA NRSA measured **river water** TN/TP, not the source-specific municipal discharge fraction. Only aquatic **benthic larval sample** relative composition was analysed; no adult insect emergence, aquatic→terrestrial subsidy, pollination visits, butterfly abundance or plant reproduction was measured. All regressions and sensitivities are observational; the additional strict-flag subsamples were inspected after original effects were known.

## Source provenance and pre-outcome freeze
- Official EPA 2018–19 original benthic, site and water chemistry CSV SHA256 pinned in [source contract](AQUATIC_TERRESTRIAL_FUNCTIONAL_COUPLING_SOURCE_CONTRACT_20261111.md).
- Original site+VISIT_NO join source gate **STOP** on two duplicate keys. [Native EPA UID correction before opening any nutrient-vs-taxon association](NRSA1819_NATIVE_UID_SOURCE_KEY_CORRECTION_BEFORE_OUTCOMES_20261111.md) is fully disclosed; each duplicate composite key contains two **different valid native UIDs**.
- Correct UID source manifest **1,221 visits / 1,125 unique physical sites / 48 states / 9 ecoregions / 729 HUC8**, original file SHA256 and frozen manifest SHA256 **`fbf53b8379d4e6564dbb51b3705b02fdd03e4a2caacad06a2314d643c8412e91`** before opening numeric outcomes. [Source-only frozen manifest SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43034107486).
- NRSA 2018/19 benthic `BERW` source sample: 326,938 individual counts in `TOTAL300` subsamples. Measured water TN literal unit **mg N/L**, TP **µg P/L** (divided by 1000 to mg P/L before modelling), original NARS flags retained. No zero-fill or invented effluent mass.
- 1,220 visits / 1,124 physical sites after river-order + conductivity original complete case; state and ecoregion, lat/lon, year, stream order, conductivity adjusted. Deterministic 5-fold held-out groups by original HUC8. HUC8 clustered descriptive intervals. Original numeric source flags were not ignored: further stricter sensitivity below.

## Original held-out and HUC8-cluster association
Response is normalized **individual fraction**, e.g. EPT counted individual fraction (not biomass, emergence or species richness).

| Response | Original measured stream exposure | Partial slope per 10× nutrient (fraction) | Descriptive HUC8-cluster95 | Incremental HUC8-disjoint 5-fold heldout MSE gain |
|---|---|---:|---|---:|
| EPT individual fraction | TN alone | **−0.08196** | [−0.11645,−0.04748] | **+1.772%** |
| EPT individual fraction | TP alone | **−0.06930** | [−0.09790,−0.04069] | **+1.840%** |
| EPT individual fraction | TN + TP joint | **TN −0.05032**, **TP −0.04775** | TN [−0.09160,−0.00905], TP [−0.08222,−0.01327] | **+2.354%** (joint) |
| Chironomidae individual fraction | TN alone | +0.05074 | [+0.00972,+0.09177] | +0.209% |
| Aquatic-insect larval fraction | TP alone | −0.07166 | [−0.10656,−0.03676] | +2.167% |
| Odonata individual fraction | TN + TP joint | uncertain | both spans zero | −0.220% |
| FFG shredder/collector/filterer exact-single-code fractions | TN+TP joint | unstable/no robust OOS gains | — | −0.578 to −0.218% |

**Crucial:** 18.5% of total original individual-count mass has **multiple FFG assignments** (e.g. `CG,SH`). The initial exact-single-code shredder/collector estimates therefore undercount assigned functional affinities and are **not rigorous negative tests of their ecosystem services**. Do not interpret these as true functional nulls.

Stream TN and TP original log concentrations are correlated **r=+0.640**; their partial coefficients cannot be assigned to wastewater N vs wastewater P sources. Main result is a repeated, small **water-nutrient vs EPT compositional negative association** across many survey units, not evidence of reduced insect emergence or of sewage-induced ecological change.

## Completely separate original-data HUC8-cluster OLS reproduction and additional outcome-seen sensitivities
The second Python/NumPy implementation **reproduced** the main joint model EPT TN and TP coefficients to numerical precision, independent of the primary model's code. Conservative posthoc subsets:

| Fixed assessment | n | joint TN beta (95% descriptive) | joint TP beta (95% descriptive) |
|---|---:|---|---|
| All eligible | 1,220 | −0.0503 [−0.0916,−0.0090] | −0.0477 [−0.0822,−0.0133] |
| **Both NARS nutrient flags blank** | **969** | −0.0650 [−0.1137,−0.0162] | −0.0415 [−0.0815,−0.0016] |
| Count subsample **≥250** | **1,005** | −0.0330 [−0.0773,+0.0114] | −0.0595 [−0.0969,−0.0221] |
| Both conservative restrictions | **796** | −0.0397 [−0.0930,+0.0136] | −0.0701 [−0.1144,−0.0259] |
| **First visit only** / distinct `SITE_ID` | **1,124** | −0.0540 [−0.0959,−0.0120] | −0.0430 [−0.0782,−0.0078] |

**Revised evidence rating:** the **TP–EPT negative relation persists across all these posthoc source-quality restrictions**, while TN is less robust when subsamples are restricted to at least 250 counted individuals. This is not independent prospective biological replication or causal source attribution.

Original [2026-11-11 actual numeric EPT/Odonata/Chironomidae/FFG source and 5-fold model SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43035012808). [Original three-source re-download + independent EPT QA/stricter source restrictions SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43035675111). Original aggregate JSON and QA ZIP are archived in that run's artifacts.

## Ecological interpretation and bridge criteria

The result supports **possible decoupling of aquatic biological quality/functional identity from simple resource quantity**: higher nutrient concentration correlates with smaller relative EPT composition; Chironomidae alone might be elevated under TN, but that conditional signal is weak and not robust in simultaneous N/P adjustment. Total adult export, guild-specific biomass, contaminant load and consumer access are not measured. No statement that the aquatic prey subsidy increased or decreased is warranted.

A complete across-ecosystem coupling test must use these **same river site coordinates and compatible years** for riparian host/flower plants and butterfly recording, then test detection effort, insect diet/predator, source-to-receptor hydrology, **and ideally real emergence/seed set**. Do not join unrelated European EEA river samples to US wetland plots as an inferred observed chain.

**Keep Draft; PR #58 stacked only on prior PR #57, no changes to gesui/main or chocho.**
