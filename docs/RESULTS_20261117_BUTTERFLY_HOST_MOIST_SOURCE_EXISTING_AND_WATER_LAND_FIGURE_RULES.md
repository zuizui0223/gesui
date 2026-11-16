# 2026-11-17 source-verification: existing larval-host moisture F in the aquatic–terrestrial wastewater model

This audit directly addresses the attached figure containing (i) reported wastewater source load/river flow → measured river total P → benthic ecological condition, (ii) agricultural/deposition/urban nonpoint environments → plant Ellenberg nitrogen affinity → butterfly host nitrogen affinity, and **dashed new water↔land paths**. The figure suggests adding a butterfly host moisture/hydric indicator to test whether wetland-host species are specifically related to the river. Key correction: **the score already exists in the current project original 303-butterfly table**; what is not yet measured is the *source-conditioned ecological/hydrologic connection*.

## Byte-authentic source and executed results

[Source-only original 303 accepted host scores and 454 H7×period test Actions #43896037393 — SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43896037393). Original immutable butterfly roster `results/20261115_sync/GESUI_303_HOST_NUTRIENT_SCORES.csv` source SHA256 `00c1a36a9237a3d1a6d77334aa8baa5dc8b329a49925c2e31f033fa2d302ecdb`; original 454 authentic same European H7×period science panel SHA256 `1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca`, **232 H7,90 parent H6,15 countries and 2 fixed periods**. Butterfly sample identities remain 303 original speciesKeys including non-ACCEPTED, never recoded as accepted.

| Original source support | Result |
|---|---:|
| Full original butterfly keys, 5 families | **303** |
| ACCEPTED, no used-alias conflict | **290** |
| Strict accepted NUTRI-host score finite | **186** |
| Strict accepted host MOIST score finite | **184** |
| **Strict accepted BOTH at butterfly mean level** | **180** |
| Lycaenidae strict N+MOIST | **44** |
| Nymphalidae / Hesperiidae / Pieridae / Papilionidae both | 85 /23 /18 /10 |
| Descriptive butterfly host MOIST ≥7 (not a preregistered wetland diagnostic) | **14/180** |
| Among MOIST≥7 Lycaenidae | **2 species only** |
| Pearson, same 180 butterfly means of N and MOIST | **+0.30754** (descriptive, not river exposure) |

Actual examples from the original table: Lycaena dispar mean larval-host MOIST=8.90 / NUTRI=7.15, Polyommatus icarus MOIST=3.733/NUTRI=3.171; counts/coverage of actual host entries may differ per trait. Scores are **plant habitat indicator means**, NOT measured larval feeding on moist plants, uptake of sewage N, adult puddling, water-runoff hydrologic connection, host abundance, population increase, pollination or food-web services.

### Previously tested plant CWM MOIST is not butterfly host MOIST

Tichý et al. 2023 original EIV N vs MOIST source paired *plant species* accepted exact group n3,722, same source 454 H7×period primary gate; pairwise scientific [original results](https://github.com/zuizui0223/gesui/blob/main/docs/RESULTS_20261115_EIV_NUTRIENT_SPECIFICITY_FIVEPAIRS.md) show plant NUTRI retains H6-heldout butterfly host N score prediction after **plant CWM MOIST** (+33.95%) but adding plant CWM MOIST beyond N gives only +1.46% and uncertainty across zero. It is scientifically incorrect to claim this is a test of **butterfly host-specific moisture score** or river TP×host moisture interaction.

### Hard source STOP, why not immediately fit dashed arrows

The existing `host_moisture_mean` and `host_nutrient_mean` are pre-computed from the butterfly's documented hosts but may use **different subsets of host species with numeric original EIV MOIST vs NUTRI**. Existing `trait_coverage` and `n_scored_hosts` primarily reflect N-host scoring; 180 butterfly means with both numeric values does **not** prove identical individual host taxa scored for both. Before any statistically identifiable N versus MOIST comparison on butterflies, recover the frozen original butterfly×host-species EIV source and calculate both means using **exact same double-scored host species** and explicit coverage. Do not convert unknown moisture host scores to zero or split into just F≥7 vs others with 14 putatively wet species (2 Lycaenidae). Use MOIST continuous, with family/genus structural controls and taxon coverage.

The other crucial gap is sampling geography: original general EU H7×period 454 (2017–19/2020–22) butterfly **GBIF occurrence composition** and the genuinely 2022 Austrian 67 river TP+benthic source have NOT yet been proven same physical catchments, site-year observations, or actual same receiving wastewater. The earlier EU H7 exact river TP to butterfly host N or plant NUTRI did not show strong regional external prediction. These data CANNOT justify a directed river benthos→butterfly larval host mechanism.

### Revised research hypothesis / figure semantics

**Better cross-boundary hypothesis:** conditioned on independent source wastewater and nonpoint agriculture/urban/natural site covariates, **river-associated habitat moisture/connectivity** might moderate the association between measured river quality and the *composition* of butterflies using moister versus drier larval host plants. Formally treat water indicator EQR and butterfly host trait scores as potentially **parallel responses to connected river landscape**, not EQR causing the butterfly larval host score. Actual river P to vegetation depends on flooding/groundwater access; aquatic emergence to riparian spiders is a *different* trophic transfer route that does not automatically involve butterfly larval hosts.

Proposed separate test hierarchy:
1. source host double-scored overlap and physical H7×year×station identity match gate BEFORE any effect;
2. H6-group holdout 2020–22 water/terrestrial conditional association on genuine same H7 cohort, with continuous host MOIST and N from SAME larval host taxa, plant moisture/N from SAME plant taxa and effort/family/period/country-method controls;
3. site proximity and independently measured riparian inundation/floodplain connectivity as a precommitted interaction, compare with equally wet but unexposed rivers. Direct NN+aquatic ecological pathway fit prohibited if source year/receiver/hydrological path unmatched;
4. true actual local butterfly abundance/larval feeding, root/leaf elemental concentration, or adult resource access required before food-web nutrient-transfer or wastewater-dependence claims.

**Observational outcome from this QA:** numeric host MOIST source exists, **184 strictly eligible** butterfly species and **180 with both N+MOIST at butterfly level**, but empirical cross-water response, paired per-host trait coverage and wetland-specific sewage effect remain UNTESTED. The 2022 Austria67 S2 natural land near/outer technical gate failed 5/10 against frozen 6/10 and the whole-upstream RiverATLAS clipped source fails 0/67; keep these negative source gates unchanged. Do not make a figure arrow look like a supported causal path.
