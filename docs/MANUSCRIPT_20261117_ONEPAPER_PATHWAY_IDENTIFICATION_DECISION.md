# 2026-11-17 | 下水研究の一本化：水陸経路別の識別性と次の実証試験

Status: **review-only synthesis / manuscript design**, not a new fitted ecology analysis, preregistration, or successful source acquisition. Based on source-checked gesui PR #112, #114, #115, #118, #119 and explicitly identified public literature. No main/chocho merge or changes to any previous STOP decision.

## 0. 一つの問いと、まだ主張できないこと

**General ecological question:** When can human wastewater-associated resource and contaminant gradients be detected in aquatic and terrestrial community structure, and what observable physical and trophic channels govern their propagation?

**Why wastewater specifically?** Effluent is a directional anthropogenic point input that may combine (i) water/baseflow maintenance, (ii) N/P/organic resource inputs, and (iii) chemical, thermal and contaminant stress. These processes may benefit some tolerant or resource-using organisms and harm sensitive ones. A facility map, reported annual N/P mass, received river TP, modelled effluent fraction, and species habitat-affinity indices are **not interchangeable measures**.

**Existing support:** two *within-compartment* observational patterns: actual receiving-water TP predicts same-river benthic EQR; plant-community nutrient and moisture affinity predict corresponding *documented larval-host habitat affinities* of observed butterflies. Austria additionally shows hydrologic-source-normalized reported load relates to same-water EQR **after** ambient TP and local landcover, but is not chemical-specific.

**Explicit NO:** no observed wastewater-to-riparian vegetation or wastewater-to-butterfly transfer in a single physically verified outfall–receiver–bank–consumer unit. Neither the absence of extra cross-domain predictive gain nor the technical inability to observe the riparian endpoint proves ecological insulation. **Never** infer a causal WWTP→TP→vegetation→butterfly SEM by pooling different populations.

## 1. Frozen evidence table: populations and effect interpretation

| Evidence / source | Observed unit, exposure / indicator | Key numerical result | Permitted interpretation |
|---|---|---|---|
| PR #112 (EEA 2018/2022) | 2634/2475 actual benthic stations, year-specific reported upstream N/P | 2018 EQR MSE gain after upstream N/P +1.19/+0.98%; 2022 +1.40/+1.10% respectively | modest out-of-region prediction; reported mass not true delivered dose; site/method/report selection |
| PR #114 (FR73) | 73 H7/24 H6, actual French 2018/2022 same-river TP and benthic I2M2 EQR; GBIF plants/butterflies at H7 spatial/3-year windows | logTP→EQR H6 gain +24.683%, beta −0.08455; logTP→plant NUTRI −2.979%; logTP→butterfly HNS after plants −3.899%; plant NUTRI→butterfly HNS after TP +30.209% | strong aquatic signal and terrestrial trait-composition matching; **not sewage-specific**: only 4 physical WWTP-source positive river-station pairs /2 H7/2 H6 |
| PR #118 (Austria67) | original 67 actual 2022 BMLUK TP + EEA EQR stations, EEA reported upstream annual P/N, hydrologic Q, WorldCover2021 1 km | WorldCover→EQR +26.582%; ambient TP after WorldCover→EQR +23.518%; upstream P/Q after both +8.388%, 5/5 holdout regions; upstream N/Q +12.199%; P/Q after N/Q −1.317%; N/Q vs P/Q Spearman 0.9453 | *composite source/flow/landscape predictor*, not phosphorus-specific or a causal load effect; only 5 geographical blocks, all posthoc |
| PR #118 terrestrial source smoke | same fixed 67 river monitoring stations, 2022 cloud-masked S2 natural land 250m vs 250–500m | 67 scene metadata source PASS; prelocked 10-site true pixel pilot **5/10 vs required 6/10** (STOP) | no fitted NDVI, no terrestrial vegetation/biodiversity/biomass response |
| PR #119 strict matched HOST | 447 EU H7×period /89 H6; 162 butterflies from **same larval host species** with pair-scored EIV N/F; plant EIV scores | plant MOIST→butterfly host MOIST conditional H6 MSE +22.005% [8.68,34.19]; EQR after terrestrial traits +0.586% [−2.41,3.09]; static W −0.876% [−3.59,1.76] | terrestrial community/habitat-affinity covariation, shared EIV reference and uncertain observational support; not actual feeding, abundance or aquatic-to-land transport |
| PR #115 near-outfall geometry | 16 original WWTP outlets, directed 1.5km up/down river; GBIF plants | geometry 16/16 PASS; 250m both banks source opportunity 0/16 in either 3-year window; 1000m 2020–22 10/16 but 2017–19 4/16 | no stable repeated local terrestrial biological response panel; 1km observations may not be riparian |

Metric comparison rule: an MSE gain from one holdout scheme/cohort is **not** the same estimand as a gain from another; do not compare these percentage magnitudes as a cross-taxonomic meta-analysis. H6-year/period panels are not the original Austrian physical monitoring points.

## 2. Biological pathways: choose one actual observable transfer route, not an attractive diagram

### A. Aquatic exposure to aquatic community status — observed observationally
Reported annual upstream N/P, Q, local land cover, measured ambient TP, benthic EQR. Source load divided by Q may encode several variables (N, P, co-contaminants, chronic upstream pressures, river morphology). Observed prediction is not proof of direct N/P causality.

### B. Hydrological seepage/floodplain → riparian plants — **not observed**
Requires bank-to-water hydro connectivity, riparian soil/porewater chemical uptake, true current vegetation (composition, production, foliar N/P), directed up/down/control site-year. Modelled sewage W, 2021 WorldCover *landcover covariate*, and 2022 S2 scene metadata do not measure this. Preserve all STOP gates.

### C. Emergent aquatic insects → riparian spiders — established in *external literature*, **not measured here**
An identifiable food-web channel should quantify aquatic emergence biomass/taxon, spider diet or tracer/chemistry at matched banks, and source-normalized load, potentially modified by insect emergence and dispersal traits. Do not claim an EQR value measures emergence biomass, nor substitute butterflies' host-affinity for spider diet.

### D. Terrestrial vegetation composition ↔ butterfly **larval** host habitat affinity — observational support
Consistent nutrient/moisture habitat matching across H7, subject to habitat geography, plant/butterfly trait-reference overlap, identification and GBIF observation bias. No proof of current host use or actual trophic uptake. This is an important *terrestrial association* but not a direct aquatic food-web edge.

### E. Butterfly **adult** puddling / carrion / mineral resource — separate, currently untested hypothesis
Adult Lepidoptera can take up sodium/protein-rich substrate resources; butterfly families vary in puddling preferences. Experimental and review studies make Lycaenidae a plausible *trait-based prediction*, not a sewage-sensitive guild established by gesui. Fish carcasses, dung, wet sediment and sewage-influenced riverbanks are different substrate exposures. Test Na+, conductivity, dissolved N and organic materials at **actual feeding points**; quantify observed adult visitation per search effort vs control and sex/family, excluding larval host score from adult-intake measurement. Do not convert GBIF proximity to feeding or sewage dependence.

## 3. Reconcile scales, years, and causality

**Only valid common-grain claims:** (i) the Austria67 physical water-monitor point 2022 aquatic source→water chemistry/EQR comparison; (ii) the France73 H7 2018/2022 receiving-water/landscape-composition descriptive comparison with river monitoring stations nested within H7; (iii) the broader EU447 H7×period terrestrial plant/butterfly association. These are three different populations even when geographical labels are superficially similar. France GBIF 2017–19/2020–22 vs river 2018/2022 is a *window*, not same-day; source W static is not a repeated natural experiment.

**Design principle for a new actual cross-boundary cohort:** one WWTP outfall or matched river reach is the clustering unit; within the same year/season collect upstream, downstream and independently verified no-interposed-source control reaches, each with river chemistry, Q, physical riverbank connectivity, standardized terrestrial responses and effort/detection. Define one *primary terrestrial response* in advance (prefer actual emergent insect/spider resource transfer if objective is physical aquatic→land food-web transfer; alternatively vegetation/leaf N/P if objective is hydraulic nutrient transport). Exclude invalid/overlapping controls and avoid response-based site selection.

**Falsifiers before ecology:** 1) upstream reported WWTP presence and loads must not be reconstructed by downstream placebo; 2) a Q-only model or landcover/upstream-catchment model must not fully explain any incremental source term; 3) N/Q vs P/Q cannot be separately interpreted when collinear; 4) spatial support needs at least two truly independent regions plus cluster-level validation; 5) true near and far terrestrial sampling opportunity must pass before inspecting biological endpoints. Technical STOP is not a biological zero.

**Interaction of interest (future, NOT fitted):** terrestrial biological response ~ source load/Q + physically observed bank/connectivity + (source load/Q × connectivity) + source-independent covariates; if emergence channel, replace 'connectivity' with actual insect-emergence flux/dispersal. Report all main effects, heldout geographic predictive gains, blocked bootstrap on outfalls, negative controls and exact uncertainty. No mediation product until exposure, mediator and response are measured in the same observational unit, and causal assumptions are assessed.

## 4. What can be *one article* today, and what still cannot

**Conditional one-paper framing:** “Pathway-specific ecological footprints of wastewater-associated gradients across aquatic and terrestrial communities” is a *research question*, NOT a demonstrated causal finding. A paper may transparently report a layered observational evidence comparison, with primary verified aquatic source–condition results and separate terrestrial trait association plus predeclared cross-domain null prediction. Its central contribution would be **the difference between within-compartment predictability and failed transportability across unmatched proxies** — not a demonstration that effluent stops at the bank. The manuscript must keep populations/years/estimands separate rather than calling these a meta-analysis.

**Review-stage publication decision:** robust aquatic **associational** manuscript is presently easier to defend: source-reported N/P relative to Q and ambient TP/landcover explain *some* additional variation in river benthic EQR, chemical specificity remains unidentified. Cross-compartment causal wastewater–ecosystem-service article is **NO-GO** under current evidence. Do not claim beta+ or MSE+ as beneficial effects on biodiversity, biomass, fitness or services; EQR higher/lower is ecological status.

**Tentative four figures**:
1. Explicit geography/time/population and channel-based conceptual diagram (observed association solid, hypothesized transfer dashed, technically STOP separate).
2. EEA 2018/22 aquatic load→EQR and Austria67 actual TP, reported N/Q P/Q, landcover, 5-region heldout evidence.
3. France73 true river TP→benthic versus TP→terrestrial traits and plant→butterfly host N (separate estimands clearly labeled).
4. Strict same-host N/F: plant→butterfly spatial heldout and no incremental EQR/W; transparent nulls and no causal bridge.

## 5. Literature novelty audit (not analyses from gesui)
- Hamdhani, Eppehimer & Bogan (2020), FWB 65:1657–1670, https://doi.org/10.1111/fwb.13519: 147 effluent-fed-stream quantitative studies; water-quality emphasis; ecological harm and environmental-flow benefit context.
- Pereda et al. (2020), Environmental Pollution 258:113719, https://doi.org/10.1016/j.envpol.2019.113719: ecosystem manipulation shows complex subsidizing and inhibitory processes.
- Diesbourg et al. (2025), Environmental Pollution 372:125948, https://doi.org/10.1016/j.envpol.2025.125948: aquatic–riparian food-web microbiomes across effluent.
- Behrens & Bernhardt (2025), Freshwater Science 44, https://doi.org/10.1086/735929: emergence biomass and metals vary with community structure and wastewater-influenced reach.
- Mianecki et al. (2025), ACS Environmental Au, https://doi.org/10.1021/acsenvironau.5c00021: neonicotinoid transfer to riparian spiders.
- Wardlaw et al. (2026), Environmental Pollution 410:129116, https://doi.org/10.1016/j.envpol.2026.129116: 8 WWTPs, some taxon-specific up/down microplastic transfer via emergence, weak magnification evidence.
- Beck et al. (1999), tropical butterfly puddling, https://pubmed.ncbi.nlm.nih.gov/28308154/; Butterflies often use adult protein/sodium resources with family differences; not a sewage result.
- Puddling in butterflies: current knowledge and new directions (2025), https://doi.org/10.1093/aesa/saaf007: adult mineral and N foraging is a distinct behavioral link.
- Landscape aquatic–terrestrial spider fatty acid study (2025), https://doi.org/10.1002/ecm.70025: traits and landscape shape consumer subsidy connectivity.

**Novelty boundary:** “Wastewater changes water/land food webs” is already known at local field scale. A new **multi-region, same-source, directional, transportability/trait×channel** test would be novel only after adequate observation support, not by re-labeling separate existing results as a pathway.

## 6. Next action, frozen without changing failed source gates
1. Manuscript evidence check: independent table of numerator/denominator, units, timing, selection and ecological interpretation for each proposed figure; retain source URLs, branch SHA and signed Actions.
2. For future causal physical transfer, prioritize a *new*, non-outcome-selected same-outfall bank sampling frame with actual water chemistry and either emergence/spider or actual local plants; do not relax the fixed AT67 S2 near/outer 5/10 technical failure after inspecting outcomes.
3. If no source-confirmed land response, produce the aquatic observational article first, preserve PR119's terrestrial result as a distinct test and discuss hypotheses, not observed sewage transfer.

Linked internal evidence: [main integrated PR112](https://github.com/zuizui0223/gesui/pull/112), [FR73 PR114](https://github.com/zuizui0223/gesui/pull/114), [near-outfall STOP PR115](https://github.com/zuizui0223/gesui/pull/115), [Austria67 PR118](https://github.com/zuizui0223/gesui/pull/118), [strict-host PR119](https://github.com/zuizui0223/gesui/pull/119).
