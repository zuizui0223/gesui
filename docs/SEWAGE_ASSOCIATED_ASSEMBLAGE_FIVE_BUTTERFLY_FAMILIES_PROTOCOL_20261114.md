# Sewage-associated assemblages: butterfly-family resource-axis comparison (14 Nov 2026)

Status: **exploratory POSTHOC protocol, not a pre-registered finding**. This is a **new five-family question** after prior Lycaenidae-only results; all prior nulls, false starts, source missingness and sign uncertainty retained.

## Scientific question
Does the spatial geography of registered/modelled municipal wastewater outfalls associate with differences in recorded butterfly-family composition, after controlling river-basin land cover and sampling opportunities? Why do Lycaenidae appear around wet, organic-rich outfalls? This is **not** an attempt to prove butterflies depend on sewage. Hypothesized processes are distinct:

1. **Adult direct use:** sodium, water, nitrogen-/protein-associated puddling; manure/carrion are independent decaying-resource controls. An adult visiting a wet outfall is not evidence of feeding, reproduction, pollination, or wastewater dependence.
2. **Plant-mediated bottom-up:** dry-season water subsidy, nettle/legume/other host turnover, flowering resources, larval recruitment. Species-specific host use and ant attendance cannot be represented by family-level GBIF counts alone.
3. **Aquatic-to-terrestrial:** aquatic insect emergence, scavenger/detritivore resources, spiders and other predators; a real trophic link requires emergence, diet, stable isotope or repeated same-site feeding observations.
4. **Stress:** oxygen depletion, contaminants, modified river bed, salinity and pathogen hazards; declines in sensitive lineages can coexist with aggregations of tolerant/resource-seeking taxa.

## Why Lycaenidae is a testable focus, not a privileged positive
Comparative puddling studies (Beck et al. 1999; Boggs & Dau 2004; 2025 puddling review) suggest Lycaenidae / some Nymphalidae / Hesperiidae often use protein/organic resources, while Papilionidae and Pieridae often use sodium; variation is large within families and latitudes. Lycaenid larvae are often plant-feeding and many are ant-associated, but Pieridae also favour nitrogen-rich host plants (Pellissier et al. 2012). **Nymphalidae** tests whether purported organic affinity is Lycaenidae-specific. **Pieridae** is a sodium-/host-plant comparison. **Papilionidae** is another sodium-puddling reference, **Hesperiidae** a broad adult comparison. A raw increase of Lycaenidae observations near urban wastewater is insufficient evidence of any mechanism.

References (not outcome data):
- Beck et al., *Mud-puddling behavior in tropical butterflies: in search of proteins or minerals?* https://pubmed.ncbi.nlm.nih.gov/28308154/
- Boggs & Dau 2004, *Resource Specialization in Puddling Lepidoptera* https://doi.org/10.1603/0046-225X-33.4.1020
- *Puddling in butterflies: current knowledge and new directions* (2025) https://doi.org/10.1093/aesa/saaf007
- Pellissier et al. 2012, *High host-plant nitrogen content: a prerequisite for evolution of ant–caterpillar mutualism?* https://doi.org/10.1111/j.1420-9101.2012.02555.x

## Exact preexisting source and inference grain (no new GBIF pulls)
**454 source-frozen European HydroBASINS L7 catchment × 2017–19 / 2020–22 rows** in 232 L7 basins; these are two periods of the same catchments, not 454 independent basins. Original family counts separately queried from **234 originally sampled catchments × 2 source periods × 5 families** under fixed GBIF checklist and same flags; source shard run #43204891739. The EEA + butterfly combined original cohort was saved by GitHub Actions run #43215754868. True Level6 containing parent geography (90 H6) is from original valid-only + explicit repair audit #43338311684. HydroATLAS upstream crop/urban/flow from original run #43225986249. All input provenance verified before model; no repull, no sign filtering.

### Primary biological comparisons (frozen before inspecting new five-family effects)
- **Lycaenidae / Pieridae record log-ratio:** log((Lyc + 0.5)/(Pier + 0.5)), all 454 valid originally selected rows.
- **Lycaenidae / Nymphalidae record log-ratio:** log((Lyc + 0.5)/(Nymph + 0.5)), same rows.

These are **compositional observation-index outcomes**, not true relative abundance, actual adult puddling or mating performance. Numerator and denominator may each respond; absolute increases cannot be deduced from ratios.

### Secondary outcomes (ALL, never cherry-picked)
- For each of Lycaenidae, Nymphalidae, Pieridae, Papilionidae, Hesperiidae: log((family count +0.5)/(sum(other four) +0.5)).
- Lycaenidae share as previously documented is retained as the pre-existing null reference. Butterfly counts have source recording effort and possible taxon-specific GBIF bias.

### Exposure and baseline
Exposure: log1p(HydroWASTE modelled upstream-to-*basin outlet* annual/day source flow / (upstream area+1)) from prior 454 original CSV. It is largely static, NOT exact receiver-specific wastewater delivered dose or wastewater N/P/Na/protein content.

Baseline: source country, river benthic method, period, area, EEA station count, river flow, GBIF all-animal / plant / butterfly opportunity, upstream cropland and urban % and natural flow. Static WWTP added to baseline. **This is a geographic association, NOT source-specific chemical mechanism.** Study outcome and covariate cohort is kept fixed across every family, and source errors are NA (never zero).

### Generalization checks and claims
- Primary evaluation: fivefold **true H6-parent-disjoint** heldout assessment, H6-equal MSE, plus H6-parent bootstrap for uncertainty. Both records of any L7 fall into the same heldout H6 group.
- Independent country-left-out transfer with no country fixed effects, equal-country diagnostics limited to countries with >=10 period-rows and >=5 H7; France-dominant study cannot be labelled globally representative.
- Report **signed** MSE changes, including worsening; no optimizing predictors or splitting families conditional on positive results.
- Direction: adjusted within-sample coefficient for source gradient is descriptive and cannot establish association in new countries if holdout fails; 2 primary comparisons get Holm multiplicity control if formal p-values computed (current script reports descriptive bootstrap, not formal tests).
- Falsifier: Lycaenidae/Pieridae fails H6 external or country transfer; Lycaenidae-Nymphalidae does not diverge; source effects disappear after upstream crop and urban controls; all five families shift together; 2017–19 and 2020–22 disagree. Such outcomes are scientific results, not permission to subset for significance.
- Data types: butterflies are GBIF counts, EQR is *river ecological quality*, AFDM is *biofilm standing organic matter*, richness is separate. Never infer a meta-analytic pooled effect, true population biomass, ecological dependence, trophic cascade, fish-mortality pathway, plant host consumption or actual pollination from this model.

## Next evidence layer if butterfly family contrast succeeds
Carry out **paired local field observations at independently verified WWTP outfalls and background wet sites**. Minimum data: same date/effort, butterfly species and sex, per-individual behaviour (puddling vs flower vs carrion vs resting), host plant and larvae, ants, water temperature/flow, Na+, NH4+, conductivity, DO, dissolved organic nitrogen, flowers, exposed moist area. Compare riparian untreated streams, natural mineral springs and carrion/humidity separately: **wastewater is not chemically equivalent to dead fish**. Do not use biological sewage or fish carcasses as experimental baits without safety approval. Measure larval recruitment and spatial survival before claiming habitat improvement. More general ecosystem interactions require aquatic insect emergence and predator/diet data.

## Meta-analysis bridge, distinct from this macro analysis
Independent published upstream-downstream, before-after, and paired-control studies will be extracted as study/site-specific effects with **hierarchical study/river random effects**: (i) butterfly puddling frequency (ii) larval host usage (iii) pollination visitation (iv) ant attendance (v) aquatic emergence (vi) predator feeding, (vii) richness and biomass separately. Do not pool coefficients, MSE improvements or data from overlapping sites as independent effect sizes.

Original draft PR evidence #63–83 is retained; do not auto-merge main or modify chocho, no eBird or TinyFish.