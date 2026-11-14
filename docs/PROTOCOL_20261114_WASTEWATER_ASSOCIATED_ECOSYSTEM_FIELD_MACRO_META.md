# Testing a *wastewater-associated ecosystem* rather than wastewater dependence

**Field–macro–meta multi-level test, 14 November 2026. Design note only, NOT executed or approved for field sampling.**

## A. Working biological definition (operational, falsifiable)
A **wastewater-associated assemblage** is a replicated location-specific assemblage of plants, aquatic biota, visitors/consumers, scavengers and predators whose **functional interaction opportunities and observed resource-use frequencies** consistently differ near **verified wastewater discharges** from both upstream and independently matched *non-wastewater, equally wet/organic-rich* reference sites. Source identity, distance from discharge, effluent fraction, time/season, and site-specific resource substrates must be recorded.

This is *not* proof of evolutionary adaptation, obligate dependence, adult survival or successful pollination. No one taxon, especially Lycaenidae, is required to show a uniform positive abundance effect. A location can have more adult puddling but less larval recruitment and lower EPT biodiversity.

### Four distinct channels to measure, do not force a single causal pathway
1. **Adult direct mineral/organic intake:** wet surface, Na+, ammonium, dissolved organic nitrogen, protein-rich decomposition. Observe butterfly species/sex and unambiguous proboscis-to-substrate feeding; group into mud/wet soil, flower, carrion, dung, resting. **Fish death is not evidence sewage caused death**, and presence near fish is not confirmed ingestion.
2. **Plant-mediated:** persistent moisture, nutrient preference and flowers, plant host cover and leaf traits, actual larval occupancy/feeding/survival; floral visitation is not confirmed pollen deposition or seed set. Ant-tending Lycaenidae are a *species-level* secondary interaction, not universal.
3. **Aquatic to riparian:** EPT/chironomid larvae, adult emergence **dry mass per unit trap effort**, web-building spider diet, aquatic prey fraction, birds/bats where observed; compare food-resource subsidy with exported contaminants.
4. **Stress/selection:** measured DO, salinity, pharmaceuticals/metals/microplastics if available, aquatic taxa tolerance and habitat conditions; these may oppose the potential resource subsidy.

## B. Comparative butterfly question
A family-level distribution is too coarse to measure puddling. All five previously audited GBIF family categories are **comparators**, not species-level ecological syndromes. Include Lycaenidae, Nymphalidae, Pieridae, Papilionidae and Hesperiidae in the species×behavior observations; consider Riodinidae in regions where present. Record family *and species*, sex, behavioural substrate, flower resources, morning/afternoon, weather, observer, exact effort and available wetted surface area.

Predeclare interaction estimand rather than simple family prevalence:

`P(verified puddling | butterfly present, site, family/species, sex, substrate, survey opportunity)`

Ask whether `Wastewater exposure × measured Na/N/organic resource × family (or taxon)` modifies *behavior-specific visit rate* relative to natural wet controls. If species abundance differs, use **two-part** models for (i) any butterfly/visitor detection and (ii) feeding fraction conditional on detection. Species nested in family and geographic site, photographer/survey day as observation effects.

An adult puddling resource hotspot can be an **ecological attractor** without being a high-quality breeding habitat, so separate (i) adult direct feeding, (ii) floral nectar visit, (iii) larval host and abundance, (iv) ant attendance, (v) realized reproduction.

## C. Nested same-scale measurements, instead of mixing incomparable units
**Patch 1–10 m** (resource substrate and individual behavior) → **site/reach 50–500 m** (water quality, flower/host availability, aquatic emergence) → **independent discharge-receiving catchment** (site/source identity, landscape/urban/crop/drought) → **H7/H6 and country** (macro geography). Never assign basin-outlet modelled discharge as a measured dose to a specific butterfly feeding site.

A feasible *pilot* is 4–6 independently identified discharge systems, each observed over repeated weather-matched days at three focal microhabitats: discharge-associated wet bank, upstream background bank, and separately matched natural wet/non-sewage source. This is a pilot feasibility suggestion **not power-justified replication**; full sample sizes must be set after measuring baseline variance and independent-source count. A stormwater drain is a **different source class**, not untreated sewage. Multiple patches under one outfall are non-independent biological subsamples.

Exposure verification: source from operator or agency records, flow orientation and discharge season, water electrical conductivity, Na, NH4, NO3, PO4, DO, temperature, DOC/DON, direct source fraction/dilution if measured. Compare to controls with similar humidity, river size, shade, flowers, host plants, temperature and urban setting. A source reporting no wastewater cannot be assumed absent. Do not deliberately introduce fish carcasses or handle wastewater without field safety approval and appropriate PPE.

## D. Counterfactual contrasts
- Verified WWTP downstream minus upstream *in the same catchment*, and independently matched natural wet-bank down–up contrast to separate general moisture/water effects.
- High vs low effluent fraction, **only with verified flow direction and dilution**; measure distance gradient from the actual outfall, not only radius from a plant facility point.
- Naturally occurring decomposition and carrion vs clean wet soil at similar moisture/Na separately; no artificial exposure to pathogens.
- Dry vs rainy season (explicitly preserve sampling dates and treatment operation), urban vs rural receiving water, high vs low biofilm and host abundance.

## E. Analyses
### Field interaction model
`log E[N_feed_patch] = log(survey_minutes * accessible_wet_area) + β0 + β1 source_class + β2 log(Na) + β3 log(DON) + β4 family + β5 source_class:family + β6 season + u_outfall + u_day + u_species`

Use a defensible hurdle/negative-binomial generalized mixed model, or direct visits per actual observer effort with detection if visitor visibility differs; never infer female/male ratios from raw opportunistic photos. Test family-specific effects with partial pooling and multiplicity, report uncertainty and all five families. If adult feeding counts rise but eggs/larvae/ant attendance don't, label **adult aggregation only**.

### Cross-layer network linkage
At the same reach×survey interval, observed edges should count actual `plant→herbivore`, `plant→flower visitor`, `ant↔lycaenid larva`, `aquatic-emergent insect→spider` and other available prey interactions. Co-occurrence alone constitutes a **potential** link; confirmed feeding/visitation requires video, predation record, emergence-to-web sampling, diet and ideally chemical/isotopic source tracing for actual material transport. Assess context-specific edge density, beta diversity of interactions, and relative use of aquatic/terrestrial resources. To infer mediation, all edges must be measured in compatible site/time/source units; don't synthesize a chain from independent regional regressions.

### Independent-study meta analysis, different estimand
Include paired up/down, before/after, control-impact and BACI primary studies. Build one row per **study×independent discharge×contrast×response×time**, with shared-site dependence tracked. Pool **within comparable response families only**: ln response ratio for positive abundance/biomass/density, appropriate count-rate/odds effect for feeding/visitation, standardized contrast for biochemical/isotopic or EQR outcomes (ideally keep EQR method stratification). Study random intercept and river/outfall hierarchical effects, cluster-robust uncertainty, correlated contrasts, treatment type, effluent dilution, climate and ecosystem stage as moderators. Literature about nutrient preference *without actual wastewater sources* belongs to mechanisms, **not** direct wastewater effect pool.

## F. Independent verified literature anchors / eligibility roles
| Source | Data design | What it actually shows | Meta eligibility |
|---|---|---|---|
| Beck et al. 1999, DOI 10.1007/s004420050770 | Borneo two environments, protein vs salt bait and butterfly families | Taxonomic differences in adult nutrient bait attraction; not sewage | Mechanism only |
| Boggs & Dau 2004, DOI 10.1603/0046-225X-33.4.1020 | Butterfly species×puddling substrate | Na/substrate preferences vary by species | Mechanism only |
| Papilio polytes ammonium ^15N work, PMID 22842938 | Stable isotope/ammonium feeding experiment | Exogenous NH4 N incorporated into reproductive/somatic components in swallowtail | Mechanism only; excludes Lycaenidae-exclusive nitrogen mechanism |
| Whole-stream effluent addition, DOI 10.1007/s10750-025-06103-w | Treated effluent BACI/paired field experiment | Biofilm standing stock ~2.93×, shift towards green dietary path, food chain length unchanged | Direct wastewater stream effect, one independent experiment |
| Municipal WWTP aquatic–riparian microbiomes, DOI 10.1016/j.envpol.2025.125948 | Aquatic insects and riparian spiders upstream/downstream 3 WWTPs with isotopes | MWWE-associated food web transfer / microbiome alterations | Direct but different outcome; extract site-level contrasts |
| Neonicotinoids stream→spiders, DOI 10.1021/acsenvironau.5c00021 | Effluent-dominated Iowa stream and riparian spiders | Neonicotinoid transfer to predators; strong chemical tracers | Direct occurrence/transfer; not density/biomass meta |
| Microplastics insect→riparian spiders, DOI 10.1016/j.envpol.2026.129116 | WWTP-associated aquatic and riparian consumers | Microplastic transport across water–land interface | Direct contaminant transfer; not a resource-subsidy biomass outcome |

Do not conflate installation of treatment plant **reducing raw wastewater** with adding tertiary-treated discharge to a pristine stream (opposite interventions). Neither water nutrients nor host/pollinator records should be mislabelled as true plant pollination services or butterfly dependencies.

**Study goal:** describe *which specific resource/stress channels reorganize biotic interactions* around wastewater, including partial positive/negative signals, not insist on universal increased diversity or "sewage-dependent butterflies".