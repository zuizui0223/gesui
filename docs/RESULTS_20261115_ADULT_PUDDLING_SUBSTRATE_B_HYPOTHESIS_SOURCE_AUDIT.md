# Adult resource use hypothesis B — before equating iNaturalist annotations with real butterfly feeding

Date 2026-11-15 JST. This is a separate, **exploratory** adult-resource evidence strand, not an amendment or a rescue of the frozen 2022 source-specific HNS test. All results of H7 PR #88–#95 are left unchanged. No new PR, no main/chocho merge.

## Biological mechanism, not species-HNS conflation

A genuine adult-resource mechanism concerns uptake of **mineral salts, amino acid/nitrogen-rich fluids or other solutes** from wet soil, human/vertebrate excretions, decomposing carcasses, or water (puddling). It is conceptually different from larval host plant habitat nutrient affinity (HNS), adult abundance, ant mutualisms and directly measured wastewater nitrogen in a receiving river.

The supported prior shows taxonomic and diet heterogeneity:
- Experimental adult puddling choice (Beck et al., 1999/2017 abstract indexed PubMed PMID 28308154): lycaenids preferentially attracted to *albumin* relative to NaCl, whereas papilionids/pierids more often visited sodium solutions; nearly all experimentally attracted butterflies were males; protein bait responses are not measured sewage nitrogen or wastewater dependency.
- Boggs & Dau 2004, *Environmental Entomology* DOI 10.1603/0046-225X-33.4.1020: species differed in preference among mud, herbivore dung and carnivore dung, sodium differs across natural substrates.
- Lamie et al. 2025 review *Annals of the Entomological Society of America*, DOI 10.1093/aesa/saaf007: reviewed puddling on mud, dung, carrion and sweat, with substrate, taxonomic and sex heterogeneity.

A dying fish near an urban river is not automatically a wastewater nutrient subsidy, and a butterfly standing nearby does not prove it is feeding. Natural carrion/stormwater, water temperature, oxygen, fish mortality and observer selection are alternative explanations.

## New original iNaturalist source fields (user-generated)

The public [iNaturalist observation field 500](https://www.inaturalist.org/observation_fields/500) is named **Butterfly Puddling Substrate**, with official description "Substrate on which butterfly is puddling". Its free-text value can contain mud, sand, gravel, dung or other names. Values are NOT standardized, not individually photo-verified, and do not quantify chemical composition or assimilated nutrient.

The more widely used [field 17428](https://www.inaturalist.org/observation_fields/17428) currently called **Lepidopterans aggregating/Roosting/mud puddling** admits non-ingestion aggregation, roosting and sometimes seeking nectar or water. It is **not** a positive puddling behavior field, nor is a blank field evidence of nonpuddling. The original field 16376 categorizes scavenging/decomposition vs nectar/pollen vs parasitism but applies to many taxa, and field 18609 "Behavior: Consuming carrion" is not butterfly-only.

The authoritative iNaturalist help states observation fields can be user-created and are not closely regulated; multiple fields capture overlapping behaviors. Original field-name/description audit from 2026-11-15 [GitHub Actions #43702388311 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43702388311) surveyed 97 unique fields and found 69 candidate keyword matched titles/descriptions. Keyword lists include unrelated dung beetles, bird roosts, prey comments and plant wounds; this is a **metadata search**, not 69 valid butterfly assays. Within the original field 500 a pilot API request reported **746 total field-tagged observations**, and the earliest 80 distinct returned records had 19 Lycaenidae, 27 Nymphalidae, 13 Papilionidae, 10 Pieridae, 7 Hesperiidae and 4 other/unclassified. These **80 are an ID-ordered biased subset**, not an independent random sample or family-level selection ratio. One observation owner accounted for 33.75% of this subset. The sample recorded 45 mud/soil/sand, 9 rock/gravel, 3 dung, 1 water-only, 11 ambiguous multiple categories and 11 unclassified text labels. No source-water linkage was tested.

## Reproducible rules and immediate stop tests

1. Original iNat API source field500 exact source IDs/values and standard taxonomic family ancestry; confirm JSON source receipt and all eligible records. This is a self-selected platform annotation subset, not an unbiased number of biological encounters.
2. Per-annotated-record candidate substrate categories **do not assume ingestion or chemistry**:
   - mud / soil / sand / clay / silt (substrate only; NOT Na measured);
   - gravel / rock / stone (substrate only);
   - dung / feces / animal droppings (organically enriched *possible* substrate, no measured N);
   - carrion / decomposing animals / fish carcasses (possible adult access, no measured protein intake);
   - urine / sweat / explicitly saline (other and chemically distinct);
   - water-only, multiple, free-text unknown.
   Any mixed/free-text ambiguous value is reported as ambiguous, not forced to favorable class. The category regex must remain inspectable and excluded from inferential biology until manual substrate/photo QA.
3. Before ecological claims, validate each putative positive's actual photograph/video for **adult proboscis making contact with identifiable substrate** (not merely resting), classifier confidence and life stage, and co-occurring behaviors. If photo cannot resolve substrate or ingestion, retain "annotated, unverified".
4. **No population propensity inference** from (field500-recorded Lycaenidae)/(field500-recorded butterflies) because there is no standardized survey denominator; 17428 and all iNat butterfly records do not provide comparable non-feeding controls. Differences may reflect photographer, project, time, geography and platform field usage. Observer concentration must be reported. Differences in documented type are at most **descriptive composition conditional on somebody applying an annotation**.
5. Do not merge adult puddling statistics with 2022 Europe HNS coefficients. To test source specific behavior later, need site-year verified receiving-water discharge/Na/nutrients and **same-standardized sampling effort** at naturally wet, wastewater-source, stormwater and no-registered-source sites, with source-to-site hydrologic direction and matched within-reach controls. Unreported wastewater source is not zero exposure.

## Separate testable ecological predictions

- **B1, substrate niche (adult):** Among verified active adult resource-feeding events (if enough exist), Lycaenidae may show a greater propensity for protein/organic-rich bait than NaCl bait compared with Papilionidae/Pieridae, as a *species/region contingent* hypothesis, not a universal species diagnosis.
- **B2, wastewater association (adult):** Verified resource-feeding use (per independently standardized observation opportunity) may vary by *measured* downstream resource concentration and wet reference substrate, not by crude WWTP proximity alone. Need nutrient/salt chemistry or source-exposure contrast.
- **B3, trophic transfer:** A fish carcass resource may attract adults directly regardless of wastewater; this cross-subsidy requires documented fish mortality/resources and feeding, not merely river fish, faeces or raw iNat photograph proximity.

**Status:** Original field metadata established. Additional source-tagged observation audit is nonconfirmatory and must honestly report sample/census limits; no biological sewage effect, no field protocol equivalence, no claims of adult nutritional dependency. No eBird or TinyFish.


## FULL source census achieved — actual 746 iNaturalist field500 annotations

**New original source result** — [GitHub Actions original-source census #43703011984 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43703011984); [SHA-checked source-result archive #43703509062 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43703509062). Reproducible [observation / taxonomic-family annotation reader](../analysis/20261115_adult_resource_behavior/audit_inat_puddling_substrate_original_obs.py), [exact original-field 97-name definition audit](../results/20261115_adult_resource_behavior/INAT_ADULT_FIELD_DEFINITIONS_SOURCE_ONLY.json), and [original observed 746 field500 source census with family×substrate, observer-concentration and HTTP receipts](../results/20261115_adult_resource_behavior/INAT_ORIGINAL_SUBSTRATE500_FIELD_SOURCE_AUDIT.json).

The iNaturalist public API at analysis time reported **746 observations whose original user-generated annotation field500 was populated**; **all 746 distinct original observation IDs were retrieved and source field-ID content individually matched**, not just the first 80. No photograph content, observation owner's personal identity, precise GPS or HNS source values are republished. Standard scientific family membership was derived from original iNaturalist taxon authority ancestor IDs. This is the **universe of that field's available public annotations**, NOT the universe of global butterfly feeding, puddling, or exposed habitats.

### Source field500 free-text substrate labels (all 746 annotations)

| Provisional raw-text substrate bin, NOT photo-validated behavior | Original annotations | % of 746 |
|---|---:|---:|
| Mud, soil, sand, clay or silt | **426** | 57.10% |
| Rock, gravel, pebbles | **73** | 9.79% |
| Dung / animal faeces or droppings | **46** | 6.17% |
| Carrion/dead-organism regex labels | **3** | 0.40% |
| Urine / sweat / explicit saline | **5** | 0.67% |
| Water only | **18** | 2.41% |
| Text clearly mixing two substrate classes | **41** | 5.50% |
| Unclassifiable free text | **134** | 17.96% |
| **Total** | **746** | **100.00%** |

There were **204 distinct unstandardized free-text values**. Many records include words like "road", "tree sap", "stream", "fruit" or "animal remains", which may not be captured by the conservative pre-written regex, so these bins are neither validated substrate types nor proof of actual nutrient absorption. Specifically, "animal remains" without a regex-matched carrion keyword can fall into unclassified; the n=3 in carrion is **not a biological prevalence estimate or a proven upper bound**. The **175 mixed/unclassified labels** represent 23.46% of source annotations. Photograph-based curation must be implemented before comparing actual organic-resource use.

### Family × provisional substrate from all original 746 observations

| iNat taxon authority butterfly family | All source-tagged | Mud/soil/sand | Rock/gravel | Dung/faeces | Carrion regex | Mixed/unclassified | Other annotated |
|---|---:|---:|---:|---:|---:|---:|---:|
| **Lycaenidae** | **200** | 137 | 16 | 9 | 0 | 33 | 5 |
| Nymphalidae | 300 | 143 | 33 | 24 | 2 | 88 | 10 |
| Papilionidae | 115 | 61 | 13 | 6 | 0 | 30 | 5 |
| Pieridae | 64 | 51 | 1 | 0 | 0 | 11 | 1 |
| Hesperiidae | 48 | 22 | 8 | 5 | 1 | 10 | 2 |
| Riodinidae | 6 | 4 | 1 | 0 | 0 | 1 | 0 |
| Other/unclassified | 13 | 8 | 1 | 2 | 0 | 2 | 0 |
| **Total** | **746** | **426** | **73** | **46** | **3** | **175** | **23** |

**Scientific limitation: DO NOT interpret dung 9/200 for Lycaenidae versus 24/300 for Nymphalidae as biological family substrate preference, nutrient requirements, or absence of carrion feeding.** The original iNat field is independently user-filled and incomplete, and most butterfly observations do NOT receive these fields. The current family denominator is a set of **annotated public sightings**, not standardized samples with substrate choice and available sites. The observation owner's original public IDs were used internally only for concentration QA; they are not republished. One observer accounts for **20.64%** and the **five most active observers 50.27%** of the entire 746-annotation dataset (**181 distinct observation-owner IDs**). Counts may reflect a handful of photograph-focused contributors and their locality/season.

### Two other fields are NOT independent confirmation

- Original field17428 ("Lepidopterans aggregating/Roosting/mud puddling") currently has **665 public matching observations** by API search (a stale website scrape showed fewer). Its name and definition explicitly allow resting, aggregation and even actively seeking nectar; **affirmative is not proboscis intake**. Only the first 30 older observations were checked for taxonomy; do not treat their 25 Pieridae etc as family totals.
- Field16376 covers all insects/organisms and has **24,479** source observations; only **20** oldest were inspected; this cannot establish puddling in butterflies without taxonomic restriction plus human photo verification. Field18609 "Behavior: Consuming carrion" had only **six** source records, none confirmed as butterflies by iNat current family ancestry. A different field on carcass use may exist; this is not a global absence statement.

### Result and next admissible ecological inference

**Original data establish a documented adult substrate-niche axis, not a wastewater subsidy effect.** Public annotations include several natural wet and organic substrates and a classifiable multi-family source cohort. But the expected Lycaenidae-specific organic subsidy from tropical bait choice experiments is **NOT replicated or falsified** by this self-selected annotation set, because platform users did not record all nonpuddling butterflies, substrate availability, sex, independent individuals, source nutrient concentrations or unbiased site visits. Evidence of feeding by a fish carcass also cannot be inferred from the field text alone.

A defensible next independent macro analysis would:
1. Obtain photo/video verification for original adult proboscis feeding and explicit substrate, including negative controls where adults are photographed but not feeding; annotate substrate as mud/sand, feces, verified carcass, urine/salt, fruit/sap etc; independent two-reader coding and uncertainty required.
2. Calculate **within-observer, season, location and availability strata** a difference in verified substrate choice among butterfly families *only where sampling opportunity for alternative substrates is documented*; do not report global puddling propensity from absence of field500.
3. Separately test whether observed **source-specific receiving-water nutrient/Na/organic-solute** concentrations and outfall-direction geography predict substrate-use within matched riparian sites; H7 geographical W index alone is not nutrient dose, and original 2022 fine-site butterfly study remains source-gated.
4. Keep larval HNS positive macro association and adult annotated resource use as **distinct testable paths**, not a causal path analysis without joint measurement.

This full census is explicitly **exploratory original-data annotation/source-quality work**, not a new confirmatory ecological test. All previous HNS and other null results remain untouched, and there is no sewage/biomass/ecosystem-services causal conclusion. Keep Draft PR #95; no new PR/branch, no merge to main/chocho, no eBird/TinyFish.


## Independent sparse v2 re-audit, same 746 original records — observer matching and photo eligibility

**2026-11-15 original data SUCCESS:** [iNaturalist v2 sparse 746-source audit](https://github.com/zuizui0223/gesui/actions/runs/43704447819), [exact SHA-verified source archive GH Actions](https://github.com/zuizui0223/gesui/actions/runs/43704784922). Original full source 746 distinct observations were re-fetched through the *separate* iNaturalist v2 sparse-field interface using only original field500 values, accepted taxon ancestor IDs, owner IDs in transient memory, observed year and native photo URL. The **previous v1 complete 746 family and original regex category counts were EXACTLY reproduced**, so this is genuine source validation, not 746 new observations.

The [complete inspectable v2-source QA JSON](../results/20261115_adult_resource_behavior/INAT_FIELD500_V2_OBSERVER_PHOTO_REVIEW_QA.json) records per-family revised field-text classes, source count and SHA receipts, observer-concentration sensitivities; [the **78-row publicly linked original-photo REVIEW QUEUE**](../results/20261115_adult_resource_behavior/INAT_FIELD500_ORGANIC_SUBSTRATE_PHOTO_REVIEW_QUEUE.csv) contains direct links back to each original iNat record/photo plus explicitly empty **human substrate**, **proboscis contact** and **independent second reader** columns. No images, geolocations, owner names or hidden original observation data were redistributed. **All 78 are NOT human reviewed.**

### Corrected free-text vocabulary (not biological confirmations)

The original "carrion=3" regex category had **false classifications**: labels such as "rotten fruit" and "decaying vegetation" were assigned to carrion or ambiguous mixed animal/mud classes, while "animal remains" (two source records) was missed. A refined, biologically explicit text distinction now separates:

| Revised source text category | Field-tagged observations | Interpretation |
|---|---:|---|
| Soil/mud/sand | **443** | Substrate; not chemical salt concentration |
| Rock/gravel | **70** | Substrate only |
| Dung/feces only | **46** | Possible animal-derived resource; no intake confirmed |
| **Animal-remains/carrion candidate** | **3** | **Two "animal remains" + one "carcass of amphibian"**; no verified adult feeding and *no explicit fish carcass in this set* |
| Fruit / sap / decaying plant material | **23** | Plant sugars/exudates/debris, **NOT animal carrion** |
| Mixed feces+soil | 3 | Both mentioned; do not force one source |
| Mixed fruit/sap+soil | 2 | Mixed substrate |
| Mixed fruit/sap+sweat/urine | 1 | Mixed resource, not safe single assignment |
| Mixed soil+gravel | 37 | Physical substrate mixture |
| Sweat/urine/explicitly saline only | 4 | Label only, no solute measurements |
| Water only | 18 | No chemical composition |
| Unknown/unclassified | 96 | Free text insufficient |
| **Total** | **746** | Original annotated records, not ecological prevalence |

Note: the revised number "3 carrion candidates" happens to equal the old "3 carrion" *by coincidence*; the **identities and classification logic are different**. "Animal remains" can be an ambiguous visual label; photographs must establish actual carcass and proboscis contact before promoting it. The full raw-to-revised source vocabulary is retained in the JSON. Revisions to 182 source-record classifications are an **exploratory text recoding**, not 182 newly discovered behavior events.

### An observer-conditioned limit on family differences

Across all 746 source-tagged public observations, **181 distinct public observation owner IDs** were found in transient memory; top one contributed **20.64%** and top five **50.27%**. **Only 18 owners recorded *both* Lycaenidae and Nymphalidae in the same field500 census, and only 8 owners had both family names plus at least one "organic candidate" and one mineral substrate tag**. Even those eight do not necessarily have a within-owner, within-site, within-season substrate choice experiment, so **a matched preference coefficient would be severely unsupported**. Do not interpret the source cross-table as taxon-specific organic nutrient uptake.

Descriptive sensitivity, **conditional on voluntary annotations only**:

| Family | Any revised organic candidate label, all tagged records | Exclude single most frequent tagger | Exclude five most frequent taggers |
|---|---:|---:|---:|
| Lycaenidae | **11/200 (5.5%)** | **9/141 (6.4%)** | **7/70 (10.0%)** |
| Nymphalidae | **51/300 (17.0%)** | **32/227 (14.1%)** | **25/138 (18.1%)** |
| Hesperiidae | **6/48 (12.5%)** | 6/47 (12.8%) | 2/36 (5.6%) |
| Papilionidae | **8/115 (7.0%)** | 8/109 (7.3%) | 6/75 (8.0%) |
| Pieridae | **0/64 (0%)** | 0/53 (0%) | 0/41 (0%) |

**These percentages are not ecological proportions of puddling, sex-specific feeding propensity, nor substrate-choice preference.** They are a diagnostic for contributor imbalance within source-annotated observations only. Excluding top contributors is also posthoc and does not create representative negative controls. No p-values, family biological effect estimates, nutrient assimilation or WWTP association are claimed.

### Photo verification gate, not completed in this step

All **78** amended organic candidate annotation records return at least one public original photo URL, including **3 ambiguous carrion candidates, 46 dung-only labels, 23 plant sap/fruit, and six mixed organic labels**. A human or visually verified two-reader review must separate:
1. substrate plainly visible / uncertain / off-screen;
2. butterfly alive/adult and taxon identification credible;
3. **proboscis tip visible and touching named substrate** (versus mere perching/resting);
4. source class actually animal scat, animal remains, plants, soil, or ambiguous;
5. relevant duplicate photographs, one individual reported more than once, shared observer/location;
6. site exposure only if independently measured/verified source chemistry and spatial direction become available.

The 78 original iNat links are a **review queue, not 78 confirmed feeding observations**. Do not download and republish users' copyrighted photos in the public GitHub repository. No original sewage water quality, source-to-river incidence, host plant quality, or biomass measurement is present. In particular **zero text-labelled fish carcass candidates does not establish absence of fish-carcass feeding** in nature.

The former 2022 eight-arm Lycaenidae-HNS primary remains separately STOP-gated, original H7 macro positive stays observational, and no claim of adult wastewater resource dependency has been established.
