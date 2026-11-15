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


## Original source photo-ID reuse and observer-day audit (2026-11-15, no photos viewed)

**New independent code:** [audit_field500_observer_day_independence.py](../analysis/20261115_adult_resource_behavior/audit_field500_observer_day_independence.py). Verified [source run #43705857353 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43705857353) and immutable byte-checked [archive #43706046893 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43706046893), with full aggregate [JSON](../results/20261115_adult_resource_behavior/INAT_FIELD500_OBSERVER_DAY_PHOTO_DUPLICATE_SOURCE_QA.json). The 746 publicly annotated records were fetched independently from the sparse v2 original endpoint; original 746 family and 78 photo-review-queue IDs were source-verified. We did **not** fetch, inspect or redistribute any original photo **pixels**, exact coordinates or observer account names. The existence of a photo URL does **not** establish actual feeding.

### What an independent sample can and cannot mean

| Source metric | Actual result | Interpretation limit |
|---|---:|---|
| iNat observations with field500 | **746** | A voluntary, annotation-selected source frame, NOT a sampled population of butterflies |
| Original distinct observation owners | **181** | Not 181 controlled feeding surveys |
| Effective number of observation owners (inverse Simpson of owner record shares) | **13.965** | Concentration diagnostic **only**, not biological sample size or n_eff of independent carcasses |
| Owner × observation-date groups across 746 | **473** | Without source site and animal ID, cannot certify independent butterfly feeding events |
| Organic-labeled 78 original observations | **78** | Includes dung, animal remains, plant sap/fruit and mixed notes, not confirmed intake |
| Distinct owner × date groups among organic-labeled observations | **62** | A proxy for repeated reporting; same owner-day may include separate local encounters |
| Total photo ID tokens across all 746 (including secondary photos) | **1,642 distinct IDs** | Photos have not been inspected |
| Photo IDs referenced in more than one original observation | **3** | Clear **photo reuse** (not necessarily same insect/site) |
| Original observations affected by at least one reused photo ID | **6** | Do not treat every tagged observation as an independent visual behavior record |

The pre-existing 78-image first-photo review list contains one same-first-photo pair (rotting banana); the new all-photo metadata audit finds **3 distinct shared photo IDs across 6 different original observations**. This explains why counting observations rather than independent photographic events can overstate independent support. It does NOT prove only three biological events are duplicated; separate images can still document the same event, and image reuse can have nonbiological reasons.

### The three "animal-remains" annotated observations are not three confirmed carrion-feeding events

| Public original observation | Original user-entered field500 value | Source owner-day comparison | Photo/proboscis validation |
|---|---|---|---|
| [#150764924](https://www.inaturalist.org/observations/150764924) | "animal remains" | **2015, same original owner and same day as #150764928** | 1 photo ID, **pixels not viewed** |
| [#150764928](https://www.inaturalist.org/observations/150764928) | "animal remains" | **2015, same original owner and same day as #150764924** | 2 photo IDs, **pixels not viewed** |
| [#367297297](https://www.inaturalist.org/observations/367297297) | "carcass of amphibian" | **2026, another owner-day** | 1 photo ID, **pixels not viewed** |

**Proven:** three text annotations, two distinct owner-day clusters, two source observers, no identical photo ID among those three. **Not proven:** whether the first two are the same carcass, whether any photographed animal remains are genuine carrion, whether a visible adult butterfly extended its proboscis, whether any nutrient was swallowed, whether there is fish carcass use in this source, and whether any site had wastewater input. Two owner-days cannot be labelled "two independent feeding events."

### Within-owner family annotation comparison — do not promote to biological choice

Exactly **18 of the 181** original annotating owners contributed both Lycaenidae and Nymphalidae records, **8** of whom also had both an organic-labelled and a mineral-labelled record somewhere in their original source annotations. Across the 18, the *unweighted mean owner-specific difference* in the **share of voluntarily annotated organic candidates** (Lycaenidae minus Nymphalidae) was **−0.09789**. Leave-one-owner-out range was **−0.10365 to −0.07424**: **8** owners negative, **10** zero, **0** positive in their own annotated source. These are *only* conditional labels within a highly selected, overlapping observation-source sample; they are NOT measured comparative ecological substrate preferences. Observers may have sampled different habitats, sites, seasons, species and substrate availability for each family, even under the same account. No independent non-puddling denominator, source chemistry, sex, matched local substrate offers, or controlled feeding outcome is present.

The data do **not** independently corroborate a generalized special Lycaenidae propensity for organic/carrion intake, nor do they refute experimental protein-bait choices. A selected observational label distribution is not a standardized preference experiment. **Do not calculate a confirmatory p value or wastewater-mediated path from these data.**

### Decisions for the research record

- Stage HNS (larval host trait): unchanged 2022 riparian one-shot, still STOP due source/geometry/observer opportunity.
- Stage B1 (adult substrate labels): original full 746 / 78 organic review-ready records, photo-ID and observer-day dependence audited. **No human-reviewed proboscis contact** yet.
- Stage B2 (wastewater-specific adult resource): still lacking measured site-specific organic/Na solutes, exact source→receiver exposure, habitat opportunity and balanced butterfly behavior observation.
- Stage trophic pathway (fish carcass): original metadata in this selected field have no explicit fish-carcass annotation and only three broader animal-remains texts in two owner-days; **does not rule out fish-carcass feeding documented elsewhere**. Original local fish death and receiving-water chemistry must be independently distinguished before an ecosystem-service link.

All revised original metadata stored on existing Draft PR #95. No images were redistributed, no new analysis PR or merge, and no purported proof of sewage-dependent butterflies.


## Original CC BY-NC photograph scene overlap and same-site check

**2026-11-15, further genuine source audit, not visual feeding confirmation.** The [three original iNaturalist records' photo metadata](https://github.com/zuizui0223/gesui/actions/runs/43707346799) were re-fetched from the original iNaturalist API and exact photo-ID/source-field text reverified. All four photos (original IDs 260028346, 260028363, 260029204, 670748495) are **CC BY-NC** under the author's original photo licence; no photo was assumed freely redistributable. Three original observations correspond to two distinct owner/date units, including the 2015-08-13 original records #150764924 and #150764928 from the same observer, and the 2026-06-01 record #367297297 labelled "carcass of amphibian".

To check **photo-scene overlap, not the organism or its feeding behavior**, the [noncommercial transient-source image similarity script](../analysis/20261115_adult_resource_behavior/audit_three_carrion_photo_similarity.py) fetched only the four source CC BY-NC **medium-resolution photos into the GitHub Actions runner's process memory**. The input images were discarded; the code did not commit original photos, upload copyrighted image bytes, publish their pixels, or perform manual biological image classification. The source-only numeric workflow [#43707869730](https://github.com/zuizui0223/gesui/actions/runs/43707869730) succeeded; its SHA-verified derived [JSON feature result](../results/20261115_adult_resource_behavior/INAT_3_CARRION_ORIGINAL_TRANSIENT_PHOTO_FEATURE_QA.json) is archived separately.

**New source findings:**
- Original public **2015 observation coordinates coincide exactly** in the original API (computed public-coordinate Haversine difference **0.0000 km**). This is a same *published coordinate*, not proof the butterflies or their resource occupied literally the identical point (coordinate georeferencing/rounding can collapse nearby encounters).
- In the original photo #260028346 (iNat observation #150764924) and **secondary** original photo #260029204 (observation #150764928), SIFT-free ORB descriptor matching found **77 Lowe-ratio0.7 candidate feature pairs**, of which **75 are geometrically consistent within a 5-pixel RANSAC homography (97.4% inliers)**. That is strong feature-based support for **shared visual scene content**, potentially different frames/crops/angles. No source photo was visually presented for proboscis inspection.
- Original photo #260028346 versus **primary** #260028363 produced only **one** candidate ORB match (no homography); the two photos of #150764928 also had zero descriptor matches in that algorithm. Nonmatches do **not** prove separate ecological encounters; the camera angle, zoom and subject can differ.
- Original **three text-based animal remains annotations** are therefore best treated as **two independent observer-day clusters**, with two 2015 records likely sharing a scene. The number of independent **carcasses**, **butterfly individuals**, and **actual feeding acts** is still unknown.
- The 2026 original source explicitly says **amphibian carcass**, not fish. Neither the 2015 labels nor the 2026 original source establish actual fish-carcass puddling, source effluent, or nutrient uptake. A registered WWTP upstream cannot be inferred from a photo or an original user-entered field.

**Interpretation discipline:** a homography confirms geometric alignment in image features, **not** visible proboscis contact or ingestion, and cannot identify a specific carcass without actual visual substrate review. The photo source has not undergone independent human substrate/proboscis scoring. Report the two 2015 **source observations** (with an evidence-supported scene link) separately from the ecological denominator of independently sampled resources. Retain CC BY-NC attribution and direct original observation links, not any copied image material, in all public artifacts.

**Original links:** [#150764924](https://www.inaturalist.org/observations/150764924), [#150764928](https://www.inaturalist.org/observations/150764928), and [#367297297](https://www.inaturalist.org/observations/367297297). All three remain **feeding-unverified**. The 2022 fixed same-H7 outfall→HNS test is unchanged, as are larval-HNS and negative-control conclusions; this is an independent adult-resource evidence-quality audit only.


## Targeted dung 46 source audit — full source record/photograph independence, no feeding confirmation

**2026-11-15 next adult-resource test; no change to locked 2022 H7 source-specific HNS.** Re-fetched **all 746** original iNaturalist field500 annotations (the same source population) through the sparse iNaturalist v2 original API. The frozen revised substrate mapping yields **46 source records labelled unambiguously dung/faeces**, plus **3 mixed dung+mud labels** deliberately excluded from this strict 46. An exact source-only unit audit counted **46 observations / 33 public observation owners / 44 owner × observation-date units**. GitHub Actions [#43711396724 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43711396724), with byte-verified derived original numerical output [#43711586394 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43711586394).

### How much actual independent biological evidence does it support?

| Original iNat field500 strict dung metadata | Observations | Original owners | Owner × day units | Distinct original source photo IDs |
|---|---:|---:|---:|---:|
| Lycaenidae | **9** | **7** | **8** | **14** |
| Nymphalidae | **24** | **20** | **24** | **62** |
| Papilionidae | **6** | **5** | **5** | **26** |
| Hesperiidae | **5** | **4** | **5** | **9** |
| Pieridae | **0** | **0** | **0** | **0** |
| Other/unclassified | **2** | not individually listed | not individually listed | not individually listed |
| **Total** | **46** | **33** | **44** | **113** |

Total per-family owners and days cannot simply be summed as independent entities without checking whether an owner annotated different butterfly families. Notably **only TWO owners** tagged dung observations in both Lycaenidae and Nymphalidae, reducing the potential for a source-matched family comparison. The inverse-Simpson effective original observation-owner number for the 46 tags is **23.0** (concentration diagnostic only, not ecological sample size), and the top five owners account for **36.96%** of the dung-tagged source frame.

Across **all 113 original photo IDs** attached to the 46 dung-tagged observations, **no identical photo ID was reused across distinct observations**. This does **not** mean all photographs represent different individuals or dung piles. There were **two observer-day clusters with two separate observations each** (4 observations within duplicate owner-day groups). Of the original 46 records, only **16** had a description, and **6** had a verbatim description containing a feeding-related keyword: neither source field description nor keyword establishes **proboscis contact**.

Photo-specific original source licence metadata were **19 CC BY, 62 CC BY-NC, and 32 without an explicit re-use licence** across 113 photo IDs. Do not retrieve/publish unlicensed originals; a published observation does not mean the photograph is licensed for reuse. The original [46-row observation review queue](../results/20261115_adult_resource_behavior/INAT_DUNG46_POTENTIAL_PHOTO_REVIEW_QUEUE.csv) contains public source observation links, source family, original substrate text and photo licence checks, with every human proboscis and substrate verification field explicitly **NO/NOT_REVIEWED**. The original [numerical independence ledger](../results/20261115_adult_resource_behavior/INAT_DUNG46_ORIGINAL_INDEPENDENCE_QA.json) retains total source, counts, licence classes and source HTTP receipt hashes, without raw photo pixel redistribution or owner identifiers.

### Original two same-owner/day cow-manure clusters show why owner-day alone is not an individual feeding unit

The two original Lycaenidae source observations [#16195221](https://www.inaturalist.org/observations/16195221) and [#16195223](https://www.inaturalist.org/observations/16195223) both have original field500 value "cattle manure" and the same photographer/observation date, but **their public source coordinates are 6.12855km apart**. Their two original photograph licences are missing, so **no image pixels from either were fetched**. These observations may represent distinct site-level encounters, though a reported coordinate can be generalized and butterfly intake remains unverified.

The original Papilionidae source observations [#18670197](https://www.inaturalist.org/observations/18670197) and [#18670234](https://www.inaturalist.org/observations/18670234) likewise share observer/day and "cattle manure" label, but have **exactly coinciding public coordinates**. Across their **7+14=21** source photograph assets all report CC BY-NC, so a noncommercial transient in-memory photo feature comparison was allowed. Original different-observation photos had **a maximum 63 ratio-0.7 ORB feature matches and 63 geometrically consistent 5px homography RANSAC inliers**, with **16 cross-observation photo pairs** having ≥8 inliers. Strong source evidence of **shared photographic scene content**. This supports not treating the two observation records as independent substrate-use events by default, but not the biological identity of a specific manure pile or butterfly, and NOT feeding itself.

[Original source-only 2018 family-pair geometry and licensed transient photo feature run #43711773520 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43711773520). [Original derived SHA-verified pair summary](../results/20261115_adult_resource_behavior/INAT_2018_LYCA_PAPILIO_COWMANURE_PAIR_SOURCE_QA.json) was archived with neither original image pixels nor restricted photographer metadata. The combined 46-case owner-day and source-coordinate audit demonstrates that **same-day repeated tags can refer to distant encounters or shared scenes**. Simply dividing all raw 46 by owner-days, unique photo IDs or unique field values would be a biased proxy for independent feeding opportunities.

### Scientific conclusion for hypothesis B, separate from hypothesis A

- **Supported:** voluntary iNat annotations explicitly place several butterfly families on faecal/mud/organic substrates; original metadata can be independently deduplicated, and the strict dung-tagged source frame is more geographically and observer clustered than the record count suggests. It includes original 9 Lycaenidae and 24 Nymphalidae dung-labelled observations, but only 2 observers ever annotated dung for both.
- **Not demonstrated:** photographed adult proboscis contacting dung, drinking or nutrient uptake, field-independent substrate choice, family-specific dung preference, sewage-N or Na resource subsidy, larger butterfly biomass, or ecosystem services. **0 verified photo-feeding observations** in this code workflow means **none was visually reviewed**, not a proven absence of real feeding.
- **Selection:** original field500 is a user-generated, user-filled annotation; the vast majority of ordinary butterfly sightings have no field500, and they cannot be treated as nonfeeding controls. Even a robust comparison of tagged-substrate composition cannot estimate family relative attraction without knowing available competing substrates and standardized observations.
- **Next admissible measurement:** licence-aware *human* photo/field observation verification of substrate plainly visible, adult proboscis **tip** visibly contacting substrate, possible different individuals and same manure patch, sex when possible and independent observers. A confirmed photograph can establish **contact** at best, not sodium/nitrogen assimilation. Separate matched hydrologically valid source versus unaffected wet-reference conditions with source and water chemistry are still needed to test wastewater dependence.

The original 2022 H7 Lycaenidae larval HNS study remains separately frozen and STOP-gated. Existing PR #95 is Draft; no new PR, no main/chocho merge, no eBird or TinyFish.
