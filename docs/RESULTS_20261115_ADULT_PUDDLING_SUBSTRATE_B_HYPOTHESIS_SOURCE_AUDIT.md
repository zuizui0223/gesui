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
