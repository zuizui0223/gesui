# 2026-11-10 Denmark eBMS: first original standardized Lycaenidae individual-count test near reported WWTP outfalls

**Scope:** Source-backed and independently replayed **2023–2024 standardized survey counts**, NOT iNaturalist photographic reporting. The test was conceived after earlier DE/AT positive WWTP–photo exploratory associations: **POSTHOC** ecological hypothesis, not a new independent causal confirmation. The 2022 EEA reported N/P outfall location and annual reported mass are static source proxies, not documented real 2023–24 water chemistry or butterfly exposures. All findings stay on the separate Draft branch.

## Independent public source and frozen gate sequence

1. Original licensed Denmark eBMS v1.4 DwC-A source ZIP, SHA256 **2cc64ec352754223db06a92eef5b01803968eac94bae34057c17f1abdb6da79b**, 15,606 Event core rows and 21,617 Occurrence extension rows. Survey event dates 2023 (3,100) and 2024 (4,461), official CC-BY-NC 4.0. Original EEA DK2022 original one-FW-outfall both-NP-reported measured source **203 outlets** across Danish country bounding region; previous "69" belonged only to an old Sweden-oriented geographic crop, an audited pre-butterfly correction.
2. First **source-only** stage: exact original source, same 5-decimal-coordinate site identity, <=500m reported geo uncertainty, >=2 distinct visit dates in each 2023 and 2024. **107 repeated original public survey-site coordinate rows**. 2022 WWTP source distance => 23 public site rows <=5km, 38 >=10km, 3 original 1° source-and-remote shared cells. At 2km: **0** near; no 2km biological effect calculation. The original 5km source-availability gate passed before reading any species outcomes. Corrected 5-decimal raw-coordinate exact equality implementation bug and history are retained.
3. Independent spatial pseudo-replication source-only audit: **23 near site rows represent 4 distinct actual WWTP outfalls** (site cluster sizes 7,7,5,4); 38 far site rows form 8 components at 0.5km connected distance. **Effective independent near source count is four, NOT 23.**
4. Original source-verified GBIF backbone scientificName taxonomy frozen before original individualCounts were opened: 112 distinct scientificName labels, **21 confidently matched Lycaenidae species, 7 Pieridae, 1 Papilionidae**, secondary Nymphalidae 28 and Hesperiidae 8. Other unassigned names include dragonflies and moths not part of the three-family primary endpoint, and ambiguous genus; none was opportunistically assigned to the desired family. Original original source-name archive and licensed attribution retained.
5. Original eBMS DwC-A Event and Occurrence linked with **100% successful original coreid join across 21,617 occurrence rows**; same-source original 107 fixed sites × 2023/24 = **214 site-year count panel**, 4,731 eligible original occurrence rows, 2 missing/invalid/zero counts excluded, 4,226 confident accepted butterfly-family occurrence rows. Original published occurrenceStatus was "present" for all 4,731 site-year source rows. Independent count replay SUCCESS. No observer identities, personal recording IDs or original licensed raw records republished.

## Observed family compositions: all original source-qualified repeated sites

The original survey counts (individualCount), NOT the number of photographs:

| Standardized survey year | 5km near registered 2022 WWTP: Lycaenidae / (Lycaenidae + Pieridae + Papilionidae) | >=10km remote | Near and remote original distinct survey days |
|---|---:|---:|---:|
| **2023** | **74/(74+128+0) = 36.63%** (11/23 sites with any focal 3 families) | **387/(387+355+0) = 52.16%** (28/38 sites) | 102 / 378 |
| **2024** | **34/(34+151+0) = 18.38%** (18/23 sites) | **494/(494+481+0) = 50.67%** (31/38 sites) | 219 / 617 |

Original **Lycaenidae counted individuals per documented survey day** are 2023 near **0.7255** vs remote **1.0238**, and 2024 near **0.1553** vs remote **0.8006**, also lower for near in these pooled observational site sets. This does not establish a population or within-source effect: the 23 near site coordinate rows represent four source clusters, between-site and source distances are strongly grouped, and survey dates/effort differ.

**Crucial:** papilionid species were present in the frozen taxonomic dictionary but **no counted Papilionidae individuals were recorded on the original 107×2 source cohort**. The primary three-family proportion thus reduces to Lycaenidae versus Pieridae in these sample years; no evidence is claimed about Papilionidae response.

## Geographic-matched 1° comparison—where the original sampling supports estimability

Original 5km near and >=10km far public repeated visit sites share exactly **three** 1° source cells: **55|10, 55|12, 56|9**. The primary descriptive endpoint is the equal-site mean log odds [log((L+0.5)/(Pieridae+Papilionidae+0.5))] among sites with >=1 count of the three focal families, near minus remote. Original site-year with zero focal 3-family counts is **undefined for composition**, not arbitrary zero-logit evidence.

| Original 1° block | 2023 near-vs-far site-equal log-odds | 2024 near-vs-far site-equal log-odds |
|---|---:|---:|
| 55|10 | not estimable (near: 0 Lycaenidae, 0 Pieridae, 0 Papilionidae) | **+0.622** (but near Lycaenidae still 0; continuity correction strongly matters) |
| 55|12 | not estimable (near: 0 of focal 3 families) | **−0.474** |
| 56|9 | **−1.038** | **−1.658** |

- **2023:** only 1 of the 3 shared source cells had focal butterfly family counts on **both** sides, so **no 3-region or two-year confirmatory result can be estimated**.
- **2024:** all 3 shared cells had observed target butterflies on both sides. Site-equal three-region mean contrast **−0.503** log-odds; 2 of 3 regional point estimates negative and one positive. The leave-one-region-out mean range is **[−1.066,+0.074]**, crossing zero. The pooled-count regional log-odds version is **−0.817** equal-region average, also not a causal test. The apparent positive 55|10 site-equal coefficient is unstable to denominator continuity correction and occurs with **0 near Lycaenidae individuals**.
- There are **only 3 shared geographic cells and 4 registered source-outfall clusters (one in the fourth near region with no far control)**. Repeating site records in 2023/2024 does not create more independent WWTP interventions. No meaningful confirmatory p-values or high-confidence causal effects are claimed.

## Scientific interpretation and decisions

**New observational result:** The previously positive source-P geography association with **opportunistic Germany/Austria iNaturalist photographic reporting composition** does not directly generalize as a positive effect in a different country's **standardized 2023/2024 family individualCounts**. In this Danish sample, the broad pooled contrasts are **negative in both years**; the geographically paired regional comparison supports **no stable positive signal**, and 2023 is inestimable in two of three matched regions.

This does **NOT** show WWTPs harm Lycaenidae. 2022 WWTP reports indicate facilities, not actual source-to-plot hydraulic connectivity, aquatic chemistry, 2023/24 effluent changes, sodium/protein puddling, larval host quality or ant mutualism. A 5km outfall radius may lie wholly outside a relevant river corridor. Source statuses are not randomized; land use, vegetation, phenology, participant geography, repeated survey sections and within-site detection confounding remain. Original remote >=10km means no nearby **qualified registered WWTP** in the strict source file, not no sewage. Thus **stop claiming source-associated Lycaenidae ecological benefit** until independently measured exposure/habitat and genuine source-level replication exist.

**Recommendation for future macroecology:** do not rescue the result by choosing new radii or a favorable year. Use this as a genuinely different, source-audited standardized observational validation panel, with clear negative/heterogeneous directions. Direct test requires larger independently exposed WWTP units, physical river direction and verified substrate/canopy/host ecology and genuine repeated butterfly monitoring.

## Reproducibility and source precedence

- [Original source-only Event × EEA WWTP source gating SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42757540447)
- [Original source outfall/near-site independence QA SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42757948325)
- [Frozen blind GBIF scientificName family taxonomy SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42758342905)
- [Original licensed family count Event-Occurrence joins and quality SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42758924673)
- [Source-near vs remote three-cell annual contrasts SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42759299290)
- [Independent arithmetic/counts/rate replay SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42759571075)
- [Frozen site-year standardized counts and full source QA](https://github.com/zuizui0223/gesui/blob/analysis/ebms-denmark-wwtp-source-feasibility-20261110/results/ebms_denmark_2023_2024_WWTP_source_only/DENMARK_EBMS_2023_2024_FAMILY_SITEYEAR_COUNTS_QC_BEFORE_WWTP_CONTRAST.json)
- [Full 2023/24 region-year positive/negative and nonestimable results JSON](https://github.com/zuizui0223/gesui/blob/analysis/ebms-denmark-wwtp-source-feasibility-20261110/results/ebms_denmark_2023_2024_WWTP_source_only/DENMARK_EBMS_2023_2024_LYCAENIDAE_SOURCE5KM_VS_REMOTE10KM_POSTHOC_REGION_RESULTS.json)
- [Independent count replay](https://github.com/zuizui0223/gesui/blob/analysis/ebms-denmark-wwtp-source-feasibility-20261110/results/ebms_denmark_2023_2024_WWTP_source_only/DENMARK_LYCAENIDAE_2023_2024_INDEPENDENT_COUNTS_SOURCE5KM_QA.json)

**Reporting rules:** Preserve raw licensed source attribution/metadata; no commercial redistribution of original CC-BY-NC raw occurrence rows here; no original personal observer IDs stored. This result and all scientific caveats remain on an isolated Draft analysis branch; do not modify the independent chocho introduction-host manuscript or gesui main, and do not overwrite pre-existing positive source/photo evidence from Germany/Austria.
