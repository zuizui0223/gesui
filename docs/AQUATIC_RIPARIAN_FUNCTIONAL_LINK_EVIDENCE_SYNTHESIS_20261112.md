# Aquatic–riparian evidence bridge: same-site river nutrients, EPT and vegetation; separating layers from causal chains

**Checked 2026-11-12 | exploratory, noncausal | stacked unmerged analysis**

## The decisive rule

Do **not** multiply the separately estimated associations from European WWTP-EQR, USA river TP-EPT, USA wetland plant diversity and island/port butterfly photographic records and report them as a connected wastewater→insect→plant→pollinator path. These observations have different geography, outcome types, dates, effort processes and sampling units. A multigroup SEM on those separate summaries would estimate an ecological mechanism that has **never been observed on the same units**.

The scientifically coherent question instead is:

> **When is freshwater biological integrity coupled to riparian habitat structure, when is aquatic-to-terrestrial biological resource supply maintained despite taxonomic replacement, and where is the coupling broken by nutrient and human disturbance gradients?**

"Water quality", "aquatic biodiversity", "potential aquatic insect subsidy", "plant host/flower resources", "butterfly guild", and "measured ecosystem service" must remain different outcomes.

## Recomputed same-visit bridge: US EPA NRSA2018–19, source locked

**Original sources**: 2018–19 EPA NARS native UID metadata + original **BERW benthic counts**, original measured river total N/TN and total P/TP, and **PHab field-measured** `XGveg` ground vegetation, `XPCAN` canopy reach presence, `XFC_OHV` stream overhanging cover, `W1_HALL` human disturbance. Same visit, no GBIF record joins. Pre-outcome UID manifests SHA256 `fbf53b...` (1,221 original) and `23cf9546...` (riparian 1,220). Original source file SHA pinned, code written before reading new riparian magnitude.

### Prespecified main analysis

- n **1,184 native UID visits / 1,091 physical sites / 48 states / 9 eco9 / 722 HUC8**, after all same-visit background variables and original assumed 0–1 limit for `XGVEG` (35 original values above one were excluded). Five **HUC8-disjoint** spatial holdout folds, Ridge alpha 10, identical cases, train-only standardization; adjusted for state+eco, lat/lon, original stream order, conductivity, year and field-observed `W1_HALL`.
- Model B = those controls; N = B+log10(TN)+log10(TP); H = B+XGVEG+XPCAN+XFC_OHV; NH = both; NHI = NH + single predeclared TP×(XGVEG−0.5) interaction.

| Prespecified heldout contrast, site-equal MSE | Improvement | HUC8 cluster resample descriptive 95% |
|---|---:|---|
| **B→N** nutrients add predictive information | **+2.354%** | [+0.568, +4.091]% |
| **B→H** riparian structure added alone | **+1.535%** | [−0.317, +3.281]% |
| **N→NH** riparian structure after chemistry | **+1.284%** | [−0.453, +2.919]% |
| **H→NH** nutrient chemistry after riparian structure | **+2.105%** | [+0.393, +3.793]% |
| **NH→NHI** extra TP×vegetation term | **+0.315%** | [−0.504, +1.165]% |

Adjusted exploratory OLS on same original n1,184 NH model, HUC8 cluster descriptive 95%:
- river **TP coefficient −0.04659** [−0.08172,−0.01147] in EPT *individual fraction per tenfold TP concentration*;
- river TN −0.04196 [−0.08501,+0.00109], interval crosses zero;
- ground vegetation XGVEG coefficient −0.11544 [−0.18130,−0.04957];
- canopy presence XPCAN +0.05284 [+0.00853,+0.09715]; overhanging vegetation XFC_OHV −0.01022 [−0.10999,+0.08954].
- NHI TP × (XGVEG−0.5) interaction **+0.13393** [0.04806,0.21981] by descriptive clustered OLS, **but** adds only +0.315% heldout fit, uncertainty spans zero. No demonstrated vegetation causal moderation, ecological restoration or emergence subsidies.

### Critical source-value-range audit (NEW, **strictly posthoc**)

The original EPA physical habitat `XGveg` is an additive metric, not a simple [0,1] mutually exclusive percentage: in the authentic source file, `XGVEG = XGH + XGW` to machine precision for the **35 visits with XGVEG>1** (max 1.22727). The EPA source metadata does **not** prescribe a max of 1. The main model's prespecified upper-bound rule therefore excluded 35 source-defined observations without good ecological grounds. It is **retained transparently** as the primary original specification, not overwritten or advertised as a strict independent confirmation.

The independent Python original four-source analysis re-implemented OLS slopes and heldout fits and verified the original n1,184 interaction *exactly*. New posthoc tests preserving **all authentic original XGVEG values**, n **1,219** visits / **1,123** physical sites:
- riparian addition after N: **+1.482%** HUC8 heldout gain;
- nutrient + habitat OLS TP slope **−0.04915** (descriptive 95 [−0.08365,−0.01465]);
- TP×ground interaction **+0.13361** [0.04919,0.21804], but **incremental OOS gain just +0.308%**.
- Capping original XGVEG >1 at 1 gave +1.493% and interaction gain +0.286% (posthoc sensitivity).
- With *both* clean original nutrient NARS flags and ≥250 counted benthic individuals, n **795**, vegetation-after-N heldout gain **+2.343%**, but **interaction gain reverses to −0.210%** and descriptive 95% for interaction [−0.01025,+0.25214] spans zero. TP slope in NH model stays negative **−0.07242** [−0.11690,−0.02794].

**Final evidence rating:** same-site environmental co-variation is demonstrated. The *nutrient TP–EPT negative direction is more robust* than the independent value of riparian predictors or the vegetation×nutrient interaction. The positive interaction can guide a future functional buffering test, but is currently a fragile posthoc tendency rather than a proven service-mediated buffering process.

[Original same-visit EPT-water-riparian result SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43050507703) · [Independent original-source ground-cover range and interaction sensitivity SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43050963972).

## Evidence-link matrix (do not misstate join type)

| Proposed edge | Actual source & grain | Evidence strength | Can claim causal cascade? |
|---|---|---|---|
| **EU measured WWTP N/P → river benthic EQR** | EEA European 2018 n2,634 and 2022 n2,475; upstream reported WWTP source, downstream controls; independent original years but 1,167 shared sites | small reproducible **predictive** gains P+0.98/+1.10%, N+1.19/+1.40% | **No**: lacks delivered dose, agriculture/urban controls, source-biased EQR |
| **US actual river TN/TP ↔ benthic EPT count fraction** | 2018–19 NRSA n1,220, many states, identical native UID | strong reproducible **observational negative TP association**, out-of-HUC8 prediction | **No**: cannot apportion TP to WWTP vs agricultural sources |
| **US riparian vegetation ↔ benthic EPT/nutrient co-organization** | same n1,184 valid originally / posthoc n1,219 as above | some independent spatial predictive signal, wider block CI includes zero; interaction OLS positive but OOS weak | **No**: field structural cover, simultaneous habitat and chemistry |
| **Benthic insects → winged insect emergence → predators** | no emergence flux nor diet on identical NRSA sites | literature only; no `gesui` matched endpoint | **No** |
| **Port wastewater/WWTP ↔ plant/host catalogues** | 1,055 global ports / 42 countries; opportunistic GBIF records | raw records/species positive but rarefied/Animalia/Fungi adjusted weak | **No** |
| **Candidate host plants → butterflies** | `chocho` island+port resource geography; `gesui` DE/AT photo and Danish eBMS | candidate resource availability positive; actual local Lycaenidae positives fail within-region and effort controls | **No** |
| **Actual butterflies → pollination, crop yield, plant seed set** | no standardized same-site interactions or plant reproductive success | unmeasured direct service | **No** |

**Critical same-site butterfly source gate:** correcting invalid iNaturalist ID `47923` to GBIF Lycaenidae `5473` produces 324/360 original requests successful, but only **9** fixed 5 km 2017–2021 river sites support Plantae≥50 **AND** Lycaenidae≥20; preregistered ≥100-site gate FAILS. No linked butterfly regression was run. [GBIF ID and source correction SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43048245515).

## Recommended manuscript architecture

**Headline**: *Aquatic biological integrity and terrestrial functional potential can be decoupled across human-modified watersheds* — this must be described as a hypothesis / working title, not a demonstrated across-ecosystem service finding.

**Main empirical result**: the US true river-nutrient × EPT association and same-visit riparian habitat information, with source QA and HUC8 holdout. Europe adds independent reported source–EQR *predictive context*, not a pooled mediation coefficient.

**Mechanism figure**: two incoming environmental branches (WWTP/other nutrient pressure→river nutrients; terrestrial habitat→river) and two *potential*, dashed outgoing routes (aquatic insect emergence→predators and riparian plants→butterflies). Mark observed EPT/field habitat co-organization, and show distinct endpoint types. Do **not** assert arrows to pollination or pest control.

**Strongest next data requirement**: measured aquatic **emergence biomass / taxon composition** and flower/host plant or insect consumers from the same watershed/period. For a bridge with real flux endpoints, public [2022 H.J. Andrews LTER adult aquatic insect emergence observations](https://andrewsforest.oregonstate.edu/data/datacatalog/AS015) cover 17 streams with measured emerging insect density/biomass, and public [James River estuary emergence data](https://doi.org/10.5061/dryad.41ns1rnk1) include biomass per m² per day. Both are *separate source programmes and habitats*; neither should be artificially joined to NRSA national sites or used to claim national transfer without a new source-first measurement study. No eBird.

**No model deletion, no main merge, no automatic changes to `chocho`.** All results remain Draft and exploratory.
