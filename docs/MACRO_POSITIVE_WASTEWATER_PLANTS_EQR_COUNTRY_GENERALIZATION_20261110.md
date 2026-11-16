# Macroscale wastewater ecology: verified positive associations, independent-country limits and research priority

**2026-11-10 | draft research synthesis | no causal wastewater-to-biodiversity effects claimed**

## Central macro ecological question

Do modeled wastewater production and officially reported wastewater nutrient loads covary with plant-community recording composition and aquatic biological condition *across countries* after controlling coarse geography, urban/rural habitat and observation/reporting opportunity? These are TWO different exposure measurements and TWO different biological response datasets, not one common N/P mechanism.

Stop breaking the question into Lycaenidae taxa and a handful of nearby outlets before establishing large-scale positive co-geography, independent replication and scope. Geographic scale and response type must be reported beside every coefficient.

## A. The leading positive geographic association: 2015 modeled wastewater production and port indexed plant records

**Original same-sample fully reproduced:** 1,596 original WPI ports across 142 countries with original HydroWASTE 58,502 WWTPs, independent PANGAEA Jones2015 5-arcminute modeled wastewater production WWp (million m³/year) and ESA WorldCover2021 built/crop pixels. Original 1,055 ports across 42 countries passed same country-fixed-effect/WorldCover/plant-indexed source QC. Exact original regression:

`log((alien candidate butterfly host plant species recorded+.5)/(native candidate host plant species recorded+.5)) ~ log1p(modeled WWp_2015 within 25km) + log1p(registered WWTP count within 25km) + WorldCover built/crop + absolute latitude + periodic longitude + log1p(total GBIF-indexed candidate hosts) + country fixed effects`.

The original positive regression beta **+0.0586977993**, 42-country-cluster 95% **[+0.0131035655,+0.1042920331]**, nominal p **0.01163**, R² 0.647753 on full controls. This reproduces the archived original positive **bit-level numerical estimate within floating point tolerance**, not a new ecological mechanism or all-flora result. 2015 WWp is modeled **WASTEWATER GENERATED**, not measured river effluent, stream/soil N/P load, or direct vegetation water supply. Original plants are candidate butterfly larval hosts, **not ALL Plantae** nor observed introductions at each port.

### NEW large-scale heldout tests on exactly the original 1,055-port/42-country cohort

| Heldout design (all chosen before new evaluation) | Sample | Relative MSE change after adding WWp | Supported interpretation |
|---|---:|---:|---|
| Original country FE, no heldout (historic) | 1,055 ports; 42 countries | WWp beta +0.0587, CI fully positive | Original positive source-geographic association |
| Delete each country and re-fit SAME country-FE model | 42 independent country deletions | β positive **42/42**, min **+0.05062**, max **+0.08061** | Positive coefficient not driven by any ONE country; does not test new countries |
| Same COUNTRY, different geographic 2° blocks in 5 folds | **996 ports, 37 countries, 446 source 2° blocks** | **+1.277% country-equal** / +0.542% port-pooled | Modest cross-region forecast improvement within already represented countries, 23/37 country strata improve |
| Entirely unseen country (42 country holdouts), NEVER country FE | **1,055 ports, 42 countries** | **−1.225% country-equal** / −0.619% port-pooled | Adding WWp does NOT improve across-country extrapolation; country-cluster fixed-sample bootstrap 95% **[−3.641%,+1.196%]** |

Within 30 source-supported countries with ≥8 ports and ≥3 WWp/plant-response unique ranks, 20 have positive original country-specific Spearman and 10 negative. Thus pooled within-country beta positive does **not** imply a sign-uniform country law.

The separate original candidate-recording controls are both positive: WWp β on log(1+alien candidate indexed spp) **+0.22862** [0.15272,0.30451] and WWp β on log(1+native candidate indexed spp) **+0.17539** [0.09137,0.25941], after geography/WWTP/landcover and country FE. That demonstrates WWp correlates with recorded richness of **both** components; human settlement, port activity and botanical indexing opportunity are plausible confounders. The relative alien/native positivity does not independently establish plant colonization or nutrient enrichment.

**Primary finding A:** a real positive global pooled observational association exists and does not depend on one country. It is not yet a transferable universal relation across *new* countries, and is not a causal wastewater-driven plant enrichment effect. Its positive result should be retained, not hidden by later local failures.

## B. Independent positive macroscale wastewater load information in river benthic invertebrates (not a positive biological abundance slope)

Official EEA UWWTD **reported measured N/P annual source mass**, reported to official FW outfalls and joined to modeled HydroRIVERS directed source, predicts measured official EEA benthic-invertebrate EQR modestly beyond original geography/method/river baseline, and this positive **prediction increment** was replicated in 2018 and 2022:

| Source report year | Measured benthic EQR stations | Country/region support | Added upstream measured P mass MSE gain | Added N mass MSE gain |
|---|---:|---|---:|---:|
| **2018** | **2,634** | 15 countries / 52 regional blocks, 202 source-positive | **+0.982%** | **+1.190%** |
| **2022** | **2,475** | 13 countries / 45 regional blocks, 147 source-positive | **+1.103%** | **+1.397%** |

Prior 2018 and 2022 EQR×registered source mass rank correlations were **negative** (2018 P −0.131, N −0.132; 2022 approx −0.072). The positive is **forecast-information gain**, NOT increased invertebrate abundance, richness, ecological health, or a beneficial nutrient effect. Benthic EQR index is based on multiple country/method reference schemes. Strong source/reporting missingness: 2022 807 Polish benthic stations lack numeric EQR; missingness is WWTP-source-associated. Phytobenthos shows little or inconsistent generalization. N and P mass are strongly correlated, so do not infer separate P vs N chemistry causality.

Reproducibility: [2018 original source and independent EQR time replication PR #17](https://github.com/zuizui0223/gesui/pull/17), [2022 original EQR PR #14](https://github.com/zuizui0223/gesui/pull/14), and [2022 source-bias audit PR #16](https://github.com/zuizui0223/gesui/pull/16). Their original executed analyses and separate independent fixed-fold checks were successfully archived.

## C. Direct water-chemistry positive associations are local, exploratory and not yet macro generalization

- Six **independent German WWTP-source rivers**: 2,500 genuine TP monitoring samples, 836 paired months, modeled wastewater water fraction vs downstream/upstream TP concentration-gradient Spearman **ρ +0.8286**, n=6; ordinary asymptotic p=0.0416, but **EXACT two-sided 6! rank-permutation p=42/720=0.0583** (not below 5%). The TP-to-plant species-record gradient is null (~−0.029). [PR #7](https://github.com/zuizui0223/gesui/pull/7).
- Germany 2022 EEA reported-source-interposed vs registered-sourcezero real WISE6 four-arm same-month TP gradient **+0.1686 log mg(P)/L** with descriptive 95% [−0.0770,+0.4428], 7 independent quartets. A different 2023 source-constrained comparison gave +0.2252 in 5 quartets but exact small-n reference non-significant. Exact site×calendar month four-quartet 2022 +0.2962 and 2023 +0.2398, no year-change evidence. [PR #33](https://github.com/zuizui0223/gesui/pull/33), [PR #34](https://github.com/zuizui0223/gesui/pull/34). These are limited direct physical endpoints and **not** evidence of globally increased butterfly or plant populations.

## D. Butterfly positive photo-share belongs to a separate noncausal and denominator-unstable panel

2022 DE/AT original EEA P source vs Lycaenidae fraction of photographed Lycaenidae/Pieridae/Papilionidae counts was positive, country/WorldCover adjusted posthoc 2° block heldout MSE gain **+13.67%** on 43 cells, descriptive interval [−14.18,+33.04]. But 2023 photo history, geography-matched 2° block shuffles, 2025 Denmark independent standard transects, 3-family vs 5-family denominator, and actual 2021 land cover do not give stable source-positive Lycaenidae recorded abundance. [PR #25](https://github.com/zuizui0223/gesui/pull/25), [PR #40](https://github.com/zuizui0223/gesui/pull/40), [PR #41](https://github.com/zuizui0223/gesui/pull/41). Keeping only favorable family denominators would be outcome selection. **Do not make this the macro principal endpoint.**

## Research line and immediate next measurement

**Main hypothesis:** more human wastewater generation/exposure is associated with shifts in the **relative composition of flora and water-linked biological assemblages**, with effects (if any) gated by river routing, landcover, dilution, recording and regional history. Currently this is a question, not a proven WWTP subsidy.

**First independent larger-scale target must be ALL PLANT SPECIES, not just larval-host candidate plants.** Use the same 1,596 original source-fixed port locations (and 1,055 country/urban matched cohort) to construct species-ID-accepted **all-Plantae** public GBIF inventories, with fixed 5km spatial footprint and explicit original observation year/effort. Never infer native/alien from host-catalogue status alone; source a native status registry before an alien/native ALL-flora ratio can be claimed. All-recorded-plant species richness, trait diversity and site-frequency controls are distinct estimands. The available original archived candidate-host table cannot simply be renamed 'all flora'. Unregistered GBIF records are missing observations, not absence. All future new sample/gates must be frozen before new plant results.

Second macro line: independently measured multiple-year (2018+2022) benthic EQR and source-reported **mass** versus modeled WWp *volume* (clearly separate), with region/habitat quality and independent negative-controls. Do not paste the two sample tables as if they measured the same organisms/outcome.

**Methodological stop:** do not continue splitting this question by butterfly species, single outlets or site buffers for now; those detailed local Draft PRs stay as secondary falsifiers. Keep all source positive findings alongside their country-generalization failures. New macro claims require independent all-Plantae response or measured receiving-water ecology, not only a reregression of a selectively chosen favorable host/photo denominator.

## This new macro cross-country reproducibility

- [Original source-verified country/2°/WWp multi-country fully detailed results](https://github.com/zuizui0223/gesui/blob/analysis/global-2015-WWp-positive-country-transfer-20261110/results/global_WWp2015_country_transfer_20261110/GLOBAL_1055_PORT_42_COUNTRY_WWp_OUTOFCOUNTRY_AND_GEO_HOLDOUT.json) — [source recovery and exact old β + geographic/country heldout run SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42791826892).
- [Independent original file SHA, every heldout country/port MSE, geographic source block allocations and 2000 whole-country resampling audit](https://github.com/zuizui0223/gesui/blob/analysis/global-2015-WWp-positive-country-transfer-20261110/results/global_WWp2015_country_transfer_20261110/GLOBAL_1055PORT_42COUNTRY_WWp_INDEPENDENT_OFFLINE_HOLDOUT_QA.json) — [independent QA run SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42792225053).
- [Exact source-frozen posthoc country-generalization protocol, including original 1,596 port, 2015 Jones, 2021 WorldCover and original chocho model](https://github.com/zuizui0223/gesui/blob/analysis/global-2015-WWp-positive-country-transfer-20261110/docs/GLOBAL_WWp2015_1055PORT_42COUNTRY_POSITIVE_TRANSFER_PROTOCOL_20261110.md).

**KEEP Draft / no `gesui/main` merge; original `chocho` manuscript and archived earlier positive records are unchanged.**