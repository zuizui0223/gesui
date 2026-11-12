# Empirical aquatic–terrestrial ecological-function synthesis: quality, quantity, and timing are distinct (2026-11-12)

**Status:** new HJA AS015 genuine original-data *posthoc* seasonal reanalysis and independent arithmetic audit, plus official Dryad public-source metadata inventory for two later field panels. Draft / NOT merged / no invented ecosystem service mediation.

## New measured result: seasonal emergence is pulsed; EPT and midges do not usually compensate in time

**Original source**: H.J. Andrews Forest LTER AS015 v2, summer 2022, DOI `10.6073/pasta/6e9b83e2547e3a82d92dab068476f1d7`. Original exact public CSV MD5 `f3e6e400990da883ae0d341731431749`, SHA256 `c8ae00e91e7b210fbc86f5b08c6047ceaeade30bb2d410c3748f495fea0de4b2`. Source-only pre-outcome **17-stream ×23-period×4-order = 1564** keys frozen in original SHA256 `51354ce5cc3ab3acc63aecc28dd93534438ab3ffaf09b86ab2b0994b51874923` before adult numeric magnitudes were accessed. Original source verification [SUCCESS #43125631204](https://github.com/zuizui0223/gesui/actions/runs/43125631204).

**Current analysis** chose seasonal covariance and pulses *after* the earlier EPT share versus whole emergence biomass sign had been viewed. Thus new outputs below are exploratory, not independently preregistered confirmation.

Observation unit: 17 first/second-order Oregon streams **all in the same Lookout Creek basin**, nine intermittent and eight perennial, 23 3–4-day sampling intervals in 2022. Four source groups E+P+T and Diptera (Chironomidae); total emergence dry biomass measured/estimated per area (g/m² per sampled interval). These are *potential aquatic-to-terrestrial prey resource* metrics, not actual predator assimilation, all-insect emergence or ecosystem services.

| Original site-level source-authentic 2022 metric | Newly computed result |
|---|---|
| Portion of each stream's 23-event summed biomass contained in its **largest five intervals** | **median 0.49988 (50.0%)**, site range [0.37093,0.70367] |
| Baseline equal-duration fraction for five of 23 intervals | **5/23=0.21739**, therefore top-5 biomass share **2.30×** this time fraction at the median site |
| Temporal CV of 23 interval total biomass | **median 0.87092**, range [0.48403,1.37837] |
| EPT adult **biomass** fraction pooled over 23 events | **median 0.98812**, site range [0.93033,0.99806] |
| EPT adult **abundance** fraction pooled over 23 events | **median 0.93116**, range [0.61176,0.98630] |
| Difference between per-stream EPT biomass and abundance fractions | **median +0.06084**, showing counts vs dry mass weight taxa differently |
| Between-group covariance of EPT and Diptera biomass fluctuations across 23 periods | **16/17 positive**, **1/17 negative**; no temporal antiphase stabilization across most streams |
| Two-component temporal compensation index `C=−2 cov(EPT,D)/(Var(EPT)+Var(D))` | **median −0.01807**, mean −0.02988; 10k stream bootstrap descriptive95 median [−0.04612,−0.00323] |
| Synchrony index `Var(EPT+D)/(SD(EPT)+SD(D))²` | **median 0.97482**; caution: EPT overwhelming mass dominance influences this value |
| In original 23-period circular-shift comparison, observed variance below shifted-series median | **only 1/17 streams**, descriptive timing check; seasonal shift not randomized intervention |
| Intermittent nine-stream mean top-5 fraction | **0.53541** (versus perennial eight-stream 0.47563); **no causal regime contrast established** |

**Interpretation:** these source-authentic adult insects form a **strongly time-concentrated subsidy pool**, and relative adult abundance does not translate directly to biomass or temporal continuity. Diptera (chironomids) were not an observed strong temporal stabilizer for EPT biomass in this *one basin*, partly because EPT dominate measured dry mass. No evidence that wastewater created pulses or that larger pulses improve predation/pollination is available.

[New original measured 17-stream temporal source-SHA run SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43142773724). [Separate pure-Python independent check of all 17 Var identities (max residual 1.73×10^-18), EPT fraction/5-of-23 totals and source SHAs SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43143044572). This arithmetic audit is independent code, **not a new field sample**.

### Relationship to previous original aquatic result (PR #60)

In identical 17 streams, original stream-season EPT adult *individual count fraction* vs total emerging dry biomass Spearman **ρ−0.16176** (20k shuffled labels reference fraction 0.5318; 10k stream bootstrap range [−0.708,+0.457]). This is not a proof of no effect. Varying the EPT share definition flips numeric signs (count-pooled +0.054, biomass-pooled +0.387), so only a *measurement limitation* is supported. That supports treating **quantity, composition, body size, and temporal supply** as distinct axes, not collapsing all into EQR or an absolute subsidy effect.

## Newly verified source availability of real riparian consumers

**Panel A, South Nation watershed (Ontario, 29 physical river/ditch locations):**
- Original [Dryad dataset](https://doi.org/10.5061/dryad.sbcc2frfk) published 2025-01-28, public metadata version6 / 9 original source files in official API version342267. Original list verified through official API, including `benthic_classifications.csv`, `benthic.csv`, `bird.csv`, `emergence_seasonality.csv`, `emergence.csv`, `environment.csv`, `sites_bird.csv`, `sites_full.csv` and `README.md`. All original SHA digests and bytes are in [official public version filelist audit run](https://github.com/zuizui0223/gesui/actions/runs/43143134955).
- Field data within same named study sites cover insect emergence + bird point counts in **2022**, but water chemistry and benthic insect material in **Sept 2023**. **Do NOT call the full N/P→benthos→emergence→birds chain same-year or estimate a causal serial mediation**.
- Original authors [Rideout et al. 2025, doi:10.3389/fsufs.2024.1484377](https://doi.org/10.3389/fsufs.2024.1484377) reported that high drainage-ditch emergence biomass did not predict higher insectivorous riparian bird fractions, while aquatic and bird assemblages covaried (reported Procrustes r=0.44, p=0.023), subject to landscape and single bird point-count detection controls. **Published conclusion only, not a newly reproduced result.**

**Panel B, Quebec two urban rivers /20 field sampling sites:**
- Original [Dryad dataset](https://doi.org/10.5061/dryad.v6wwpzh8z) published 2026-03-10, public metadata version7 / 2 source files including `File_3.zip` and `README.md`. File metadata and original SHA digest pinned in same [official public filelist audit](https://github.com/zuizui0223/gesui/actions/runs/43143134955).
- Original investigators [Gagnon et al. 2026, doi:10.1002/oik.12080](https://doi.org/10.1002/oik.12080) reported emergence body-size downsizing with urban imperviousness, without strong total biomass change; isotope mixing-model Tetragnatha aquatic diet fraction **mean 52%, range 18–85%**. This is *an existing published result*, not independent reanalysis, and 20 stations are not 20 independent catchments.
- Field data from **May–June 2024** include emerging insect size/abundance, habitat, conductivity, spider family counts and spider/prey isotope ratios and MixSIAR outputs (per public README). Same-site alignment still requires original CSV join auditing.

**Data access limitation:** Public Dryad API v2 datasets, versions, and file-name+checksum metadata returned HTTP200 for both originals; the official public binary file routes `/api/v2/files/<id>/download` returned **401 Unauthorized**, despite CC0 metadata. The web downloadable routes appear authorization or human-access gated in 2026. We **did not** bypass authentication or claim that original raw consumer coefficients were independently verified. No new Dryad actual data/site join or ecological model was performed. The source-only gate is currently **STOP_PENDING_AUTHORIZED_ORIGINAL_BYTES**. Preserve original source-file lists and retrievable identifiers, and do not fabricate values.

[Source-only dataset+version+file list and 401 audit](https://github.com/zuizui0223/gesui/actions/runs/43142032918) · [Complete source-only file metadata official API audit](https://github.com/zuizui0223/gesui/actions/runs/43143134955).

## Exact evidence architecture for a future paper

| Layer | Original dataset / unit | Current inference |
|---|---|---|
| Registered WWTP reported N/P → benthic quality | EEA Europe 2018/2022, regional holdouts, distinct years with some repeated stations | tiny reproducible predictive increment; NOT confirmed delivered dose or ecosystem service |
| Actual stream TN/TP ↔ EPT count fraction | EPA NRSA US 2018–19, n~1,220 visits | positive robustness of **negative TP–EPT fraction association**; NOT wastewater-attributed |
| Riparian vegetation structure ↔ EPT fraction | EPA NRSA same-UID field structural cover, n~1,184 original/1,219 raw complete cases | observational co-organization, independent vegetation predictive gain and TP interaction uncertain |
| Adult EPT abundance fraction ↔ total emerging dry mass | HJA AS015 summer2022, n17 Oregon headwaters one catchment | no reliable transport equation, wide uncertainty; not same sites as EPA NRSA |
| Emerging adult biomass **temporal continuity** | **NEW same HJA AS015**, 17×23 actual intervals | half of total emerges in top five intervals (median), mostly EPT biomass, weak compensation by midges |
| Emergence body-size/quantity ↔ terrestrial predators' realized diet | *Published* 20 Quebec urban-stream sites; original Dryad bytes inaccessible to this runner | size & diet changes reported by authors, not independently validated here |
| Emergence biomass ↔ field bird community | *Published* 29 Ontario agricultural locations, data temporally misaligned on water chemistry | community co-response but no clear quantity→bird insectivory relation as reported by authors |
| Butterfly larval hosts/flowers ↔ realized pollination service | Chocho port-island resources and sparse GBIF; same-river source threshold only 9 supported butterfly sites | NO same-site pollination, pollen transfer or seed-set evidence |

**Single unifying question:** *When does environmental stress decouple aquatic taxonomic integrity from the quantity, timing and quality of subsidies available to terrestrial consumers?*

For **global/macro ecological inference**, 17 HJA streams and 20/29 sites in other watersheds are still far too few independent regions; the current three published field designs can inform an empirical **measurement framework**, not a general global wastewater-service causal estimate. Construct a new global `water–land subsidy` source registry keyed by original study, river, site, year, area, trapping protocol, emergence biomass vs abundance vs trait distribution, consumer evidence, source chemical attribution and spatial dependence. Avoid pooling mixed-unit coefficients or selling a global SEM on 17+20+29 as 66 independent river ecosystems.

Keep draft and unmerged, no changes to original `gesui/main`, `chocho`, PR #59 or #60. Avoid eBird.
