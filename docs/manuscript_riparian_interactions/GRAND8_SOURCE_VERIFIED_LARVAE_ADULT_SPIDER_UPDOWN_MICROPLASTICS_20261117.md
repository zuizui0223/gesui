# 2026-11-17｜下水処理場8地点：水生幼虫→羽化成虫→河岸クモで実測したマイクロプラスチック（原著再解析）

**Verified original data, exploratory reanalysis of Wardlaw et al. 2026, NOT a new independent publication, experiment, or source-specific N/P flux measurement.** Site/time-level specimen source gates were checked before numerical outcome magnitudes. This is external mechanistic context for the ongoing `gesui` macro paper rather than evidence that the France122 reported-effluent-negative sample contains a WWTP positive.

## Real original public source and site/time support

- [Original Grand River 2026 paper](https://doi.org/10.1016/j.envpol.2026.129116), primary published source [Borealis doi:10.5683/SP4/U0MN9S](https://doi.org/10.5683/SP4/U0MN9S), original file `22.23.MP.Invert.Summary.tab`, source file ID **1213261**, SHA256 **1b69eb88fe38f601b11ff466d6a9b3177546a9b489a74c61eaf2a64d770dec66**.
- **906 true insect / riparian spider original material-measurement rows**, `Site`, `Location` upstream/downstream, `Time` in 1,2,3, `Taxa`, `Life_Stage` with normalised metric `MP_per_mg` (microplastic particles / wet tissue mg).
- **8 original distinct WWTP identifiers**: Ayr, Caledonia, Guelph, Hespeler, Kitchener, Paris, Preston, Waterloo, all within **one Grand River network** rather than eight independent regional catchments.
- The same two Borealis original DOIs contain independent genuine effluent (**77**), biofilm (**240**), sediment (**240**), source aquatic 2022 insects (**426**), cross-aquatic-adult-spider specimens (**906**), stable isotopes (**918**) with **all 8 facility names shared**. Source-only original six-file SHA/countercheck: [GitHub Actions #44103189023 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/44103189023), original artifact **14337271487**. The `Time` codes of `906` and the effluent DOI were NOT validated as the same calendar-year/date, so **do not fit the six-file effluent mass→larvae→adult→spider mediation**. Complete ISO year field spans 2022/23 whereas original `906` table does not have explicit calendar year.

## Strict source-first paired response

Pre-read unit of eligibility: original site×`Time`×taxon with **at least 2 actual sample rows both upstream and downstream**, and **at least 2 matching visits for ≥6 of the 8 WWTPs**. One original facility is the independent spatial unit for describing up–down differences; within-facility repeated 3 visits are averaged equally before the 8-facility summary. No imputation or zero replacement for missed assays, no signed outcome-based target filtering, no pooling different taxon stages as equal independent recipients.

| Fixed source taxon/stage | Original site×visit pairs with genuine numeric `MP_per_mg` on both banks | Facilities with ≥2 valid visits (out of 8) | Source decision |
|---|---:|---:|---|
| Hydropsychidae larvae | **24** | **8** | GO exploratory |
| Hydropsychidae adult | **23** | **8** | GO exploratory |
| Tetragnathidae riparian spider | **24** | **8** | GO exploratory |
| Heptageniidae larvae | **24** | **8** | GO exploratory |
| Chironomidae larvae | **22** | **8** | GO exploratory |
| Chironomidae adults | **15** | **5** | **STOP source coverage**; may report descriptive observation, no full-stage pathway |
| Heptageniidae adults | **0** with ≥2 source samples per side | **0** | **STOP source coverage**, not biological zero |

Full original site/visit QA, exact source SHA and effect table: [GitHub Actions **#44104639473 SUCCESS**](https://github.com/zuizui0223/gesui/actions/runs/44104639473), artifact **14337581735**, reproducible [source-locked Python](https://github.com/zuizui0223/gesui/blob/docs/wastewater-riparian-interactions-outline/analysis/20261117_emergence_macro_source/grand8_true_WWTP_updown_larval_adult_spider_original_MP_source_gated.py). These are source-authentic; no third-party raw data files are redistributed.

### Actual 8-facility results for measured microplastic particle concentrations

For each original `Site×Time×Taxa`, compute mean upstream/downstream `MP_per_mg` across eligible actual original specimens; then average original visits equally within each facility, and the eight facility contrasts equally. No pseudocount added to make log ratios. Ratios below are **the median of the original eight facility-specific down/up ratios**, not ratio of pooled specimen counts.

| Taxon/stage | Facility-mean original upstream particles/mg | Downstream particles/mg | Facility-equal downstream−upstream particles/mg | Facilities with positive difference | Median of facility down/up ratios |
|---|---:|---:|---:|---:|---:|
| **Hydropsychidae larvae** | 0.02743 | 0.04296 | **+0.01553** | **8/8** | **1.62×** |
| Hydropsychidae adults | 0.02718 | 0.03236 | +0.00518 | 6/8 | 1.40× |
| **Tetragnathidae spiders** | 0.01087 | 0.01954 | **+0.00867** | **8/8** | **1.69×** |
| Heptageniidae larvae | 0.01324 | 0.02570 | +0.01246 | **8/8** | 1.73× |
| Chironomidae larvae | 0.06616 | 0.05290 | −0.01325 | 3/8 | 0.90× |
| Chironomidae adults | 0.03990 | 0.07346 | +0.03356 | 4/5 | 1.82× (5 facilities only, source STOP) |

All *eight* original Tetragnathidae site-level concentration differences downstream−upstream were >0, source facilities examples: Ayr +0.00310, Hespeler +0.01797, Waterloo +0.01508 particles/mg. Tetragnathidae ratios ranged from ~1.40× to ~2.76× among eight facilities; neither the absolute difference nor the ratio proves actual individual prey consumption. In contrast, Hydropsychidae adult changes were positive for only six original facilities; chironomid larvae changed in mixed directions. Hence **the identity of the biological channel matters**, even though several groups are part of the same aquatic→riparian network.

The facility-level results are *reanalyses of published author measurements* and are consistent in broad direction with the original authors' reported downstream spider elevation. The observational upstream/downstream design, repeated facilities within a river network, uneven measurement support, and non-observed physical food ingestion prevent a novel independent causal or global generality claim. **8/8 is a real, bounded descriptive pattern, not a valid global sign-test p-value from eight independent rivers**.

## Independent nutrient-assimilation contrast: Bow River 2022 5 sites

The separately published 2025 Bow River microbe/isotope study [source](https://github.com/Karen-Kidd-Lab-Ecotoxicology/Bow-River-Microbiome-2022) is a different river, pollutant/resource, and sample grain; actual isotopes **184** source organisms, **155** author-selected subset; five sites (3 upstream unexposed, 2 exposed), nested beneath **one** Bow River and three study treatment outfalls. The author isotope data show Tetragnathidae site-equal `δ13C` exposed−unexposed about **+3.14‰**, and a site-baseline-corrected `Tetragnathidae δ13C−local terrestrial leaves δ13C` difference ~**+2.61‰**, but δ15N after local leaf subtraction differs in sign (**−2.19‰**). 2023 isotope CSV `Date` is laboratory analysis date, not river sampling year. [Bow River original-data isotope audit](https://github.com/zuizui0223/gesui/actions/runs/44103188904). This supports isotope baseline heterogeneity; the small n5 and local geography do not yield five independent treatment WWTP effects.

**Do not pool these isotope differences (‰) with Grand River particle/mg or with France `TP→NDVI` into one numerical meta-analysis.** They are *complementary modes of identifying a cross-ecosystem channel*, not interchangeable measures of subsidy or one continuous wastewater→flower→butterfly trophic chain.

## Macro synthesis and next scientific gate

The France 2018→22 physically matched **122** river sites measured true summer TP and near/far 10m vegetation NDVI but original **2018 reported upstream WWTP P-positive=0/122** and no true transferred recipient. This Grand River dataset is an authentic **source-specific point-input and animal recipient** dataset with eight source locations, real up/down and three visits but all sources in one watershed. Together these demonstrate a **measurement-domain and biological channel difference**: plant greenness not detected along ambient TP changes is not a falsification of particle/specialized prey→spider transfer from known WWTPs.

Next defensible macro move, without selecting taxon based on the positive sign, is to check **additional independently sourced treatment catchments** for the *same unit* of downstream−upstream stage-specific material exposure, explicit source year, and ≥2 observed time pairs per facility. Need true river-catchment IDs to avoid assigning eight nested WWTPs to eight independent basins. No status from other basins is assumed.
