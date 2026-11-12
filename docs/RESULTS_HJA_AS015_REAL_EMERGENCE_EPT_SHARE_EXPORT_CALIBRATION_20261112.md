# Direct measurement of aquatic-insect emergence: does EPT composition proxy total aquatic→terrestrial insect biomass?

**2026-11-12 | SOURCE-LOCKED numerical test + posthoc independent arithmetic sensitivity | DRAFT, not causation.**

## Source, study, and process

Authentic original **H.J. Andrews LTER AS015 v2** (Latorre-B & Arismendi 2026), DOI [10.6073/pasta/6e9b83e2547e3a82d92dab068476f1d7](https://doi.org/10.6073/pasta/6e9b83e2547e3a82d92dab068476f1d7), CC BY 4.0; EML [AS015_2.xml](https://andrewsforest.oregonstate.edu/sites/default/files/lter/fsdbdata/eml_archive/AS015_2.xml), official original file as specified in EML [AS01501_v1.csv](https://andrewsforest.oregonstate.edu/fsdbdata/cache/AS01501_v1.csv). Source original MD5 **`f3e6e400990da883ae0d341731431749`**, SHA256 **`c8ae00e91e7b210fbc86f5b08c6047ceaeade30bb2d410c3748f495fea0de4b2`**, EML SHA256 `fb4a6eb41374e7d682bb74b77a0841bd802198f8285951ca7c47a86786fc962e`.

The EDI PASTA entity-names API returned HTTP **403** from the GitHub runner, so source discovery switched **before any actual AS015 numerical outcomes** to the *publisher's original EML*; exactly the same original file MD5 was verified (not an alternative study or random website). All URLs/receipts archived in source audit run [43125631204](https://github.com/zuizui0223/gesui/actions/runs/43125631204).

The original pre-outcome verified native source had **17 distinct first/second-order headwater streams in one Oregon Lookout Creek basin**, 9 intermittent/8 perennial, **23 time intervals per stream**, four aquatic-insect groups **E+P+T and Chironomidae (Diptera)** each event, **391 complete site-period observations and 1564 complete taxon rows**, zero missing/count-invalid abundance or biomass. Source site/order/period/date ID manifest SHA256 **`51354ce5cc3ab3acc63aecc28dd93534438ab3ffaf09b86ab2b0994b51874923`** was frozen **before reading any adult emergence amounts**.

Data are actual measured adult **emergence abundance (per m²) and estimated dry biomass (g/m² per collection interval)**, from 4 emergence traps per reach. Sampling intervals were approximately 3–4 days. Do not label the published g/m² values as g/m²/day without duration correction. Chironomidae biomass is partly extrapolated from size estimates, per original metadata.

## First source-frozen adult emergence proxy experiment

**Primary spatial independent unit = one stream, n=17**. Each site mean uses 23 calendar collection periods. Do not treat 391 repeated time events as 391 independent streams. Average EPT share at each stream was the equal-interval arithmetic mean of `(E+P+T abundance)/(E+P+T+Chironomidae abundance)`, only on periods with positive total adult abundance; **18 of 391 site-periods had zero total adult abundance**, excluded from ratio (not set to zero). The response was site-mean total emergence dry biomass in reported g/m² per collection interval. Also compared adult count and flow regime.

**Original result [source-frozen action SUCCESS 43126089434](https://github.com/zuizui0223/gesui/actions/runs/43126089434):**

| 17-stream rank relationship or regime contrast | Result |
|---|---:|
| Adult EPT abundance fraction vs **total emergent biomass** | Spearman **ρ=−0.16176**, Monte-Carlo 20,000 rank-permutation two-sided reference fraction **0.5318** |
| 10,000 whole-stream resampled descriptive 95% on same rho | **[−0.7083, +0.4566]** |
| Adult EPT abundance fraction vs total emergent **abundance** | ρ **−0.17892**, descriptive 95% [−0.7372,+0.4801] |
| Site-season mean EPT adult abundance fraction | across streams mean **0.92490** |
| Within-stream seasonal count-fraction–biomass association | **12/17** site-level signs negative; median site rho **−0.3250**; these repeated observations are not independent regional replications |
| Stream mean EPT fraction I vs P | **0.9061 vs 0.9460**, exact exhaustive 9/8 label-permutation two-sided fraction 0.2569 |
| Stream mean adult biomass I vs P | **0.06977 vs 0.04029 g/m² per interval**, difference +0.02947, exhaustive label fraction **0.05471** |
| Stream mean adult abundance I vs P | **22.93 vs 16.98 individuals/m² per interval**, difference +5.95, label fraction 0.2660 |

**Not a proof of no relation.** Original EPT fractions are already high (ceiling effects) and the 17 sites share a catchment. Inferences are far less independent than 391 samples and dry biomass is size-derived. Both primary CI and permutation reference fraction are fully consistent with positive or negative underlying correlations.

## Newly independent pure-Python aggregate QA: numerator/denominator sensitivity

[Independent SHA-pinned site-level arithmetic audit SUCCESS 43126473521](https://github.com/zuizui0223/gesui/actions/runs/43126473521) used the original saved 17-stream numeric outcomes and reimplemented rank correlations without importing the first script or querying the field source a second time. Reproduced original Spearman −0.16176 exactly.

**The numeric association changes when the EPT 'proportion' is defined differently**, even on identical original stream sites:

| EPT composition definition vs SAME stream mean total biomass | Spearman rho |
|---|---:|
| Mean across time intervals of EPT **adult abundance proportion** (prespecified primary) | **−0.162** |
| Proportion after pooling adult **individual counts** across all intervals | **+0.054** |
| Mean per-interval EPT **biomass proportion** | **+0.157** |
| Proportion after pooling measured/estimated **biomass** across all intervals | **+0.387** |

Therefore **do not claim the sign of an EPT-composition→resource-export slope is stable across definitions**. Correlations with total biomass were small under several reasonable abundance-based definitions, but some biomass-weighted definitions yield positive descriptive rank associations (partly mathematical coupling of numerator/denominator). None provides a transferable quantitative prediction from benthic EPT larvae to emergence.

Other original-data QA:
- Across 17 streams, adult EPT mean per-interval abundance share range **0.682–0.990**, median **0.936** (high and range restricted). Mean total emergence biomass range **0.01450–0.14338 g/m² per interval**, a **9.89×** spread. So similar high EPT dominance is compatible with very different total exports in the observed sample.
- Excluding any one stream: Spearman primary range **−0.3471 to +0.00588** (16/17 negative, one weakly positive); this is influence testing **within the same Oregon basin**, not 17 replicated ecoregions.
- Sites with at least 20 periods showing positive adult abundance: **16** streams, Spearman **−0.1265**; one stream with 19 positive sampling periods is excluded only as a documented **posthoc** sensitivity.
- Within hydrological regimes, site-level rho **+0.417** for intermittent (n9), **−0.524** for perennial (n8); incompatible direction **at small n** is hypothesis-generating, cannot establish interaction or regime effect.
- Median-split (posthoc) opposite categories: **4** lower-EPT/higher-biomass streams and **5** higher-EPT/lower-biomass streams; classification depends on arbitrary sample medians.

## What this actually adds to the wastewater / ecosystem-service study

**IMPORTANT: This is cross-study *triangulation*, not observed serial mediation.**

| Step | Data / result | What is known | Missing connection |
|---|---|---|---|
| Reported WWTP point loads ↔ river biological integrity | EEA 2018/2022 real reported WWTP N/P + benthic EQR (separate Europe sites), modest positive heldout predictive gain | Source *geography* contains predictive information | No confirmed effective river concentration or causal exposure change |
| Stream measured TP ↔ larval benthic EPT fraction | EPA NRSA2018–19, original 1,220 same native UID visits, many states; stronger negative TP slope than TN | Within-site aquatic functional *composition* covaries negatively with TP | The measured TP is not attributed to WWTP, and benthic samples aren't adult emergence |
| Riparian structure ↔ aquatic EPT composition | EPA NRSA same original 1,184 visits, HUC8 heldout habitat addition after TP ~+1.28%, bootstrap includes 0 | Vegetation structure may add geographic information | No documented mediation / external trophic subsidy |
| Adult EPT share ↔ actual emergence biomass | **Independent AS015 17 Oregon streams 2022** | Direct adult resources measured; simple adult composition fraction does **not** reliably determine total emergence biomass in this small cohort | None of these 17 are the NRSA stream sites; no direct benthic larvae→adult emergence calibration under N/P |
| Plant hosts/flower resources → butterflies → seed set | GBIF record source limited; original same-site Lycaenidae quality gate n9 not 100 | Candidate host distribution in port/island data; not local pollination | Actual floral traits, butterfly visits and pollination processes unmeasured |

The new AS015 result is a **necessary measurement caution**: an observed change in **relative** aquatic EPT proportions cannot be converted into "X% less cross-boundary insect resource" or "X% less pollination" with no total emergence measurement. Equally, a weak correlation at n17 **does not prove resource supply remains stable** after nutrient pollution.

The macro ecological question should be framed:

> **When do changes in aquatic community integrity translate into changes in resource quantity and realized terrestrial functions?**

The first step toward a defensible connected causal study would be **same-watershed, same-year benthic larval counts + emerging adult biomass + measured nutrients + receptor diets/plant visits**, at multiple independent catchments (or ideally a standardized multi-site programme). Avoid joining separate data products' site means by ecological analogy or claiming pathway coefficients/SEM without shared units and time.

**Preserve existing positive results, scripts, and PR #59, #58; retain a separate Draft and do not merge to gesui/main or chocho.** No eBird.
