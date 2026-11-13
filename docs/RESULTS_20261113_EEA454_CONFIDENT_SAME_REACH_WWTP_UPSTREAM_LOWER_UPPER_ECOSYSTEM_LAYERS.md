# Source direction specificity after resolving within-reach WWTP and EEA monitoring-site positions

**13 November 2026 · Original-source-locked 232 European native HYBAS7 /454 native H7×2017–19/2020–22 records · Posthoc observational macroecology · Draft/unmerged.**

## Why this additional audit was necessary

PR #68 showed that the older EQR heldout MSE gain from **registered HydroWASTE wastewater flux accumulated to the H7 BASIN OUTLET (+11.129%)** shrank markedly when discharge was routed only from source river reaches geometrically upstream of the **actual EEA monitoring reach (+1.434%, uncertain)**. A greater additional gain was also observed from the **DOWNSTREAM** facility geography (+5.782%, uncertain), which cannot directly cause pollution flowing upstream by ordinary gravity-fed river transport. The downstream comparator is not balanced for source density and is not a randomized placebo. The original geometry source still classified **471 of 3,525 original station-year records** as having same-reach facility (within native `HYRIV_ID`) but unresolved source/station ordering. Address them through original river geometry, not by expanding radius or changing source H7 cohort.

## Original source-only geometry quality gate, before inspecting new biological results

- Original 2018 and 2022 EEA native river monitoring source CSVs and geographic station coordinates (Git blobs `0f30bfbc3dac13f6a182c8a910e9899c04778322` and `325f199e985d9862ffee5b0edb205c61a7be3833`), unchanged original 454 H7×periods /232 H7. Original EEA original site N=3,525; original 2,382 distinct monitored HYRIV reach IDs.
- Official original [HydroRIVERS Europe/Middle-East v1](https://www.hydrosheds.org/products/hydrorivers), **938,544 reach attributes and actual polyline geometry**, full publisher ZIP SHA256 `500da7d36ceee0aa4c82dd625f5c53d80ba2f7ac14ad022cd3a6601bf55068e5`. `NEXT_DOWN` links official line reaches; official technical documentation does NOT assert that a stored shape vertex order is reliably upstream-to-downstream. Infer orientation by which reach endpoint joins the **original next-downstream** neighboring polyline endpoint in **locally metric projected coordinates**, not the first or last shape coordinate itself.
- Genuine HydroWASTE publisher original ZIP SHA256 `6e19d8f31d95d533b2b29056b823eeb915d933fa22ea090f234dc38f6427d044`, original 58,502 registered WWTP with native `HYRIV_ID`, `LAT_OUT`, `LON_OUT`, modelled/mixed-provenance `WASTE_DIS` m³/day; STATUS unreported for most facilities.
- Source-only fixed quality rules **before new biological coefficients**: single-part river line, length >250m; correct upstream/downstream endpoint determined using `NEXT_DOWN` neighbor with connection distance ≤500m and ≥150m advantage to one endpoint; original EEA monitoring site and outfall must each snap within 500m of that exact original polyline; along-river source–station separation ≥500m and ≥2 times their combined snapping offsets. Otherwise SOURCE INDETERMINATE, not a genuine no-emission case.
- **284 original river reaches** in this same-source candidate set; 278 can be oriented from next downstream geometry and 6 cannot. **471 station-year** source cases correspond to **700 original station×registered WWTP pairs**, among which **352 are confidently UPSTREAM**, **222 confidently DOWNSTREAM**, and **126 remain unresolved** because direction/length/snap-distance does not justify a definite order. Original geometry gate passed with 464/471 station-year source cases position eligible, 306 distinct original monitoring IDs, 12 EEA source countries and 352 new confident upstream pairs. The 352 **plant×site pairs** may repeat the same facility or station across years, so are NOT 352 independent plants, independent field stations or causal treatment replicates.
- Original 152 H7×period cells have at least 1 same-reach registered facility confidently upstream; **151** have positive original estimated `WASTE_DIS` m³/day from those confidently ordered upstream sources. Rejected unknown source can only enter an explicit possible-UPPER bound, not a selected positive upstream lower-bound.
- Source-only original 471-station-derived artifact SHA256 **`9bf89e7c372cfd294113f913955f48dc1b6262d9435b8918b208758e9005d302`**.
- [Genuine original source-only river vector/next-downstream endpoint orientation and 700-source pair quality gate SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43364964322).

## Same exact 454 H7×period three-biodiversity-layer comparison with position uncertainty bounds

All earlier original source data, 232 H7 basins, 90 authentic HYBAS6 parent geographic holdout units, EEA country/BQE methods, true HydroATLAS upstream crops/urban/natural flow, original modeled site discharge, native original GBIF Animalia/Plantae/five-family source counts, 2018/22 original EEA EQR and original plant speciesKey S50 (50-rarefied) / Lycaenidae five-family source log odds unchanged.

For each original EEA monitoring site-year in a fixed H7×period, add newly source-position-confirmed **same-reach upstream WASTE_DIS /(site modeled river discharge ×86400)** to the previously validated STRICTLY upstream DISTINCT-reach source. Average once across original **all** EEA field monitoring stations within the H7-period, not only stations with a same-reach plant; take log1p of the resulting site-equal mean. Compute a **possible UPPER** bound by also adding *unresolved* same-reach source, without claiming it is upstream. Same-reach confidently downstream source is added only to the original downstream geographic **comparison**, not the upstream numerator. All source measures are static modeled/mixed-source facility geography, NOT observed annual nutrient concentrations, verified operations or an intervention.

Fixed train-only scaling, Ridge α10, parent HYBAS6 90-unit 5fold group CV. Additional registered source uses same source-only thresholds and identical 454 rows. Original country-out source validation is conducted separately with country fixed effects removed (therefore not an identical model specification).

### Signed H6-parent-heldout model results (all same 454 original rows)

| WWTP source pressure added to identical landcover/observation B2 baseline | EEA actual benthic EQR OOS MSE gain | GBIF indexed plant S50 | GBIF Lycaenidae 5-family log odds |
|---|---:|---:|---:|
| **Previous registered modeled flow to H7 BASIN OUTLET** | **+11.129%** | −0.196% | −0.762% |
| Previous STRICT distinct-reach UPSTREAM to monitoring reach | **+1.434%** | −1.297% | +0.449% |
| **NEW strict upstream + confidently upstream SAME-reach registered facilities (LOWER bound)** | **+1.874%** | −1.598% | +0.462% |
| NEW strict upstream + confirmed + all unresolved same-reach possibly upstream (**UPPER source bound**) | +2.422% | −1.832% | +0.057% |
| Downstream 50km **registered source geography only** | +5.782% | +0.614% | +1.977% |
| Downstream 50km + confirmed downstream SAME-reach source | +5.848% | +0.617% | +1.987% |

Original 90-H6-parent 2,000 resample descriptive 95%:
- Prior strict upstream EQR +1.434% **[−0.580%,+3.577%]**;
- **NEW confidently ordered same-reach lower bound EQR +1.874% [−0.661%,+4.624%]**, still crosses zero;
- NEW possible upper bound EQR +2.422% **[−0.647%,+5.728%]**, still crosses zero;
- Adding confident same-reach source to previous strict receiver model itself yields a modest **+0.447%** incremental EQR gain **[−0.186%,+1.163%]**, still crosses zero;
- After adding confidence-ranked same-reach upstream source, including the same old downstream 50km geographic index further improves original EQR OOS MSE +6.584% **[−0.482%,+13.792%]**, still not a balanced randomized negative control.
- Old modeled outlet +11.129% H6 interval [+3.944%,+18.267%] is **not** a true site-specific source exposure effect. New N/P numerical biological coefficients are not implied.

**Country-heldout transfer** (11 stable original H7-country groups): EQR strict source to receiver +1.544% as before; NEW confident same-reach upstream lower bound **+1.547%**; possible upper +2.615%; modeled basin outlet +0.953%; downstream registered 50km +11.432%, and downstream-with-geometry +11.521%. This does NOT resolve independent transnational source-response transfer; country FE was omitted for unseen-country methods.

**Empirical conclusion:** resolving the previously uncertain within-reach WWTP/source direction increases the genuine position-confirmed upward pressure predictive gain only slightly (**+1.43%→+1.87%**) and the 90 hydrologic-parent bootstrap interval still spans zero. The earlier strong **+11.1% predicted EQR pattern is mostly a basin-outlet regional/source-geography signal**, not evidence for a similarly strong discharge directly reaching the sampled EEA benthic biological community. Registered facility geometry, source-date/operation reporting, river size, downstream urban settlement and country method differences remain confounds. Plant S50 and Lycaenidae indexed reports show no consistent confirmed-upstream positive effect. **No evidence for a causal water→plants→butterflies pollination/service cascade.**

[Original 454 native H7×period three-layer actual river geometry model and H6/country heldout results SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43365344082).

## Independent second program numerical QA

A second **scikit-learn Ridge/SVD** implementation using the same exact original source-frozen EEA 2018/22 station modelled river flows and confident 471 same-reach source entries, original 454 H7-period source and existing environment recomputes the original manual NumPy regression **3 original biological response families ×17 WWTP exposure model specifications ×2 geographic evaluation designs =102 OOS MSE values**, all identical to numerical tolerance. Original biological source population is unchanged and not independently field-replicated.

[Independent source-locked 102-model original geometry/river pressure QA SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43365570643).

## What is missing for broader, genuinely ecological-service macro analysis

The immediate defensible next step is a **source observation-process audit**: original facility-level active status (mostly absent), actually reported *effluent nutrient concentration* matched to original outfalls/year, source reporting completeness by country and intervention-date, and direct receiving-water TN/TP when available. Another promising physical check is equal source-opportunity matching of upstream/downstream river distance, landcover/river size, reporting intensity; the present downstream control is geographically heterogeneous and should not be called a source-specific causal falsification.

Globally unifying the existing 1219 HYBAS7 source ×2 time bins for aquatic EQR/EPT, plant species index and Lycaenidae indexed butterfly composition was achieved in PR #63, but lacks true shared biological **function** measurements: adult aquatic insect emergence/prey biomass, riparian consumer diet, plant larval host changes, standardized butterfly visitation/pollen delivery and seed set. Existing nonoverlapping Oregon/Canada studies cannot be synthetically chained into independent 454-EU watershed rows.

**Preserve** earlier PR #68/67/66/65 and genuine positive/negative outcomes. No merged main/chocho, no eBird, no raw licensed HydroRIVERS source redistribution. All new original source/model records remain explicitly **POSTHOC exploratory observational**.
