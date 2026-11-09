# Germany 2022 actual river TP source-interposed vs source-zero matched upstream/downstream monitoring

**November 9, 2026. Completed source-locked original EEA2026 WISE6 TP direct original-sample analysis.**

**Scientific status:** new actual same-river environmental monitoring **observational** comparison. Original source-conditioned monitoring quartet selection was committed **BEFORE** looking at any new 2022 TP numeric site-year values for these matched quartets. However the overall research direction follows an earlier 2022 WWTP–Lycaenidae photographic-positive analysis and is POSTHOC overall. The field data are not randomized effluent inputs. This study tests river Total-P **only**; no butterfly outcome, adult puddling protein/Na, host-plant fitness or sewage-driven biology.

## 1. Source-only 16 river monitoring pairs (8 exposure, 8 registered-source-zero) locked before chemistry

Original nationwide WISE6 2022 river TP date-qualified public DE/AT monitoring metadata and SHA-locked HydroRIVERS Europe v1 third-party converted directed reach topology identified **115 possible** same-country same-directed-river monitoring pairs within 50km, **69** with at least one qualifying original EEA 2022 both-positive measured N/P single direct freshwater (FW) outfall on the strictly intervening upstream-to-downstream river/tributary flow graph, **40** with no qualifying original registered intervening source, **6** ambiguous/unresolved. These are river measurement-location source classes, not actual station date-aligned effluent chemistry.

Without opening numeric TP, the source-only original design in [Draft PR #30](https://github.com/zuizui0223/gesui/pull/30) locked **8 matched quartets** (each one source-interposed upstream/downstream TP monitoring pair and one no-qualifying-registered-intervening-source upstream/downstream pair) and **32 different original WISE6 2022 monitoring station IDs**, no station reused across quartets. Source matching favored geographic block and similar river path length; original sources and controls 2022 same year. Six distinct source 2° geographic blocks among eight quartets; exactly four initially in the same country×2° upstream block as their controls. WWTP source candidates include tributaries entering the river path, not just the main channel; their P mass annual reports (t/year) do NOT establish the dose measured at the downstream station.

Original source-only dataset:
- [8+8 matched 2022 German river-monitor quartets without TP values](https://github.com/zuizui0223/gesui/blob/analysis/wise6-2022-national-DEAT-river-TP-pairs-20261109/results/wise6_2022_national_DEAT_TP_sourceonly/GERMANY_2022_WISE6_TP_8SOURCE_8ZERO_SAME_RIVER_MONITOR_QUARTETS_SOURCE_ONLY_LOCK.json).
- Physical source flow graph is original source-snap on SHA-locked EU HydroRIVERS **third-party converted** GeoParquet and EEA official 2022 report codes; NOT independently compared to the official HydroRIVERS shapefile. This limitation remains.

## 2. Independent original numeric 2022 Total-P site QC (source still fixed)

Original EEA2026 WISE6 `DisaggregatedData` publicly retrieved original full archive ZIP **2,024,284,753 bytes** SHA256 `2a2ef2a43a16ec068b7974a80ae5fca49e44697660a46ce3b91fb1813abe62d2`; all **102,569,057 original water-quality rows** streamed. Original river water Total phosphorus `CAS_7723-14-0` in **mg{P}/L**; DE 2022 original sampled records 3,690. Of the 32 source-frozen stations, **19** had already been independently recorded in previous (photo-region selected) numeric WISE6 quality source; **13 additional original national stations** were now directly read under exactly the same original-source unit, numeric and LOQ handling.

Frozen before numeric opening:
- 2022 WISE6 TP site annual median after median of all valid same-day replicate assays;
- below-LOQ records excluded, with **no zero or half-LOQ substitution**;
- ≥4 distinct quantified sample days in ≥3 original months, belowLOQ fraction ≤0.50;
- paired upstream and downstream stations must both pass and have ≥**2 quantified calendar months in common**;
- ≥**6 complete quartets** across ≥**4 original country×2° blocks**, or STOP without source/no-source TP gradient calculation.

**All 32 source-selected TP monitoring sites themselves met numeric LOQ/day/month sample QC.** However, the registered-source-zero *control* pair of original quartet **#1** had **zero overlapping valid sample months** between its upstream and downstream monitoring stations, so it was excluded without replacement. Exactly **7 complete quartets / 28 distinct original WISE6 monitoring stations / 5 source-defined country×2° blocks** passed the month-overlap and source-selection gate. Remaining source-positive-between station pairs have original EEA annual reported intervening P mass **0.50–6.04304 tonnes/year**; 13.573 t/year original source path quartet #1 was excluded for source-independent control-month mismatch, not for measured TP effect.

[Original 32-site EEA numerical TP and 7-quartet quality gate workflow SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42587838833). [Source-first 32-site numeric QC and complete missingness](https://github.com/zuizui0223/gesui/blob/analysis/wise6-2022-DE-quartet-32-station-numericTP-20261109/results/wise6_2022_DE_8source8zero_quartet/2022_ORIGINAL_WISE6_32STATION_NUMERIC_QUALIFICATION_BEFORE_TP_DIFFERENCES.json).

## 3. Actual same-river source-versus-control 2022 TP gradient contrast

On each of the 7 matched quartets, compute:

- Original WWTP-interposed river monitoring **source gradient** = `log(TP_downstream + 0.001 mg P/L) − log(TP_upstream + 0.001 mg P/L)`.
- Registered-source-zero paired river monitoring **control gradient** = same upstream→downstream log concentration change.
- Original quartet contribution = **source gradient − matched registered-zero control gradient**.

The +0.001 mg P/L additive constant, seven-quartet support and month-overlap thresholds were fixed before viewing these new TP concentration differences. Each complete quartet is one matched observation, rather than treating the 28 original monitoring stations as independent source-treatment units.

| Observational matched 2022 TP gradient result | Value |
| --- | ---: |
| **Complete independent quartets** | **7** |
| Distinct original monitoring stations | 28 |
| Geographic 2° source blocks | 5 |
| **Mean matched source-minus-control log-TP gradient** | **+0.09045** |
| Median matched source-minus-control log-TP gradient | +0.20103 |
| Quartets favoring source relative increase / decrease | **5 positive, 2 negative** |
| 9,999-quartet resampling 95% descriptive mean interval | **[−0.18161, +0.36327]** |
| Exact two-sided matched sign-flip reference fraction, 128 configurations | **0.546875** (NOT causal/randomized p-value) |
| Median raw source-pair upstream/downstream logTP gradient | +0.11248 |
| Median raw registered-source-zero-pair logTP gradient | +0.18590 |

**Key influence check (posthoc, transparent):** leave one of seven original selected quartets out at a time; the mean source-minus-control gradient ranges **−0.00810 to +0.17538**. Excluding original quartet **#6** changes the mean sign to **−0.00810**, hence the apparent positive direction is **fragile to one original monitoring quartet**. Source pair #4 has a positive matched contrast even though both its river arms have *lower downstream TP*; the control pair's decline is stronger. Thus positive matched contrast is not necessarily actual downstream eutrophication.

**Conclusion:** The preselected 7 German matched source+control quartets show a weak positive average logTP *relative-gradient* contrast, but the distribution spans zero and the sign is sensitive to individual original pairs. There is **no robust identifiable increase** in TP caused by registered WWTP outfall interposition. A negative contrast would also not disprove actual physical effluent contributions; these are not randomized controls and original substrate/flow dates were not matched day by day.

[Actual original source-quartet observational TP calculation SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42588536120). [Every original 2022 measured river TP upstream/downstream pair and matched quartet difference](https://github.com/zuizui0223/gesui/blob/analysis/wise6-2022-DE-quartet-32-station-numericTP-20261109/results/wise6_2022_DE_8source8zero_quartet/2022_WISE6_REAL_TP_MATCHED_SOURCE_MINUS_SOURCEZERO_UPSTREAM_DOWNSTREAM_GRADIENT.json).

## 4. Evidence hierarchy for wastewater and butterfly ecology

1. **Completed source-specific environmental negatives:** The independent 2022 168-river-reach WISE6 TP study on 2021–2023 photo-source eligible regions had **strictly upstream registered 2022 WWTP P** negative out-of-region incremental TP prediction **−1.88%** MSE gain, and downstream source-placebo **−2.43%**. Posthoc source-label permutations did not support a unique geography-specific WWTP source indicator. [Different PR #31](https://github.com/zuizui0223/gesui/pull/31).
2. **This new national German within-river contrast:** independently identifies 16 hydrologically connected source-labeled TP monitoring pairs and 7 complete WWTP-interposed/control quartets under strict source-first matching. The mean source-v-control TP contrast is **+0.09 log units**, but 95% descriptive interval crosses zero and sign flips on excluding one quartet. This neither rescues the earlier null nor proves a treatment effect; it provides a distinct, physical-scale, source-specific comparison.
3. **Earlier DE/AT country×1° reported WWTP P vs Lycaenidae iNaturalist photographic reporting** initially had +17.8% posthoc out-of-region MSE gain without the original 2021/2022 butterfly photographic histories, but historical 2021–2023 photo controls, ESA WorldCover2021 land, and true 2022 WISE6 river TP did not identify a causal mechanism; recorded photo composition is **not population size**.
4. **Unchanged separate positive ecological observations:** `chocho` original 22/25 Lycaenidae species expanded *potential introduced-host larval resource geography*, not WWTP effects; `gesui` aquatic 2018/2022 source-aligned macroinvertebrate EQR predictive gains remain a different river response.

## 5. Limits to fix before scientific submission

- **Hydrology verification:** original HydroRIVERS converted third-party v1 GeoParquet not independent verified against original official release; source WWTP reporter annual P t/year and exact within-reach order unknown; the original source interposition includes contributory tributaries.
- **Time alignment:** annual 2022 WISE6 upstream/downstream monitoring site medians may represent **different observed sample dates** despite ≥2 common months. Same-day upstream/downstream TP sampling would provide a cleaner natural experiment.
- **Water-quality QC:** original WISE6 extra expert `resultObservationStatus`, flag-based exclusion and Eionet LOQ codes beyond numeric LOQ status not independently audited, may alter acceptance.
- **Small n and confounding:** 7 matched quartets in only Germany, some original 2° blocks repeat; one extreme #6 determines positive average sign. Natural tributary TP, road/crop/built/industrial runoff and river dilution unresolved. Exact source-positive original WWTP is not randomized, original source-zero means only no *qualifying registered* WWTP in the given flow path, not water free of human pollutants.
- **No Lycaenidae observation on these sites:** no adult puddling substrate, Na/protein, ant association, larval plant nutrient quality, or standardized butterfly individual counts. Do not treat measured water TP as nutrient delivered to butterflies.

**Scientific take-home:** Source-specific real water-quality cohort and paired comparisons are now reproducible and substantially more directly hydrologic than broad photo/source co-location. Yet the data do **not** establish a repeatable wastewater TP enrichment pathway or source-specific butterfly ecology effect. All original PRs remain Draft and original positive/negative estimates are retained.
