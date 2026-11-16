# Original EEA upstream vs downstream N/P source-report negative-control audit (2026-11-15)

**Status: exploratory, SOURCE-LOCKED, noncausal.** Stacked on PR #91 (same-site Water Quality ICM 2026 source feasibility) and PR #90 (annual reported N/P HNS/EQR/plant analysis). Source-only direction audit GitHub Actions run [#43679382449](https://github.com/zuizui0223/gesui/actions/runs/43679382449) SUCCESS. Local 7 tests PASS, independent numerical reimplementation checked.

## Source population & critical units

- 454 original EEA H7×periods, 232 original H7 basin identifiers, 90 true H6 parent basins, both periods nested within H7; original EEA biology EQR monitor n=3,525 station×year.
- Upstream reported WWTP N/P source positive in **132 station×years; 56 H7×periods**. Downstream annual-source-report placebo is **249 station×years; 61 H7×periods**; both upstream and downstream positive within **35 H7×periods**.
- Original frozen 2018 EEA source CSV Git blob `0f30bfbc3dac13f6a182c8a910e9899c04778322`; 2022 blob `325f199e985d9862ffee5b0edb205c61a7be3833`. Exactly reuses same native EEA source sites and HydroRIVERS direction labels as earlier EEA reports.
- Every upstream report-positive count and original per-H7 source N/P magnitude matches the frozen original upstream-only source SHA256 `012ad79feabc8c11e0d74f803730b5b7913b96a64a6af26fb0c97263d9aa7c26` to numerical tolerance 1e-9.
- `downstream` fields are *topologically classified annual reported EEA WWTP source loads* computed in existing original EEA source result, **not a Waterbase-observed river concentration**, field experiment, or truly randomized placebo treatment.
- Report-matching index 0 = no matched positive reported source, **not** absence of treatment plants, actual N/P emissions, or baseline upstream chemistry. H7-wide GBIF host resources do not necessarily occur at the EEA monitored reach.

## H6-disjoint true-spatial-parent out-of-fold predictions

Identical original 454 H7×period rows and same 90 H6 splits for each outcome; PR90 / PR88 train-only standardization and country×EEA method + period categorical baseline, original crop, urban, flow, area, EEA/GBIF recording opportunity, static HydroWASTE W, and HNS trait coverage + five-family record ratio controls. Ridge alpha10; H6-equal out-of-fold MSE. Each one-direction model adds one predictor to the same baseline; **do not subtract MSE gains from different bases**.

| Outcome | Source add | Upstream gain % [descriptive H6 bootstrap95] | Downstream gain % [descriptive H6 bootstrap95] | Up-minus-down gain pp [95] |
|---|---|---:|---:|---:|
| EQR | report-positive indicator | +1.920 [-1.298,+5.122] | +2.170 [-0.863,+5.294] | -0.249 [-4.132,+4.003] |
| EQR | N report/flow log equiv | -0.009 [-4.338,+3.645] | +3.497 [-0.957,+8.208] | -3.506 [-9.153,+2.126] |
| EQR | P report/flow log equiv | +1.872 [-1.366,+5.486] | +2.311 [-1.007,+6.035] | -0.439 [-4.452,+3.980] |
| HNS | report-positive indicator | +3.276 [-0.893,+6.304] | -0.093 [-4.673,+4.367] | +3.369 [-3.152,+9.555] |
| HNS | N report/flow log equiv | +0.530 [-1.382,+2.054] | +0.560 [-3.915,+4.292] | -0.030 [-4.392,+4.491] |
| HNS | P report/flow log equiv | +1.308 [-0.850,+3.232] | -0.068 [-3.357,+3.894] | +1.376 [-3.167,+5.446] |
| plant_S50 | report-positive indicator | -0.339 [-0.670,+0.099] | -2.200 [-5.220,-0.306] | +1.861 [+0.033,+4.886] |
| plant_S50 | N report/flow log equiv | +0.192 [-0.834,+1.394] | +1.326 [-0.832,+4.652] | -1.134 [-4.988,+1.752] |
| plant_S50 | P report/flow log equiv | -0.159 [-1.599,+1.174] | +0.123 [-1.157,+1.980] | -0.283 [-3.293,+1.973] |

### Interpretation after the falsification check

- **EQR:** downstream presence (+2.170%) and upstream presence (+1.920%) both improve MSE descriptively; the up-versus-down difference (-0.249pp) bootstrap crosses zero. Downstream N report intensity (+3.497%) numerically exceeds upstream (-0.009%), but the direction difference -3.506pp [-9.153,+2.126] crosses zero. No demonstrated directional specificity; environmental or source-reporting geography may help both.
- **Butterfly host-N affinity:** upstream source-positive presence +3.276% while downstream -0.093%, but up-down difference +3.369pp [-3.152,+9.555] crosses zero. Upstream and downstream N intensity gains +0.530% and +0.560% are essentially indistinguishable. This H7-wide GBIF metric is especially weak for a monitoring-reach downstream neg-control.
- **Plant S50:** upstream presence -0.339% and downstream presence -2.200%, difference +1.861pp [+0.033,+4.886] in posthoc H6 fixed-error resampling. Among 27 tested models this cannot be treated as a confirmatory wastewater effect or evidence of actual plant richness or resource subsidy.
- Original source has strong overlap between upstream and downstream reported positive (35/56 upstream H7 periods and 35/61 downstream periods). Directional source labels and source reporting intensity are not randomized or independent. Source downstream dilution and geographic source density can be correlated. These are predictive **stress tests**, not causal negative controls.

## Numerical integrity and open source gate

- Local source unit tests: **7/7 PASS**. Independent sklearn StandardScaler + OneHotEncoder + Ridge across **27 exact model specifications** agrees with manual weighted H6 ridge; maximum absolute H6 MSE-gain discrepancy 1.79e-12 percentage points. These use the **same original biological records**, not independent ecological field replication.
- Source original GH run [#43679382449](https://github.com/zuizui0223/gesui/actions/runs/43679382449) ensures native station/period original count, exact upstream source equality and original Git blob source SHA. The local output is separately auditable and must be rerun on GitHub before claiming GitHub numerical QA completion.
- Concurrent PR91 source preflight: official 2026 EEA WISE6 same-station year-2018/2022 measured TN/TP has NOT been found in the finite tested pilot; original historical Water Quality records mainly pre-2013 for the 7 inspected actual reported source stations. A full 2026 Part2 publisher download has NOT been obtained and its full population cannot be called missing.
- No causal wastewater -> water chemistry -> algae/benthos/plant/host/larva/pooled butterfly -> predator pathway is estimated. No measured butterfly adult puddling or reproductive outcomes. Preserve previous positive and null claims without optimistic sign selection.

## Why this matters for next research step

Prioritize exact observed 2018/2022 receiver chemistry and verified temporal source availability; obtain reach-specific effluent flow fraction and upstream/downstream plant operation. Until then, a structural causal SEM cannot distinguish resource subsidies from sediment stress, urban land use, source reporting, or observation bias. The direction placebo is a diagnostic, not a substitute for true water chemical exposure.
