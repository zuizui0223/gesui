# Actual 2018/2022 reported upstream WWTP N/P sources: fixed watershed×period integrity gate

**2026-11-12 JST. SOURCE AUDIT BEFORE NEW ORIGINAL REPORTED-MASS × CROSS-TAXON EFFECT MODELS.** The recent 454 basin-period static HydroWASTE spatial EQR negative signal is already known. This analysis is explicitly **posthoc** and must not be presented as independent prospective corroboration.

## Why the new source is different

Original European EEA 2018 and 2022 **reported treatment-plant nutrient mass in tons/year**, already linked to its original monitored downstream benthic EQR station through HydroRIVERS upstream direction and ≤50 km source topology in PR #17/#14:
- 2018 original station CSV: `2018_measured_upstream_P_t_per_year`, `2018_measured_upstream_N_t_per_year`, `2018_upstream_qualifying_WWTP_count`, `DIS_AV_CMS`, original station key and BQE procedure.
- 2022 original station CSV: `upstream_2022_EEA_P_measured_t_per_yr`, `upstream_2022_EEA_N_measured_t_per_yr`, `DIS_AV_CMS`. No dedicated source-count field; **do not invent a comparable source count** for 2022.
- Previously original source inventory showed many stations with **no positive matched reported sources**. A numeric zero here indicates no *matched reported upstream load in the source inventory*, NOT zero real sewage or N/P concentration.
- These original EEA reported masses are not measured chemical TN/TP concentration at the receiver, and use source-year-specific national reporting practices. N and P highly correlated and reporting geography may change between years; no causal WWTP treatment interpretation.

Analytical geography fixed: same original **454 HydroBASINS level7 × period** rows / 232 unique basins, 227 early and 227 late; genuine source-fixed original plant S50, five-family Lycaenidae and EQR rows in [Actions #43215754868](https://github.com/zuizui0223/gesui/actions/runs/43215754868). Additional original site→basin manifest from [Actions #43165326244](https://github.com/zuizui0223/gesui/actions/runs/43165326244), SHA256 `a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89`. Do not change sample by outcome or convenience.

## Source-only Stage 0

Re-read real EEA 2018/2022 original source files, pin source Git blob SHA, decode original CSV, **do NOT open the EQR numeric columns** during this audit. For each of the original 454 H7-period units from the saved previous original three-outcome row artifact, gather *only* native original station IDs per proper period; audit:
- station counts and native H7 assignment one-to-one;
- distribution of matched reported-source positive stations and basins (`P>0 OR N>0`);
- within-station original modeled river DIS_AV_CMS valid >0, and reported P/N valid ≥0;
- per-basin year mean of site-level reported upstream annual P and N t/year; label values source-reported-only and **do not** treat 0 as actual concentration or true nonemission;
- source-pressure per monitoring station scaled by modeled river discharge: `(source_t_per_year / river_DIS_AV_CMS) * 10^6 / (365.25*24*3600)` yields **mg/L-equivalent from full upstream annual reported source mass assuming full delivery/no retention**, NOT actual observed river nutrient concentration; use site-equal mean after unit conversion and check raw P/N source. The factor is fixed independent of source values and cannot be changed to chase signs.
- source-positive row count by both years, countries, and number of native H7 basins. No EQR/plant/butterfly numerical association or sign is computed at this stage.
- Original EEA upstream matching is bounded to reported sources within 50km and source missing for many stations; a 0 means **no positive matched reported load**, which is not a genuine controls group without observation-reporting coverage adjustment.

**Predeclared research support gate for any inferential original 454 unit study:** at least **30 H7-period rows with positive matched actual reported upstream P or N loads**, at least **20 unique source-positive Level7 basins**, source-positive coverage both time bins with ≥10 H7 period rows each, and at least **3 countries** with positive original H7 load. If fail, STOP biological response modelling; still preserve source descriptive report. This gate is fixed before looking at *this new source distribution*; prior broad station source sparsity is acknowledged. Do not rescue by aggregating extra river countries, redefining time or switching response thresholds posthoc.

## Stage 1 only if source gate passes

On EXACT same original 454 basin-period rows, compare model B (country + original BQE method, period, log area, site observation n, log modeled river flow and GBIF opportunity covariates) with B+P (reported upstream P river-discharge-scaled), B+N (reported upstream N scaled), separately for original measured EEA EQR, GBIF plant S50, and five-family Lycaenidae log odds. Do not interpret an N vs P difference as mechanistic independent nutrient contribution: high source collinearity.

Use original H7-basin-hash 5fold Ridge10 with both periods per basin same fold; descriptive H7 basin-resample 95 intervals and leave-France-out/predeclared country-weight sensitivity. Independently validate numerical coefficients if successful; otherwise preserve untested status. No direct pollination, insect emergence or predator assimilation on H7 sites.

Source count/reporting confounds matter: a positive source is more likely in countries with good WWTP reporting, urbanization, agricultural inputs and hydrological station coverage. An improved heldout EQR prediction would be *additional reported source geography*, not identified sewage harm or ecosystem services.

Keep this as a new unmerged Draft. Do not alter original PR #65, PR #64, PR #63, `gesui/main`, `chocho`, or use eBird.
