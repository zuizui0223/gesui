# One-scale macroecology: watershed × period master analytical unit (frozen 2026-11-12)

**Replace serial literature case studies** with a source-first **global macroecological sampling frame** for reported WWTP emissions, actual in-river biology, flora/plant-host resources and butterflies.

## Main geography: HydroBASINS Level 7 by Pfafstetter basin ID

Source: [HydroBASINS v1c](https://www.hydrosheds.org/products/hydrobasins) / [HydroATLAS](https://www.hydrosheds.org/hydroatlas), global public standardized basin polygons, 12 nested levels and upstream/downstream topology; use **Level 7 sub-basins** as **primary biology geography** and Level 6 as source-first **predeclared sensitivity**. Different basin areas remain different: adjust drainage area, river length, sampled site counts and non-water extent. A basin is a geographical analytic unit, **not** a hydrologically downstream exposed reach without actual directed routing.

**H3 r5 and r4/r6 are ONLY the first source-coverage audit**, not substitute for hydrological connectivity. r5 cells are approximately area-balanced but not exactly equal area; include actual H3 area, never declare every cell identical. r4 (coarser) and r6 (finer) are source coverage diagnostics, not choose-the-positive-effect model search.

## Time bins, a priori

Two nonoverlapping periods: **P1 = 2017–2019** and **P2 = 2020–2022**.

- Original EEA WWTP-source × EQR n2,634 **2018** → P1; 2022 n2,475 → P2. Exact monitored site IDs permit repeated-station audits; national assessment differences prevent blindly pooling EQR across countries.
- Original EPA NRSA **2018–19** n1,221 source-qualified native UIDs / 1,125 physical stations → **P1 only**, with actual TN/TP and benthic EPT fractions and field riparian structure. No fabricated P2 USA observations.
- WPI original **1,055 ports**, 42 countries, **2015 modeled wastewater generation** and GBIF photos 2011–24 form a DIFFERENT COASTAL LOCATION SAMPLE FRAME. Do NOT treat WPI generated wastewater as 2018/2022 measured WWTP effluent or the 2011–24 facets as P1/P2 records. The current spatial overlap gate includes only original port **coordinates** to quantify geography incompatibility; no biological model, no sewage effect from those ports.
- Future genuine point source: globally georeferenced [HydroWASTE](https://www.hydrosheds.org/products/hydrowaste), **58,502 registered treatment plants** with modelled location, outfall, effluent flow and treatment. HydroWASTE point characteristics are mixed measured/imputed and essentially a **static pressure index**, not directly N/P tonnes per year. Europe's 2018/2022 actual reporting can be retained as separate stronger source measurement. Follow directed reach/HYDRORIVERS to calculate upstream exposed load/dilution; don't sum arbitrary WWTPs within 25 km as delivered dose.
- Future biological layers per **same HydroBASINS 7 × period**: GBIF raw **accepted plant species occurrences**, candidate butterfly host plant incidence, Lycaenidae and other butterfly occurrences in same polygon/period, and standardized monitoring where available (eBMS). **Reissue original year+polygon species queries**; source port 5km species facets cannot be assigned to an arbitrary basin or years. Rarefied richness is interpretable only after source effort/coverage gates. GBIF recording is not true plant abundance or pollination service.
- The 17 Oregon LTER adult emergence streams / Canadian 8 WWTP plastic→spider / Ontario bird sampling are **external functional calibration/validation**, not 17+8+29 invented national macro basin observations.

## Row-level master table — one row per (HYBAS_ID, PERIOD), no unmatched-join inflation

| Field | Grain / meaning |
|---|---|
| `HYBAS_ID`, `PERIOD`, `country`, `area_km2`, `river_km`, `hydroATLAS_env` | Full coverage of basin and period, static attributes marked as static |
| `WWTP_static_effluent_flow`, `N/P_upstream_tpy` | Separate estimated plant discharge pressure from measured upstream N/P loads; missing explicit |
| `river_EQR_n`, `river_EQR_site_mean`, `river_EQR_method`, `river_EQR_source_country` | Europe only; methods harmonized within country/assessment method; not directly equivalent to USA EPT |
| `EPT_frac_n`, `EPT_frac_site_mean`, `river_TN_TP`, `riparian_structure` | US only in P1, each site/UID weighted once as defined, visit repeats clustered |
| `plant_occ_record_opportunity`, `plant_rarefied_species`, `host_occ_rarefied` | Independent basin × period GBIF raw occurrences, same source completeness rule for all cells |
| `butterfly_obs_days_or_event_effort`, `lycaenidae_relative_detection`, `butterfly_family_mix` | Field survey if available; GBIF presence-only relative recording with controls, not butterfly population |
| `effort_quality_flags`, `source_SHA`, `source_method` | No value is filled with zero for a missing taxa-source API response |

## Comparable outcomes: common predictor and common support, NOT fake equal biological units

On identical jointly observed basin-period rows, compare *within-layer* standardized exposure slopes, using **layer-specific observation models**, not EQR≈EPT fraction≈plant richness≈Lycaenidae counts. Optionally define source-standardized within method Z, but retain response families separated and use a hierarchical multivariate model only when the SAME basin-period sample has support.

- First-stage macro estimand = multivariate environmental association with `log1p(upstream reported actual WWTP NP delivered dose / flow)` **where validated** or `log1p(static HydroWASTE effluent / streamflow)` as *separate source predictor*, adjusting cultivation, urban, flow, elevation, river length, recording effort, country and year.
- Benthic quality: country/method-specific ordinal standardized EQR or local EPT trait fraction as **different responses**; no direct slope comparison without link model.
- Plants: effort-rarefied host/flowering-resource richness and composition *not actual flowering or pollen transfer*. Butterfly family detection/composition with observer effort and non-Lepidoptera background.
- `WWTP × habitat × taxon layer` interaction asks which response layer associates differently with the same environmental gradient, **not serial causal mediation**. EQR may respond negatively while host resources/recording opportunity appear positive. No requirement that a positive endpoint be found.
- Service quantity (adult emergence biomass, realized predation, plant seed set) remains NA in unmeasured macro basins and MUST NOT be filled with predictions from unrelated HJA/Canada site coefficients.

## Stage 0 source geometry quality gate (the task implemented first)

Exact original **EEA 2018/2022 site CSVs** and **US EPA NRSA frozen pre-outcome UID site geography** plus original **1,055 port source point metadata** are used to compute:
- one-to-one geography assignment at H3 r5 and r4/r6, unique cells and IDs per original data source and temporal bin,
- early vs late original EEA **same river monitoring site** and **same H3 cell** overlaps,
- global port point same H3 cell as EEA/US river station **GEOGRAPHIC ONLY**; port 5km GBIF flora is NOT mapped to that cell and generated 2015 water is not routed to the station,
- per-source country/state counts, original URL/blob SHA, null geolocation and H3 area ranges,
- strict exclusion of noncoordinate and out-of-domain values; retain all original rows as source counts, mark missing instead of zero.

Stop if <95% valid location for any original source, original row counts drift, or source blob SHA differs. **No ecological coefficient calculated**, no erosion of originally positive results.

## Stage 1 (after this source audit)

Acquire actual HydroBASINS Level7 (North America + Europe first; global other continents next), original HydroWASTE point version, and compute crosswalk `original_station_ID → H7 basin`; save multi-year source coverage. For GBIF plant/butterfly basin×period requery, predeclare minimum independent taxa observations/effort on each native basin and country support before inspecting associations. Use one H7 basin-time table with all originally eligible missingness flags; report independent temporal replication and region-block holdout. Do not increase basin level or search radius in order to rescue a favorable slope.

**All old PRs, gesui/main, chocho unchanged. No eBird. No new causal/ecosystem-service claim until genuine same-unit flow and functional data align.**
