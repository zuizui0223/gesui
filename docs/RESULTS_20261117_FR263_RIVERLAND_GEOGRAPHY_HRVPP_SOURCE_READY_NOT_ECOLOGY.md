# 2026-11-17: FR263 water-to-riparian macro source recovery + raster pilot readiness

**State: successful source geometry join; computational local synthetic pixel test; actual HR-VPP data not acquired / no new terrestrial ecology coefficient.** Technical and ecological STOP/GO are separate and old tests are unaltered.

## 1. Input lineage and new independently checked source geometry
- Original 2018+2022 genuine measured TP and physically identical EEA I2M2 benthic sites: n263 station keys, n73 H7 native basins, n24 H6 parents, source Actions #43373639391 and month-matched #43378917089 (signed original artifacts).
- Original repo EEA2022 2,475 rows: `results/eea2022_biology_bqe_20261108/BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv`, Git blob SHA-1 `325f199e985d9862ffee5b0edb205c61a7be3833`.
- Exactly the 263 original EEA native site IDs reselected by exact string (NEVER nearest alternative or favorable ecology value); **263/263 geographical matches**, 263 unique coordinate pairs, **262 distinct original HYRIV_ID reaches**, all finite in France mainland geographic rectangle (latitude 45.02031..50.11889, longitude -3.71043..8.08586). Saved: `results/20261117_FR263_FIXED_ORIGINAL_EEA2022_STATION_COORDINATES_SOURCE_ONLY.csv`.
- Across 34,453 original station pairs, approximate haversine distances: pairs <=1 km **0**, <=2 km **5** (10 involved stations), <=5 km **21** (41 involved), <=10 km **81** (108 involved); nearest-station min 1.15085 km, median 10.80325 km, max 61.64413 km. One true same reach ID `20503419` contains two distinct stations (1.21607 km straight-line apart). These are NOT independent river basins or hydrological along-network distances. Saved: `results/20261117_FR263_GEO_DISTANCE_PRE_RASTER_QA.json`.

## 2. Actual external HR-VPP products independently verified (NO download completed)
- CLMS 10m annual Total Productivity, 2017-present, 2018 and 2022 within listed availability: https://land.copernicus.eu/en/products/vegetation/high-resolution-total-productivity ; Quality Flag: https://land.copernicus.eu/en/products/vegetation/high-resolution-vegetation-phenology-parameters-quality-flag .
- Official manual issue 2.5: https://land.copernicus.eu/en/technical-library/product-user-manual-of-seasonal-trajectories/@@download/file . **VPP-QFLAG values 10 high, 7–9 medium**; exclude <7 for primary. `TPROD` has digital-to-physical factor 0.1 (`PPI × day`); exclude nodata 65535. Product is integrated greenness / phenological index, NOT observed biomass production in carbon units.
- Prefer 10m EPSG:3035 LAEA product for consistent France/2018/2022 zone (official EEA also offers `CLMS_HRVPP_VPP_TPROD_SEASON1_10M_LAEA` and analogous QFLAG). Authenticate/or supply WEkEO or CLMS product data; official CLMS downloads require login/token or alternative provider endpoint. No claim of downloaded rasters, valid pixel counts or successful n263 terrestrial source gate.
- Static candidate cover mask: 2021 ESA WorldCover tree 10, shrub 20, grass 30, wetland 90; NOT pre-2018 land cover and could itself induce change-selection bias. All-vegetation sensitivity must be separate. Copernicus Riparian Zones 2018 only selected rivers, must pass coverage gate if used.

## 3. Fixed geometry / response extraction
- Use authentic original `HYRIV_ID` line geometry, not nearest river search and not a catchment polygon surrogate. Project physical field station and matching reach into EPSG:3035, source gate station-line <=250m (original source matching WQ/EEA sites <=250m; not guaranteed line fidelity). Slice the SAME actual reach at ±500m along-channel from projected station and require >=900m of genuine line before allowing pixel sampling; do not automatically route through confluences.
- Both banks, true same-reach, 10m pixel CENTER-in-polygon mask. Predefined perpendicular river-buffer bands 20–100, 100–250, 250–500, 500–1000m; exclude river itself, farms and built-up in natural-cover primary, and invalid quality/season pixels.
- Station/year/band mean TPROD and valid/source-available pixel counts; primary close-versus-far contrast `mean(TPROD_20-100) - mean(TPROD_500-1000)` for **each** year. Only then compare with measured river TP and separate I2M2 EQR without pretending mediation.
- New proposed data gate (not ecological finding): both 2018+2022 years; >=100 genuine filtered natural-vegetation pixels per band/year/site; >=120 independent station pairs, >=40 H7 and >=15 H6. Report all exclusions; KEEP previous S2 5/10<6/10 STOP as previously failed independent investigation.

## 4. Tested implementation
A stand-alone Python extraction script `gesui_FR263_HRVPP_extract.py` was generated for this conversation and synthetically tested, **not executed on real HR-VPP tiles**. Inputs: 263-source geo CSV, original signed FR263 ZIP, authentic 262-reach vector subset, 2018 and 2022 LAEA TPROD + QFLAG COGs, WorldCover 2021 COG. Source/pixel/QA checks are fail-closed.

Local mathematical smoke: a 10-m EPSG:3035 synthetic raster and exact 1600m synthetic line returned 2560/4800/8000/16000 candidate pixels for four bands, mean 10 PPI-day from raw TPROD 100/10 for each band. Replacing half QFLAG with 'filled'/low reduced valid near-band pixels to 1280/2560; quality selection validated. These are **synthetic code tests, not ecological field or satellite evidence**.

## 5. Next physical-source GO gates, fixed before looking at any satellite outcomes
1. Acquire exact genuine HR-VPP 2018/2022 TPROD and matched QFLAG for all 263 original points and same adjacent bank; capture product versions, nodata, checksum and year-specific data shape.
2. Acquire genuine HydroRIVERS original 262 reach geometries in the same recognized river topology + verify line-to-site gap, river reach segment length and valid natural-vegetation pixel opportunity, all outcomes hidden.
3. Apply the fixed >=120/40/15 and >=100-pixel gates and preserve every failure. If STOP, do not relax on the basis of TP/vegetation signs. If pass, calculate the river-localized TPROD contrast vs real TP/EQR and whole-H6 holdout scores, preserving 2018→22 delta null as a comparison.
4. Austria67 may be only a SEPARATE year-2022 external check with actual upstream reported WWTP P/Q and 5 region blocks, not a combined France wastewater mediation.

**Scientific conclusion now:** Water-side original 263 repeated observation dataset is truly geolocatable down to coordinates and exact river IDs. Whether any river-bank productivity covaries with the water signal remains **untested**. A genuine macro-ecological water–land bridge is feasible to attempt but not yet demonstrated.
