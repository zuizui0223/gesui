# 2026-11-16 Austrian BMLUK 2022 real TP ↔ EEA benthic EQR: country replication and hard wastewater-flow gate

**Exploratory same-station aquatic observational analysis.** Official original national 2020–2022 yearbook source and original EEA source; original measurements directly paired by exact government PK = EEA PKINT without nearest-site imputation. National result independently complements the France 2018/22 true river TP↔benthic EQR findings but **DOES NOT establish sewage origin or cross-ecosystem transfer**. Do not pool as an experimental wastewater treatment effect.

## 1. New original 2022 government chemical source discovered

- [BMLUK genuine Wassergüte Jahresbericht 2020–2022 / official results ZIP](https://www.bmluk.gv.at/themen/wasser/wasserqualitaet/wassergute-jahresbericht-2020-2022.html), signed exact [Actions #43872894930 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43872894930). Published original ZIP size **601,117 bytes**, SHA256 **6a0db0d0b8360c11f580333e5eabc3479c68fe8a5d623b87ef61f0185e5bd31a**. Original file FW1_Fließgewässer_Chemie_2020_2022.xlsx 116,777 bytes, SHA256 **17de22544e9d043f83a928154847fb3fcd938894ce7007cfd06c000444b6eeb2**.
- Actual sheet FW_Chemie_Mittelwerte_2022 (A1:W103): 100 original unique FW government river site IDs with 2022 annual means and percentile statistics; genuine **Gesamtphosphor-P (unfiltriert) [mg/l] - Mittelwert** field index 18 (Excel R). All 100 numeric TP annual means; reported 2022 nutrient concentration is annual stream P, not actual upstream WWTP discharge apportioned phosphorus. The workbook does NOT contain actual total nitrogen TN. Other variables include nitrate-N 90th percentile, ammonium-N mean, orthophosphate-P 90th percentile, which must **NOT** be equated to actual total N or actual total P annual mean.
- FW2 in same 2020–22 original source is operational monitoring list **2019–2021**, not same-year 2022 benthic data. Historical main EEA2022 Austrian true benthic site source Git blob **325f199e985d9862ffee5b0edb205c61a7be3833**, 68 original Austrian station IDs and all MI-RW-AT original classification, normalized EQR source measured 2022. Exactly **68/68 official original BMLUK FW1.1 station IDs** match those same EEA2022 site PKINT IDs; no inferred coordinate or proxy join for this Austria ecological comparison. 68 points in original 5 region_3deg blocks (15_3:9,15_4:15,15_5:21,16_4:13,16_5:10).

## 2. Same Austrian 68 actual stations: independent, spatially heldout river ecology

Analysis method frozen in [protocol](https://github.com/zuizui0223/gesui/blob/analysis/gesui-20261116-wisa-gzuev-water-benthos-2022-source-access/docs/PROTOCOL_20261116_BMLUK_FW11_2022_TP_EEA2022_BENTHIC_AT68_GEOBLOCK.md) before inspecting all 68 numeric TP→EQR pairs. Primary predictor ln(2022 TP mg P/L annual mean), response 2022 normalized benthic EQR. Source controls ln(1+mean discharge), ln(1+upstream catchment area), latitude and longitude. All 68 original observation sites retained; 2022 water concentration mean is not assured to be sampled on the same day as benthic EQR field evaluation. Raw originally published TP range **0.002992–0.431667 mg P/L**, EQR **0.23–0.94**. Exploratory data-driven source discovery, no confirmatory multiple-testing inference.

| Statistic | Result |
|---|---:|
| Original station N; genuine original spatial regions | **68**, **5** |
| Spearman total-P vs EQR raw ρ | **−0.586214** |
| Pearson ln(total-P) vs EQR raw r | **−0.627653** |
| Full model adjusted EQR β per 1 source SD lnTP | **−0.099483** |
| 5 original region blocks heldout, equal block MSE baseline | **0.044610** |
| Same heldout with TP add, equal block MSE | **0.029489** |
| **Geographic-region-out MSE predictive improvement** | **+33.8964%** |
| Regions with positive improvement | **5/5** |

Each original complete region held out entirely, same controls and train-only source scaling with Ridge alpha10 (unpenalized intercept); original block sizes listed in section 1. By source region:

| Heldout native 3-degree block | genuine site count | extra TP MSE gain |
|---|---:|---:|
| 15_3 | 9 | **+29.4003%** |
| 15_4 | 15 | **+25.1984%** |
| 15_5 | 21 | **+44.0245%** |
| 16_4 | 13 | **+39.3786%** |
| 16_5 | 10 | **+36.8692%** |

Full model adjusted lnTP beta stays negative after excluding each independent original region (range **−0.111865 to −0.079804**). Numeric source transfer QA against full genuine 2022 EEA2022 68 records independently reproduced original sum(EQR)=42.74, weighted source-index EQR sum=1506.29, source flow sum=5132.796, area sum=251118.7 and coordinates sums; no manually altered EQR rows. Original heldout 5-fold ridge gain independently reproduced via augmented least squares rather than normal equations (MSE discrepancy < 3e−16). These are *numerical checks on same observations*, not independent field studies.

**Important uncertainty boundary:** There are only FIVE spatial groups. Leave-region-out improvement and negative sign are robust **for this observational source**, but classical cluster robust 95% CIs/p-values would be fragile with five clusters, and the national dataset/source discovery itself followed repeated exploratory analyses. The 68 sites are from one country, one monitoring year, one method. This is independent ecological source evidence relative to French 73 watershed ×2 TP/EQR but not a separate wastewater intervention.

## 3. Strict sewage-source geography: TP+BQE actually downstream of original known measured-load outlets?

Within *original 225 government GZÜV sites* geographically within 2km of *original 117 true annual reported 2022 N+P WWTP outfalls*, exactly **3** are present in original published 2022 FW1 TP annual means:
- Original FW54110117 (Saalach): has EEA2022 benthic field site ID, BUT government-vs-EEA source coordinates disagree **954.614m**.
- Original FW73190907 (Ötztaler Ache): has 2022 TP annual mean but NOT among the 68 original same-year EEA benthic site IDs.
- Original FW73290907 (Ziller): has 2022 TP and EEA benthic original site ID with coordinate difference about 0.5m.

For the two benthic source sites near original WWTP outfalls, authenticated full original HydroRIVERS Europe v1 downstream river graph (source SHA256 5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f) and original reported 2022 outlet true directed river paths verify **BOTH benthic reaches are actually UPSTREAM of corresponding original wastewater receiving outlet reaches**. [Signed true flow graph source audit Actions #43874361627 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43874361627): EEA2022 benthic source at FW54110117 is 2404m geodesic from its outlet and location conflicted; FW73290907 is 570m from its outlet and physically on an upstream reach. Thus **0** genuine source WWTP outlet→downstream actual same-year TP + benthic paired sites in this source-selected frame. Upstream sites can serve as historical reference opportunities but are not actual receiving effluent biological responses, and unknown multiple other inflows must be checked.

This refutes the *interpretation of geographic closeness as exposure*, not the ecological hypothesis that WWTP changes river organisms. National FW1 2022 annual means and original EEA BQE identities cannot substantiate causal wastewater N/P reaching plant leaves, downstream benthos, butterfly feeding, or greenhouse/ecosystem services. Original 2022 positive wastewater outlet P tonnes per year are NOT the measured river TP value.

## 4. Implications for the ONE-paper choice

- **Strong cross-country within-water empirical result:** original France73 2018/2022 genuine river TP→real I2M2 EQR source 263 physical sites, 24 original H6 and outside-parent +24.6% fixed-source Ridge MSE improvement; now independently Austria68 2022 genuine river TP→normalized benthic EQR +33.9% outside geographic block. Different source units, controls and blocks, so *do not average 24.6% and33.9% into a meta-effect*.
- **Strong terrestrial within-layer association (separate population):** original European 451 H7×period plant species-equal EIV NUTRI↔butterfly reference larval-host nutrient-habitat HNS positive, but same-basin measured river TP→plant NUTRI or butterfly HNS is not robust. Both scores share a plant-habitat nutrient source and species-versus-genus composition structure.
- **Water→land sewage mechanism is STILL unproven:** the original outlet cohort has 0 verified downstream 2022 actual TP+benthic at same source, local riparian GZÜV observations and satellite field phenology not sampled simultaneously with original source river TNTP; EEA Part1+Part2 original 225 station2022 exact water rows 0. Austrian national WISA/H2O public actual station-year source needs explicit access or official data steward release.
- **Defensible one-paper now:** basin trait-filtering and compartment-dependent predictability **without falsely identifying wastewater as the common cause** could be written as an exploratory matched-geography synthesis, but a source-specific "wastewater subsidizes aquatic and terrestrial biomass/food webs" causal manuscript is NO-GO. An original-source, single-outfall near-field coordinated water/benthos/riparian measured-data design remains necessary.

## 5. Source and reproducibility provenance

Original source public annual workbook temporarily downloaded in Actions #43872894930 (not committed or raw-file rehosted). Numerical work here done in current-turn scientific Python using artifact_tool to import exact published workbook, stdlib CSV for original BQE site data, NumPy/SciPy for all 68 numeric pairs and true geographic-heldout test; second Ridge solver independently cross-reproduced original coefficients. The derived **aggregate statistics, per-region losses and source integrity** are stored in a separate inspectable JSON result, not all 100 original workbook rows. Public publisher page and official EEA BQE original source are cited in this document. Keep PR117 in draft and do not merge into main or touch chocho.
