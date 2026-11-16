# 2026-11-16 Austria68: upstream WWTP reported P → real 2022 river TP → benthic EQR

Status: source verified and independent numeric 4/4 QA PASS. Exploratory within-water observational association, not causal mediation and not a water–land food-web response.

## Authentic sources and fixed geography

Original Austrian BMLUK 2022 national FW1 annual TP, genuinely unfiltered Gesamtphosphor-P in mg P/L; original FW1 XLSX SHA256 17de22544e9d043f83a928154847fb3fcd938894ce7007cfd06c000444b6eeb2. EEA2022 observed benthic EQR official source gitblob SHA1 325f199e985d9862ffee5b0edb205c61a7be3833, all 68 Austria same official station IDs with 2022 stream TP; true five spatial 3-degree blocks with n=9/15/21/13/10. Existing directed river 50km source mass: 53/68 WWTP annual reported P/N positive UPSTREAM; 54/68 DOWNSTREAM reported P positive. Both-positive reporters only; zero means no matched eligible reported source, NOT no sewage. Upstream river connectivity comes from version-pinned third-party HydroRIVERS.

## Three main observational relationships

First original signed main source run: https://github.com/zuizui0223/gesui/actions/runs/43882137861 . Five genuine source-region leave-one-out Ridge10, background log river Q, log area, lat/lon; train-only numeric standardization.

| Predictor → original response | Five region-equal MSE improvement | Region count improving |
|---|---:|---:|
| UPSTREAM 2022 reported P tonnes/yr → real stream lnTP | **+5.2917%** | 4/5 |
| DOWNSTREAM reported P tonnes/yr placebo → stream lnTP | +1.0646% | 4/5 |
| UPSTREAM annual reported P normalized by same original Q → stream lnTP | **+16.3021%** | 5/5 |
| real 2022 stream lnTP → same-site benthic EQR | **+33.8964%** | 5/5 |
| UPSTREAM raw reported annual P → EQR | +0.4671% | 3/5 |
| DOWNSTREAM raw reported annual P placebo → EQR | **+13.2950%** | 4/5 |
| real stream lnTP ALREADY controlled, then UPSTREAM raw annual P → EQR | −2.0963% | 1/5 |
| UPSTREAM P/Q potential concentration → EQR | **+19.2691%** | 5/5 |
| real stream lnTP controlled, then UPSTREAM P/Q → EQR | +6.7216% | 4/5 |

Raw Spearman annual upstream reported P versus observed TP rho+0.20986, downstream reported P versus TP rho−0.06191. Real stream TP versus EQR rho−0.58621. Full covariate-adjusted coefficient for source log1p upstream P predicting actual stream lnTP β+0.59929 per one source exposure SD, not a physical concentration unit. Separate France same 2022 river TP/I2M2 association is an independent water-replicating population but not a WWTP experiment.

## Symmetric flow-adjusted downstream negative-control adversarial source test

Original [Actions #43883032442 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43883032442): 4 core results independently reproduced via augmented-lstsq Ridge instead of normal-equation solver. Source P potential concentration proxy factor = 1e9/(365*24*3600*1000) ≈ 0.03170979 mg P/L for 1 tonne P/year per 1 m³/s, **FULL delivery** assumed and no actual measured P from WWTP separated from other sources. Apply identical annual P/Q normalization to UPSTREAM and DOWNSTREAM reported source.

| Predictor (both normalized by source Q) | Heldout stream lnTP gain | Heldout EQR gain |
|---|---:|---:|
| UPSTREAM P/Q source | **+16.3021%** (5/5) | **+19.2691%** (5/5) |
| DOWNSTREAM P/Q negative-control | **−45.9936%** (1/5) | **−16.4225%** (1/5) |

The downstream-source proxy is unstable: lnTP prediction on block 16_5 worsened ~236%, so its large negative gain is NOT proof of causal upstream source identity. Adding UPSTREAM P/Q after actual stream lnTP improves EQR prediction +6.7216% (4/5), but simultaneously adjusting upstream and downstream normalized source components makes models unstable. Fixed five-block holdout-error exhaustive resampling intervals (NOT independent sampling CI): UP P/Q→stream lnTP +7.61 to+27.59%, UP P/Q→EQR +14.47 to+24.40%. Trimming original source P/Q top3 leverage sites by predictor alone (posthoc, n65) gives UP stream TP +26.14%, EQR +18.85%, DOWN stream TP −58.31%, EQR −30.34%. Do not replace n68 primary with trimmed sensitivity.

## Biological and causal stop boundaries

Actual national government 2022 annual stream TP is ALL-SOURCE aquatic nutrient concentration, not quantified wastewater-derived P. River Q annual average does not measure seasonal dilution or unique discharge event. Reported annual source loads only cover select EEA UWWTD reporters; monitoring location, other agricultural and urban P, habitat/sediment/temperature, multi-outfall urban network effects and source selection confound ecological identification. Within-water same-site chemical and EQR correlations are robust geographically, but the five blocks are insufficient for ordinary cluster-based confirmatory inference. Source 2022 annual chemistry and EQR same year are not same-day samples, and EQR is ecological quality, not biodiversity or invertebrate biomass.

The original Austrian reported-positive outfall-selected cohort of 117 facilities near national stations and FW1 2022 yielded only 3 original near-outlet TP sites; of 2 with same-2022 benthic observation, BOTH river reaches are hydrologically UPSTREAM of their corresponding outfall, so no local causally attributable downstream test exists in that separate source population. Original 68 monitored sites have 53 regional upstream reported sources by modeled 50km network but NOT true effluent nutrient delivery at their specific plant/butterfly sites.

**Single paper conclusion:** possible aquatic leg of a watershed nutrient and biological-response mechanism now has actual colocated wastewater source proxy + measured water TP + benthic condition in the same 2022 year. An aquatic-to-terrestrial WWTP ecology causal paper is still NO-GO without same physically localized 2022 riparian plant biomass/greenness and actual receiving nutrient source/clean reference. Do not pool separate European 451 H7 plant NUTRI–butterfly HNS affinity with 68 Austria river stations as if an observed causal trophic chain.

Source code:
- analysis/20261116_austria_wisa_official_source/AT68_2022_real_WWTP_P_to_ambient_TP_to_benthic_EQR.py
- analysis/20261116_austria_wisa_official_source/AT68_equal_Q_updown_mass_adversarial_qa.py

Both source verified Actions artifacts retained. main/chocho unchanged; Draft PR #117 not merged.
