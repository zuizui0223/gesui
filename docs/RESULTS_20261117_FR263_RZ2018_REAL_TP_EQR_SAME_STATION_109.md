# 2026-11-17 | Genuine SAME 109 French river stations: 2018 measured TP and benthic I2M2 × external 2018 Copernicus riverbank LCLU

**Exploratory, posthoc source-checked numerical report. No causal sewage→aquatic→terrestrial transfer, no 2018–22 bank vegetation productivity, no butterfly feeding or biomass measured.**

## Identity and pre-outcome source selection
- From exactly 263 original **same physical 2018/2022 French WQ TP plus EEA I2M2 EQR station IDs**, 73 H7/24 H6; 2018 and 2022 observed TP source [GitHub Actions #43373639391](https://github.com/zuizui0223/gesui/actions/runs/43373639391). Original exact 263 2018/2022 dates source SHA256 `2190ac905326f4f2f4f278eb0debb4a7bfda7c4c71347d4941b54cc19d3de264`.
- Independent 2018 **real mapped** Copernicus Riparian Zones Level1 terrestrial class area at the very same 263 EEA station IDs, from source-original 2018 RZ full source [#44004357868](https://github.com/zuizui0223/gesui/actions/runs/44004357868). Actual source RZ covers both 0–250m circle and 250–500m surrounding annulus at >=70% **109/263**. Both rings have mapped non-water terrestrial area >=10% at **109/263**. Exactly these *source-selected* sites contain **49 original H7, 20 original H6, 8 original 3-degree geographic blocks**. All source-only gates fixed *without looking at numeric TP/EQR*: >=80 sites, >=30 H7, >=12 H6, >=5 geographic blocks, passed.
- **Only 1/109** of these genuine same-site France river points has original eligible *reported upstream 2018 WWTP P mass-positive*. Therefore **NO sewage-specific exposure or upstream vs clean control inference** is possible within this panel.
- Numeric originally signed source authentic input CSV `FR263_SAME_PHYSICAL_STATION_2018_2022_REAL_TP_AND_EEA_EQR.csv` SHA256 `1824b0d9b84d35a460b4d1d0513659f93a16c59d7f5666a2d8707be1c6da9942`. Same prior original H7 geography controls original CSV SHA256 `f2ff4e56f245e1aa8f84f5ac099d68ec1c7e27a786c1067439db388b186dd90b`.
- Cohort loss matters: selected 109 original 2018 TP median=0.083 mg/L vs unselected 154=0.0975; EQR median selected 0.72658 vs unselected 0.70235. These are **descriptive selection diagnostics**; coverage is not random.

## Outcomes and fixed method
Mapped *terrestrial* natural-cover fraction per ring = 2018 RZ `CODE_1_18` 3 woodland, 4 grassland, 5 scrub/heath, 7 wetland, divided by mapped **non-water** land. Primary terrestrial outcome `DeltaNat = nearNatural_0-250m − farNatural_250-500m`. On the 109 field points near mean=**0.57989**, far mean=**0.58233**, near−far mean **−0.00244**. This is not greenness, foliage nitrogen, productivity, or a near-actual-riverline buffer (circle is centered on station).
Original genuine **2018 measured ambient river TP** log mg/L and **2018 measured I2M2 EQR**; 2022 TP only exploratory stability check, NOT 2022 bank response.

All models **on exactly same 109 stations**; original 5 source-fold IDs are **0..4**, each H6 in exactly one fold. H7-equal weight per station at fitting (1 / source stations in H7), Ridge α=10, scaling from training only. Baseline predictors held fixed: original site log river modeled Q; source H7 log area, outlet modeled flow, upstream cropland and urban; same-site far 250–500m natural, urban, cropland land composition. Group-level out-of-fold error aggregation equal H6 (n20); 2,000 fixed-H6 error resamples are *descriptive* 95% ranges, **not independent trial p-values**.

| Test, one addition to the same prespecified baseline | Original parent-H6 equal MSE improvement | Fixed 20-H6 OOF error resampling descriptive 95 |
|---|---:|---:|
| add genuine same-year 2018 logTP to predict **DeltaNat** | **−0.776%** | [−5.686,+3.110]% |
| add **DeltaNat** to predict 2018 river I2M2 EQR, controlling 2018 real TP | **−4.181%** | [−11.284,+0.707]% |
| add **DeltaNat** to predict 2018 real logTP, controlling original upstream landscape | **−1.020%** | [−6.001,+4.225]% |
| temporal stability only: add 2022 measured logTP to predict the **same static 2018 DeltaNat** | **+2.112%** | [−3.133,+6.365]% |

Raw 109-site Spearman: **2018 logTP–EQR ρ=−0.5651**, confirming water-side co-gradient; **2018 logTP–DeltaNat ρ=+0.0601**; same fixed habitat vs 2022 logTP ρ=+0.0620. These are observational; the 2018/2022 comparisons on same sites are not independent years of habitat measurement.

## Interpretation
A genuine same-physical-station water chemistry–benthic condition *spatial* association exists within selected France 2018 cohort. But this one-time riverbank *mapped class fraction* does **not** add transferable predictive signal for local water TP, bottom fauna EQR, or vice versa at the fixed observation scale. This is a **positive source-integration result and negative local terrestrial class-composition diagnostic**; it does not prove absence of hydrologic water→land subsidy. Avoid switching 250/500m, natural-cover codes, country, year, H6 fold, or site-selection criteria to chase an affirmative coefficient.

**Next distinct physically observable endpoint**: fixed 2018+2022 original true satellite near–far vegetation productivity/NDVI; prior 24-H6 STAC metadata-only pilot now source PASS 24/24, but original 10m clear land pixels and authentic riverline-based distance not yet validated. Original official HydroRIVERS ZIP returned HTTP403 (Actions #44008867522); an independently source-SHA locked **third-party converted HydroRIVERS v1 GeoParquet** was used earlier in PR115 and can be separately verified without masking the failed official ZIP download.

## Audit links
- New local real-source SHA-locked script `fr263_RZ2018_realTP_EQR_same_station.py` and JSON `gesui_FR263_RZ2018_TP_EQR_same_station_20261117.json` generated in 2026-11-17 conversation; SHA of all input bytes checked before fitting. PR #120 remains Draft; this report is not new underlying satellite collection, not a submitted paper, and not `main`.
- Original 2018 RZ source PR #123, [#44004357868](https://github.com/zuizui0223/gesui/actions/runs/44004357868). Original FR263 actual water source PR #71, [#43373639391](https://github.com/zuizui0223/gesui/actions/runs/43373639391).
