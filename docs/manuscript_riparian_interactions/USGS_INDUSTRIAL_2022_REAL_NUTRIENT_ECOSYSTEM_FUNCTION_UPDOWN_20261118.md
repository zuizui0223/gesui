# 2026-11-18｜独立した産業排水源：原本2022年の栄養塩・有機物分解・酸素変化率を同一河川上下流で実測比較

**New original-source numeric reanalysis. Exploratory, one industrial effluent-receiving stream, NOT an independent causal sewage–ecosystem-service trial, NOT a municipal WWTP study, and NOT a demonstrated stream→riparian nutrient flux.** Original publication Kotalik et al. 2025 has already examined this source and 2026 Perrotta et al. reused the same outfall.

## Original verified public source, actual observation grain

Original CC0 USGS data: [DOI 10.5066/P13BSPG4](https://doi.org/10.5066/P13BSPG4), ScienceBase item `68af0ddcd4be0253bb2c7c6b`, 16 original supplementary files. The original analysis is [Kotalik et al. 2025](https://doi.org/10.1021/acs.est.5c04867) and the associated related 2026 microbiome paper [Perrotta et al.](https://doi.org/10.1021/acsesttox.6c00051). **Do not double-count the same food-processing outfall as two independent studies.**

[Original source all files #44220860908 PASS](https://github.com/zuizui0223/gesui/actions/runs/44220860908) with 2021 PFAS assay n1,268, 2022 PFAS assay n2,839, plus stream-function original tables. The original PFAS source names `Site,Sample_datetime,Sample_type,Lifestage,Parameter,CAS,Result_value,Remark,Reporting_level,Unit,Report_level_code,Method` were tab-delimited **.txt**, checked in the corrected source program. Most chemical assay rows include multiple analytes per biological source sample and many `BDL`/nondetect entries (2022 `Result_value` empty in 2,267/2,839 lines): do not count rows as true independent animal specimens, and do not code BDL as observed zero.

[Original source unit/month/class fixed #44243893723 PASS](https://github.com/zuizui0223/gesui/actions/runs/44243893723) on 2021+2022 PFAS source and 2022 nutrients/source ecosystem functions. Original source table SHA256:
- `SoyoilFoodWeb_Table13_Nutrients.txt`: `eafe63566c2f071dc42d21eac438ee56f1887f7f10f1900f8356e0661f3e95f8` (81 rows = three sites × three dates × three named analytes × three original replicate slots, some missing).
- `SoyoilFoodWeb_Table14_Decomposition.txt`: `058fa757dab0e2b55d41f197595d6cc8dfd0b29fd05ca82c50d68ef1411cfa8e` (12 original individual bags, upstream/downstream × two dates × 3 bags).
- `SoyoilFoodWeb_Table15_Metabolism.txt`: `83aac8818151ce265ce1330ee391038d6832f99a28d8787afb5ebb2343448395` (36 individual oxygen assay rows = upstream/downstream × two dates × 3 bags × 3 technical readings).
- Field data `SoyoilFoodWeb_Table11_Fieldparameters.txt`: 18 rows, original same 2021/22 survey visits. Original source source catalogue is one actual receiving river system; three sampling locations `Upstream, Outfall, Downstream` are **not three independent streams or municipal WWTP replicates**.

## Genuine original results: three same-year summer visits in 2022

[Original numeric reproducible GitHub Actions #44244436328 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/44244436328). [Reproducible code](https://github.com/zuizui0223/gesui/blob/docs/wastewater-riparian-interactions-outline/analysis/20261118_industrial_PFAS_source/kotalik2022_original_three_nutrients_decomposition_metabolism_descriptive.py) uses genuine original-source checksum-locked values and keeps location/day/assay units and original technical replication. No invented treatment replication; no confirmatory p-values or CIs.

All three nutrient classes **NH4-N, NO3-N, SRP** shown, without posthoc selecting just positive endpoint.

| Source measured analyte (original unit) | July 5 downstream−upstream | July 19 downstream−upstream | Aug 2 downstream−upstream |
|---|---:|---:|---:|
| NH4-N (micrograms N/L) | **−34.00** | **+33.77** | **+40.25** |
| NO3-N (milligrams N/L) | **+0.900** | **+0.133** | **+0.600** |
| SRP (micrograms P/L) | **+56.80** | **−1.67** | **−15.68** |

Real original **outfall** NO3-N: 22.5, 7.4333, 4.0 mg N/L on same three days, versus upstream 2.0667, 1.6667, 0.50 and downstream 2.9667, 1.80, 1.10. Outfall SRP: 1,906.73, 1,308.0, 580.8 micrograms P/L; this was much higher than either up/down source river level all days. Thus observed outfall N/P concentration increases do **not** imply the downstream river's SRP concentration is higher on all three days; background, flow, transformations and other inputs are not identified.

### Original decomposition assay: observed fraction of initial bag mass lost

| Original date | Upstream (3 genuine bags) | Downstream (3 genuine bags) | Down−up, fraction and percentage points |
|---|---:|---:|---:|
| 2022-07-19 | **0.46009** (46.0%) | **0.42403** (42.4%) | **−0.03606 (−3.61 pp)** |
| 2022-08-02 | **0.72725** (72.7%) | **0.51764** (51.8%) | **−0.20962 (−20.96 pp)** |

Two survey dates in one source stream only. Original input `Initial_mass,Final_mass` ratio; we did not compute annual litter production, annual ecosystem decomposition services, animal resource export or true river-to-land transfer. Original potential variability among 3 bags/date remains.

### Original signed oxygen-rate measurement (NOT direct ecosystem productivity)

Predefined original reported field `Oxygen_time_mass` averaged over **three technical readings within each of 3 independent bag IDs**, then bag means averaged per date×site. Negative values retained exactly as reported.

| Original date | Upstream original oxygen-time/mass | Downstream | Signed downstream−upstream |
|---|---:|---:|---:|
| 2022-07-14 | **−2.07956** | **−1.51800** | **+0.56156** |
| 2022-08-02 | **−1.36733** | **−0.84922** | **+0.51811** |

The signed oxygen-time variable is a specific original method-dependent assay; **not** a direct stream gross primary production (GPP), net primary production, carbon mass/energy flux, or pollutant export. A positive downstream-minus-upstream value is a *less negative signed rate*, not proof that aquatic biomass production improved.

## Scientific interpretation and next gate

1. **Nitrogen versus phosphorus differ at the receiver.** NH4-N down−up changes sign, NO3-N is positive for all 3 surveys, SRP only positive for 1/3 even though the discharge source is strongly enriched. Even source-confirmed chemical concentration is not a universal same-direction receiver response.
2. **Chemical loads and function are different response axes.** The same source study has downstream-lower bag mass-loss fraction on both 2022 paired dates, while reported normalized oxygen-rate field also differs. That is a localized *descriptive* contrast, not a causal estimate of wastewater impact or positive ecosystem service subsidy.
3. **There is real aquatic→riparian material exposure in the published 2025 source** (PFAS in adult insects and spiders), but the true chemical PFAS stages in 2022 were predominantly sampled **Aug 1–2** on the up/down river, and function visits were Jul 14/19 and Aug 2. A single late-summer near-date does not support a 2021–22 replicated fully synchronized material→function→trophic transfer SEM. The 2021 source has only one matched spider/insect 2021 field visit.
4. **One actual industrial effluent source is one spatial treatment.** Even if multiple dates/specimens have the same sign, do not treat them as independent watersheds in the user’s Europe-wide sewage-impact macro project. Multiple observational studies of this identical outfall including Perrotta (2026) do not create new independent sources.

The next source-only test, if scientifically needed, is whether original `PFOS` and `6:2 FTS` are actually detected above their original reporting-level/censor codes in the **same original 2022 Aug1–2 event** for effluent/river water, aquatic Trichoptera larvae/adults and riparian Tetragnathidae spiders. Maintain original `ng/L` versus original `ng/g wet or dry` without mixing units. Do not infer direct physical ingestion from concentration ratios and do not fit a new food-web hypothesis on one outfall after inspecting signs. The effect of N/P discharge on broader terrestrial vegetation or butterflies remains unidentifiable from this industrial-source case.
