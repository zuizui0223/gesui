# 2026-11-18｜USGS Farmington 5河川：水生幼虫・羽化成虫・河岸クモのPFAS原本に本当に共通する観測単位

**SOURCE-ONLY authentic raw file ID/year audit. No PFAS, isotope or wastewater causal coefficient fitted.** Current manuscript Draft PR #124 stays review-only, no `main` merge.

## 1. Independent published original dataset

USGS 2026 CC0 public data release: [10.5066/P139C2SN](https://doi.org/10.5066/P139C2SN). Original ScienceBase item ID `69b80121b66b010d85582e8b`, original metadata SHA256 `a1ffa182ccc51c5f0018cda1cca99bda1a2407346f3d03fd9b716ae20f832524`. Source five rivers are **five sampling streams inside one Farmington River watershed, NOT five independent wastewater treatment systems**. Original publication [Campbell et al. 2026, Environmental Science & Technology](https://doi.org/10.1021/acs.est.5c16942) has *already* published an ecosystem-scale model of these data and reports compound-dependent bioaccumulation and mostly decreasing `Σ12PFAS` adult-insect→spider concentrations. Do not claim a newly discovered water→spider migration.

GitHub Actions [true source files #44215442112 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/44215442112) verifies seven original files and 4 genuine tabular scientific sources. [Four-source field-stage verification #44216178144 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/44216178144) matches the published raw SHA256 for original chemical, body-mass, isotope and site-location files. [Exact specimen link #44216745073 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/44216745073) verifies real sample IDs with no chemical magnitude fitting.

| Original data file | Real rows | Actual grain |
|---|---:|---|
| `Farmington_PFAS_FWstudy.csv` | **3,192** | **12 chemical analyte rows per assayed sample**, includes lab/field blanks |
| `Farmington_InsectMasses.csv` | 128 | Composite wet mass / individual mass and natural life-stage metadata |
| `Farmington_PFAS_Stable_Isotopes.csv` | 76 | Stable C/N isotope source records across 2022 and 2023 |
| `Farmington_PFAS_Survey.csv` | 2,156 | Basin-wide **water** PFAS sampling; its 22 labelled locations are not 22 river food-web sites |
| `SiteLocations.csv` | 22 | Published geographical places/source site types |
| `pfas_names.csv` | 22 | Analyte nomenclature, NOT 22 uniformly sampled compounds in all original assay matrices |
| Source XML | 1 | Technical original metadata, not a biological specimen table |

Original `FWstudy` has **96 blank-control analyte rows** = 8 QA specimens, removed **before** biological sample counting. Remaining **3,096 analyte rows / 12 compounds = 258 real unique source samples**, with all 5 named stream sites. Their measured `sample_type` and `units` show 5 sites each has original **water, basal biofilm, aquatic larvae, emergent adults and riparian Tetragnathidae**. Chemistry water units **ng/mL**, biota **ng/g wet weight**, so mixing raw values or multiplying ratios by biomass is not an acceptable metric. Some chemical biological `collected` dates were blank: **date missing ≠ same-year assumed**.

## 2. Fixed source 5 physical streams: site and stage support

| Real stream/year of water | Water PFAS source samples | Larval PFAS specimens | Emergent adult PFAS specimens | Riparian spider PFAS specimens | Biofilm PFAS samples |
|---|---:|---:|---:|---:|---:|
| Burr Pond Brook / 2023 | 6 | 16 | 8 | 4 | 6 |
| Hop Brook / 2022 | 6 | 20 | 13 | 6 | 6 |
| Pequabuck River / 2023 | 6 | 11 | 10 | 5 | 6 |
| Ratlum Brook / 2023 | 6 | 15 | 17 | 4 | 6 |
| Russell Brook / 2022 | 6 | 10 | 6 | 4 | 6 |
| **Total** | **30** | **72** | **54** | **23** | **30** |

Original wet mass file calls one biological location **Pequabuck Brook**, whereas USGS field chemistry, isotope, `SiteLocations.csv` name it **Pequabuck River**. Explicit reference alias is source only, not a globally inferred hydrological crosswalk; chemical/isotope individual sample identity must still be independently checked.

## 3. Exact specimen/source-year integration: critical constraint

- Out of **258 true chemical samples vs 76 stable isotope rows**, exact original specimen ID matching using chemistry `field_id/cese_id` and isotope `u_conn_sample_id/usgs_sample_id` succeeds for **10 IDs** only, **all 10 are biofilm**. This is a statement about **same individual specimen** support, not a claim that independently measured trophic groups cannot be compared at a site level.
- No individual **aquatic larva, emergent adult, or spider** isotope specimen ID matched directly to the same original chemical assay specimen. Therefore do not construct false individual records with both `PFAS concentration` and `δ13C/δ15N` using site name alone.
- **3/5 real site-years** have isotope `larva + adult + spider` all present at least once (Burr Pond 2023, Hop 2022, Pequabuck 2023). Under **at least two genuine isotope records per group** it is only **2/5** (Burr Pond/Pequabuck). Hop has only **1 spider**.
- Ratlum: actual isotope spiders from **2022**, but local larval/adult/biofilm isotope records from **2023**; do not put them in one 2023 local trophic path. Russell: 2022 larva/adult isotopes have only **1** each and **no spider isotope records**. The chemical PFAS source file has spiders at all five streams but missing year in many biota `collected` fields.
- **0/5 site-years** have original specimen-level `same unique specimen chemistry + isotope` links for **all three animal stages**. The *valid published site-level aggregate transfer model* is not disproved by this stricter ID check; this is only a limitation for **direct individual coupling in our proposed new macro dataset**.

[Full source site/year/support table](evidence/Farmington5_USGS_2022_2023_true_PFAS_isotope_siteyear_stage_support.csv), verified original source-only scripts:
- `analysis/20261118_Farmington_5stream_PFAS/farmington_USGS_2020_23_real_insect_spider_source_schema.py`
- `analysis/20261118_Farmington_5stream_PFAS/farmington_USGS_original_chem_spider_mass_isotope_siteyear_gate.py`
- `analysis/20261118_Farmington_5stream_PFAS/farmington_USGS_chemical_stableisotope_true_sample_ID_year_source_audit.py`.

## 4. Decision for the gesui one-paper macro question

This authentic source completes an independent **PFAS family × five stream × actual aquatic-to-riparian biological matrix opportunity** in a second major watershed relative to published Grand River WWTP-microplastic records. It **does not** provide an independently repeated **municipal WWTP source→recipient** contrast. Numerically combine neither these PFAS concentrations (ng/g wet weight) with Grand River microplastic particle counts (particle/mg) nor Bow River isotope shifts (‰).

Main comparators remain:
- **Grand River, 8 municipal WWTPs in 1 basin:** real up/down MP at larva/adult/spider, riparian Tetragnathidae concentration **8/8 downstream-positive**, but 2022/23 original true effluent-source-year×three-stages ≥2-visit full-chain gate only **4/8 < predefined 6/8**, STOP.
- **Bow River, 1 basin, 5 sites:** 2 exposed vs3 unexposed actual isotopes, author reanalysis, not five independent sources.
- **France122 river stations:** ambient TP vs NDVI matched at site/month, but **0/122 original 2018 positively WWTP-reported P sites** and NDVI not actual energy/carbon/foodweb recipient.
- **Farmington, 5 streams in 1 basin:** actual PFAS matrix in multiple animal life stages, but source-specific municipal up/down exposure unmeasured, individual biochemical and isotope analyses often nonmatched, and year incompletely observed in chemical biota field.

**No cross-basin causal effect or global nutrient export effect is identified by these four separate datasets.** A macro paper could instead propose and falsify the explicit **measurement-domain hypothesis**: apparent spatial nutrient/host-affinity gradients need not track aquatic→riparian chemical transfer; the latter is more directly identifiable with a tracer and recipient, but exact source×time×biological channel matters. This is a study-framework result, not proof of ecological subsidy magnitude. Publication value requires independent multibasin assays, not new p-values from the same one-basin authors' published data.

**Scientific priority:** Acquire an independently sampled **municipal WWTP receiving stream in a third independent basin** with effluent/date/year, aquatic larval and adult exposures and riparian spider recipient material in identical units (or honest original-source comparison without numerical pooling). Use explicit field-source gates before model fitting.
