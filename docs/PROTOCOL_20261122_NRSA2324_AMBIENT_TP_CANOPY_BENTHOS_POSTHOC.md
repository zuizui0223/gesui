# Original EPA NRSA2023–24 ambient TP × field canopy × benthic communities — posthoc protocol (2026-11-22)

**Not sewage-specific**, and not a new independent national sampling campaign. This follows source failure of historical POTW effluent TP × canopy at PR #139 (31 sites total/26 BETB). Do not describe measured in-stream TP as wastewater-origin TP.

Original official source SHA:
- NRSA2023–24 SiteInfo f7d14cf1574d12828b3aa97bd7b08b95f22aec2115753d5b75514f8eed43d013
- WaterChem b34d3f9c4b4bdc763cb5e53a856bbab329f1646d248f223f78e9404bc4560d71
- PHab 8a6094bb87f8a57cf45a35a1272b137e6a9327cefd158f4bee88b025eaf9b2f4
- BenthicCount 5290d679f7eb444e97c8d45096beb4531b450186cef190580f41f4b9430b0805
- BenthicMMI 5ef607e8e4726ed8c3d24b9d15fe351e7300784859654b812fa595c09a4009dc

Source-only schema [Actions #45014068137](https://github.com/zuizui0223/gesui/actions/runs/45014068137) established common native UID **2139** across all 5 sources; exact method-compatible first-level source counts **BERW/WADEABLE 1307** and **BETB/BOATABLE 825**, other 7 nonmatching methods excluded. Every original response remains hidden at this selection stage.

## Estimands and no outcome-based choices

- Exposure: field-measured ambient total phosphorus `PTL_RESULT` in original `PTL_UNITS`, correctly converted ug P/L to mg/L, positive and blank original `PTL_NARS_FLAG`, log10-transformed. This is NOT effluent load or source-reached TP.
- Primary response: original standardized benthic multimetric score `MMI_BENT`. Secondary response: counted EPT individual fraction from exactly original `TOTAL300` and `ORDER` (Ephemeroptera/Plecoptera/Trichoptera). The latter is neither species richness nor aquatic emergence biomass.
- Two strata analyzed **separately**: BERW/WADEABLE using midchannel canopy `XCDENMID` (%), BETB/BOATABLE using bank canopy `XCDENBK` (%). No protocol-agnostic canopy pooling.
- Baseline explanatory terms: native state + AG_ECO9, lat/lon, original Strahler stream order, year and within-year season sine/cosine, log10 conductivity, original ambient TN, log10 wetted stream width, **both** logTP and canopy main effects. Added term in full model: `z(logTP) × z(canopy)`. All complete-case eligibility is source-variable/validity determined, not sign selected.
- Hold out entire original `HUC8` in deterministic fivefold Ridge(alpha=10) cross-validation, fit predictor transformers on training folds only. Compare the **same exact observations** and control terms between baseline and interaction. Primary comparison = equal-HUC8 held-out MSE gain, not random site split or in-sample p. Report original-site-visit-weighted MSE as a secondary comparison.
- For descriptive effect size fit standardized-response OLS with HUC8-cluster covariance. 2000 fixed-out-of-fold HUC8 bootstrap resamples give descriptive 95% for MSE gain. Do not call this a randomized or causal effect; initial question selected after prior 2018 and 2023 biofilm/habitat findings (posthoc).
- Gate method separately before fitting: at least 250 complete visits/100 independent HUC8, at least 3 HUC8 in each fixed fold, no UID/date mismatches, all canopy 0–100%, outcome MMI finite, EPT denominator >0. If a method/outcome fails report STOP, do not borrow observations from the other protocol.
- Report both endpoints and both strata, including nulls, failures and sensitivity limitations. A nonzero interaction need not imply canopy changes actual nutrient retention, sewage flow, ecosystem-service production or aquatic→terrestrial insect export. Any new mechanism would need direct light, algae, oxygen and flow.
