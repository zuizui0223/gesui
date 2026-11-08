#!/usr/bin/env python3
"""Posthoc independent cross-year EEA monitoring/site/source-overlap quality audit.

The genuine 2018 WWTP mass and EQR study was completed and its model/folds
frozen before this posthoc audit. Compare with 2022 source/year original data
by exact Waterbase monitoring site_key and fixed BQE code. No models refitted,
no outcome/site selection and no causal before/after or WWTP treatment claim.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr

FILES={
"benthic_invertebrates":{
"2018":"benthic_invertebrates_2018_WWTP_source_real_EQR.csv",
"2018sha":"2a1a72d43656aca588a832dfae13860054c01a651e76bf25f984c5a374709b02",
"2022":"BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv",
"2022sha":"5e1f3c31b389d0bd769984ab055d3e2c54b63ccaa4f5eedfb2d82bcdc83c28e6",
"n2018":2634,"n2022":2475},
"phytobenthos":{
"2018":"phytobenthos_2018_WWTP_source_real_EQR.csv",
"2018sha":"8302daf2afa3f3556debc5ca769ee667bd06816d0de2e63afd4cb6b538587ad7",
"2022":"PHYTOBENTHOS_2022_SOURCE_EQR.csv",
"2022sha":"7fa5dfff9023337ac0a62a87fe406eb71783979858d9d924ab97019c6e4e5fc3",
"n2018":2670,"n2022":2784}
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def distribution(v):
    return {str(k):int(n) for k,n in v.fillna("MISSING").value_counts().items()}
def fraction(num,denom):return num/denom if denom else None

def audit_group(year18,year22,config,group):
    if len(year18)!=config["n2018"] or len(year22)!=config["n2022"]:
        raise ValueError("The original biological result site denominators have changed")
    for d in (year18,year22):
        if d.site_key.duplicated().any() or d.site_key.isna().any():
            raise ValueError("Nonunique or missing official EEA monitoring site key")
    a=year18.copy()
    b=year22.copy()
    a["positive2018"]=pd.to_numeric(a["2018_measured_upstream_P_t_per_year"],errors="raise").gt(0)
    b["positive2022"]=pd.to_numeric(b["upstream_2022_EEA_P_measured_t_per_yr"],errors="raise").gt(0)
    a["Npositive2018"]=pd.to_numeric(a["2018_measured_upstream_N_t_per_year"],errors="raise").gt(0)
    b["Npositive2022"]=pd.to_numeric(b["upstream_2022_EEA_N_measured_t_per_yr"],errors="raise").gt(0)
    # Original study restricts all registered WWTP sources to *both* positive N/P;
    # do not call a recorded source-zero year "no sewage".
    if (a.positive2018!=a.Npositive2018).any() or (b.positive2022!=b.Npositive2022).any():
        raise ValueError("Year source status does not align with original both-positive eligibility")
    joined=a.merge(b,on="site_key",how="inner",suffixes=("_18","_22"),validate="one_to_one")
    n_overlap=len(joined)
    if not n_overlap:
        return {"group":group,"status":"NO_EXACT_MONITORING_SITE_OVERLAP"}
    q18=joined["numeric_EQR_2018"].to_numpy(float)
    q22=joined["normalised_EQR"].to_numpy(float)
    if not (np.isfinite(q18)&np.isfinite(q22)).all():
        raise ValueError("Original 2018/2022 final EQR panels include missing numeric scores")
    if not (joined.countryCode_18==joined.countryCode_22).all():
        raise ValueError("Exact monitoring site key maps to different reporting country")
    year18pos=joined.positive2018.to_numpy(bool)
    year22pos=joined.positive2022.to_numpy(bool)
    method18=joined.procedureClassificationSystem_18.fillna("MISSING").astype(str)
    method22=joined.procedureClassificationSystem_22.fillna("MISSING").astype(str)
    same_method=method18==method22
    # Exact EEA monitoring site code may have been re-geocoded and sometimes
    # a different reach. Report it and DO NOT use this cross-year overlap
    # as a before/after WWTP treatment experiment.
    same_reach=joined.HYRIV_ID_18.astype(str).eq(joined.HYRIV_ID_22.astype(str))
    same_region=joined.region_3deg_block_18.astype(str).eq(joined.region_3deg_block_22.astype(str))
    same_country_blocks=joined.countryCode_18.astype(str)+"|"+joined.region_3deg_block_18.astype(str)
    if same_method.sum()>=10 and np.unique(q18[same_method]).size>1 and np.unique(q22[same_method]).size>1:
        rhosame=float(spearmanr(q18[same_method],q22[same_method]).statistic)
    else:rhosame=None
    transition={(f"{x}_to_{y}"):int(((year18pos==bool(x))&(year22pos==bool(y))).sum())
                for x in (0,1) for y in (0,1)}
    bycountry=[]
    for country,sub in joined.groupby("countryCode_18"):
        p18=sub.positive2018.astype(bool)
        p22=sub.positive2022.astype(bool)
        bycountry.append({
            "country":str(country),
            "exact_site_overlap":len(sub),
            "2018_source_positive":int(p18.sum()),
            "2022_source_positive":int(p22.sum()),
            "same_biology_classification_method":int((
                sub.procedureClassificationSystem_18.fillna("MISSING").astype(str)==
                sub.procedureClassificationSystem_22.fillna("MISSING").astype(str)).sum()),
            "same_HydroRIVERS_reach_ID":int((
                sub.HYRIV_ID_18.astype(str)==sub.HYRIV_ID_22.astype(str)).sum()),
            "2018_to_2022_source_presence_change":int((p18!=p22).sum()),
        })
    return {
      "group":group,
      "status":"SOURCE_LOCKED_CROSS_YEAR_MONITOR_SITE_OVERLAP_QC_ONLY",
      "original_2018_monitored_site_n":len(a),
      "original_2022_monitored_site_n":len(b),
      "exact_same_site_2018_and_2022":n_overlap,
      "2018_sites_reused_in_2022_fraction":fraction(n_overlap,len(a)),
      "2022_sites_also_observed_2018_fraction":fraction(n_overlap,len(b)),
      "2018_only_sites":len(a)-n_overlap,
      "2022_only_sites":len(b)-n_overlap,
      "shared_country_n":joined.countryCode_18.nunique(),
      "2018_2022_source_reported_load_status_transition_on_same_sites":transition,
      "same_site_source_presence_changed_count":int((year18pos!=year22pos).sum()),
      "same_site_exact_biological_method_code_count":int(same_method.sum()),
      "same_site_biological_method_code_changed_count":int((~same_method).sum()),
      "same_site_HydroRIVERS_reach_ID_count":int(same_reach.sum()),
      "same_site_river_reach_ID_changed_count":int((~same_reach).sum()),
      "same_site_geographic_3degree_block_count":int(same_region.sum()),
      "same_site_geographic_3degree_block_changed_count":int((~same_region).sum()),
      "descriptive_overlapping_site_2018_EQR_2022_EQR_spearman":float(spearmanr(q18,q22).statistic),
      "descriptive_identical_BQE_method_site_year_EQR_spearman":rhosame,
      "overlap_site_2018_EQR_median":float(np.median(q18)),
      "overlap_site_2022_EQR_median":float(np.median(q22)),
      "overlap_site_2022_minus_2018_EQR_median":float(np.median(q22-q18)),
      "exact_source_station_overlap_by_country":bycountry,
      "repeated_site_year_differences_are_not_instrumented_NP_effects":True,
      "source_zero_means_unreported_qualified_2018_or_2022_WWTP_load_not_no_effluent":True,
      "source_years_have_different_reporting_firm_selection_and_source_status_definitions":True,
      "2018_2022_source_data_not_from_two_independent_biological_agencies":True,
      "analysis_predeclared_original_Ridge_models_NOT_refitted":True,
    }
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--year2018",type=Path,required=True)
    p.add_argument("--year2022",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    groups={}
    for group,config in FILES.items():
        p18=a.year2018/config["2018"]
        p22=a.year2022/config["2022"]
        if sha(p18)!=config["2018sha"] or sha(p22)!=config["2022sha"]:
            raise ValueError("Original 2018/2022 scientific EQR or source files changed: "+group)
        df18=pd.read_csv(p18,dtype={"site_key":str},low_memory=False)
        df22=pd.read_csv(p22,dtype={"site_key":str},low_memory=False)
        out=audit_group(df18,df22,config,group)
        groups[group]=out
        print("GESUI_2018_2022_SAME_MONITOR_SITE_COUNTRY_QC",
          json.dumps({key:value for key,value in out.items()
              if key not in ("exact_source_station_overlap_by_country",)},sort_keys=True),flush=True)
    result={
      "schema":"gesui_temporal_2018_2022_monitored_site_and_reported_WWTP_source_status_overlap_QC_v01",
      "date_utc":datetime.now(timezone.utc).isoformat(),
      "status":"SOURCE_SHA_VERIFIED_TWO_YEAR_SITE_OVERLAP_QC_COMPLETE",
      "original_year2018_WWTP_source_and_EQR_run":42442266935,
      "original_year2022_WWTP_source_and_EQR_run":42432384196,
      "original_site_source_SHA256":{
          f"{group}_{yr}":config[yr+"sha"]
          for group,config in FILES.items() for yr in ("2018","2022")},
      "source_years_both_independently_reported_WWTP_NP_mass":True,
      "fixed_original_model_prediction_not_recomputed_or_reselected":True,
      "cross_year_site_overlap_posthoc_descriptive":True,
      "same_EEA_Biology2026_data_release_used_both_years":True,
      "does_not_establish_2018_2022_causal_WWTP_treatment":True,
      "results":groups,
      "limits":[
          "Same monitoring site appearing twice means time replication does not imply spatially independent ecological observations",
          "Model source population and country EEA reporting differ across 2018 and 2022",
          "Source presence/absence changes can be changes in reporting and plant location quality, not treatment intervention",
          "In-river EQR year differences are not a controlled treatment contrast or nutrient concentration mass balance",
          "No agriculture or other WWTP/P/N sources were retrospectively selected to maximize EQR correlation"
      ]
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GESUI_TWO_YEAR_SOURCE_EQR_SITE_OVERLAP_QC_COMPLETE",
          json.dumps({x:{"site_overlap":v.get("exact_same_site_2018_and_2022"),
                          "original_2018_n":v["original_2018_monitored_site_n"],
                          "original_2022_n":v["original_2022_monitored_site_n"]}
             for x,v in groups.items()},sort_keys=True),flush=True)
if __name__=="__main__":main()
