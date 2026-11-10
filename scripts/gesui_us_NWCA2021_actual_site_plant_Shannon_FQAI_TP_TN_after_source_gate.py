#!/usr/bin/env python3
"""US EPA NWCA2021 standardized wetland plant Shannon × measured N/P.

STRICTLY matched original source probability sites/IDs prior to seeing chemical
or H_ALL responses. Actual water chemistry is NOT labeled WWTP effluent.
Records, flags, real measured P ug/L vs N mg/L, 9 ecoregions, original design
weights and 1st visit are audited. All results (positive/negative/null) kept.
"""
from __future__ import annotations
import argparse,hashlib,json,math,re
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import spearmanr

BASE=Path("results/us_NWCA2021_plant_TNTP_20261110")
GATE=BASE/"EPA_NWCA2021_ORIGINAL_VEG_WATERCHEM_SITEONLY_GATE.json"
FILES={"SITEINFO":"nwca21_siteinfo-data.csv","WATERCHEM":"nwca21_waterchem_wide-data.csv",
       "VEGMETRICS":"nwca21_vegmetrics-data.csv","HYDROSOURCE":"nwca21_hydrology_sources_wide-data.csv"}
Y=("HALL","HNAT","HALIEN","FQAIALL","FQAINAT")
def norm(s):return re.sub("[^A-Z0-9]","",str(s).upper())
def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for buf in iter(lambda:f.read(1024*1024),b""):h.update(buf)
    return h.hexdigest()

def read(path):
    x=pd.read_csv(path,low_memory=False,encoding="utf-8")
    columns={col:norm(col) for col in x.columns}
    if len(set(columns.values()))!=len(columns):
        raise ValueError("Original official EPA column canonicalization collision")
    x=x.rename(columns=columns)
    if "SITEID" not in x:raise ValueError("Original EPA source lacks SITE_ID")
    x["SITEID"]=x.SITEID.astype("string").str.strip()
    x["VISITNO"]=pd.to_numeric(x["VISITNO"],errors="coerce") if "VISITNO" in x else 1
    return x

def select_and_join(root,files):
    audit=json.loads((root/GATE).read_text(encoding="utf-8"))
    if not audit["precommitted_SOURCE_gate_before_measured_NTLP_and_plant_Shannon"]["may_open_2021_plant_diversity_and_measured_TNTP_results"]:
        raise ValueError("Previously frozen NWCA source opportunity failed")
    if not audit["all_original_2021_plant_H_ALL_FQAI_and_water_NTL_PTL_values_NOT_READ"]:
        raise ValueError("Original source gate contaminated with 2021 ecology outcomes")
    for key,p in files.items():
        if sha(p)!=audit["official_EPA_original_source_csv_receipts"][key]["official_source_sha256"]:
            raise ValueError("Original EPA measured chemistry/plant source bytes changed after source-only gate")
    site=read(files["SITEINFO"])
    chem=read(files["WATERCHEM"])
    veg=read(files["VEGMETRICS"])
    hydro=read(files["HYDROSOURCE"])
    site["WGTTPCORE"]=pd.to_numeric(site["WGTTPCORE"],errors="coerce")
    for key in ("LATANALYS","LONANALYS","YEAR"):
        site[key]=pd.to_numeric(site[key],errors="coerce")
    source_selected=(site.SITETYPE.astype("string").str.upper().eq("PROB").fillna(False)&
        site.WGTTPCORE.gt(0).fillna(False)&
        site.LATANALYS.between(24,50).fillna(False)&
        site.LONANALYS.between(-126,-66).fillna(False)&
        site.YEAR.eq(2021).fillna(False)&site.VISITNO.eq(1).fillna(False))
    if "STUDY" in site:source_selected&=site.STUDY.astype("string").str.contains("NWCA",case=False,na=False)
    site=site.loc[source_selected].copy()
    key=["SITEID","VISITNO"]
    if site.duplicated(key).any():raise ValueError("Duplicated source eligibility site first visit")
    qset=set(audit["source_joined_original_SITE_ID_SHA256_not_plant_or_chemistry_values"])
    if len(qset)!=audit["source_SITE_ID_and_VISIT_NO_matched_chemistry_vegmetrics_rows_n"]:
        raise ValueError("Original precommitted qualifying site receipt keys duplicated")
    site_valid_hash={hashlib.sha256((str(a)+"|"+str(int(b))).encode()).hexdigest()
                     for a,b in zip(site.SITEID,site.VISITNO)}
    if not qset.issubset(site_valid_hash):raise ValueError("Source-sample site hash differs from original site metadata")
    def unique_table(table,name):
        t=table.copy()
        if name=="WATERCHEM" and "SAMCODE" in t:
            sam=t.SAMCODE.astype("string").str.upper().str.strip()
            categories=dict(Counter(sam.dropna()))
            # The original official NWCA 'REGULAR' sample may be encoded R or REG;
            # do not arbitrarily combine laboratory duplicates with site value.
            wanted=sam.isin(("R","REG","REGULAR","1"))
            if wanted.any():t=t.loc[wanted].copy()
        else:categories={}
        t=t.loc[t.VISITNO.eq(1)].copy()
        if name=="VEGMETRICS" and "SAMPLETYPE" in t:
            ptype=t.SAMPLETYPE.astype("string").str.upper().str.strip()
            if (ptype=="PLANT").any():t=t.loc[ptype=="PLANT"].copy()
        if t.duplicated(key).any():
            raise ValueError("Unresolved duplicate original EPA source visit after only category-based filtering: "+name)
        return t,categories
    chem,chem_cat=unique_table(chem,"WATERCHEM")
    veg,_=unique_table(veg,"VEGMETRICS")
    hydro,_=unique_table(hydro,"HYDROSOURCE")
    paired=site.merge(chem,on=key,how="inner",validate="one_to_one",suffixes=("","_WC"))
    paired=paired.merge(veg,on=key,how="inner",validate="one_to_one",suffixes=("","_V"))
    if len(paired)!=len(qset):
        raise ValueError("Original source-only matched 568 source sites not retained in the measured outcome stage")
    original={hashlib.sha256((str(a)+"|"+str(int(b))).encode()).hexdigest()
              for a,b in zip(paired.SITEID,paired.VISITNO)}
    if original!=qset:
        raise ValueError("Any original site changed because of plant/nutrient outcome data")
    hydrolabels=hydro[key+[name for name in ("INFLOWPRESENT","GRDWATERPRESENT","OVBKFLOODPRESENT","OTHERSOURCEPRESENT") if name in hydro]]
    paired=paired.merge(hydrolabels,on=key,how="left",validate="one_to_one",suffixes=("","_HYDRO"))
    return paired,chem_cat,audit

def summarize_flags(d):
    out={}
    for term in ("NTLNARSFLAG","PTLNARSFLAG","NTLRESULT","PTLRESULT","HALL","HNAT","HALIEN","FQAIALL","FQAINAT"):
        if term not in d:raise ValueError("Expected independently predeclared actual EPA NWCA outcome column missing "+term)
        v=d[term]
        if term.endswith("FLAG"):out[term]=dict(Counter(v.astype("string").fillna("NA").str.strip().str.upper()))
        else:out[term]={"nonmissing_numeric":int(pd.to_numeric(v,errors="coerce").notna().sum())}
    return out

def valid_rows(d):
    out=d.copy()
    for name in (*Y,"NTLRESULT","PTLRESULT"):
        out[name]=pd.to_numeric(out[name],errors="coerce")
    if "NTLNARSFLAG" in out and "PTLNARSFLAG" in out:
        def ok(col):
            q=out[col].astype("string").str.strip().str.upper().fillna("")
            return q.isin(("","NA","L"))
        # 'L' = analyte detected above MDL but below reporting limit:
        # retain values but show QC flag; remove ND, NR, H, Q.
        good=ok("NTLNARSFLAG")&ok("PTLNARSFLAG")
    else:good=pd.Series([False]*len(out),index=out.index)
    temp=out.loc[good&
        out.NTLRESULT.gt(0).fillna(False)&out.PTLRESULT.gt(0).fillna(False)&
        out.HALL.ge(0).fillna(False)&np.isfinite(out.HALL.fillna(np.nan).astype(float))].copy()
    temp["log10_TN_mgL"]=np.log10(temp.NTLRESULT)
    temp["log10_TP_mgL"]=np.log10(temp.PTLRESULT/1000.)
    temp["orig_site_WGT_TP_CORE"]=pd.to_numeric(temp.WGTTPCORE,errors="coerce")
    temp["orig_LAT_ANALYS"]=pd.to_numeric(temp.LATANALYS,errors="coerce")
    for k in ("AGECO9","WETCLSGRP","HYDCLS","PSTLCODE"):
        temp[k]=temp[k].astype("string").str.strip()
    temp=temp.dropna(subset=["AGECO9","WETCLSGRP","HYDCLS","PSTLCODE","log10_TN_mgL","log10_TP_mgL","HALL"])
    if temp.duplicated(["SITEID","VISITNO"]).any():raise ValueError("Duplicate source site in actual 2021 plant analysis")
    return temp

def fitted(source,y,terms,weighted=False):
    d=source.loc[pd.to_numeric(source[y],errors="coerce").notna()].copy()
    if len(d)<100:return None
    formula=y+"~"+terms+"+C(AGECO9)+C(WETCLSGRP)+C(HYDCLS)+orig_LAT_ANALYS"
    if weighted:
        # Original survey weights represent population wetland area, not
        # independent treatment replication; WLS exploratory sensitivity.
        w=np.asarray(d.orig_site_WGT_TP_CORE,dtype=float)
        w=w/np.median(w)
        m=smf.wls(formula,data=d,weights=w).fit(cov_type="cluster",cov_kwds={"groups":d.PSTLCODE})
    else:m=smf.ols(formula,data=d).fit(cov_type="cluster",cov_kwds={"groups":d.PSTLCODE})
    result={"formula":formula,"N_actual_site_visits":len(d),"N_US_states":int(d.PSTLCODE.nunique()),
            "N_NWCA_ecoregions":int(d.AGECO9.nunique()),"model_kind":"NWCA_pop_area_weighted_WLS" if weighted else "site_equal_OLS",
            "source_original_state_clustered_not_WWTP_treatment_randomized":True,
            "R2":float(m.rsquared),"coefficients":{}}
    for x in ("log10_TP_mgL","log10_TN_mgL"):
        if x in m.params:
            result["coefficients"][x]={
                "beta_per_10fold_concentration_in_original_NWCA_diversity_units":float(m.params[x]),
                "US_state_cluster95":list(map(float,m.conf_int().loc[x])),
                "nominal_p":float(m.pvalues[x])
            }
    return result

def run(root,files):
    original,labtypes,gate=select_and_join(root,files)
    qc=summarize_flags(original)
    qual=valid_rows(original)
    sufficient=bool(len(qual)>=200 and qual.PSTLCODE.nunique()>=10 and qual.AGECO9.nunique()>=3)
    receipt={
      "schema":"gesui_US_EPA_NWCA2021_real_wetland_plot_diversity_TN_TP_same_site_post_source_gate_v01",
      "UTC":datetime.now(timezone.utc).isoformat(),
      "original_source_only_nwcaSITEID_plant_waterchem_SHA256":hashlib.sha256((root/GATE).read_bytes()).hexdigest(),
      "original_EPA_four_source_csv_sha256":{name:sha(p) for name,p in files.items()},
      "original_2021_probability_site_plant_water_568_sites_preserved_before_outcome":True,
      "source_site_visit_first_only":True,
      "original_EPA_waterchem_lab_sample_codes":labtypes,
      "observed_original_source_columns_QC_nonmissing_flags":qc,
      "source_568_veg_and_waterchem_joined_site_rows":len(original),
      "source_actual_water_TP_ugP_L_converted_to_mgP_L_divide1000":True,
      "source_actual_NTL_mgN_L":True,
      "original_2021_H_ALL_Shannon_is_standard_plant_community_diversity_not_GBIF_richness":True,
      "source_joint_TNTP_and_Shannon_nonmissing_quality_sites":len(qual),
      "source_joint_actual_ecoregions":int(qual.AGECO9.nunique()),
      "source_joint_actual_states":int(qual.PSTLCODE.nunique()),
      "precommitted_stage1_actual_plant_chemistry_gate":{
          "minimum_clean_H_ALL_and_TNTP_sites":200,
          "minimum_distinct_US_states":10,
          "minimum_distinct_NWCA_ecoregions":3,
          "may_test_actual_NP_plant_standardized_community_correlations":sufficient
      },
      "US2021_measured_TNTP_are_NOT_proven_WWTP_sewage_water":True,
      "waterchem_observed_only_when_water_available_site_selection_risk":True,
      "original_2021_NWCA_probability_weight_not_random_wastewater_treatment":True,
      "never_impute_missing_plant_or_NP_measurement_as_zero":True
    }
    if not sufficient:
        receipt["status"]="STOP_INSUFFICIENT_REAL_NTP_HALL_SOURCE"
        return receipt
    receipt["status"]="EXPLORATORY_STANDARDIZED_WETLAND_PLANT_CHEMISTRY_ANALYZED"
    def qtile(a):
        z=np.asarray(a,dtype=float)
        return [float(np.quantile(z,q)) for q in (.1,.25,.5,.75,.9)]
    receipt["original_qualified_waterchem_and_plant_metric_descriptive"]={
      "TP_mgP_L_quintiles":qtile(qual.PTLRESULT/1000),
      "TN_mgN_L_quintiles":qtile(qual.NTLRESULT),
      "H_ALL_Shannon_quantiles":qtile(qual.HALL),
      "NWCA_state_counts":dict(Counter(qual.PSTLCODE)),
      "NWCA_ecoregion_counts":dict(Counter(qual.AGECO9)),
      "wetland_class_counts":dict(Counter(qual.WETCLSGRP)),
      "hydrologic_class_counts":dict(Counter(qual.HYDCLS))
    }
    plots=[]
    for outcome in Y:
        if outcome=="HALL":modelspecs=[
            ("H_ALL_primary_TP_TN", "log10_TP_mgL+log10_TN_mgL"),
            ("H_ALL_only_TP","log10_TP_mgL"),
            ("H_ALL_only_TN","log10_TN_mgL")
        ]
        else:modelspecs=[(outcome+"_TP_TN","log10_TP_mgL+log10_TN_mgL")]
        for name,spec in modelspecs:
            if len(qual.loc[qual[outcome].notna()])>=100:
                plots.append((name,outcome,spec,False))
                if outcome=="HALL" and name=="H_ALL_primary_TP_TN":
                    plots.append(("H_ALL_primary_TP_TN_design_weighted",outcome,spec,True))
    receipt["actual_standardized_wetland_plant_models"]={
        name:fitted(qual,outcome,terms,weighted)
        for name,outcome,terms,weighted in plots}
    if receipt["actual_standardized_wetland_plant_models"]["H_ALL_primary_TP_TN"] is None:
        raise ValueError("Primary source-supported wetland Shannon model not estimable")
    def raw_rho(x,y):
        try:return float(spearmanr(qual[x],qual[y],nan_policy="omit").statistic)
        except:return None
    receipt["unadjusted_posthoc_spearman"]={
      "TP_vs_HALL":raw_rho("PTLRESULT","HALL"),
      "TN_vs_HALL":raw_rho("NTLRESULT","HALL"),
      "TP_vs_HNAT":raw_rho("PTLRESULT","HNAT"),
      "TP_vs_HALlEN":raw_rho("PTLRESULT","HALIEN")
    }
    return receipt

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--source-dir",type=Path,required=True)
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args()
    root=Path(__file__).resolve().parents[1]
    files={key:args.source_dir/p for key,p in FILES.items()}
    out=run(root,files)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":out["status"],
        "source_gate":out["precommitted_stage1_actual_plant_chemistry_gate"],
        "n":out["source_joint_TNTP_and_Shannon_nonmissing_quality_sites"],
        "models":out.get("actual_standardized_wetland_plant_models"),
        "rho":out.get("unadjusted_posthoc_spearman")},indent=2))

if __name__=="__main__":
    main()
