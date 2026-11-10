#!/usr/bin/env python3
"""Independent EPA NWCA2021 real nutrient/plant Shannon source+numeric audit.

Separate original EPA CSV parser and NumPy least-squares / state-cluster
sandwich, not original statsmodels formula estimator. Verify exactly 568
first-visit probability sites and 499 nutrient QC sites, then reproduce H_ALL
two-nutrient and separate-N/P coefficients+CI without asserting effluent source.
"""
from __future__ import annotations
import argparse,hashlib,json,math,re
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

D=Path("results/us_NWCA2021_plant_TNTP_20261110")
S=D/"EPA_NWCA2021_ORIGINAL_VEG_WATERCHEM_SITEONLY_GATE.json"
F=D/"EPA_NWCA2021_REAL_PLANT_SHANNON_TNTP_WETLANDS_STATE_ECOREGION_QA.json"
FILES={"SITEINFO":"nwca21_siteinfo-data.csv","WATERCHEM":"nwca21_waterchem_wide-data.csv",
       "VEGMETRICS":"nwca21_vegmetrics-data.csv","HYDROSOURCE":"nwca21_hydrology_sources_wide-data.csv"}

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()
def canon(s):return re.sub(r"[^A-Z0-9]","",str(s).upper())
def read(p):
    x=pd.read_csv(p,low_memory=False)
    x=x.rename(columns={k:canon(k) for k in x.columns})
    x["SITEID"]=x.SITEID.astype("string").str.strip()
    x["VISITNO"]=pd.to_numeric(x.VISITNO,errors="coerce")
    return x

def sample(root,dir):
    source=json.loads((root/S).read_text())
    output=json.loads((root/F).read_text())
    if output["original_source_only_nwcaSITEID_plant_waterchem_SHA256"]!=sha(root/S):
        raise ValueError("Pre-response EPA NWCA2021 source gate SHA changed")
    files={k:dir/name for k,name in FILES.items()}
    for k,p in files.items():
        if sha(p)!=source["official_EPA_original_source_csv_receipts"][k]["official_source_sha256"]:
            raise ValueError("Official original EPA NWCA2021 CSV file bytes modified")
    a=read(files["SITEINFO"])
    b=read(files["WATERCHEM"])
    c=read(files["VEGMETRICS"])
    a["WGTTPCORE"]=pd.to_numeric(a.WGTTPCORE,errors="coerce")
    valid=(a.SITETYPE.eq("PROB")&
           a.WGTTPCORE.gt(0)&pd.to_numeric(a.YEAR,errors="coerce").eq(2021)&
           pd.to_numeric(a.LATANALYS,errors="coerce").between(24,50)&
           pd.to_numeric(a.LONANALYS,errors="coerce").between(-126,-66)&
           a.VISITNO.eq(1)&a.STUDY.astype(str).str.contains("NWCA",case=False,na=False))
    a=a[valid].copy()
    b=b[b.VISITNO.eq(1)].copy()
    if "SAMCODE" in b and b.SAMCODE.astype("string").str.upper().isin(("R","REG","REGULAR","1")).any():
        b=b[b.SAMCODE.astype("string").str.upper().isin(("R","REG","REGULAR","1"))]
    c=c[c.VISITNO.eq(1)].copy()
    if "SAMPLETYPE" in c and c.SAMPLETYPE.astype("string").str.upper().eq("PLANT").any():
        c=c[c.SAMPLETYPE.astype("string").str.upper().eq("PLANT")]
    key=["SITEID","VISITNO"]
    if any(x.duplicated(key).any() for x in (a,b,c)):
        raise ValueError("Original official EPA source visit cardinality/duplicate issue")
    joined=a.merge(b,on=key,validate="one_to_one").merge(c,on=key,validate="one_to_one")
    ids={hashlib.sha256((str(x)+"|"+str(int(y))).encode()).hexdigest()
         for x,y in zip(joined.SITEID,joined.VISITNO)}
    if ids!=set(source["source_joined_original_SITE_ID_SHA256_not_plant_or_chemistry_values"]):
        raise ValueError("Independently reconstructed eligible 568 original EPA NWCA2021 joint sites differ")
    if len(joined)!=568:raise ValueError("Different original wetland plant-water site cohort")
    for col in ("NTLRESULT","PTLRESULT","HALL"):
        joined[col]=pd.to_numeric(joined[col],errors="coerce")
    qc=joined.NTLNARSFLAG.astype("string").str.upper().str.strip().fillna("").isin(("","NA","L"))&\
       joined.PTLNARSFLAG.astype("string").str.upper().str.strip().fillna("").isin(("","NA","L"))
    obs=joined.loc[qc&joined.NTLRESULT.gt(0)&joined.PTLRESULT.gt(0)&joined.HALL.ge(0)].copy()
    for c0 in ("AGECO9","WETCLSGRP","HYDCLS","PSTLCODE"):
        obs[c0]=obs[c0].astype("string").str.strip()
    obs=obs.dropna(subset=["AGECO9","WETCLSGRP","HYDCLS","PSTLCODE","HALL","NTLRESULT","PTLRESULT"])
    if len(obs)!=499 or obs.PSTLCODE.nunique()!=48 or obs.AGECO9.nunique()!=9:
        raise ValueError("Original 499/48state/9ecoregion chemistry quality gate not recreated")
    obs["logTP"]=np.log10(obs.PTLRESULT/1000.)
    obs["logTN"]=np.log10(obs.NTLRESULT)
    obs["lat"]=pd.to_numeric(obs.LATANALYS,errors="coerce")
    return obs,files,output

def independent_state_OLS(d,include_P,include_N):
    terms=["Intercept"]
    X=[np.ones(len(d))]
    if include_P:terms.append("logTP");X.append(d.logTP.to_numpy(dtype=float))
    if include_N:terms.append("logTN");X.append(d.logTN.to_numpy(dtype=float))
    for col in ("AGECO9","WETCLSGRP","HYDCLS"):
        levels=sorted(d[col].astype(str).unique())
        for level in levels[1:]:
            terms.append(col+"["+level+"]")
            X.append((d[col].astype(str)==level).to_numpy(dtype=float))
    terms.append("lat");X.append(d.lat.to_numpy(dtype=float))
    mat=np.column_stack(X)
    y=d.HALL.to_numpy(dtype=float)
    coefficients=np.linalg.lstsq(mat,y,rcond=None)[0]
    residual=y-mat@coefficients
    B=np.linalg.pinv(mat.T@mat)
    allstates=sorted(d.PSTLCODE.unique())
    meat=np.zeros_like(B)
    for state in allstates:
        mask=(d.PSTLCODE==state).to_numpy()
        v=mat[mask].T@residual[mask]
        meat+=np.outer(v,v)
    n=len(d);g=len(allstates);k=np.linalg.matrix_rank(mat)
    clusterfactor=(g/(g-1))*((n-1)/(n-k))
    covariance=B@meat@B*clusterfactor
    result={"n":n,"state_clusters":g,"matrix_columns":mat.shape[1],"design_rank":int(k),"modelled_parameters":{}}
    for name in ("logTP","logTN"):
        if name not in terms:continue
        idx=terms.index(name)
        beta=float(coefficients[idx])
        se=math.sqrt(float(covariance[idx,idx]))
        result["modelled_parameters"][name]={
            "beta_per_10fold_original_measured_nutrient":beta,
            "US_state_cluster95":[beta-1.95996398454*se,beta+1.95996398454*se],
        }
    return result

def test(root,dir):
    joined,files,reported=sample(root,dir)
    model=reported["actual_standardized_wetland_plant_models"]
    configs={
        "H_ALL_primary_TP_TN":(True,True),
        "H_ALL_only_TP":(True,False),
        "H_ALL_only_TN":(False,True)
    }
    out={}
    for name,(P,N) in configs.items():
        one=independent_state_OLS(joined,P,N)
        previous=model[name]
        assert one["n"]==previous["N_actual_site_visits"]==499
        assert one["state_clusters"]==previous["N_US_states"]==48
        for k,v in one["modelled_parameters"].items():
            canon_key="log10_TP_mgL" if k=="logTP" else "log10_TN_mgL"
            ref=previous["coefficients"][canon_key]
            if abs(v["beta_per_10fold_original_measured_nutrient"]-ref["beta_per_10fold_concentration_in_original_NWCA_diversity_units"])>1e-8:
                raise ValueError("Independent EPA wetland logTP/TN Shannon slope does not match original")
            for a,b in zip(v["US_state_cluster95"],ref["US_state_cluster95"]):
                if abs(a-b)>1e-5:
                    raise ValueError("Independent US state-cluster 95%CI differs from original H_ALL: calculated "+repr(v["US_state_cluster95"])+", original "+repr(ref["US_state_cluster95"])+", rank "+str(one["design_rank"])+", K "+str(one["matrix_columns"]))
        out[name]=one
    ranks={"TN_vs_HALL":float(spearmanr(joined.NTLRESULT,joined.HALL).statistic),
           "TP_vs_HALL":float(spearmanr(joined.PTLRESULT,joined.HALL).statistic)}
    for k,v in ranks.items():
        if abs(v-reported["unadjusted_posthoc_spearman"][k])>1e-12:
            raise ValueError("Original EPA measured nutrient / vegetation rank nonreproducible")
    return {
        "schema":"gesui_US_NWCA2021_499samewetland_2021_measured_TNTP_Shannon_independent_numpy_statecluster_QA_v01",
        "source_original_csv_sha256":{k:sha(v) for k,v in files.items()},
        "source_frozen_568site_gate_SHA256":sha(root/S),
        "source_original_first_results_SHA256":sha(root/F),
        "source_568_site_probability_firstvisit_joint_preoutcome_matched":True,
        "measured_chemistry_flag_QUALITATIVE_exclusions_S_shipping_Q_issues_no_NR_ND":True,
        "independently_recomputed_complete499_48states_9ecoregions":True,
        "three_primary_and_single_nutrient_beta_and_US_state_cluster95_exact":out,
        "unadjusted_measured_TN_and_TP_Spearman_to_actual_H_ALL_reproduced":ranks,
        "observed_standardized_wetland_Shannon_not_opportunistic_GBIF_indexing":True,
        "measured_wetland_TNTP_not_proven_to_come_from_sewage":True,
        "study_is_US_wetlands_not_global_nonwetland_plants":True
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--source-dir",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    arg=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    x=test(root,arg.source_dir)
    arg.out.parent.mkdir(parents=True,exist_ok=True)
    arg.out.write_text(json.dumps(x,indent=2,sort_keys=True)+"\n")
    print(json.dumps(x["three_primary_and_single_nutrient_beta_and_US_state_cluster95_exact"],indent=2))
if __name__=="__main__":
    main()
