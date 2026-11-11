#!/usr/bin/env python3
"""Independent USDA/EPA NWCA2021 original plot H_ALL vs 2021 DMR source QA.

Original 568 wetland source frame and 2021 exact EXO total N/P monitored
near/far populations were frozen without rereading biodiversity outcomes.
Reconstruct 75 and 177 source-selected NWCA plant plots from official EPA CSV,
then independently estimate 12 original source model coefficients and
state-clustered 95% intervals via direct weighted NumPy OLS, not statsmodels.
The DMR monitoring existence is NOT a sewage-to-wetland nutrient dose.
"""
from __future__ import annotations
import argparse,hashlib,json,math,re
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd

D=Path("results/us_NWCA2021_plant_TNTP_20261110")
SOURCE=D/"EPA_NWCA2021_ORIGINAL_VEG_WATERCHEM_SITEONLY_GATE.json"
GATE=D/"NWCA2021_CAL2021_ACTUAL_DMR_TNTP_NEARFAR_STATEECO_SOURCE_ONLY_GATE.json"
FIT=D/"NWCA2021_2021_DMR_MONITORED_SAME_EXO_FEATURE_NEARFAR_HALL_FQAI_POSTHOC.json"
RESPONSES=("HALL","HNAT","HALIEN","FQAIALL","FQAINAT")
OUTCOMES=("DMR_actual2021_TN_or_TP","DMR_actual2021_BOTH_TN_TP_same_EXO")
CI=1.959963984540054

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as fp:
        for blob in iter(lambda:fp.read(1024*1024),b""):h.update(blob)
    return h.hexdigest()

def clean(x):return re.sub(r"[^A-Z0-9]","",str(x).upper())
def read_csv(p):
    df=pd.read_csv(p,low_memory=False)
    df=df.rename(columns={c:clean(c) for c in df.columns})
    df["SITEID"]=df.SITEID.astype("string").str.strip()
    df["VISITNO"]=pd.to_numeric(df.VISITNO,errors="coerce")
    return df

def samples(root,dir):
    src=json.loads((root/SOURCE).read_text())
    geo=json.loads((root/GATE).read_text())
    fit=json.loads((root/FIT).read_text())
    if fit["source_2021_DMR_exact_feature_NP_monitoring_sites_pre_response_SHA256"]!=sha(root/GATE):
        raise ValueError("Original source-only DMR N/P site sample changed")
    if not geo["SOURCE_ONLY_no_2021_NWCA_plant_H_ALL_H_NAT_FQAI_or_measured_wetland_water_values_reopened"]:
        raise ValueError("DMR near/far strata selected using biology")
    fp={k:dir/f"nwca21_{k}-data.csv" for k in ("siteinfo","vegmetrics","waterchem_wide","hydrology_sources_wide")}
    for name,file in fp.items():
        orig={"siteinfo":"SITEINFO","vegmetrics":"VEGMETRICS","waterchem_wide":"WATERCHEM","hydrology_sources_wide":"HYDROSOURCE"}[name]
        if sha(file)!=src["official_EPA_original_source_csv_receipts"][orig]["official_source_sha256"]:
            raise ValueError("Original 2021 EPA probability sample source CSV mutated: "+name)
    site=read_csv(fp["siteinfo"])
    for c in ("WGTTPCORE","YEAR","LATANALYS","LONANALYS"):site[c]=pd.to_numeric(site[c],errors="coerce")
    eligible=(site.SITETYPE.astype(str).str.upper().eq("PROB")&
              site.WGTTPCORE.gt(0)&site.YEAR.eq(2021)&site.VISITNO.eq(1)&
              site.LATANALYS.between(24,50)&site.LONANALYS.between(-126,-66)&
              site.STUDY.astype(str).str.contains("NWCA",case=False,na=False))
    site=site.loc[eligible].copy()
    if len(site)!=920:raise ValueError("Original 920 probability wetland first visits missing")
    key=["SITEID","VISITNO"]
    water=read_csv(fp["waterchem_wide"])[key].drop_duplicates()
    veg=read_csv(fp["vegmetrics"])
    veg=veg.loc[veg.VISITNO.eq(1)].copy()
    if "SAMPLETYPE" in veg and veg.SAMPLETYPE.astype(str).str.upper().eq("PLANT").any():
        veg=veg.loc[veg.SAMPLETYPE.astype(str).str.upper().eq("PLANT")].copy()
    if veg.duplicated(key).any():raise ValueError("Duplicate species plot survey sample for NWCA site")
    plots=site.merge(water,on=key,how="inner",validate="one_to_one").merge(veg,on=key,how="inner",validate="one_to_one")
    if len(plots)!=568:raise ValueError("Original EPA 568 site plant-water joint source mismatch")
    plots["source_sha"]=[
      hashlib.sha256((str(s)+"|"+str(int(v))).encode()).hexdigest()
      for s,v in zip(plots.SITEID,plots.VISITNO)]
    if set(plots.source_sha)!=set(src["source_joined_original_SITE_ID_SHA256_not_plant_or_chemistry_values"]):
        raise ValueError("Original 568 prior source-only wetland site IDs mismatch")
    for y in RESPONSES:plots[y]=pd.to_numeric(plots[y],errors="coerce")
    plots["STATE_ECO"]=plots.PSTLCODE.astype(str).str.strip()+"|"+plots.AGECO9.astype(str).str.strip()
    outputs={}
    for key in OUTCOMES:
        source=geo["original_source_DMR_2021_ANY_NP_and_BOTH_NP_strict_same_state_ecoregion_controls"][key]
        near=set(source["original_NWCA2021_NEAR_source_site_SHA256s"])
        far=set(source["original_NWCA2021_FAR_source_site_SHA256s"])
        if not near.isdisjoint(far):raise ValueError("Original DMR near and far site sets overlap")
        q=plots.loc[plots.source_sha.isin(near|far)].copy()
        if len(q)!=len(near)+len(far):raise ValueError("DMR 2021 source site missing original measured plant record")
        q["NEAR"]=q.source_sha.isin(near).astype(int)
        if len(q.STATE_ECO.unique())!=source["state_ecoregion_cells_with_NEAR_AND_FAR"]:
            raise ValueError("Source region change after 2021 DMR near/far gate")
        outputs[key]=q
    return src,geo,fit,outputs,fp

def numpy_fit(d,y,weighted):
    q=d.loc[d[y].notna()].copy()
    states=sorted(q.PSTLCODE.astype(str).unique())
    X=[np.ones(len(q)),q.NEAR.to_numpy(float)]
    for column in ("STATE_ECO","WETCLSGRP","HYDCLS"):
        text=q[column].astype(str)
        classes=sorted(text.unique())
        for c in classes[1:]:
            X.append((text==c).to_numpy(float))
    X.append(q.LATANALYS.to_numpy(float))
    mat=np.column_stack(X)
    response=q[y].to_numpy(float)
    weights=q.WGTTPCORE.to_numpy(float) if weighted else np.ones(len(q))
    weights=weights/np.median(weights)
    w=np.sqrt(weights)
    xm=mat*w[:,None]
    ym=response*w
    b=np.linalg.lstsq(xm,ym,rcond=None)[0]
    r=ym-xm@b
    bread=np.linalg.pinv(xm.T@xm)
    meat=np.zeros_like(bread)
    for s in states:
        mask=(q.PSTLCODE.astype(str)==s).to_numpy()
        v=xm[mask].T@r[mask]
        meat+=np.outer(v,v)
    N=len(q);G=len(states);K=mat.shape[1]
    if K>=N:raise ValueError("Too many original source model regressors relative to 2021 NWCA sites")
    correction=(G/(G-1))*((N-1)/(N-K))
    se=math.sqrt(float((bread@meat@bread*correction)[1,1]))
    beta=float(b[1])
    return {
       "n":N,"states":G,
       "n_near":int(q.NEAR.sum()),"n_far":int((q.NEAR==0).sum()),
       "source_matrix_columns":K,"source_design_rank":int(np.linalg.matrix_rank(xm)),
       "independent_beta_near10_minus_far20":beta,
       "independent_state_cluster95":[beta-CI*se,beta+CI*se],
       "weighted_survey_area":weighted
    }

def compare(actual,old,name):
    if old is None or old.get("STUDY_STOP"):
        if old is None:return None
        if actual["n"]<50:return actual
        raise ValueError("Original 2021 DMR NWCA sample originally stopped but independent fit would not")
    if actual["n"]!=old["n_firstvisit_study_wetlands"]:
        raise ValueError("Original 2021 DMR monitoring plant sample n changed "+name)
    if actual["n_near"]!=old["n_near_2021_DMR_monitored_same_external_feature"] or actual["n_far"]!=old["n_control_far20km"]:
        raise ValueError("Original site group count differs "+name)
    beta=old["beta_near_DMR_reported_POTW_minus_far_in_original_response_units"]
    if abs(actual["independent_beta_near10_minus_far20"]-beta)>1e-7:
        raise ValueError("Original standardized plant recorded difference differs "+name)
    for a,b in zip(actual["independent_state_cluster95"],old["state_cluster95"]):
        if abs(a-b)>1e-5:
            raise ValueError("Independent 2021 DMR wetland state-cluster CI mismatch "+name+
                             " new "+repr(actual["independent_state_cluster95"])+
                             " original "+repr(old["state_cluster95"]))
    return actual

def run(root,dir):
    src,geo,fit,samples,files=samples_fn(root,dir)
    models={}
    for group,d in samples.items():
        previous=fit[group]
        after={}
        for resp in RESPONSES:
            val=numpy_fit(d,resp,False)
            after[resp]=compare(val,previous["original_2021_standardized_plants_vs_2021_reported_near_POTW_feature_models"][resp],group+"/"+resp)
        val=numpy_fit(d,"HALL",True)
        after["HALL_NWCA_designWeighted"]=compare(val,previous["original_2021_all_plant_Shannon_NWCA_area_weighted_sensitivity"],group+"/HALL_weighted")
        models[group]=after
    return {
       "schema":"gesui_EPA_NWCA2021_160POTW_2021DMR_2021_source75_177wetland_12_24state_independent_numpy_sourceQA_v01",
       "original_568_US_EPA_NWCA_source_only_SHA256":sha(root/SOURCE),
       "original_EPA_2021_DMR_source_monitored_nearfar_gate_SHA256":sha(root/GATE),
       "original_primary_2021_NWCA_DMR_wetland_models_SHA256":sha(root/FIT),
       "original_EPA2021_raw_site_vegetation_water_source_SHA256":{k:sha(v) for k,v in files.items()},
       "original_568_source_site_SHA256s_verified_against_preoutcome_sample":True,
       "independent_75_and_177_source_plant_HALL_HNAT_HALlEN_FQAI_and_weighted_HALL_numpy_OLS_state_cluster":models,
       "historical_reported_outfall_NP_source_is_NOT_actual_wetland_water_hydrology_or_treatment":True,
       "site_subset_analysis_POSTHOC_and_small_12USstate_clusters_in_strict_both_group":True
    }

samples_fn=samples
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-dir",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    root=Path(__file__).resolve().parents[1]
    o=run(root,a.source_dir)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({k:{r:{"beta":v["independent_beta_near10_minus_far20"],"CI":v["independent_state_cluster95"],"n":v["n"]}
             for r,v in g.items()} for k,g in o["independent_75_and_177_source_plant_HALL_HNAT_HALlEN_FQAI_and_weighted_HALL_numpy_OLS_state_cluster"].items()},indent=2))

if __name__=="__main__":main()
