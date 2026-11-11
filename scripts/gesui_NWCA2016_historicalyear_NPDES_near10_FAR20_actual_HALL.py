#!/usr/bin/env python3
"""2016 as-of-2016 permitted sewage POINT proxy × standardized NWCA H_ALL.

SOURCE-frozen 50 NEAR10 / 162 FAR20, 19 US states, 23 state/ecoregion bins.
This is POSTHOC (2016 2021-proxy effect already known). EPA ECHO Nov2026
permit effective dates do not verify actual 2016 WWTP effluent delivery.
2011 historical source n12 near fails; no historical 2011 effect fitted.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter
from pathlib import Path
from datetime import datetime,timezone
import pandas as pd
import numpy as np

from gesui_NWCA2011_2016_2021_POTW_near10_HALL_HNAT_FQAI_temporal_source_replay import model,sha,q

D=Path("results/us_NWCA2021_plant_TNTP_20261110")
SITE=D/"NWCA2011_2016_ORIGINAL_SITE_VEG_SOURCE_BEFORE_PLANT_OUTCOMES.json"
OUT=D/"NWCA2016_2026_NPDES_PERMIT_EFF_ASOF2016_SOURCE_ONLY_NEAR10_FAR20.json"
A=D/"NWCA2011_2016_EPA_ECHO2026_POTW_ACTUAL_SURVEY_YEAR_PERMIT_SOURCE_QA.json"

def hashsite(s,n):
    return hashlib.sha256((str(s)+"|"+str(int(n))).encode()).hexdigest()

def load(root,dir):
    manifest=json.loads((root/SITE).read_text())
    geo=json.loads((root/OUT).read_text())
    previous=json.loads((root/A).read_text())
    if geo["original_2016_veg_SITE_source_SHA256"]!=sha(root/SITE) or geo["original_historical_year_permit_age_audit_SHA256"]!=sha(root/A):
        raise ValueError("Original registered year2016 NWCA source site and outfall SHA mismatch")
    if not geo["original_2016_historical_source_vegmetric_plant_H_ALL_NOT_READ"] or not geo["2011_historical_permit_source_failure_acknowledged_do_not_fit_old_2011"]:
        raise ValueError("Source not blind to biological H_ALL result or failed 2011 source must be excluded")
    frame=manifest["source_years"]["2016"]
    p=dir/"2016_siteinfo.csv"
    v=dir/"2016_vegmetrics.csv"
    if sha(p)!=frame["original_official_siteinfo_SHA256"] or sha(v)!=frame["original_official_vegetation_metrics_SHA256"]:
        raise ValueError("Official original EPA 2016 files modified after site source gate")
    ids=frame["source_siteinfo_ID_columns_resolved"]
    needed=("SITE_ID","VISIT_NO","LAT","LON","SITETYPE","WEIGHT","YEAR","STATE","ECO","WETCLS","HYDCLS")
    ss=pd.read_csv(p,usecols=list({ids[k] for k in needed}),low_memory=False)
    ss=ss.rename(columns={ids[k]:k for k in needed})
    ss.SITE_ID=ss.SITE_ID.astype("string").str.strip()
    ss.VISIT_NO=pd.to_numeric(ss.VISIT_NO,errors="coerce")
    ss.YEAR=pd.to_datetime(ss.YEAR,format="mixed",errors="coerce").dt.year
    for col in ("LAT","LON","WEIGHT"):ss[col]=pd.to_numeric(ss[col],errors="coerce")
    clean=(ss.SITETYPE.astype("string").str.upper().str.strip().eq("PROB")&
            ss.YEAR.eq(2016)&ss.VISIT_NO.eq(1)&ss.WEIGHT.gt(0)&
            ss.LAT.between(24,50)&ss.LON.between(-126,-66)&
            ss.STATE.notna()&ss.ECO.notna()&ss.WETCLS.notna()&ss.HYDCLS.notna())
    ss=ss.loc[clean].copy()
    ss["SITE_HASH"]=[hashsite(s,v) for s,v in zip(ss.SITE_ID,ss.VISIT_NO)]
    original_fullset=set(frame["source_site_hashes_only"])
    ss=ss.loc[ss.SITE_HASH.isin(original_fullset)].copy()
    if len(ss)!=966 or set(ss.SITE_HASH)!=original_fullset:
        raise ValueError("Original 966 plant-survey source opportunity frame changed")
    source=geo["2016_original_source_status_site_hash_receipts"]
    bysite={r["original_2016_NWCA_site_visit_hash"]:r for r in source}
    if len(bysite)!=212 or len(source)!=212:
        raise ValueError("Historical 2016 registered outfall source matched sample must be 50+162")
    ss=ss.loc[ss.SITE_HASH.isin(bysite)].copy()
    if len(ss)!=212 or set(ss.SITE_HASH)!=set(bysite):
        raise ValueError("Historical 2016 asof-2016 source 212 original site IDs changed")
    ss["NEAR"]=[int(bysite[s]["original_2016_2026_current_permit_asof2016_NEAR10_or_FAR20"]=="near10") for s in ss.SITE_HASH]
    ss["STATE_ECO"]=ss.STATE.astype("string").str.strip()+"|"+ss.ECO.astype("string").str.strip()
    if Counter(ss.NEAR)!=Counter({1:50,0:162}) or ss.STATE_ECO.nunique()!=23 or ss.STATE.nunique()<19:
        raise ValueError("Original year2016 source near vs far sample before outcomes changed")
    wanted=["SITE_ID","VISIT_NO","H_ALL","H_NAT","H_ALIEN","FQAI_ALL","FQAI_NAT"]
    vm=pd.read_csv(v,usecols=wanted,low_memory=False)
    vm.SITE_ID=vm.SITE_ID.astype("string").str.strip()
    vm.VISIT_NO=pd.to_numeric(vm.VISIT_NO,errors="coerce")
    if vm.duplicated(["SITE_ID","VISIT_NO"]).any():
        raise ValueError("Original vegetation metric source duplicates same site visit")
    ss=ss.merge(vm,on=["SITE_ID","VISIT_NO"],how="left",validate="one_to_one")
    if len(ss)!=212:raise ValueError("Historical source vegetation site merge changed size")
    for y in ("H_ALL","H_NAT","H_ALIEN","FQAI_ALL","FQAI_NAT"):ss[y]=pd.to_numeric(ss[y],errors="coerce")
    ss=ss.rename(columns={"H_ALL":"HALL","H_NAT":"HNAT","H_ALIEN":"HALIEN",
        "FQAI_ALL":"FQAIALL","FQAI_NAT":"FQAINAT","STATE":"PSTLCODE",
        "WETCLS":"WETCLSGRP","HYDCLS":"HYDCLS","LAT":"LATANALYS","WEIGHT":"WGT",
        "ECO":"AGECO9"})
    return ss,geo

def run(root,dir):
    ss,source=load(root,dir)
    out={
        "schema":"gesui_NWCA2016_ECHO2026_asof2016_POTW_source50near162far_HALL_SENSITIVITY_POSTHOC_v01",
        "utc":datetime.now(timezone.utc).isoformat(),
        "original_site_2016_ALL_VEG_SOURCE_SHA256":sha(root/SITE),
        "original_year2016_source_only_permit2026_near_far_SHA256":sha(root/OUT),
        "NPDES2026_2016_PERMIT_EFFECTIVE_not_observed_2016_actual_discharge":True,
        "2016_source_preoutcome_near10_n":50,"2016_source_preoutcome_far20_n":162,
        "2016_NWCA_original_state_eco_strata_n":23,
        "2011_historical_outfall_date_source_insufficient_no_biological_model":True,
        "NPDES_historical_date_audit_after_previous_2016_2021proxy_beta_was_known":True,
        "NPDES_river_wetland_direction_and_2016_TNTP_discharge_unmeasured":True
    }
    outcome={}
    for name in ("HALL","HNAT","HALIEN","FQAIALL","FQAINAT"):
        outcome[name]={
            "site_equal_OLS":model(ss,name,False),
            "observed_source_HALL_n":int(ss[name].notna().sum())
        }
        if name=="HALL":
            outcome[name]["original_NWCA_design_weight_WLS"]=model(ss,name,True)
            outcome[name]["near_raw_HALL"]=q(ss.loc[ss.NEAR==1,name])
            outcome[name]["far_raw_HALL"]=q(ss.loc[ss.NEAR==0,name])
    out["original_2016_state_ecoregion_NWCA_plant_diversity_by_yearvalid_REGISTERED_POTW_source"]=outcome
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-dir",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    data=run(root,a.source_dir)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(data,sort_keys=True,indent=2)+"\n")
    print(json.dumps({"2016":data["original_2016_state_ecoregion_NWCA_plant_diversity_by_yearvalid_REGISTERED_POTW_source"]["HALL"]},indent=2))
if __name__=="__main__":main()
