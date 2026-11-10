#!/usr/bin/env python3
"""New independently SOURCE-frozen NWCA2011/2016 plant H_ALL across POTW distance.

Original source selection:
- 2011=967 and 2016=966 probability first-visit vegetation survey site IDs,
- same 2026 EPA ECHO exact source SHA and 2021-permit *proxy* (NOT historic discharge),
- independent source within state×ecoregion near10 and far20 statuses frozen
  before reading H_ALL species-level or 2011/2016 quality responses.

This is temporal comparison of DIFFERENT NWCA probability surveys, posthoc
after 2021 weighted positive and site-equal negative. No sewage causality.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path
from datetime import datetime,timezone
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from scipy.stats import spearmanr

D=Path("results/us_NWCA2021_plant_TNTP_20261110")
SITE=D/"NWCA2011_2016_ORIGINAL_SITE_VEG_SOURCE_BEFORE_PLANT_OUTCOMES.json"
GEO=D/"NWCA2011_2016_POTW2026_SOURCE_ONLY_NEAR10_FAR20_2CYCLE.json"
YLIST=("H_ALL","H_NAT","H_ALIEN","FQAI_ALL","FQAI_NAT")
MIN_A=50
MIN_B=60
MIN_STATES=10
MIN_REGIONS=20

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(2**20),b""):h.update(b)
    return h.hexdigest()
def sitehash(a,b):return hashlib.sha256((str(a)+"|"+str(int(b))).encode()).hexdigest()
def source_frame(year,root,dir,site_manifest,geo_manifest):
    meta=site_manifest["source_years"][year]
    geo=geo_manifest["source_2011_2016_years"][year]
    if meta["source_gate_STATUS"]!="SOURCE_SITE_VEG_PASS":
        return None,"ORIGINAL_SITE_VEG_METRICS_NOT_SUPPORTED"
    if not geo["source_precommitted_independent_year_gate"]["may_open_original_NWCA_H_ALL_and_fit_2011_2016_comparison"]:
        return None,"ORIGINAL_SOURCE_2011_2016_POTW_GEOGRAPHY_UNSUPPORTED"
    site=dir/(year+"_siteinfo.csv")
    veg=dir/(year+"_vegmetrics.csv")
    if sha(site)!=meta["original_official_siteinfo_SHA256"] or sha(veg)!=meta["original_official_vegetation_metrics_SHA256"]:
        raise ValueError("EPA ORIGINAL survey source changed after PRE-H_ALL site and geometry gates")
    a=meta["source_siteinfo_ID_columns_resolved"]
    needed=("SITE_ID","VISIT_NO","LAT","LON","SITETYPE","WEIGHT","YEAR","STATE","ECO","WETCLS","HYDCLS")
    raw=pd.read_csv(site,usecols=list({a[k] for k in needed}),low_memory=False)
    raw=raw.rename(columns={a[k]:k for k in needed})
    raw["SITE_ID"]=raw.SITE_ID.astype("string").str.strip()
    raw["VISIT_NO"]=pd.to_numeric(raw.VISIT_NO,errors="coerce")
    raw["YEAR"]=pd.to_datetime(raw.YEAR,format="mixed",errors="coerce").dt.year
    for k in ("WEIGHT","LAT","LON"):raw[k]=pd.to_numeric(raw[k],errors="coerce")
    mask=(raw.SITETYPE.astype("string").str.upper().str.strip().eq("PROB")&
          raw.WEIGHT.gt(0)&raw.YEAR.eq(int(year))&raw.VISIT_NO.eq(1)&
          raw.LAT.between(24,50)&raw.LON.between(-126,-66)&
          raw.STATE.notna()&raw.ECO.notna()&
          raw.WETCLS.notna()&raw.HYDCLS.notna())
    raw=raw.loc[mask].copy()
    raw["SITE_HASH"]=[sitehash(x,y) for x,y in zip(raw.SITE_ID,raw.VISIT_NO)]
    geo_rows=geo["original_historical_source_matched_site_date_weight_and_outfall_receipts"]
    mapping={r["site_SHA256"]:r for r in geo_rows}
    if len(mapping)!=len(geo_rows):
        raise ValueError("Original historic registered sewage near/far source site duplicated")
    raw=raw.loc[raw.SITE_HASH.isin(mapping)].copy()
    if len(raw)!=len(mapping) or set(raw.SITE_HASH)!=set(mapping):
        raise ValueError("Original pre-response historical source near/far site matched receipts have drifted")
    if raw.duplicated(["SITE_ID","VISIT_NO"]).any():raise ValueError("Repeated original probability wetlands")
    raw["NEAR"]=[1 if mapping[x]["source_group"]=="near10" else 0 for x in raw.SITE_HASH]
    raw["STATE_ECO"]=raw.STATE.astype("string").str.strip()+"|"+raw.ECO.astype("string").str.strip()
    if any(v not in ("near10","far20") for v in (r["source_group"] for r in geo_rows)):
        raise ValueError("Original historical source status contains unallocated intermediate wetland")
    groups=Counter(raw.NEAR)
    source_qa=geo["original_shared_state_and_ecoregion_near10_far20_source_n"]
    if groups[1]!=source_qa["near10"] or groups[0]!=source_qa["far20"]:
        raise ValueError("Historical source study original near10/far20 split differs")
    for col in ("STATE","ECO","WETCLS","HYDCLS"):
        raw[col]=raw[col].astype("string").str.strip()
    columns={k:v for k,v in meta["source_metric_ID_columns_resolved"].items() if v}
    if not meta["metric_has_2011or2016_H_ALL_header"]:
        raise ValueError("Original EPA H_ALL metric absent; cannot substitute another metric")
    wanted=list({*YLIST,*columns.values()})
    b=pd.read_csv(veg,usecols=wanted,low_memory=False)
    b=b.rename(columns={v:k for k,v in columns.items()})
    b["SITE_ID"]=b.SITE_ID.astype("string").str.strip()
    b["VISIT_NO"]=pd.to_numeric(b.VISIT_NO,errors="coerce")
    if b.duplicated(["SITE_ID","VISIT_NO"]).any():
        raise ValueError("Original source vegetation metrics duplicate firstvisit site; no choosing favorable H_ALL record")
    df=raw.merge(b,on=["SITE_ID","VISIT_NO"],how="left",validate="one_to_one")
    if len(df)!=len(raw):
        raise ValueError("H_ALL source merge expanded or shrank frozen source site cohort")
    if not df.H_ALL.notna().any():raise ValueError("Original H_ALL field all missing")
    for k in YLIST:df[k]=pd.to_numeric(df[k],errors="coerce")
    df=df.rename(columns={"LAT":"LATANALYS","WEIGHT":"WGT","STATE":"PSTLCODE",
                          "WETCLS":"WETCLSGRP","HYDCLS":"HYDCLS","ECO":"AGECO9",
                          "H_ALL":"HALL","H_NAT":"HNAT","H_ALIEN":"HALIEN",
                          "FQAI_ALL":"FQAIALL","FQAI_NAT":"FQAINAT"})
    return df,None

def model(df,outcome,weighted):
    sample=df.loc[df[outcome].notna() & np.isfinite(df[outcome])].copy()
    actual=Counter(sample.NEAR)
    if actual[1]<MIN_A or actual[0]<MIN_B or sample.PSTLCODE.nunique()<MIN_STATES:
        return {"status":"RESP_OUTCOME_NONMISSING_SUPPORT_FAIL", "n":len(sample),
                "near":actual[1],"far":actual[0],"states":int(sample.PSTLCODE.nunique())}
    # Prespecified source stratum first; a missing H_ALL can eliminate
    # one near/far side of a stratum. Drop whole stratum only if it lacks
    # any source response for both groups, not based on observed sign.
    matched=[key for key,g in sample.groupby("STATE_ECO") if g.NEAR.nunique()==2]
    dropped=len(sample)-int(sample.STATE_ECO.isin(matched).sum())
    sample=sample.loc[sample.STATE_ECO.isin(matched)].copy()
    nn=Counter(sample.NEAR)
    if nn[1]<MIN_A or nn[0]<MIN_B or len(matched)<MIN_REGIONS:
        return {"status":"RESP_BOTH_GROUP_AREA_SUPPORT_FAIL", "source_near":nn[1],"source_far":nn[0],
                "source_regions":len(matched),"dropped_n":dropped}
    formula=outcome+"~NEAR+C(STATE_ECO)+C(WETCLSGRP)+C(HYDCLS)+LATANALYS"
    if weighted:
        w=sample.WGT.to_numpy(dtype=float)
        w=w/np.median(w)
        fit=smf.wls(formula,data=sample,weights=w).fit(cov_type="cluster",cov_kwds={"groups":sample.PSTLCODE})
    else:
        fit=smf.ols(formula,data=sample).fit(cov_type="cluster",cov_kwds={"groups":sample.PSTLCODE})
    b=float(fit.params["NEAR"])
    lo,hi=list(map(float,fit.conf_int().loc["NEAR"]))
    kish=(np.sum(sample.WGT)**2)/(np.sum(sample.WGT**2))
    return {
        "status":"ESTIMATED_NONCAUSAL_SOURCE_PROXIMITY",
        "original_2021_compatible_NWCA_country_state_ecoregion_model":formula,
        "site_equal_vs_design_weighted":"design_weighted_WLS" if weighted else "site_equal_OLS",
        "n_source_reporting_plant_Hall_sites":len(sample),
        "near10_plant_source_sites":nn[1],
        "far20_plant_source_sites":nn[0],
        "original_US_state_clusters":int(sample.PSTLCODE.nunique()),
        "original_matched_STATE_ECO_strata":len(matched),
        "H_ALL_nonmissing_response_area_drops_n":dropped,
        "NWCA_survey_weight_effective_sample_size_Kish":float(kish),
        "near_minus_far_beta":b,
        "state_clustered_descriptive95":[lo,hi],
        "nominal_posthoc_p":float(fit.pvalues["NEAR"]),
        "R2":float(fit.rsquared),
        "survey_weight_is_not_randomized_2021_sewage_treatment":True
    }

def q(x):
    v=pd.to_numeric(x,errors="coerce").dropna().to_numpy(dtype=float)
    return {"n":len(v),"mean":float(v.mean()),"median":float(np.median(v)),
            "q25":float(np.quantile(v,.25)),"q75":float(np.quantile(v,.75))} if len(v) else None

def run(root,source_dir):
    site=json.loads((root/SITE).read_text())
    geo=json.loads((root/GEO).read_text())
    if geo["original_2011_2016_vegetation_source_only_manifest_SHA256"]!=sha(root/SITE):
        raise ValueError("Original independent 2011/2016 pre-plant source SHA mismatched")
    if not geo["same_original_2021_strict_sewage_source_filter_reused_without_tuning"]:
        raise ValueError("Historic registered outfall geometry differs from 2021 original source")
    result={
        "schema":"gesui_NWCA_2011_2016_vs2021_registered_POTW_2026_NEAR10_FAR20_HALL_temporal_comparison_POSTHOC_v01",
        "utc":datetime.now(timezone.utc).isoformat(),
        "original_pre_HALL_EPA_2011_2016_veg_site_source_SHA256":sha(root/SITE),
        "original_pre_HALL_2011_2016_POTW_geo_source_SHA256":sha(root/GEO),
        "original_ECHO2026_sewage_source_snapshot_SHA256":geo["original_2026_NPDES_official_ECHO_2021permit_PROXY_ZIP_SHA256"],
        "2026_outfall_is_not_verified_2011_2016_wastewater_delivery":True,
        "years":{},
        "2021_reference_only_not_reestimated":{
            "original_year":2021,"site_n":240,"near10":103,"far20":137,"states":27,
            "state_eco_strata":32,"site_equal_beta":-0.12043442865445983,
            "site_equal95":[-0.2505774852985684,0.00970862798964875],
            "NWCA_survey_weighted_beta":0.18062300867624317,
            "NWCA_survey_weighted95":[0.010378490906395199,0.35086752644609115]
        },
        "old_2021_result_already_known_before_2011_2016_source_requested":True
    }
    for y in ("2011","2016"):
        df,stop=source_frame(y,root,source_dir,site,geo)
        if stop:
            result["years"][y]={"status":"SOURCE_GATE_STOP","reason":stop}
            continue
        groups=Counter(df.NEAR)
        metrics={"source_fixed_n":len(df),"source_near":groups[1],"source_far":groups[0],
                 "source_states":int(df.PSTLCODE.nunique()),
                 "source_state_eco_strata":int(df.STATE_ECO.nunique()),
                 "actual_original_plants_HALL_reported_sites_n":int(df.HALL.notna().sum()),
                 "source_near_HALL_nonmissing":int((df.NEAR.eq(1)&df.HALL.notna()).sum()),
                 "source_far_HALL_nonmissing":int((df.NEAR.eq(0)&df.HALL.notna()).sum()),
                 "source_year_original_plant_group_na":{v:int(df[v].isna().sum()) for v in ("HALL","HNAT","HALIEN","FQAIALL","FQAINAT")},
                 "outcomes":{},
                 "matched_original_HALL_raw_descriptive":{
                    "near10":q(df.loc[df.NEAR.eq(1),"HALL"]),
                    "far20":q(df.loc[df.NEAR.eq(0),"HALL"])
                 }}
        for v in ("HALL","HNAT","HALIEN","FQAIALL","FQAINAT"):
            metrics["outcomes"][v]={
                "site_equal_OLS":model(df,v,False),
                "site_design_weighted_WLS":model(df,v,True) if v=="HALL" else None
            }
        result["years"][y]=metrics
    return result

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--source-dir",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    a=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    out=run(root,a.source_dir)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n")
    print(json.dumps({
        "year_results":{y:{
           "source_n":v.get("source_fixed_n"),"validHALL":v.get("actual_original_plants_HALL_reported_sites_n"),
           "near_far":[v.get("source_near"),v.get("source_far")],
           "site_equal":v.get("outcomes",{}).get("HALL",{}).get("site_equal_OLS"),
           "design_weighted":v.get("outcomes",{}).get("HALL",{}).get("site_design_weighted_WLS")
        } for y,v in out["years"].items()},
        "original_2021":out["2021_reference_only_not_reestimated"]},indent=2))
if __name__=="__main__":main()
