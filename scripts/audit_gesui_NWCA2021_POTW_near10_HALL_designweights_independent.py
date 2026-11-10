#!/usr/bin/env python3
"""Independent NumPy QA and design-weight sensitivity on strict WWTP source.

Never reselect original 568 NWCA plant/water sites or 2026 NPDES outfalls.
Recompute SAME 240 state+ecoregion matched source-only wetland H_ALL,
site-equal OLS and NWCA area-weight WLS cluster CI independently, then
posthoc predeclared 95/99% weight caps and 27 state-deletion leverage audit.
Neither modeled outcome is proof of actual 2021 effluent to wetland flow.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
import pandas as pd

D=Path("results/us_NWCA2021_plant_TNTP_20261110")
SOURCE=D/"EPA_NWCA2021_NPDES_SEWAGE_OUTFALL_568SITE_SOURCE_ONLY_GEOGRAPHY_GATE.json"
OLD=D/"EPA_NWCA2021_ORIGINAL_VEG_WATERCHEM_SITEONLY_GATE.json"
OUTCOME=D/"EPA_NWCA2021_NPDES_WETLAND_HALL_FQAI_AND_ACTUAL_NP_STRICT_NEAR10_FAR20_RESULTS.json"
FILES={"SITEINFO":"nwca21_siteinfo-data.csv","WATERCHEM":"nwca21_waterchem_wide-data.csv",
       "VEGMETRICS":"nwca21_vegmetrics-data.csv","HYDROSOURCE":"nwca21_hydrology_sources_wide-data.csv"}
ONE_POINT_NINE_FIVE=1.959963984540054

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for a in iter(lambda:f.read(1024*1024),b""):h.update(a)
    return h.hexdigest()

def original_site_hash(a,b):
    return hashlib.sha256((str(a)+"|"+str(int(b))).encode()).hexdigest()

def get_original_source_rows(root,source_dir):
    src=json.loads((root/SOURCE).read_text())
    gate=json.loads((root/OLD).read_text())
    old=json.loads((root/OUTCOME).read_text())
    if src["source_precommitted_GEO_eligibility_gate"]["permitted_2021_outfall_distance_association_source_support_PASS"] is not True:
        raise ValueError("Original 2021 WWTP source geographically insufficient")
    if old["original_568_sites_ECHO_2026_source_only_geography_SHA256"]!=sha(root/SOURCE):
        raise ValueError("Original 240 same-state source model from different source")
    if src["original_preoutcome_EPA_NWCA2021_568joint_SITE_ID_SHA256"]!=sha(root/OLD):
        raise ValueError("WWTP source geometry not on original preplant 568 sites")
    for k,name in FILES.items():
        if sha(source_dir/name)!=gate["official_EPA_original_source_csv_receipts"][k]["official_source_sha256"]:
            raise ValueError("Official NWCA source file SHA drift")
    geom=src["all_original_568_NWCA_source_site_to_sewage_source_distance_receipts_without_plant_or_TNTP_values"]["strict"]
    dgeo={r["NWCA2021_source_site_visit_SHA256"]:r for r in geom}
    if len(dgeo)!=568:raise ValueError("Original 568 same-source wetland geography modified")
    # Read original official SITEINFO and VEGMETRICS independently, not calling
    # the original source/response parser from primary calculation.
    site=pd.read_csv(source_dir/FILES["SITEINFO"],low_memory=False)
    veg=pd.read_csv(source_dir/FILES["VEGMETRICS"],low_memory=False)
    site=site.rename(columns=lambda s:"".join(c for c in s.upper() if c.isalnum()))
    veg=veg.rename(columns=lambda s:"".join(c for c in s.upper() if c.isalnum()))
    site["SITEID"]=site.SITEID.astype("string").str.strip()
    veg["SITEID"]=veg.SITEID.astype("string").str.strip()
    for c in ("VISITNO","YEAR","WGTTPCORE","LATANALYS","LONANALYS"):
        site[c]=pd.to_numeric(site[c],errors="coerce")
    veg["VISITNO"]=pd.to_numeric(veg["VISITNO"],errors="coerce")
    good=(site.SITETYPE.astype("string").str.upper().eq("PROB").fillna(False)&
          site.YEAR.eq(2021)&site.VISITNO.eq(1)&site.WGTTPCORE.gt(0)&
          site.LATANALYS.between(24,50)&site.LONANALYS.between(-126,-66)&
          site.STUDY.astype("string").str.contains("NWCA",case=False,na=False))
    site=site.loc[good].copy()
    if len(site)!=920:raise ValueError("NWCA source lost probability sample 920")
    site["SITE_SHA"]=[original_site_hash(x,y) for x,y in zip(site.SITEID,site.VISITNO)]
    site=site[site.SITE_SHA.isin(dgeo)].copy()
    if len(site)!=568 or set(site.SITE_SHA)!=set(dgeo):
        raise ValueError("Original source set of 568 joint vegetation and water sites changed")
    veg=veg[veg.VISITNO.eq(1)].copy()
    if "SAMPLETYPE" in veg and veg.SAMPLETYPE.astype("string").str.upper().eq("PLANT").any():
        veg=veg[veg.SAMPLETYPE.astype("string").str.upper()=="PLANT"].copy()
    if veg.duplicated(["SITEID","VISITNO"]).any():
        raise ValueError("Official 2021 site×vegetation has duplicate protocols")
    joined=site.merge(veg[["SITEID","VISITNO","HALL","HNAT","HALIEN","FQAIALL","FQAINAT"]],
                      on=["SITEID","VISITNO"],how="left",validate="one_to_one")
    if joined.HALL.isna().any():
        raise ValueError("Original source-matched 568 accepted H_ALL missing")
    assert len(joined)==568
    nearset={s for s,g in dgeo.items() if g["outfall_rows_within_10km"]>0}
    farset={s for s,g in dgeo.items() if g["nearest_outfall_geodesic_km_within25km_else_null"] is None or
                                      g["nearest_outfall_geodesic_km_within25km_else_null"]>=20}
    if len(nearset)!=128 or len(farset)!=311 or nearset&farset:
        raise ValueError("Original source-date-plausible POTW near10/far20 classification changed")
    joined["NEAR"]=joined.SITE_SHA.isin(nearset).astype(int)
    joined["FAR"]=joined.SITE_SHA.isin(farset).astype(int)
    joined["STATE_ECO"]=joined.PSTLCODE.astype("string").str.strip()+"|"+joined.AGECO9.astype("string").str.strip()
    strata=set(joined.loc[joined.NEAR.eq(1),"STATE_ECO"])&set(joined.loc[joined.FAR.eq(1),"STATE_ECO"])
    source_selected=joined.loc[joined.STATE_ECO.isin(strata)&joined[["NEAR","FAR"]].sum(axis=1).eq(1)].copy()
    if len(source_selected)!=240 or Counter(source_selected.NEAR)!=Counter({1:103,0:137}):
        raise ValueError("Original source matched 27 state/32 area 240 site sample changed")
    if source_selected.PSTLCODE.nunique()!=27 or source_selected.STATE_ECO.nunique()!=32:
        raise ValueError("Original state/eco matched site area units changed")
    if source_selected.SITE_SHA.duplicated().any():raise ValueError("Original source site ID duplicates")
    for k in ("HALL","HNAT","HALIEN","FQAIALL","FQAINAT"):
        source_selected[k]=pd.to_numeric(source_selected[k],errors="coerce")
    return source_selected,old

def design(frame):
    n=len(frame)
    pieces=[np.ones(n),frame.NEAR.to_numpy(dtype=float)]
    cols=["INTERCEPT","NEAR"]
    for field in ("STATE_ECO","WETCLSGRP","HYDCLS"):
        vals=frame[field].astype("string").astype(str)
        levels=sorted(vals.unique())
        for level in levels[1:]:
            pieces.append((vals==level).to_numpy(dtype=float))
            cols.append(field+":"+level)
    pieces.append(pd.to_numeric(frame.LATANALYS,errors="coerce").to_numpy(dtype=float))
    cols.append("LATANALYS")
    X=np.column_stack(pieces)
    if not np.isfinite(X).all():raise ValueError("Independent matched-source NWCA source/covariates not numeric")
    return X,cols

def independent_fit(frame,y,weights=None):
    n=frame[y].notna().sum()
    df=frame.loc[frame[y].notna()].copy()
    if n<80:return None
    mat,cols=design(df)
    yy=df[y].to_numpy(dtype=float)
    if weights is None:
        w=np.ones(n)
    else:
        w=np.asarray(weights.loc[df.index] if hasattr(weights,"loc") else weights,dtype=float)
        if len(w)!=n or np.any(w<=0) or not np.isfinite(w).all():
            raise ValueError("Invalid original NWCA survey area weight under same 240 sites")
        w=w/np.median(w)
    xx=mat*np.sqrt(w)[:,None]; ystar=yy*np.sqrt(w)
    beta=np.linalg.lstsq(xx,ystar,rcond=None)[0]
    resid=yy-mat@beta
    bread=np.linalg.pinv(mat.T@(w[:,None]*mat),rcond=1e-12)
    clusters=sorted(df.PSTLCODE.astype(str).unique())
    meat=np.zeros((len(cols),len(cols)))
    for state in clusters:
        subset=df.PSTLCODE.astype(str).to_numpy()==state
        influence=mat[subset].T@(w[subset]*resid[subset])
        meat+=np.outer(influence,influence)
    G=len(clusters)
    # Statsmodels' cluster small-sample correction uses the original number
    # of matrix columns (not reduced design rank), matching earlier independent
    # 2021 HALL/N/P source audits with redundant wetland-class dummies.
    K=len(cols)
    fac=(G/(G-1))*((n-1)/(n-K))
    cv=bread@meat@bread*fac
    nearbeta=float(beta[1])
    se=math.sqrt(max(0.,float(cv[1,1])))
    return {
        "n_original_source_sites":int(n),
        "original_US_state_clusters":int(G),
        "original_state_ecoregion_strata_n":int(df.STATE_ECO.nunique()),
        "design_matrix_columns":K,
        "design_rank":int(np.linalg.matrix_rank(mat)),
        "independent_numpy_near10_minus_far20_coefficient":nearbeta,
        "independent_state_cluster_sandwich95":[nearbeta-ONE_POINT_NINE_FIVE*se,nearbeta+ONE_POINT_NINE_FIVE*se],
    }

def qd(x):
    y=np.asarray(x,dtype=float)
    return {"n":len(y),"median":float(np.median(y)),
      "p90":float(np.quantile(y,.9)),"p95":float(np.quantile(y,.95)),
      "p99":float(np.quantile(y,.99)),
      "maximum":float(max(y)),"min":float(min(y)),
      "effective_N_Kish":float(y.sum()**2/np.square(y).sum()),
      "top_1_site_fraction_of_total_weight":float(max(y)/y.sum()),
      "top_5_site_fraction_of_total_weight":float(sum(sorted(y,reverse=True)[:5])/y.sum())}

def run(root,source_dir):
    df,original=get_original_source_rows(root,source_dir)
    weights=df.WGTTPCORE.astype(float)
    out={
        "schema":"gesui_EPA_NWCA2021_27state_240wetland_POTW_registered_2026_plant_Shannon_area_weight_influence_independent_numeric_QA_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "source_pre_original_568_SITE_ID_POTW_10km_geometry_sha256":sha(root/SOURCE),
        "source_original_240site_first_biology_models_sha256":sha(root/OUTCOME),
        "original_2021_NWCA_568source_2021_plant_before_NP_SHA256":sha(root/OLD),
        "exact_source_240_site_27state_32state_ecoregion_HALL_matched":True,
        "original_27_state_same_source_site_weight_distribution":{
            "all240":qd(weights),
            "near103":qd(df.loc[df.NEAR.eq(1),"WGTTPCORE"].to_numpy()),
            "far137":qd(df.loc[df.NEAR.eq(0),"WGTTPCORE"].to_numpy())},
        "posthoc_sensitivity_after_seeing_source_weighted_beta_plus0p18_and_unweighted_minus0p12":True,
        "source_permit2021_2026location_NOT_actual_2021_wastewater_delivery":True,
        "independent_main_original_240_US_site_equal_and_area_weighted_models":{},
        "fixed_95_99_percentage_clipped_NWCA_area_weight_sensitivity":{},
        "source_state_leave_one_out_OLS_and_WLS_HALL_influence":{}
    }
    outcome=original["original_wetland_standardized_plant_WWTP_geodesic_near10_far20_country_model"]
    for response in ("HALL","HNAT","HALIEN","FQAIALL","FQAINAT"):
        result=independent_fit(df,response,None)
        old=outcome[response]["source_outfall_near10_beta_and_state_cluster95"]
        if abs(result["independent_numpy_near10_minus_far20_coefficient"]-old["beta_near10_minus_far20_in_reported_response_units"])>1e-7:
            raise ValueError("Independent source original 2021 "+response+" β differs")
        for a,b in zip(result["independent_state_cluster_sandwich95"],old["original_US_state_cluster95"]):
            if abs(a-b)>1e-5:raise ValueError("Independent source original "+response+" clustered CI differs: "+repr([a,b]))
        out["independent_main_original_240_US_site_equal_and_area_weighted_models"][response+"_site_equal"]=result
    wr=independent_fit(df,"HALL",weights)
    old=outcome["HALL_design_weighted_sensitivity"]["source_outfall_near10_beta_and_state_cluster95"]
    if abs(wr["independent_numpy_near10_minus_far20_coefficient"]-old["beta_near10_minus_far20_in_reported_response_units"])>1e-7:
        raise ValueError("Independent area-weighted source Shannon β changed")
    for a,b in zip(wr["independent_state_cluster_sandwich95"],old["original_US_state_cluster95"]):
        if abs(a-b)>1e-5:raise ValueError("Independent area-weighted source 95CI changed")
    out["independent_main_original_240_US_site_equal_and_area_weighted_models"]["HALL_original_survey_area_weighted"]=wr
    for p in (0.95,0.99):
        cap=float(np.quantile(weights,p))
        clipped=weights.clip(upper=cap)
        result=independent_fit(df,"HALL",clipped)
        out["fixed_95_99_percentage_clipped_NWCA_area_weight_sensitivity"][f"clip_{int(p*100)}pct"]= {
            "source_raw_weight_clip":cap,
            "sites_originally_above_clip":int((weights>cap).sum()),
            "source_WWTP_near_minus_far_HALL_weighted_coefficient_and_state_cluster95":result,
            "posthoc_not_original_survey_2021_national_design_estimate":True
        }
    leave={}
    for state in sorted(df.PSTLCODE.astype(str).unique()):
        sub=df.loc[df.PSTLCODE.astype(str)!=state]
        if sub.NEAR.nunique()!=2:raise ValueError("State-deletion source lost near/far")
        l=independent_fit(sub,"HALL",None)
        w=independent_fit(sub,"HALL",sub.WGTTPCORE.astype(float))
        leave[str(state)]={
            "source_original_site_rows":len(sub),
            "site_equal_beta":l["independent_numpy_near10_minus_far20_coefficient"],
            "area_weighted_beta":w["independent_numpy_near10_minus_far20_coefficient"]
        }
    out["source_state_leave_one_out_OLS_and_WLS_HALL_influence"]={
        "number_US_states_deleted":len(leave),
        "site_equal_beta_positive_n":sum(x["site_equal_beta"]>0 for x in leave.values()),
        "area_weighted_beta_positive_n":sum(x["area_weighted_beta"]>0 for x in leave.values()),
        "site_equal_beta_minmax":[min(x["site_equal_beta"] for x in leave.values()),max(x["site_equal_beta"] for x in leave.values())],
        "area_weighted_beta_minmax":[min(x["area_weighted_beta"] for x in leave.values()),max(x["area_weighted_beta"] for x in leave.values())],
        "whole_state_omission_original_source_results":leave,
        "country_state_influence_not_outofstate_validation":True
    }
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-dir",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    out=run(root,a.source_dir)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({
        "weights":out["original_27_state_same_source_site_weight_distribution"],
        "verified":out["independent_main_original_240_US_site_equal_and_area_weighted_models"],
        "clipping":out["fixed_95_99_percentage_clipped_NWCA_area_weight_sensitivity"],
        "state_delete":{k:v for k,v in out["source_state_leave_one_out_OLS_and_WLS_HALL_influence"].items() if k!="whole_state_omission_original_source_results"}
    },indent=2))
if __name__=="__main__":
    main()
