#!/usr/bin/env python3
"""Exploratory region-mean EEA2016 benthic EQR: water pointN vs soil N inputs.

Uses ONLY the original already-verified 2016 region-mean EQR archive and
the fixed 2016 source-first folds; no new sites or alternative years.
Original 2018 water x 2017-19 land 2x2 STOP is not reopened.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from pathlib import Path
from collections import defaultdict
from datetime import datetime,timezone

import geopandas as gpd
import numpy as np
import pandas as pd

BASE=Path("results/eea2016_thirdyear_source_20261111")
RESULT=BASE/"2016_EEA_RIVER_BENTHIC_EQR_JRC_GREEN_MODEL_POINTNP_NUTS2_HELDOUT_RESULTS.json"
FROZEN=BASE/"2016_EEA_BENTHIC_JRC_NUTS2_FROZEN_5FOLDS_NO_EQR.json"
ORIGINAL=BASE/"2016_EEA_BENTHIC_GISCO_NUTS2_GREEN_POINT_NP_SOURCE_GATE_BEFORE_EQR.json"
GPKG_SHA="7a233d3d683dd48241b4a9e5a9ef6c5f3e294008934c5bc7f1a465eb5cbf9c0b"
ALPHA=10.
NUMERIC_B0=("lat","lon","logarea","log_source_n")
D_VARS=("atmos","mineral","organic")
BOOT=3999
SEED=20261119

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()

def safe(v):
    if isinstance(v,dict):return {k:safe(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [safe(x) for x in v]
    if isinstance(v,np.integer):return int(v)
    if isinstance(v,(float,np.floating)):return float(v) if math.isfinite(v) else None
    return v

def read_source(source_csv,source_report,gate_json,gisco):
    j=json.loads(Path(source_report).read_text())
    gate=json.loads(Path(gate_json).read_text())
    source_gate=gate.get("results",{})
    if not source_gate.get("source_only_GO_for_aquatic_only_next_stage"):
        raise ValueError("Independent EEA2016 water-only source gate is STOP")
    if source_gate.get("source_complete_NUTS2_shared_n")!=73 or source_gate.get("source_shared_country_n")!=12:
        raise ValueError("2016 source cohort n=73 and 12 country lock changed")
    if sha(source_csv)!=j.get("csv_sha256") or not j.get("NO_EEA_EQR_EIV_GBIF_BIOLOGICAL_VALUES_READ"):
        raise ValueError("JRC 269 source original digest/source-only state changed")
    if sha(gisco)!=GPKG_SHA:raise ValueError("Official GISCO 2016 GPKG digest differs from original")
    src=pd.read_csv(source_csv,dtype={"NUTS2_country_id":str})
    all_fields=["NUTS2_country_id","water_point_N_t","land_atmos_N_mg_m2",
                "land_mineral_N_kg_ha","land_organic_N_kg_ha"]
    if not set(all_fields).issubset(src.columns):
        raise ValueError("Source CSV missing essential source fields")
    src=src[all_fields].copy()
    for k in all_fields[1:]:
        src[k]=pd.to_numeric(src[k],errors="coerce")
    if src.NUTS2_country_id.duplicated().any() or src[all_fields[1:]].isna().any().any():
        # The full input includes any source incomplete rows; strictly filter earlier
        src=src.dropna(subset=all_fields[1:]).copy()
    if len(src)<269:raise ValueError("Fixed 269 six-source records incomplete")
    geo=gpd.read_file(gisco)
    if "LEVL_CODE" in geo:
        geo=geo.loc[geo.LEVL_CODE.astype(str).eq("2")]
    geo=geo.to_crs(3035)
    geo["source_id"]=geo.CNTR_CODE.astype(str).str.upper()+"|"+geo.NUTS_ID.astype(str).str.upper()
    if geo.source_id.duplicated().any():raise ValueError("NUTS2 2016 polygon duplicate")
    cent=gpd.GeoSeries(geo.geometry.centroid,crs=3035).to_crs(4326)
    geo["lat"]=cent.y.to_numpy()
    geo["lon"]=cent.x.to_numpy()
    geodict=geo.set_index("source_id")[["lat","lon"]]
    return src.set_index("NUTS2_country_id"),geodict,source_gate

def frozen_region_data(source,geodict,source_gate,compatibility):
    original=json.loads(RESULT.read_text())
    frozen=json.loads(FROZEN.read_text())
    if original["original_original5fold_SOURCE_before_EQR_SHA256"]!=sha(FROZEN):
        raise ValueError("Earlier 2016 EQR numeric analysis not linked to true frozen folds")
    if original["original_2016_sites_GISCO2016_JRC_modelled_SOURCE_before_EQR_SHA256"]!=sha(ORIGINAL):
        raise ValueError("Earlier 2016 numeric outcomes not linked to fixed GISCO SOURCE")
    c=original["results"]["benthic_invertebrates"]
    if c["precommitted_numeric_EQR_support_gate_PASS"] is not True:
        raise ValueError("Old numeric original valid-source gate was STOP")
    raw=c["original_NUTS2_disjoint_holdout_Ridge10"]["source_JRC_green_region_sampled_EEA_EQR_summary"]
    bas=frozen["source_biological_classes"]["benthic_invertebrates"]
    pre={z["original_GISCO_NUTS2"]:z for z in bas["original_NUTS2_region_source_inventory"]}
    regions=[]
    for z in raw:
        region=str(z["nuts2"])
        country=str(z["country"])
        k=country+"|"+region
        if region not in pre:raise ValueError("Original frozen 2016 BQE geographic key changed")
        if k not in source.index:continue
        if k not in geodict.index:raise ValueError("Official 2016 polygon centroid missing")
        s=source.loc[k]
        source_pointN=float(pre[region]["source_NUTS2_GREEN_pointN_t"])
        if abs(float(s.water_point_N_t)-source_pointN)>1e-8*max(1,source_pointN):
            raise ValueError("Two independent JRC model source pointN registers differ")
        fold=int(bas["fixed_region_to_fivefold_source_only"][region])
        if fold!=int(z["source_frozen_test_fold"]):raise ValueError("Original region fold shifted")
        area=float(pre[region]["source_NUTS2_area_km2"])
        n=int(pre[region]["source_original_num_2016_EEA_sites"])
        y=float(z["mean_EQR_original_numeric"])
        if n<1 or area<=0 or not 0<=y<=1 or not math.isfinite(y):
            raise ValueError("Original 2016 region mean EQR or preoutcome region source invalid")
        regions.append({
            "country":country,"nuts2":region,"source_id":k,
            "fold":fold,"n_source":n,
            "log_source_n":math.log1p(n),"logarea":math.log1p(area),
            "lat":float(geodict.at[k,"lat"]),"lon":float(geodict.at[k,"lon"]),
            "pointN":math.log1p(source_pointN/area),
            "atmos":math.log1p(float(s.land_atmos_N_mg_m2)),
            "mineral":math.log1p(float(s.land_mineral_N_kg_ha)),
            "organic":math.log1p(float(s.land_organic_N_kg_ha)),
            "EQR":y,
            "original_valid_site_n":int(z["n_original_2016_BQE_sites"])
        })
    d=pd.DataFrame(regions)
    if d.nuts2.duplicated().any() or len(d)<55 or d.country.nunique()<8:
        raise ValueError("New regional mean EQR quality gate has fewer independent regions/countries than frozen source eligibility")
    check=compatibility.get("results",{})
    if not check.get("new_aquatic_only_four_model_numeric_GO",False):
        raise ValueError("Original valid EQR x soilN compatibility source/quality gate STOP")
    original_key_set=set(check["region_ids"])
    if set(d.source_id)!=original_key_set or len(d)!=check["overlap_region_n"]:
        raise ValueError("Cohort changed after pre-new-model numeric-availability ID gate")
    if d.country.nunique()!=check["overlap_country_n"]:
        raise ValueError("Country support changed after numeric-availability ID gate")
    return d,original

def build_design(train,test,vars):
    fields=list(NUMERIC_B0)+list(vars)
    a=train[fields].to_numpy(float)
    b=test[fields].to_numpy(float)
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("Original source numeric values missing")
    mu=a.mean(axis=0)
    sd=a.std(axis=0,ddof=0)
    sd[sd<1e-12]=1
    a=(a-mu)/sd
    b=(b-mu)/sd
    cats=sorted(train.country.unique())
    enc={c:i for i,c in enumerate(cats)}
    ua=np.zeros((len(train),len(cats)))
    ub=np.zeros((len(test),len(cats)))
    for i,c in enumerate(train.country):ua[i,enc[c]]=1
    for i,c in enumerate(test.country):
        if c in enc:ub[i,enc[c]]=1
    return np.column_stack([a,ua]),np.column_stack([b,ub])

def fit_predict(train,test,extra):
    a,b=build_design(train,test,extra)
    y=train.EQR.to_numpy(float)
    if len(train)<40:raise ValueError("Insufficient regional training observations")
    xbar=a.mean(axis=0); ybar=y.mean()
    # Equivalent to Ridge(alpha=10, fit_intercept=True) on the supplied matrix.
    beta=np.linalg.solve((a-xbar).T@(a-xbar)+ALPHA*np.eye(a.shape[1]),
                         (a-xbar).T@(y-ybar))
    return ybar+(b-xbar)@beta

def cross_validate(d):
    models={"B0":(),"BP":("pointN",),"BD":D_VARS,"BPD":("pointN",*D_VARS)}
    preds={m:np.full(len(d),np.nan) for m in models}
    for f in range(1,6):
        tr=d.loc[d.fold!=f].copy()
        te=d.loc[d.fold==f].copy()
        if te.empty or set(tr.nuts2)&set(te.nuts2):
            raise ValueError("Original 2016 NUTS2 fold source data leakage or empty fold")
        for m,extra in models.items():
            preds[m][te.index.to_numpy()]=fit_predict(tr,te,extra)
    if any(not np.isfinite(p).all() for p in preds.values()):
        raise ValueError("Not all original 2016 regional mean EQR source rows out-of-fold predicted")
    out=d.copy()
    for m,p in preds.items():
        out[m+"_OOF_predicted_mean_EQR"]=p
        out[m+"_OOF_squared_error"]=(p-out.EQR.to_numpy())**2
    return out

def metrics(out):
    keys=("B0","BP","BD","BPD")
    bycountry={}
    for c,g in out.groupby("country"):
        bycountry[c]={k:float(g[k+"_OOF_squared_error"].mean()) for k in keys}
    byregion={k:float(out[k+"_OOF_squared_error"].mean()) for k in keys}
    bysite={k:float(np.average(out[k+"_OOF_squared_error"],weights=out.original_valid_site_n)) for k in keys}
    bynational={k:float(np.mean([v[k] for v in bycountry.values()])) for k in keys}
    def summarize(v):
        P=100*(v["BD"]-v["BPD"])/v["BD"]
        D=100*(v["BP"]-v["BPD"])/v["BP"]
        return {"unique_pointN_gain_pct":P,"unique_landN_block_gain_pct":D,
                "point_minus_land_unique_gain_pct":P-D,
                "point_alone_vs_baseline_gain_pct":100*(v["B0"]-v["BP"])/v["B0"],
                "land_alone_vs_baseline_gain_pct":100*(v["B0"]-v["BD"])/v["B0"]}
    main=summarize(bynational)
    countries=sorted(bycountry)
    vals=np.array([[bycountry[c][k] for k in keys] for c in countries])
    rng=np.random.default_rng(SEED)
    draws=rng.integers(0,len(countries),size=(BOOT,len(countries)))
    averages=vals[draws].mean(axis=1)
    gainP=100*(averages[:,2]-averages[:,3])/averages[:,2]
    gainD=100*(averages[:,1]-averages[:,3])/averages[:,1]
    C=gainP-gainD
    if not np.isfinite([gainP,gainD,C]).all():raise ValueError("Invalid fixed error bootstrap")
    cis={"unique_pointN_gain_pct":np.quantile(gainP,[.025,.975]).tolist(),
         "unique_landN_block_gain_pct":np.quantile(gainD,[.025,.975]).tolist(),
         "point_minus_land_unique_gain_pct":np.quantile(C,[.025,.975]).tolist()}
    leave={}
    for c in countries:
        other=[bycountry[q] for q in countries if q!=c]
        data={k:float(np.mean([v[k] for v in other])) for k in keys}
        leave[c]=summarize(data)
    folds={}
    for f,g in out.groupby("fold"):
        err={k:float(g[k+"_OOF_squared_error"].mean()) for k in keys}
        folds[str(f)]={"NUTS2":len(g),"country_n":g.country.nunique(),"contrasts_region_equal":summarize(err),
                       "MSE":err}
    return {
        "NUTS2_region_equal_MSE":byregion,"site_count_weighted_MSE_sensitivity":bysite,
        "country_equal_MSE":bynational,
        "source_unique_predictions_country_equal":main,
        "fixed_OOF_country_bootstrap_3999_descriptive95":cis,
        "country_mean_OOF_squared_error_by_country":bycountry,
        "leave_one_country_out_fixed_OOF_error_reaggregation_not_refit":leave,
        "five_frozen_geographic_fold_diagnostics":folds,
        "countries_with_positive_heldout_unique_pointN_gain_in_full_country_equal":sum(leave[c]["unique_pointN_gain_pct"]>0 for c in countries),
    }

def source_vif(d):
    x=d.pointN.to_numpy()
    dm=d[list(D_VARS)].to_numpy()
    X=pd.get_dummies(d.country,dtype=float).to_numpy()
    X=np.column_stack([X,d[list(NUMERIC_B0)].to_numpy()])
    def resid(y,z):
        return y-z@np.linalg.lstsq(z,y,rcond=None)[0]
    xp=resid(x,X)
    D=np.column_stack([resid(dm[:,j],X) for j in range(dm.shape[1])])
    if xp@xp<1e-10:return None
    xr=resid(xp,D)
    r2=float(np.clip(1-(xr@xr)/(xp@xp),0,1))
    return float(1/(1-r2)) if r2<1 else None

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--source-csv",type=Path,required=True)
    a.add_argument("--source-report",type=Path,required=True)
    a.add_argument("--source-gate",type=Path,required=True)
    a.add_argument("--compatibility-gate",type=Path,required=True)
    a.add_argument("--gisco-gpkg",type=Path,required=True)
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    source,geo,gate=read_source(args.source_csv,args.source_report,args.source_gate,args.gisco_gpkg)
    compatibility=json.loads(args.compatibility_gate.read_text())
    d,archive=frozen_region_data(source,geo,gate,compatibility)
    predictions=cross_validate(d)
    results=metrics(predictions)
    vif=source_vif(d)
    gain=results["source_unique_predictions_country_equal"]
    ci=results["fixed_OOF_country_bootstrap_3999_descriptive95"]
    supported=(vif is not None and vif<5 and gain["unique_pointN_gain_pct"]>=1
               and ci["unique_pointN_gain_pct"][0]>0
               and gain["point_minus_land_unique_gain_pct"]>0
               and ci["point_minus_land_unique_gain_pct"][0]>0)
    results["pointN_collinearity_VIF_adjusted_for_background_and_land"]=vif
    results["fixed_strong_pointN_predictive_dominance_criterion_met"]=bool(supported)
    results["classification"]="POINTN_PREDICTIVE_DOMINANCE_OBSERVATIONAL" if supported else "NO_STRONG_POINTN_DOMINANCE"
    out={
        "schema":"gesui_2016_aquatic_only_73nuts2_original_EEA_benthic_meanEQR_point_landN_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "2016_benthic_region_n":len(d),
        "2016_country_n":d.country.nunique(),
        "original_2016_frozen_numeric_EQR_summary_SHA256":sha(RESULT),
        "original_2016_source_geometry_SHA256":sha(ORIGINAL),
        "original_2016_preoutcome_fold_SHA256":sha(FROZEN),
        "source_only_2016_water_side_gate_SHA256":sha(args.source_gate),
        "independently_pre_model_numeric_availability_gate_SHA256":sha(args.compatibility_gate),
        "JRC269_soil_inputs_SOURCE_sha256":sha(args.source_csv),
        "GISCO2016_01M_sha256":sha(args.gisco_gpkg),
        "original_outcome_is_region_mean_benthic_normalized_EQR_not_tolerance_trait":True,
        "fixed_cv_alpha":ALPHA,"fixed_region_heldout_folds":5,
        "source_geographical_unit":"NUTS2 region; not individual plant or hydrological delivery reach",
        "original_2018_2017to19_water_land_2x2_source_STOP_remains":True,
        "posthoc_after_2016_point_EQR_known":True,
        "do_not_claim_2016_actual_WWTP_discharge_causal_effect":True,
        "metrics":results,
    }
    (args.out/"EEA2016_POINT_VS_LAND_N_EQR_REGION_4MODELS.json").write_text(
        json.dumps(safe(out),ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    predictions.to_csv(args.out/"EEA2016_NUTS2_REGION_MEANEQR_HELDOUT_PREDICTIONS.csv",index=False,float_format="%.12g")
    print("GESUI_EEA2016_POINT_VS_SOIL_N_4MODEL",json.dumps(safe({
        "NUTS2":len(d),"countries":d.country.nunique(),
        "primary_MSE":results["country_equal_MSE"],
        "unique_gains":gain,
        "95_descriptive":ci,
        "vif":vif,"classification":results["classification"]
    }),ensure_ascii=False),flush=True)

if __name__=="__main__":main()
