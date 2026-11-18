#!/usr/bin/env python3
"""Original published CELLDEX worldwide river vs adjacent riparian cellulose
decomposition, paired EXACT same field site. Source-pinned reanalysis, no
WWTP exposure, no aquatic→terrestrial resource flux measured.

Scientific record: Tiegs et al Science Advances 2019 DOI
10.1126/sciadv.aav0486; 2024 Science geospatial reuses same original data.
Main question: beyond biome & latitude, does measured adjacent riparian
decomposition predict heldout *stream* decomposition in geographically
independent original partner clusters?

Prespecified:
* 514 original positive river K rates with matched RIP; 540 source station IDs
* true source 131 original research partner clusters as heldout groups
* fixed 5fold GroupKFold, ridge alpha=10, train-only transformer
* source biome and abs latitude, latitude and log1p deployment days baseline
* add ONE predictor log(k_rip) to predict log(k_river), no sign fishing
* separate original degree-day rates subset n314 as sensitivity, not
  a third independent experiment
* 2000 resamples of fixed heldout per-partner squared errors are
  descriptive, not confirmatory cross-watershed experimental p-values.

This study measures paired aquatic and terrestrial PROCESSING POTENTIAL,
NOT wastewater-derived P/N treatment, actual litter mass flux, biomass or
river→riparian insects.
"""
from __future__ import annotations
import argparse,collections,csv,datetime,hashlib,io,json,math,random,statistics,time,urllib.request,urllib.error
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler,OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.linear_model import Ridge
import pandas as pd

REPO="https://raw.githubusercontent.com/dmcostello/CELLDEX2018/edddf5a9477585cbbcb8280e686fc8a6000a29e8/"
GEO="https://raw.githubusercontent.com/dmcostello/CELLDEX_geospatial/e659f61e7eebd6a4b243ddc102c320c31c6c7dbb/"
INPUTS={
 "river":(REPO+"str_k.csv","6e9aa57c7d06eea2cb9952c08c8b3c9f8dceb9d2",514),
 "riparian":(REPO+"rip_k.csv","fd171ae671be552148da8580fc29e4ba341a5338",533),
 "field_site":(REPO+"CELLDEX_SITE_DATA.csv","58e14962c9faa01d4ea4c66b91aea1345c6f20f3",540),
 "river_degree_day":(REPO+"str_k_dd.csv","4e05d1b56c9d8eabf23b8c477c9e5dba97800405",354),
 "riparian_degree_day":(REPO+"rip_k_dd.csv","66c4655509f43f0ac280ba99b7b3db80bd64419d",343),
 "deployment":(GEO+"CELLDEX_deploy_dates.csv","92b5362f5e0d363c6b6040bf625d4f2d56f59c2d",540)
}
SCHEMA="GESUI_20261118_CELLDEX_514_GENUINE_SAME_SITE_RIVER_RIPARIAN_FUNCTION_PARTNER_HOLDOUT_V1"

def git_blob(b):return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
def get(url):
    if not(url.startswith(REPO) or url.startswith(GEO)):
        raise RuntimeError("Only exact original pinned public GitHub repository ref")
    err=None
    for i in range(3):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"gesui-CELLDEX-source-pinned-field-paired-function-v1"})
            with urllib.request.urlopen(req,timeout=70) as r:
                b=r.read(1_000_001)
                if len(b)>1_000_000:raise RuntimeError("Unexpected original source file >1MB")
                return b
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as e:
            err=e
            if isinstance(e,urllib.error.HTTPError) and e.code not in (429,500,502,503,504):break
            time.sleep(1+2*i)
    raise RuntimeError(f"Published original GitHub data fetch failed: {str(err)[:240]}")
def load():
    out={}
    for k,(url,sha,expn) in INPUTS.items():
        b=get(url)
        if git_blob(b)!=sha:raise RuntimeError("Original author git blob data SHA drift "+k)
        rr=list(csv.DictReader(io.StringIO(b.decode("utf-8-sig"))))
        if len(rr)!=expn:raise RuntimeError(f"Original {k} source n mismatch {len(rr)} vs {expn}")
        out[k]=rr
        print("ORIGINAL_CELLDEX_SOURCE_VERIFIED="+json.dumps({"table":k,"rows":len(rr),"git_blob_sha1":sha}),flush=True)
    return out
def ids(t,table):
    if t in ("river","riparian","river_degree_day","riparian_degree_day"):
        return table["part.str"].strip()
    return table["partnerid"].strip()+" "+table["stream"].strip()
def mp(rows,key):
    d={}
    for row in rows:
        id=ids(key,row)
        if id in d:raise RuntimeError("True original published station identity duplicated "+key)
        d[id]=row
    return d
def num(x):
    try:v=float(str(x).strip())
    except (ValueError,TypeError):return None
    return v if math.isfinite(v) else None
def primary(o):
    r=mp(o["river"],"river");p=mp(o["riparian"],"riparian")
    sites=mp(o["field_site"],"field_site");deployment=mp(o["deployment"],"deployment")
    if len(r)!=514 or len(p)!=533 or len(sites)!=540 or len(deployment)!=540:
        raise RuntimeError("Original source unique station identity changed")
    if not set(r).issubset(p) or not set(r).issubset(sites) or not set(r).issubset(deployment):
        raise RuntimeError("Source actual paired river/riparian/date support not complete")
    items=[]
    for key in sorted(r,key=lambda k:(int(k.split()[0]),int(k.split()[1]))):
        sr,sp,s,dt=r[key],p[key],sites[key],deployment[key]
        k_r=num(sr["k"]);k_p=num(sp["k"]);lat=num(s["latitude"]);lon=num(s["longitude"]);days=num(s["deploy_time"])
        if (None in [k_r,k_p,lat,lon,days] or k_r<=0 or k_p<=0 or not(-90<=lat<=90 and -180<=lon<=180)
          or days<=0):
            raise RuntimeError("Original river or riparian processing rate/geocoordinate invalid "+key)
        date=dt["deploy_date"].strip()
        if date:
            try:year=datetime.datetime.strptime(date,"%m/%d/%y").year
            except ValueError:raise RuntimeError("True original field deployment date invalid "+key)
        else:year=None
        items.append({
           "original_exact_site_key":key,
           "partnerid":str(s["partnerid"]), "stream":int(s["stream"]),
           "latitude":lat,"longitude":lon,"abs_latitude":abs(lat),
           "biome_short":str(s["biome_short"]).strip(),
           "original_deploy_duration_days":days,
           "source_original_year":year,
           "source_original_deploy_date":date,
           "original_river_k_per_day":k_r,"original_riparian_k_per_day":k_p,
           "ln_river_rate":math.log(k_r),"ln_riparian_rate":math.log(k_p),
           "ln_river_over_riparian":math.log(k_r/k_p)})
    if len(items)!=514 or len(set(x["partnerid"] for x in items))!=131:
        raise RuntimeError("Original exact main source 514 river-riparian pairs /131 partner clusters failed")
    if len({x["original_exact_site_key"] for x in items})!=514:
        raise RuntimeError("Original duplicated actual river-riparian pair identity")
    return items
def dd(o):
    u=mp(o["river_degree_day"],"river_degree_day")
    v=mp(o["riparian_degree_day"],"riparian_degree_day")
    if len(u)!=354 or len(v)!=343:raise RuntimeError("Original degree-day source n drift")
    full=[]
    for key in sorted(set(u)&set(v)):
        x=num(u[key]["k"]);y=num(v[key]["k"])
        if x is not None and y is not None and x>0 and y>0:
            full.append((key,x,y))
    if len(full)!=314:raise RuntimeError("Source true degree-day paired support differs")
    return full

def grouped_cv(rows):
    frame=pd.DataFrame(rows)
    target=frame["ln_river_rate"].to_numpy(float)
    group=frame["partnerid"].to_numpy()
    preds={}
    for name,nums in (
      ("original_environment_only",["latitude","abs_latitude","original_deploy_duration_days"]),
      ("original_environment_PLUS_same_site_riparian_rate",["latitude","abs_latitude","original_deploy_duration_days","ln_riparian_rate"])):
        cats=["biome_short"]
        p=np.full(len(frame),np.nan,float)
        fold_members=[]
        for tr,te in GroupKFold(n_splits=5).split(frame,target,group):
            numeric=Pipeline([("scale",StandardScaler())])
            cat=OneHotEncoder(handle_unknown="ignore")
            prep=ColumnTransformer([("numeric",numeric,nums),("biome",cat,cats)])
            fit=Pipeline([("preprocess",prep),("ridge",Ridge(alpha=10.0))])
            fit.fit(frame.iloc[tr],target[tr])
            p[te]=fit.predict(frame.iloc[te])
            fold_members.append({"n_heldout_sites":int(len(te)),
              "n_heldout_independent_partner_groups":int(len(set(group[te]))),
              "partner_overlap_train_test_n":int(len(set(group[tr])&set(group[te])))})
        if not np.all(np.isfinite(p)):raise RuntimeError("Original heldout predictions incomplete")
        preds[name]={"predicted_log_k":p,"fold_source":fold_members}
    errs={}
    for name,obj in preds.items():
        ee=(target-obj["predicted_log_k"])**2
        errs[name]={g:float(ee[group==g].mean()) for g in sorted(set(group))}
    b=errs["original_environment_only"];a=errs["original_environment_PLUS_same_site_riparian_rate"]
    groups=sorted(b)
    bm=float(np.mean([b[k] for k in groups]));am=float(np.mean([a[k] for k in groups]))
    raw=100*(1-am/bm)
    rng=np.random.default_rng(20261118)
    lo=[]
    for it in range(2000):
        selected=rng.choice(groups,size=len(groups),replace=True)
        bb=float(np.mean([b[k] for k in selected]));aa=float(np.mean([a[k] for k in selected]))
        lo.append(100*(1-aa/bb))
    return {"original_5_grouped_folds":preds["original_environment_only"]["fold_source"],
            "original_partner_groups_n":len(groups),
            "MSE_baseline_environment_partner_equal":bm,
            "MSE_add_measured_riparian_processing_partner_equal":am,
            "original_partner_equal_prediction_MSE_gain_pct":raw,
            "fixed_oof_2000_partner_group_resampling_descriptive95":
                [float(np.quantile(lo,0.025)),float(np.quantile(lo,0.975))],
            "source_site_weighted_MSE_gain_pct":100*(1-np.mean((target-preds["original_environment_PLUS_same_site_riparian_rate"]["predicted_log_k"])**2)/
                                                    np.mean((target-preds["original_environment_only"]["predicted_log_k"])**2)),
            "original_heldout_error_resampling_is_descriptive_not_independent_causal_p":True}

def run(args):
    data=load()
    rows=primary(data)
    paired_dd=dd(data)
    a=np.array([x["original_river_k_per_day"] for x in rows],float)
    b=np.array([x["original_riparian_k_per_day"] for x in rows],float)
    ratio=np.array([x["ln_river_over_riparian"] for x in rows],float)
    dates=collections.Counter(x["source_original_year"] for x in rows)
    year_counts={str(k):v for k,v in dates.items()}
    stats={
       "genuine_exact_n_paired_river_riparian_sites":len(rows),
       "genuine_source_2019_river_k_rows":514,
       "genuine_source_2019_riparian_k_rows":533,
       "source_coordinates_unique_sites":540,
       "original_biogeographic_partner_clusters":131,
       "true_original_field_deployment_year_counts":year_counts,
       "original_source_Spearman_river_vs_riparian_k":float(spearmanr(a,b).statistic),
       "original_source_Pearson_log_rates":float(np.corrcoef(np.log(a),np.log(b))[0,1]),
       "river_faster_than_adjacent_riparian_n":int(np.sum(a>b)),
       "median_original_k_river_per_day":float(np.median(a)),
       "median_original_k_riparian_per_day":float(np.median(b)),
       "median_river_over_riparian_rate_ratio":float(np.exp(np.median(ratio))),
       "source_degree_day_adjusted_complete_pairs_n":len(paired_dd),
       "original_degree_day_Spearman_river_vs_riparian":float(spearmanr(
           [x[1] for x in paired_dd],[x[2] for x in paired_dd]).statistic),
       "source_genuine_same_method_cotton_cellulose_processing_both_habitats":True}
    out={
      "schema":SCHEMA,
      "status":"SOURCE_VERIFIED_PUBLISHED_GLOBAL_514_MATCHED_RIVER_RIPARIAN_DECOMPOSITION_PARTNER_HOLDOUT_ANALYSIS",
      "publication_2019_original_DOI":"10.1126/sciadv.aav0486",
      "publication_2024_same_original_experiment_reused_DOI":"10.1126/science.adn1262",
      "input_original_Git_blob_sha1":{k:v[1] for k,v in INPUTS.items()},
      "source_2019_original_CELLDEx_commit":"edddf5a9477585cbbcb8280e686fc8a6000a29e8",
      "source_2024_geospatial_deploy_commit":"e659f61e7eebd6a4b243ddc102c320c31c6c7dbb",
      "numeric_original_scientific_source_descriptives":stats,
      "original_131_partner_hard_group_heldout_prediction":grouped_cv(rows),
      "original_2018_2022_EEA_UWWTD_emission_source_NOT_joined":True,
      "original_field_2015_2016_vs_EEA2018_time_mismatch_likely":True,
      "n_source_specific_municipal_WWTP_exposure_identified":0,
      "paired_rate_correlation_not_actual_direct_transboundary_mass_transport":True,
      "original_2_habitats_are_true_same_field_site_different_decomposition_microenvironments":True,
      "current_reanalysis_not_independent_original_2019_or_2024_scientific_discovery":True,
      "site_k_vs_nutrient_export_and_sewage_not_fitted_in_this_run":True}
    args.out.mkdir(parents=True,exist_ok=True)
    with (args.out/"CELLDEX_ORIGINAL_GENUINE_514_RIVER_RIPARIAN_COMMON_STATION_COTTON_RATE.csv").open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (args.out/"CELLDEX_514_ORIGINAL_GENUINE_RIVER_RIPARIAN_FUNCTION_GROUPOUT_QA.json").write_text(
      json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("CELLDEX_514_TRUE_RIVER_RIPARIAN_MACRO="+json.dumps({
      "status":out["status"],"source_n":len(rows),"unique_partners":131,
      "original_year_counts":year_counts,
      "Spearman":stats["original_source_Spearman_river_vs_riparian_k"],
      "degree_day_n":len(paired_dd),
      "degree_day_rho":stats["original_degree_day_Spearman_river_vs_riparian"],
      "OOS_source_added_riparian_mse_gain_pct":out["original_131_partner_hard_group_heldout_prediction"]["original_partner_equal_prediction_MSE_gain_pct"],
      "OOS_partner_descriptive95":out["original_131_partner_hard_group_heldout_prediction"]["fixed_oof_2000_partner_group_resampling_descriptive95"]},ensure_ascii=False),flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    run(ap.parse_args())
