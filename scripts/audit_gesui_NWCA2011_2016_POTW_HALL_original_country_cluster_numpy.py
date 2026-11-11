#!/usr/bin/env python3
"""Independent linear algebra for historical NWCA2011/16 probability plots.

Use SHA-frozen original source-site joins from independent preoutcome
SITE+POTW gates, then independently calculate 2011 and 2016 (not 2021)
state×ecoregion adjusted H_ALL, WLS and OLS using NumPy matrix and CR1
state-cluster sandwich, verifying earlier statsmodels year outputs.
All weights are actual original design values; sources are not true
discharge/river pathway. Explicitly posthoc 3-survey sign sensitivity.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from gesui_NWCA2011_2016_2021_POTW_near10_HALL_HNAT_FQAI_temporal_source_replay import source_frame
D=Path("results/us_NWCA2021_plant_TNTP_20261110")
F=D/"NWCA2011_2016_2021_POTW_NEAR10_FAR20_HALL_THREE_SURVEY_YEARS.json"
S=D/"NWCA2011_2016_ORIGINAL_SITE_VEG_SOURCE_BEFORE_PLANT_OUTCOMES.json"
G=D/"NWCA2011_2016_POTW2026_SOURCE_ONLY_NEAR10_FAR20_2CYCLE.json"
Z=1.959963984540054

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def design(df):
    categories=("STATE_ECO","WETCLSGRP","HYDCLS")
    pieces=[np.ones(len(df)),df.NEAR.to_numpy(dtype=float)]
    fields=["constant","NEAR"]
    for col in categories:
        levels=sorted(df[col].astype(str).unique())
        for v in levels[1:]:
            pieces.append((df[col].astype(str)==v).to_numpy(dtype=float))
            fields.append(col+"="+v)
    pieces.append(df.LATANALYS.to_numpy(dtype=float))
    fields.append("LATANALYS")
    return np.column_stack(pieces),fields

def one(df,weighted):
    X,labels=design(df)
    y=df.HALL.to_numpy(dtype=float)
    N=len(y);K=X.shape[1]
    w=df.WGT.to_numpy(dtype=float)
    if weighted:w=w/np.median(w)
    else:w=np.ones(N)
    x=X*np.sqrt(w[:,None]);z=y*np.sqrt(w)
    beta=np.linalg.lstsq(x,z,rcond=None)[0]
    residual=z-x@beta
    gram=np.linalg.pinv(x.T@x)
    meat=np.zeros_like(gram)
    groups=sorted(df.PSTLCODE.astype(str).unique())
    for name in groups:
        mask=(df.PSTLCODE.astype(str)==name).to_numpy()
        score=x[mask].T@residual[mask]
        meat+=np.outer(score,score)
    G=len(groups)
    if G<3 or N<=K:raise ValueError("Too few source states or observations for 2011/2016 covariance")
    variance=gram@meat@gram*(G/(G-1))*((N-1)/(N-K))
    if variance[1,1]<0:raise ValueError("Invalid cluster variance")
    se=math.sqrt(float(variance[1,1]))
    effect=float(beta[1])
    return {
      "n_source_matched_plots":N,"states":G,
      "state_ecoregion_matched_strata":int(df.STATE_ECO.nunique()),
      "matrix_regressor_columns":K,"matrix_rank":int(np.linalg.matrix_rank(x)),
      "beta_near_minus_far":effect,
      "state_cluster95":[effect-Z*se,effect+Z*se],
      "site_equal" : not weighted,
      "historical_NWCA_weighted_WLS_not_population_acreage_claim":weighted,
    }

def main_run(root,data_dir):
    obj=json.loads((root/F).read_text())
    site=json.loads((root/S).read_text())
    geo=json.loads((root/G).read_text())
    if obj["original_pre_HALL_EPA_2011_2016_veg_site_source_SHA256"]!=sha(root/S) or obj["original_pre_HALL_2011_2016_POTW_geo_source_SHA256"]!=sha(root/G):
        raise ValueError("Source SHA no longer matches originally reported plant-year effects")
    data={}
    for year in ("2011","2016"):
        df,stop=source_frame(year,root,data_dir,site,geo)
        if stop or len(df)!=obj["years"][year]["source_fixed_n"]:
            raise ValueError("Original PRE-HALL year source site population drifted")
        if df.HALL.isna().any():
            raise ValueError("Original 2011/2016 site-matched vegetation H_ALL incomplete; different model sample")
        df=df.sort_values(["STATE_ECO","SITE_ID"]).copy()
        result={}
        for key,weighted in (("site_equal_OLS",False),("site_design_weighted_WLS",True)):
            v=one(df,weighted)
            original=obj["years"][year]["outcomes"]["HALL"][key]
            if abs(v["beta_near_minus_far"]-original["near_minus_far_beta"])>1e-8:
                raise ValueError(f"{year} {key} independent near β mismatches original {v['beta_near_minus_far']} vs {original['near_minus_far_beta']}")
            for a,b in zip(v["state_cluster95"],original["state_clustered_descriptive95"]):
                if abs(a-b)>2e-5:
                    raise ValueError(f"{year} {key} state cluster CI mismatches {v['state_cluster95']} vs {original['state_clustered_descriptive95']}; regressor_columns={v['matrix_regressor_columns']}")
            result[key]=v
        # Main geographic influence: drop one whole US state, re-solve
        # original source model for weighted and unweighted, no significance.
        drops={}
        for state in sorted(df.PSTLCODE.astype(str).unique()):
            sub=df[df.PSTLCODE.astype(str)!=state]
            drops[state]={
               "site_equal":one(sub,False)["beta_near_minus_far"],
               "survey_weighted":one(sub,True)["beta_near_minus_far"],
               "source_site_n":len(sub)
            }
        summary={
          "US_states_left_out_n":len(drops),
          "all_state_omissions":drops,
          "source_equal_beta_positive_omission_n":sum(v["site_equal"]>0 for v in drops.values()),
          "source_survey_weighted_beta_positive_omission_n":sum(v["survey_weighted"]>0 for v in drops.values()),
          "site_equal_beta_minmax":[min(v["site_equal"] for v in drops.values()),max(v["site_equal"] for v in drops.values())],
          "weighted_beta_minmax":[min(v["survey_weighted"] for v in drops.values()),max(v["survey_weighted"] for v in drops.values())]
        }
        data[year]={"independent_NWCA_2011_2016_actual_HALL":result,"leave_one_US_state_out_influence":summary}
    source2021=obj["2021_reference_only_not_reestimated"]
    return {
      "schema":"gesui_NWCA_2011_2016_2021_HALL_3cycle_weighted_positive_site_equal_nonrobust_independent_numpy_QA_v01",
      "original_2011_2016_vegetation_SITE_only_SHA256":sha(root/S),
      "original_2011_2016_ECHO2026_outfall_source_only_SHA256":sha(root/G),
      "original_2011_2016_2021_HALL_temporal_main_model_SHA256":sha(root/F),
      "original_EPA_2011_2016_SITE_values_precommitted_before_HALL":True,
      "2021_original_reference_reused_not_reestimated":source2021,
      "both_independent_older_year_beta_and_state_cluster95_reproduced":True,
      "2011_2016_independent_source_years":data,
      "fixed2026_permit2021_plausible_not_verified_2011_or2016_WWTP_discharge":True,
      "even_3cycle_survey_area_weighted_positive_not_causal_effluent_or_true_temporal_trend":True,
      "no_plant_taxa_selected_by_favorable_WWTP_association":True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-dir",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    result=main_run(root,a.source_dir)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print(json.dumps({k:{"independent":v["independent_NWCA_2011_2016_actual_HALL"],
        "loo":{x:y for x,y in v["leave_one_US_state_out_influence"].items() if x!="all_state_omissions"}}
         for k,v in result["2011_2016_independent_source_years"].items()},indent=2))
if __name__=="__main__":main()
