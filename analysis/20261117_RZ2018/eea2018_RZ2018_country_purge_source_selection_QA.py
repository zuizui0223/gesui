#!/usr/bin/env python3
"""Original 2018 source-locked RZ landcover × real benthic EQR country transfer.

New strictly POSTHOC adversarial CV after the fivefold PR123 results were known:
leave ONE entire country out; original 2639 site source and 2634 numeric EQR;
same 2018 source RZ polygon sample, same traits, no new WWTP causal inference.
"""
from __future__ import annotations
import argparse,collections,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler,OneHotEncoder
from sklearn.linear_model import Ridge
from eea2018_RZ2018_benthic_EQR_geoblock_model_if_fullsource import (
  source_frame,source_eqr,ORIG_SOURCE,ORIG_EQR,
)

SCHEMA="GESUI_20261117_RZ2018_2018_BENTHIC_CROSS_COUNTRY_AND_REPORT_COVERAGE_ADVERSARIAL_V1"
NUM_BASE=["lat","lon","LogQ","LogArea"]
NUM_LAND=NUM_BASE+["NaturalFar","UrbanFar","CropFar"]
CAT=["countryCode","countryMethod"]
COMPARE={
  "HydroGeo_B":NUM_BASE,
  "Landcover_L":NUM_LAND,
  "Landcover_plus_UP_reportedP":NUM_LAND+["LogUpP"],
  "Landcover_plus_NearFar":NUM_LAND+["NaturalNearMinusFar"],
  "Landcover_plus_UP_reportedN":NUM_LAND+["LogUpN"],
}
SEED=20261117
NBOOT=2000

def h(v):
    return hashlib.sha256(str(v).encode()).hexdigest()

def check_source(path):
    a=json.loads(path.read_text())
    if (a.get("status")!="PASS_2018_RZ2018_ALL_ORIGINAL_BENTHIC_SITES_SOURCE_ONLY_NO_ECOLOGY" or
        a.get("source_original_Github_2018_fixed_data_sha256")!=ORIG_SOURCE or
        len(a.get("original_source_pilot_true_station_classification",[]))!=2639):
        raise RuntimeError("Original real RZ source-only class census STOP")
    r=a["original_source_pilot_true_station_classification"]
    sourceplus=sum(bool(x["source_2018_upstream_P_positive_report"]) for x in r)
    if sourceplus!=204:raise RuntimeError("Original source non-biological report-positive 204 denominator changed")
    eligible=[x for x in r if x["source_near_outer_0p70_eligible"]]
    plus=sum(bool(x["source_2018_upstream_P_positive_report"]) for x in eligible)
    if (len(eligible),plus)!=(986,175):
        raise RuntimeError("Original 986 RZ and 175 original reported-source-positive coverage changed")
    total=len(r);npos=204;nneg=total-npos
    x={"report_positive_original_n":npos,"report_nonpositive_original_n":nneg,
       "report_positive_RZ_both_rings_70pct_n":plus,
       "report_nonpositive_RZ_both_rings_70pct_n":len(eligible)-plus,
       "original_RZ_coverage_positive_pct":100*plus/npos,
       "original_RZ_coverage_nonpositive_pct":100*(len(eligible)-plus)/nneg,
       "original_report_positive_vs_nonpositive_RZ_coverage_odds_ratio":
         (plus*(nneg-len(eligible)+plus))/((npos-plus)*(len(eligible)-plus)),
       "source_nonpositive_does_not_mean_no_actual_sewage":True,
       "source_selected_sample_not_random_of_original_river_census":True}
    return x

def frames(rz):
    source_cover=check_source(rz)
    r,q=source_frame(rz)
    E,s=source_eqr()
    d=E.merge(r,on="site_key",how="inner",validate="one_to_one")
    if (len(E)!=2634 or len(d)!=976 or
        d.countryCode.nunique()!=13 or d.region_3deg_block.nunique()!=49 or
        (d.countryCode!=d.RZOriginalCountry).any() or
        (d.region_3deg_block!=d.RZOriginal3DegreeBlock).any()):
        raise RuntimeError("Source original site/EQR/country support not exactly verified 976 /13/49")
    if sum(d["2018_measured_upstream_P_t_per_year"]>0)!=169:
        raise RuntimeError("Actual original 2018 report-positive 169 source drift")
    d["LogQ"]=np.log1p(d.DIS_AV_CMS.to_numpy(float))
    d["LogArea"]=np.log1p(d.UPLAND_SKM.to_numpy(float))
    d["LogUpP"]=np.log1p(d["2018_measured_upstream_P_t_per_year"].to_numpy(float))
    d["LogUpN"]=np.log1p(d["2018_measured_upstream_N_t_per_year"].to_numpy(float))
    d["countryMethod"]=d.countryCode.astype(str)+"::"+d.procedureClassificationSystem.astype(str)
    d=d.sort_values(["countryCode","site_key"],kind="stable").reset_index(drop=True)
    return d,source_cover,q

def calc(d):
    countries=sorted(d.countryCode.unique().tolist())
    group=d.countryCode.to_numpy(str)
    y=d.numeric_EQR_2018.to_numpy(float)
    losses={}; country_errors={}
    for name,features in COMPARE.items():
        X=d[features+CAT].copy()
        X[CAT]=X[CAT].astype(str)
        pred=np.full(len(d),np.nan)
        for c in countries:
            test=group==c;train=~test
            if not test.any():raise RuntimeError("Source country purged test empty")
            pipeline=Pipeline([
               ("preprocess",ColumnTransformer([
                 ("numeric",StandardScaler(),features),
                 ("categorical",OneHotEncoder(handle_unknown="ignore"),CAT)
               ])),
               ("ridge",Ridge(alpha=10.0))
            ])
            pipeline.fit(X.iloc[train],y[train])
            pred[test]=pipeline.predict(X.iloc[test])
        if not np.isfinite(pred).all():
            raise RuntimeError("Original source country-purged prediction incomplete")
        e=(y-pred)**2
        losses[name]=e
        country_errors[name]=np.array([e[group==c].mean() for c in countries])
    out={}
    for base,new in (("HydroGeo_B","Landcover_L"),
                     ("Landcover_L","Landcover_plus_UP_reportedP"),
                     ("Landcover_L","Landcover_plus_NearFar"),
                     ("Landcover_L","Landcover_plus_UP_reportedN")):
        a=country_errors[base];b=country_errors[new]
        gain=100*(1-b.mean()/a.mean())
        rng=np.random.default_rng(SEED)
        picks=rng.integers(0,len(countries),size=(NBOOT,len(countries)))
        bs=100*(1-b[picks].mean(axis=1)/a[picks].mean(axis=1))
        out[new+"_AFTER_"+base]={
            "whole_country_holdout_siteweighted_MSEgain_pct":
              float(100*(1-losses[new].mean()/losses[base].mean())),
            "whole_country_holdout_country_equal_MSEgain_pct":float(gain),
            "descriptive_fixed_country_errors_resample95_pct":
              np.quantile(bs,[.025,.975]).tolist(),
            "number_countries_showing_positive_gain":int(sum(b<a)),
            "country_count":len(countries),
            "source_country_by_country_gain_pct":{c:float(100*(1-b[i]/a[i]))
                                                  for i,c in enumerate(countries)},
            "not_a_new_independent_country_replication_or_causal_test":True,
        }
    return out,{c:int(sum(group==c)) for c in countries}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-json",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    d,sel,q=frames(a.source_json)
    out,groups=calc(d)
    result={
       "schema":SCHEMA,
       "status":"PASS_SOURCE_VERIFIED_POSTHOC_COUNTRY_PURGED_OBSERVATIONAL_DIAGNOSTIC",
       "source_original_2018_EEA_preoutcome_site_SHA256":ORIG_SOURCE,
       "source_original_2018_EEA_numeric_EQR_SHA256":ORIG_EQR,
       "original_2018_RZ_source_status":"PASS_2639_SOURCE_ONLY",
       "source_RZ_selection_by_original_upstream_report_status":sel,
       "new_source_min10pct_terrestrial_cover_exclusion":q["new_additional_preoutcome_terrestrial_land_area_exclusion"],
       "real_original_2018_numeric_EQR_river_sites":len(d),
       "source_eligible_country_sizes":groups,
       "original_complete_country_separated_predictive_stress":out,
       "the_no_report_label_does_not_prove_no_sewage":True,
       "selected_riparian_mapped_sites_not_a_random_sample_of_rivers":True,
       "mapped_vegetated_class_includes_managed_forest_and_grassland":True,
       "no_actual_plant_biomass_or_crosswater_trophic_path":True,
       "posthoc_after_existing_PR123_fivefold_results":True,
       "no_claim_of_confirmatory_or_causal_significance":True
    }
    (a.out/"SOURCE_VERIFIED_2018_RIPARIAN_EQR_REPORT_SELECTION_COUNTRY_HOLDOUT.json").write_text(
       json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("GESUI_2018_RZ2018_COUNTRY_QA="+json.dumps({
        "status":result["status"],"n":len(d),"countries":groups,
        "source_RZ_selection":sel,
        "country_equal_Holdout_Predictive_gains":{k:round(v["whole_country_holdout_country_equal_MSEgain_pct"],4)
                                                  for k,v in out.items()}},
       ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":main()
