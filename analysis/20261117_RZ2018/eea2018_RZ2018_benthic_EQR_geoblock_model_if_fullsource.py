#!/usr/bin/env python3
"""Source-gated exploratory 2018 EEA river benthic EQR × original RZ2018 terrestrial LCLU.

This posthoc study uses original same-site/year WWTP-reported mass, genuine
2018 benthic EQR and independently mapped riparian land COVER. No plant
biomass, nutrient transfer, butterflies, causality or ecosystem services.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,collections
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder,StandardScaler
from sklearn.linear_model import Ridge

RZ_SCHEMA="GESUI_2018_ORIGINAL_EEA_UWWTD_BENTHIC_RZ2018_250M_500M_SOURCE_ONLY_V1"
EQRCODE="EEA_13-01-4"
ORIG_SOURCE="5a074f5182711c4ab74ef9bdacfa3918e9fcb434dd0b75c67072a0f06459c3e1"
ORIG_EQR="2a1a72d43656aca588a832dfae13860054c01a651e76bf25f984c5a374709b02"
EQRCsv=Path("results/eea2018_independent_temporal_bqe_20261108/benthic_invertebrates_2018_WWTP_source_real_EQR.csv")
SOURCECSV=Path("results/eea2018_independent_source_frozen_20261108/2018_SOURCE_ONLY_BQE_WWTP_NP_UPSTREAM_AND_NEGATIVE_CONTROL.csv")
SCHEMA="GESUI_2018_RIVER_BENTHIC_RZ2018_LCLU_COMMON_SITE_SPATIAL_HOLDOUT_POSTHOC_V1"
NBOOT=2000
SEED=20261117
HABITAT={"3","4","5","7"}
WATER_CODE="8"
FAR="NaturalFar"
Y="numeric_EQR_2018"
YLAND="NaturalNearMinusFar"

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def self_test():
    d={"1":.2,"2":.1,"3":.35,"4":.15,"7":.05,"8":.15}
    fr=hindex(d)
    assert abs(fr["natural"]-(.35+.15+.05)/.85)<1e-12
    assert abs(fr["urban"]-.2/.85)<1e-12
    assert 0<=fr["cropland"]<=1
    return {"status":"PASS_SYNTHETIC_SOURCE_CLASS_FRACTION_AND_NONWATER_DENOMINATOR"}
def hindex(code_area):
    a={str(k):float(v) for k,v in code_area.items()}
    if not a or any(v<0 or not np.isfinite(v) for v in a.values()):
        raise ValueError("Original mapped class code area invalid")
    coverage=sum(a.values())
    terrestrial=coverage-a.get(WATER_CODE,0)
    if not (0.1<=terrestrial<=1.0001 and 0<=coverage<=1.0001):
        raise ValueError("Unmapped or all-water source, nonterrestrial support STOP")
    nat=sum(a.get(c,0) for c in HABITAT)/terrestrial
    urban=a.get("1",0)/terrestrial
    crop=a.get("2",0)/terrestrial
    if not all(-1e-8<=x<=1.0001 for x in (nat,urban,crop)):
        raise ValueError("Original official group class area invalid")
    return {"natural":nat,"urban":urban,"cropland":crop,
        "classified_terrestrial_area_fraction_of_full_circle":terrestrial}

def source_frame(path):
    data=json.loads(path.read_text())
    if data.get("schema")!=RZ_SCHEMA or data.get("status")!="PASS_2018_RZ2018_ALL_ORIGINAL_BENTHIC_SITES_SOURCE_ONLY_NO_ECOLOGY":
        raise RuntimeError("Original full source census class QA and fixed gate did not PASS")
    population=data.get("source_original_EEA2018_2018_benthic_population",{})
    if (population.get("original_benthic_all_2018_source_n")!=2639 or
        not population.get("full_census_all_original_2018_benthic_source_sites") or
        population.get("original_source_positive_P_2018_n")!=204 or
        data.get("source_original_Github_2018_fixed_data_sha256")!=ORIG_SOURCE):
        raise RuntimeError("Original census source denominator and reporting years drift")
    if (data.get("source_pilot_eligible_n",0)<200 or
        data.get("source_pilot_eligible_countries",0)<6 or
        data.get("source_pilot_eligible_3degree_regions",0)<10 or
        data.get("source_pilot_eligible_upstream_P_positive_sites",0)<20):
        raise RuntimeError("Original new full 2018 RZ source sample gate FAIL")
    r=data["original_source_pilot_true_station_classification"]
    if len(r)!=2639 or len({x["site_key"] for x in r})!=2639:
        raise RuntimeError("Original full physical 2018 census source did not return all stations")
    valid=[x for x in r if x["source_near_outer_0p70_eligible"]]
    rows=[]
    for x in valid:
        n=hindex(x["source_Level1_class_area_frac_of_entire_near250m_circle"])
        o=hindex(x["source_Level1_class_area_frac_of_entire_outer250to500m_ring"])
        rows.append({"site_key":x["site_key"],
                     "RZOriginalCountry":x["country"],
                     "RZOriginal3DegreeBlock":x["source_original_region_3deg"],
                     "OriginalUpPReported":bool(x["source_2018_upstream_P_positive_report"]),
                     "NaturalNear":n["natural"],"NaturalFar":o["natural"],
                     "NaturalNearMinusFar":n["natural"]-o["natural"],
                     "UrbanFar":o["urban"],"CropFar":o["cropland"],
                     "NaturalNearAndFarMappedTerrestrialAreaPass":True})
    return pd.DataFrame(rows),data

def source_eqr():
    if sha(EQRCsv)!=ORIG_EQR or sha(SOURCECSV)!=ORIG_SOURCE:
        raise RuntimeError("Frozen 2018 EEA site source or original EQR SHA has changed")
    t=pd.read_csv(EQRCsv,dtype={"site_key":str})
    if (len(t)!=2634 or t.site_key.duplicated().any() or
        set(t.observedPropertyDeterminandBiologyEQRCode)!={EQRCODE}):
        raise RuntimeError("Original genuine 2018 benthic numeric EQR cohort drift")
    if (t[Y].isna().any() or not np.isfinite(t[Y].to_numpy(float)).all() or
        (t[Y]<0).any()):
        raise RuntimeError("Original 2018 EQR missing numeric values; never impute")
    src=pd.read_csv(SOURCECSV,dtype={"site_key":str})
    src=src.loc[src.observedPropertyDeterminandBiologyEQRCode==EQRCODE].copy()
    if len(src)!=2639 or src.site_key.duplicated().any():
        raise RuntimeError("Original source pre-EQR study denominator missing")
    if not set(t.site_key).issubset(set(src.site_key)):
        raise RuntimeError("A real 2018 EQR site not in original prospective source")
    return t,src

NUM_BASE=["lat","lon","LogQ","LogArea"]
CAT=["countryCode","countryMethod"]
NUM_LAND=NUM_BASE+["NaturalFar","UrbanFar","CropFar"]
MODELS={
 "Original_HydroGeographic_B":NUM_BASE,
 "LandContext_L":NUM_LAND,
 "L_plus_NearMinusFar_Natural":NUM_LAND+[YLAND],
 "L_plus_UP_2018_reported_P":NUM_LAND+["LogUpP"],
 "L_plus_UP_2018_reported_P_and_NearMinusFar_Natural":NUM_LAND+["LogUpP",YLAND],
 "L_plus_DOWNSTREAM_P_placebo":NUM_LAND+["LogDownP"],
 "L_plus_UP_report_P_and_NearNatural_INTERACTION":NUM_LAND+["LogUpP",YLAND,"UpP_x_DeltaNat"],
 "L_plus_UP_2018_reported_N":NUM_LAND+["LogUpN"],
}
COMPARE=[
 ("LandContext_L","Original_HydroGeographic_B"),
 ("L_plus_NearMinusFar_Natural","LandContext_L"),
 ("L_plus_UP_2018_reported_P","LandContext_L"),
 ("L_plus_UP_2018_reported_P_and_NearMinusFar_Natural","L_plus_UP_2018_reported_P"),
 ("L_plus_UP_2018_reported_P_and_NearMinusFar_Natural","L_plus_NearMinusFar_Natural"),
 ("L_plus_DOWNSTREAM_P_placebo","LandContext_L"),
 ("L_plus_UP_report_P_and_NearNatural_INTERACTION","L_plus_UP_2018_reported_P_and_NearMinusFar_Natural"),
 ("L_plus_UP_2018_reported_N","LandContext_L"),
]
ALPHA=10.
def matrix(d,nums):
    if not np.isfinite(d[nums].to_numpy(float)).all():
        raise ValueError("Missing source original 2018 and new RZ landcover input")
    X=d[nums+CAT].copy()
    X[CAT]=X[CAT].astype(str)
    return X

def cv(frame,outcome,models=MODELS):
    if not set(frame["spatial_fold"].astype(int))=={1,2,3,4,5}:
        raise RuntimeError("Original pre-frozen 2018 five geographical folds unavailable")
    n=len(frame)
    losses={}
    fold_metrics={}
    for name,nums in models.items():
        X=matrix(frame,nums)
        y=frame[outcome].to_numpy(float)
        pred=np.full(n,np.nan)
        fold_metrics[name]={}
        for fold in range(1,6):
            train=np.flatnonzero(frame.spatial_fold.to_numpy(int)!=fold)
            test=np.flatnonzero(frame.spatial_fold.to_numpy(int)==fold)
            if len(test)<10 or len(train)<50:
                raise RuntimeError("Original one of five source folds insufficient after independent RZ coverage")
            model=Pipeline([
                ("feature",ColumnTransformer([
                    ("numeric",StandardScaler(),nums),
                    ("country_method",OneHotEncoder(handle_unknown="ignore"),CAT)
                ])),
                ("ridge",Ridge(alpha=ALPHA))])
            model.fit(X.iloc[train],y[train])
            pred[test]=model.predict(X.iloc[test])
            fold_metrics[name][str(fold)]={"n_test":int(len(test)),
                "MSE":float(np.mean((y[test]-pred[test])**2))}
        if not np.isfinite(pred).all():
            raise RuntimeError("Original source fold missing predictions")
        losses[name]=(y-pred)**2
    results={}
    blocks=frame.region_3deg_block.astype(str).to_numpy()
    unique=sorted(set(blocks))
    by_err={m:np.array([np.mean(e[blocks==b]) for b in unique]) for m,e in losses.items()}
    rng=np.random.default_rng(SEED)
    draws=rng.integers(0,len(unique),size=(NBOOT,len(unique)))
    for add,base in COMPARE:
        if add not in models or base not in models: continue
        base_site=float(np.mean(losses[base]));add_site=float(np.mean(losses[add]))
        a=by_err[base];b=by_err[add]
        if base_site<=0 or a.mean()<=0:
            raise RuntimeError("Crossvalidation source baseline EQR variance unavailable")
        boots=100*(1-b[draws].mean(axis=1)/a[draws].mean(axis=1))
        results[add+"_ADDED_AFTER_"+base]={
         "original_test_population_site_weighted_MSE_gain_percent":
             float(100*(1-add_site/base_site)),
         "original_3deg_group_equal_MSE_gain_percent":
             float(100*(1-b.mean()/a.mean())),
         "original_fixed_oof_region_error_resampling_descriptive95":
             np.quantile(boots,[.025,.975]).tolist(),
         "source_positive_3deg_regions":int(np.sum(b<a)),
         "source_region_count":len(unique),
         "source_no_WWTP_to_land_causality":True}
    return {"spatial_holdout_n_sites":len(frame),"original_source_fold_results":fold_metrics,
            "original_source_spatial_year_is_2018":True,
            "fixed_2018_Ridge_alpha":ALPHA,
            "new_model_comparison":results}

def run(a):
    land,qa=source_frame(a.source_json)
    E,source=source_eqr()
    joined=E.merge(land,on="site_key",how="inner",validate="one_to_one")
    if len(joined)<200 or joined.countryCode.nunique()<6 or joined.region_3deg_block.nunique()<10:
        raise RuntimeError("Source valid + original numeric EQR joined cohort below minimum, STOP")
    if ((joined.countryCode!=joined.RZOriginalCountry).any() or
        (joined.region_3deg_block!=joined.RZOriginal3DegreeBlock).any()):
        raise RuntimeError("Real 2018 EQR and polygon source country/geoblock mismatch")
    exposed=joined["2018_measured_upstream_P_t_per_year"].astype(float)>0
    if len(joined.loc[exposed])<20:
        raise RuntimeError("Exact source-reported upstream N/P exposed support <20")
    joined["LogQ"]=np.log1p(joined.DIS_AV_CMS.to_numpy(float))
    joined["LogArea"]=np.log1p(joined.UPLAND_SKM.to_numpy(float))
    joined["LogUpP"]=np.log1p(joined["2018_measured_upstream_P_t_per_year"].to_numpy(float))
    joined["LogDownP"]=np.log1p(joined["2018_measured_downstream_P_t_per_year_placebo"].to_numpy(float))
    joined["LogUpN"]=np.log1p(joined["2018_measured_upstream_N_t_per_year"].to_numpy(float))
    joined["UpP_x_DeltaNat"]=joined["LogUpP"]*joined["NaturalNearMinusFar"]
    joined["countryMethod"]=joined.countryCode.astype(str)+"::"+joined.procedureClassificationSystem.astype(str)
    joined["spatial_fold"]=joined.pre_2018_source_only_spatial_fold.astype(int)
    if not np.isfinite(joined[[Y,YLAND]+NUM_LAND+["LogUpP","LogDownP","LogUpN"]].to_numpy(float)).all():
        raise RuntimeError("Ineligible original EQR or map source numeric data")
    result={"schema":SCHEMA,"status":"PASS_POSTHOC_SAME_2018_RIVER_BENTHIC_VS_RIPARIAN_LCLU_ASSOCIATION_NO_CAUSALITY",
            "original_EEA2018_EQR_SHA256":ORIG_EQR,
            "original_EEA2018_preoutcome_source_SHA256":ORIG_SOURCE,
            "external_2018_RZ_source_status":qa["status"],
            "original_river_benthic_real_numeric_EQR_population":len(E),
            "original_full_2018_RZ_EEA_source_only_2639":2639,
            "eligible_LCLU_with_2018_same_site_real_EQR":len(joined),
            "countries":int(joined.countryCode.nunique()),
            "independent_3degree_regions":int(joined.region_3deg_block.nunique()),
            "upstream_report_P_positive_2018_sites":int(exposed.sum()),
            "original_all_source_positive_report_not_clean_control":True,
            "mapped_habitat_index_not_actual_vegetation_biomass":True,
            "source_site_join_not_independent_WWTP_causal_path":True,
            "2018_exact_RZ_on_original_year_not_change_in_vegetation":True,
            "original_wastewater_P_N_reports_not_inriver_delivered_NP":True,
            "the_original_EQR_2018_source_associations_were_already_known":True,
            "posthoc_not_independently_confirmatory":True,
            "new_fixed_2018_spatial_cv_benthic_EQR":cv(joined,Y),
            "new_fixed_2018_spatial_cv_mapped_nearfar_terrestrial_class_contrast":
                cv(joined,YLAND,{
                "LandContext_L":NUM_LAND,
                "L_plus_UP_2018_reported_P":NUM_LAND+["LogUpP"],
                "L_plus_DOWNSTREAM_P_placebo":NUM_LAND+["LogDownP"]})
            }
    # Keep inspectable derived aggregate with no external third-party raw polygons.
    result["new_habitat_source_descriptives"]={
       k:{"min":float(joined[k].min()),"mean":float(joined[k].mean()),
          "max":float(joined[k].max())} for k in ["NaturalNear","NaturalFar",YLAND]
    }
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-json",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    if a.self_test:
        print(json.dumps(self_test()),flush=True);return
    a.output.mkdir(parents=True,exist_ok=True)
    try:result=run(a)
    except Exception as e:
        # Stop rather than returning false ecological null or claiming success.
        result={"schema":SCHEMA,"status":"STOP_SOURCE_OR_DATA_QUALITY_NO_ECOLOGY_CLAIM",
            "error_type":type(e).__name__,"error":str(e)[:450]}
    (a.output/"2018_EEA_SOURCE_BENTHIC_RZ2018_RIVER_NEARFAR_NATURAL_LCLU_CV_RESULTS.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    fmt=lambda d: {k:(v if not isinstance(v,dict) else {
        "n":v.get("spatial_holdout_n_sites"),
        "gains":{name:round(q["original_test_population_site_weighted_MSE_gain_percent"],5)
          for name,q in v.get("new_model_comparison",{}).items()}})
        for k,v in d.items() if k in ("status","eligible_LCLU_with_2018_same_site_real_EQR",
          "countries","independent_3degree_regions",
          "upstream_report_P_positive_2018_sites",
          "new_fixed_2018_spatial_cv_benthic_EQR",
          "new_fixed_2018_spatial_cv_mapped_nearfar_terrestrial_class_contrast",
          "error_type","error")}
    print("GESUI_2018_SAME_STATION_RIPARIAN_BENTHIC_RESULT="+json.dumps(fmt(result),sort_keys=True),flush=True)
if __name__=="__main__":main()
