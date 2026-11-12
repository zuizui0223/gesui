#!/usr/bin/env python3
"""Independent (separate code) source-locked 454 H7 European repeated
basin-period ecological environmental confounding and FR/era sensitivity.
No new original field observations, all subgroup diagnostics posthoc.
"""
import csv,hashlib,json,os,pathlib,collections,math
import numpy as np
from scipy.stats import spearmanr

root=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_H7_454_LANDUSE_INDEPENDENT_QA"
root.mkdir(parents=True,exist_ok=True)
SHAS={
 "atlas":"ce6780e01d2b98ed8e852574adb827cb3125be792a164e8be9e2ebb917b3a686",
 "reported":"012ad79feabc8c11e0d74f803730b5b7913b96a64a6af26fb0c97263d9aa7c26"}
atlas=pathlib.Path(os.environ["ORIGINAL_SOURCE_1219_HYDROATLAS"])
mass=pathlib.Path(os.environ["ORIGINAL_SOURCE_454_REPORTED_NP"])
for name,path in [("atlas",atlas),("reported",mass)]:
    if hashlib.sha256(path.read_bytes()).hexdigest()!=SHAS[name]:
        raise RuntimeError("Original source SHA mismatch "+name)
a=json.loads(atlas.read_text())["original_1219_H7_covariate_by_real_ID"]
p={(x["HYBAS_ID_L07"],x["PERIOD"]):x for x in json.loads(mass.read_text())["source_H7_period_reported_N_P_upstream_to_true_original_station"]}
with open(os.environ["ORIGINAL_454_BIOLOGICAL_ROWS"],encoding="utf-8-sig",newline="") as f:
    R=list(csv.DictReader(f))
main=json.loads(pathlib.Path(os.environ["ORIGINAL_H7_LAND_ADJUSTED_MAIN_RESULT"]).read_text())
if main["schema"]!="GESUI_EU454_H7_ORIGINAL_THREE_RESPONSE_REPORTED_NP_LANDCOVER_ADJUSTED_v1":
    raise RuntimeError("Original 3 layer land use adjusted source result changed")
if len(R)!=454 or len(p)!=454 or len(a)!=1219:
    raise RuntimeError("Original source unit sizes changed")
for x in R:
    h=x["HYBAS_ID_L07"];pr=x["PERIOD"]
    if h not in a or (h,pr) not in p:raise RuntimeError("Original H7 or period source keys missing")
    fields=["crp_pc_use","urb_pc_use","dis_m3_pyr","crp_pc_sse","urb_pc_sse"]
    if any(a[h][f] is None for f in fields):raise RuntimeError("HydroATLAS genuine source missing original covariate")
    x["crop_up"]=float(a[h]["crp_pc_use"])
    x["urban_up"]=float(a[h]["urb_pc_use"])
    x["naturalout"]=math.log1p(float(a[h]["dis_m3_pyr"]))
    x["crop_local"]=float(a[h]["crp_pc_sse"])
    x["urban_local"]=float(a[h]["urb_pc_sse"])
    x["reportN"]=math.log1p(float(p[(h,pr)]["source_year_reported_N_mgL_equivalent_site_equal_mean"]))
    x["positiveNP"]=int(p[(h,pr)]["source_station_annual_upstream_NP_positive_n"]>0)
if sum(x["positiveNP"] for x in R)!=56:raise RuntimeError("Original 56 source-positive H7×period values not consistent")
if len({x["HYBAS_ID_L07"] for x in R})!=232:raise RuntimeError("Original 232 source basin frame drift")

def model(sub):
    N=len(sub)
    basins=[r["HYBAS_ID_L07"] for r in sub]
    if N<80 or len(set(basins))<35:return {"status":"insufficient","n":N,"H7":len(set(basins))}
    countries=sorted({r["country"] for r in sub})
    methods=sorted({r["original_majority_BQE_method"] for r in sub})
    X0=np.column_stack([
      np.ones(N),
      np.array([[float(r[k]) for k in (
        "log_local_H7_area_km2","log1p_original_EQR_site_n",
        "log1p_GBIF_Animalia_source_records","log1p_GBIF_Plantae_source_records",
        "log1p_GBIF_5family_butterfly_source_records","log1p_original_EEA_site_DIS_AV_CMS")] for r in sub]),
      np.array([[float(r["country"]==c) for c in countries[1:]] for r in sub]),
      np.array([[float(r["original_majority_BQE_method"]==m) for m in methods[1:]] for r in sub]),
      np.array([float(r["PERIOD"]=="2020-2022") for r in sub])
    ])
    w=np.array([float(r["W_static_upstream_MODEL_flow_per_UpArea"]) for r in sub])
    E=np.array([[r["crop_up"],r["urban_up"],r["naturalout"]] for r in sub])
    Enn=np.array([r["reportN"] for r in sub])
    X={
      "B0":X0,
      "B1":np.column_stack([X0,w]),
      "B2":np.column_stack([X0,E]),
      "B3":np.column_stack([X0,E,w]),
      "B4":np.column_stack([X0,E,w,Enn])
    }
    fold=np.array([int(hashlib.sha256(h.encode()).hexdigest()[:8],16)%5 for h in basins])
    if min(np.bincount(fold,minlength=5))<4:return {"status":"fold support too low"}
    def fit(A,Y):
        pred=np.empty(N)
        for f in range(5):
            test=fold==f;train=~test
            Z=A.copy()
            mu=Z[train,1:].mean(axis=0);sd=Z[train,1:].std(axis=0)
            sd=np.where(sd>1e-12,sd,1.)
            Z[:,1:]=(Z[:,1:]-mu)/sd
            a=Z[train];b=Z[test]
            reg=np.eye(a.shape[1])*10.;reg[0,0]=0.
            coef=np.linalg.solve(a.T@a+reg,a.T@Y[train])
            pred[test]=b@coef
        return float(np.mean((Y-pred)**2))
    output={"n_source_period":N,"n_distinct_H7":len(set(basins)),
       "source_positive_NP_periods":sum(r["positiveNP"] for r in sub),
       "source_countries":len(countries),
       "FR_periods":sum(r["country"]=="FR" for r in sub),"outcomes":{}}
    for outcome,col in [
      ("EQR","EQR_original_site_equal_basin_period_mean"),
      ("S50","plant_S50_GBIF_speciesKey"),
      ("Lycaenidae","lycaenidae_fivefamily_index_share_logodds")
    ]:
        y=np.array([float(r[col]) for r in sub])
        mse={name:fit(A,y) for name,A in X.items()}
        output["outcomes"][outcome]={
         "five_model_original_MSE":mse,
         "H7_hash_5fold_gain_B2_upstream_environment_to_B3_plus_static_W_pct":100*(mse["B2"]-mse["B3"])/mse["B2"],
         "H7_hash_5fold_gain_B3_to_B4_plus_annual_reported_N_pct":100*(mse["B3"]-mse["B4"])/mse["B3"]}
    return output
out={"schema":"GESUI_H7_454_HYDROATLAS_SOURCE_INDEPENDENT_FR_PERIOD_NUMERIC_REPLAY_v1",
     "original_source_landuse_SHA":SHAS["atlas"],
     "original_source_annual_NP_SHA":SHAS["reported"],
     "all_results_share_one_original_source_population":True,
     "not_independent_new_ecological_data":True,
     "subgroups":{}}
subsets={
 "all_original_454":R,
 "P1_2017_2019":[r for r in R if r["PERIOD"]=="2017-2019"],
 "P2_2020_2022":[r for r in R if r["PERIOD"]=="2020-2022"],
 "FR_only":[r for r in R if r["country"]=="FR"],
 "outside_FR":[r for r in R if r["country"]!="FR"],
 "no_cross_country_or_mixed_method":[r for r in R if int(r["original_country_count_in_basin_year"])==1 and int(r["original_method_count_in_basin_year"])==1]
}
for label,data in subsets.items():out["subgroups"][label]=model(data)
# Strong independent reproduction gate for all three responses and 5 exact prior source models
main_res=main["separate_layer_result_models"]
keymap={"EQR":"EEA_original_benthic_EQR",
        "S50":"GBIF_original_plant_speciesKey_S50",
        "Lycaenidae":"GBIF_original_Lycaenidae_fivefamily_logodds"}
mkey={"B0":"B0_source_recording_country_method_flow",
       "B1":"B1_plus_STATIC_HydroWASTE",
       "B2":"B2_plus_UPSTREAM_crop_urban_natural_flow",
       "B3":"B3_UPSTREAM_env_plus_STATIC_HydroWASTE",
       "B4":"B4_UPSTREAM_env_STATIC_plus_ANNUAL_reported_N"}
allcalc=out["subgroups"]["all_original_454"]["outcomes"]
for label,k in keymap.items():
    for short,long in mkey.items():
        source=main_res[k]["all_11_fixed_model_fivefold_MSE"][long]
        new=allcalc[label]["five_model_original_MSE"][short]
        if not math.isclose(source,new,rel_tol=1e-11,abs_tol=1e-12):
            raise RuntimeError("Independent original data EQR/S50/Lycaenidae environmental source CV mismatch "+label+":"+short)
out["independently_recomputed_15_original_full_numeric_CV_MSE_values_exact"]=True
for indicator,key in [("cropland", "crop_up"),("urban","urban_up")]:
    v=[r[key] for r in R]
    w=[float(r["W_static_upstream_MODEL_flow_per_UpArea"]) for r in R]
    out["source_static_WWTP_index_raw_rho_"+indicator]=float(spearmanr(v,w).statistic)
out["posthoc_country_and_year_subgroups_are_not_new_independent_field_replications"]=True
out["no_actual_measured_receiver_TNTP_or_services"]=True
dest=root/"H7_454_HYDROATLAS_LAND_WWTP_NP_INDEPENDENT_COUNTRY_ERAS_QA.json"
dest.write_text(json.dumps(out,sort_keys=True,indent=2,allow_nan=False))
print("H7_454_HYDROATLAS_INDEPENDENT_QA_START")
print(json.dumps(out,sort_keys=True,allow_nan=False))
print("H7_454_HYDROATLAS_INDEPENDENT_QA_END")
