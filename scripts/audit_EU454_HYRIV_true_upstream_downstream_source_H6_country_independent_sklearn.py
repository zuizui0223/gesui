#!/usr/bin/env python3
"""Independent sklearn Ridge verification of genuine source-locked receiver
upstream vs downstream wastewater candidate on the original n454 H7 periods.

Reopen exact public-derived source manifests, *not* original licensed raw
HydroRIVERS ZIP. Separate implementation from NumPy Ridge original code,
same 90 H6 parent groups, 11 stable country groups, all 3 original outcome
families. Diagnose non-varying operational-only index and source asymmetry.
"""
import collections,csv,hashlib,json,math,os,pathlib,statistics
import numpy as np
from sklearn.linear_model import Ridge
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_HYRIV_receiver_independent"
ROOT.mkdir(parents=True,exist_ok=True)
SHAS={
 "receivers":"d50eba0f934278cb28eb403af7b84d4da5a19af209a6d05be312a2ee44ad6d6d",
 "h6":"372b800b85abe2f395241987cea2000ff901cf5aff31bf1d1fb825304f2495fc",
 "atlas":"ce6780e01d2b98ed8e852574adb827cb3125be792a164e8be9e2ebb917b3a686",
 "np":"012ad79feabc8c11e0d74f803730b5b7913b96a64a6af26fb0c97263d9aa7c26"
}
def load(k,sha=None):
    b=pathlib.Path(os.environ[k]).read_bytes()
    if sha is not None and hashlib.sha256(b).hexdigest()!=sha:
        raise RuntimeError("Original public source SHA mismatch "+k)
    return json.loads(b)
rec=load("RECEPTOR_SOURCE_454",SHAS["receivers"])
p=load("H6_PARENT_SOURCE_MAP",SHAS["h6"])
a=load("H7_ATLAS_SOURCE",SHAS["atlas"])["original_1219_H7_covariate_by_real_ID"]
m=load("REPORTED_NP_SOURCE_454",SHAS["np"])
orig_model=load("ORIGINAL_TRUE_RECEIVER_3LAYER_MODEL_RESULTS")
if orig_model["schema"]!="GESUI_TRUE_EEA_RECEPTOR_HYRIV_UPSTREAM_50KM_VS_BASIN_OUTLET_3BIOTA_MACRO_POSTHOC_20261113_v1":
    raise RuntimeError("Original source-station receptor study model schema unknown")
with open(os.environ["ORIGINAL_H7_EEA_BUTTERFLY_PLANT_ROWS"],encoding="utf-8-sig",newline="") as f:
    rows=list(csv.DictReader(f))
if len(rows)!=454 or len({r["HYBAS_ID_L07"] for r in rows})!=232:
    raise RuntimeError("Original primary real 454/232 source response population changed")
recs={(r["HYBAS_ID_L07"],r["PERIOD"]):r for r in rec["H7_year_original_derived_receiver_site_mean_exposures"]}
np_by={(r["HYBAS_ID_L07"],r["PERIOD"]):r for r in m["source_H7_period_reported_N_P_upstream_to_true_original_station"]}
h6={r["HYBAS_ID_L07"]:r["official_parent_HYBAS_ID_L06"] for r in p["original_454_H7_study_basin_geographic_parent_mapping"]}
if len(recs)!=454 or len(np_by)!=454 or len(h6)!=232:raise RuntimeError("Original native receptor geography not one-to-one")
country=sorted({r["country"] for r in rows})
method=sorted({r["original_majority_BQE_method"] for r in rows})
base_num=[
 "log_local_H7_area_km2","log1p_original_EQR_site_n",
 "log1p_GBIF_Animalia_source_records","log1p_GBIF_Plantae_source_records",
 "log1p_GBIF_5family_butterfly_source_records","log1p_original_EEA_site_DIS_AV_CMS"]
def fixed(k):
    return [float(r[k]) for r in rows]
c=np.array([[float(r["country"]==x) for x in country[1:]] for r in rows])
me=np.array([[float(r["original_majority_BQE_method"]==x) for x in method[1:]] for r in rows])
period=np.array([float(r["PERIOD"]=="2020-2022") for r in rows])
common=np.column_stack([np.array([[float(r[k]) for k in base_num] for r in rows]),me,period])
atlas_crop=np.array([float(a[r["HYBAS_ID_L07"]]["crp_pc_use"]) for r in rows])
atlas_urban=np.array([float(a[r["HYBAS_ID_L07"]]["urb_pc_use"]) for r in rows])
atlas_flow=np.log1p(np.array([float(a[r["HYBAS_ID_L07"]]["dis_m3_pyr"]) for r in rows]))
env=np.column_stack([atlas_crop,atlas_urban,atlas_flow])
# Country columns in the same order as previous manual implementation
full_base=np.column_stack([np.array([[float(r[k]) for k in base_num] for r in rows]),c,me,period,env])
country_out_base=np.column_stack([common,env])
outlet=np.array([float(r["W_static_upstream_MODEL_flow_per_UpArea"]) for r in rows])
u=[];u50=[];dn=[];ann=[];op=[];qual=[];same=[];valid_count=[]
for r in rows:
    z=recs[(r["HYBAS_ID_L07"],r["PERIOD"])]
    assert z["network_source_valid"] and z["original_usable_network_receiver_stations"]>=1
    def v(k):
        q=float(z["station_equal_mean_"+k])
        if q<0 or not math.isfinite(q):raise RuntimeError("Original source invalid receiver derived pressure "+k)
        return q
    u.append(math.log1p(v("strict_upstream_WWTP_fraction_of_original_model_river_discharge")))
    u50.append(math.log1p(v("strict_upstream_50km_WWTP_fraction_of_original_model_river_discharge")))
    dn.append(math.log1p(v("downstream_50km_WWTP_fraction_direction_control")))
    same.append(v("source_original_same_reach_facilities_unresolved_COUNT"))
    f=math.expm1(float(r["log1p_original_EEA_site_DIS_AV_CMS"]))
    if f<=0:raise RuntimeError("Source EEA original river flow invalid")
    op.append(math.log1p(v("upstream_50km_original_STATUS_operational_model_flow_m3_day")/(f*86400)))
    qual.append(math.log1p(v("upstream_50km_original_QUAL_WASTE_1_model_flow_m3_day")/(f*86400)))
    ann.append(math.log1p(float(np_by[(r["HYBAS_ID_L07"],r["PERIOD"])]["source_year_reported_N_mgL_equivalent_site_equal_mean"])))
    valid_count.append(z["original_usable_network_receiver_stations"])
u=np.array(u);u50=np.array(u50);dn=np.array(dn);ann=np.array(ann);op=np.array(op);qual=np.array(qual)
assert len(u)==454
H6=np.array([h6[r["HYBAS_ID_L07"]] for r in rows])
H6fold=np.array([int(hashlib.sha256(x.encode()).hexdigest()[:8],16)%5 for x in H6])
cc=collections.defaultdict(collections.Counter)
for r in rows:cc[r["HYBAS_ID_L07"]][r["country"]]+=1
countries={h:sorted(vals,key=lambda x:(-vals[x],x))[0] for h,vals in cc.items()}
country_fold=np.array([countries[r["HYBAS_ID_L07"]] for r in rows])
def xall(B):
    d={
      "B2_environment_only":B,
      "B3_original_H7_basin_OUTLET_W":np.column_stack([B,outlet]),
      "B4_original_H7_basin_OUTLET_W_plus_EEA_annual_N":np.column_stack([B,outlet,ann]),
      "B2_plus_receiver_STRICT_all_upstream":np.column_stack([B,u]),
      "B2_plus_receiver_STRICT_50km":np.column_stack([B,u50]),
      "B2_plus_DOWNSTREAM_50km_direction_control":np.column_stack([B,dn]),
      "B2_plus_receiver_all_AND_old_basin_outlet":np.column_stack([B,u,outlet]),
      "B2_plus_receiver_all_AND_downstream_control":np.column_stack([B,u,dn]),
      "B2_plus_receiver_all_AND_annual_reported_N":np.column_stack([B,u,ann]),
      "B2_plus_receiver50_AND_operational_only50":np.column_stack([B,u50,op]),
      "B2_plus_receiver50_AND_QUAL_WASTE1_reported_only50":np.column_stack([B,u50,qual]),
      "B2_plus_old_basin_W_AND_DOWNSTREAM_control":np.column_stack([B,outlet,dn])
    }
    return d
# This implementation does not construct an explicit intercept column; sklearn
# handles the unpenalized intercept. Original code handled same intercept in
# its normal equations.
models_geo=xall(full_base)
models_country=xall(country_out_base)
def sklearn_cv(X,y,groups):
    yh=np.empty(len(y))
    for group in sorted(set(groups)):
        test=groups==group;train=~test
        q=X[train];q2=X[test]
        mean=q.mean(axis=0);sd=q.std(axis=0)
        sd=np.where(sd>1e-12,sd,1.)
        A=(q-mean)/sd;B=(q2-mean)/sd
        reg=Ridge(alpha=10.,fit_intercept=True,solver="svd")
        reg.fit(A,y[train])
        yh[test]=reg.predict(B)
    return float(np.mean((y-yh)**2))
outcomes=[
 ("original_EEA_benthic_EQR","EQR_original_site_equal_basin_period_mean"),
 ("GBIF_index_plant_S50","plant_S50_GBIF_speciesKey"),
 ("GBIF_index_Lycaenidae_five_family_logodds","lycaenidae_fivefamily_index_share_logodds")
]
audit={}
for label,col in outcomes:
    y=np.array([float(r[col]) for r in rows])
    row={}
    for source_label,models,fold,priorfield in [
        ("true_H6_block",models_geo,H6fold,"H6_parent90_heldout_all_models_MSE"),
        ("source_country_block",models_country,country_fold,"country_leftout_all_models_MSE")]:
        r={}
        old=orig_model["original_source_crosslayer_results"][label][priorfield]
        for k,X in models.items():
            m=sklearn_cv(X,y,fold)
            if not math.isclose(m,float(old[k]),rel_tol=1e-8,abs_tol=1e-10):
                raise RuntimeError("Independent sklearn outflow source MSE mismatch "+
                                   label+"/"+source_label+"/"+k+": "+str(m)+"/"+str(old[k]))
            r[k]=m
        row[source_label]={"all_12_original_source_models_exactly_reproduced":True,
             "independent_sklearn_MSE_by_model":r}
    audit[label]=row
summary={
 "schema":"GESUI_TRUE_HYRIV_454_SOURCE_RECEPTOR_FLOW_H6_COUNTRY_SKLEARN_INDEPENDENT_QA_20261113_v1",
 "original_receiver_source_SHA256":SHAS["receivers"],
 "original_H6_parent_geosource_SHA256":SHAS["h6"],
 "original_1219_HydroATLAS_attributes_SHA256":SHAS["atlas"],
 "original_454_year_specific_EEA_reported_N_SHA256":SHAS["np"],
 "original_H7_basins":len(h6),"same_original_basin_period_rows":len(rows),
 "actual_2018_2022_EEA_sites_original_study_n":sum(valid_count),
 "original_source_H6_parent_group_count":len(set(H6)),
 "original_stable_country_group_count":len(set(country_fold)),
 "source_true_UPSTREAM_positive_454_source_periods":sum(u>0),
 "source_true_UPSTREAM_50km_positive_454_source_periods":sum(u50>0),
 "source_DOWNSTREAM_50km_positive_454_source_periods":sum(dn>0),
 "source_SAME_REACH_ambiguity_nonzero_454_source_periods":sum(x>0 for x in same),
 "source_original_explicitly_operational_50km_flow_nonzero_454_periods":int(sum(op>0)),
 "source_original_QUAL_WASTE1_50km_flow_nonzero_454_periods":int(sum(qual>0)),
 "source_receiver_UPSTREAM_logflow_logDOWNSTREAM_pearson_r":float(np.corrcoef(u,dn)[0,1]) if np.std(u)>0 and np.std(dn)>0 else None,
 "source_receiver_UPSTREAM_logflow_old_basin_OUTLET_logflow_pearson_r":float(np.corrcoef(u,outlet)[0,1]) if np.std(u)>0 and np.std(outlet)>0 else None,
 "all_3_outcomes_12_models_H6_and_country_original_MSE_reproduced_with_independent_sklearn":True,
 "source_original_crossresponse_numeric_QA":audit,
 "note":"Independent code and saved original source only, NOT new field sampling or a causal randomized downstream negative control. Operational-only 0 indexed source is not evidence of no functioning treatment plants, because STATUS mostly unreported."
}
path=ROOT/"EU454_ORIGINAL_HYRIV_RECEIVER_WWTP_3BIOTA_H6_COUNTRY_INDEPENDENT_SKLEARN_QA.json"
path.write_text(json.dumps(summary,indent=2,sort_keys=True,allow_nan=False))
brief={k:v for k,v in summary.items() if k not in ("source_original_crossresponse_numeric_QA","note")}
print("HYRIV_RECEPTOR_SKLEARN_INDEPENDENT_QA_START")
print(json.dumps(brief,sort_keys=True,allow_nan=False),flush=True)
print("HYRIV_RECEPTOR_SKLEARN_INDEPENDENT_QA_END")
