#!/usr/bin/env python3
"""Independent, no original API, exactly source-saved 169 H7 EU 2018→22
co-change arithmetic and hash basin-disjoint 5fold 4-model Ridge10.
All outcomes were computed in main original-data job; this audit cannot be
advertised as a new independent field replication.
"""
import collections,csv,hashlib,json,os,pathlib,math
import numpy as np
from scipy.stats import spearmanr

reportpath=pathlib.Path(os.environ["H7_REAL_EQR_PLANT_BUTTERFLY_SOURCE_REPORT"])
rowpath=pathlib.Path(os.environ["H7_REAL_EQR_PLANT_BUTTERFLY_SOURCE_ROWS"])
result=json.loads(reportpath.read_text())
if result["schema"]!="GESUI_H7_170_EU_ORIGINAL_FIXED_SOURCE_EQR_PLANT_S50_BUTTERFLY5_COHANGE_v1":
    raise RuntimeError("Original n169 study source schema changed")
with rowpath.open(encoding="utf8",newline="") as f:rows=list(csv.DictReader(f))
if len(rows)!=169 or len({x["HYBAS_ID_L07"] for x in rows})!=169:
    raise RuntimeError("Original true EQR×plant×butterfly source unit must be 169 H7 basins")
if set(x["country"] for x in rows)!=set(result["sample"]["countries"]):
    raise RuntimeError("Original country taxonomic strata changed")
Y=np.array([float(z["delta_EQR_2022_minus_2018"]) for z in rows],dtype=float)
P=np.array([float(z["delta_plant_speciesKey_rarefied_S50"]) for z in rows],dtype=float)
L=np.array([float(z["delta_lycaenidae_fivefamily_index_log_odds"]) for z in rows],dtype=float)
if int(sum(int(z["n_original_same_site_both_year_EQR"]) for z in rows))!=1023:
    raise RuntimeError("Original source-station pair counts changed")
spears={
 "delta_EQR_vs_delta_S50":float(spearmanr(Y,P).statistic),
 "delta_EQR_vs_delta_lycaenidae_logodds":float(spearmanr(Y,L).statistic),
 "delta_S50_vs_delta_lycaenidae_logodds":float(spearmanr(P,L).statistic)
}
for key,value in spears.items():
    expected=result["unadjusted_within_basin_co_response_Spearman"][key]["spearman_r"]
    if abs(value-expected)>1e-12:
        raise RuntimeError("Independent original source rank correlation not reproduced "+key)
base=result["original_predictor_covariates"]["B_recording_country_area_static_flow"]
country=sorted({z["country"] for z in rows})
def source_X(name):
    predictors=result["original_predictor_covariates"][name]
    return np.column_stack([np.ones(len(rows)),
       np.array([[float(z[k]) for k in predictors] for z in rows]),
       np.array([[float(z["country"]==c) for c in country[1:]] for z in rows])])
folds=np.array([int(hashlib.sha256(z["HYBAS_ID_L07"].encode()).hexdigest()[:8],16)%5 for z in rows])
def oos(x):
    out=np.empty(len(rows))
    for fold in range(5):
        train=folds!=fold;test=~train
        v=x.copy();mu=v[train,1:].mean(axis=0);sd=v[train,1:].std(axis=0)
        sd=np.where(sd>1e-12,sd,1.)
        v[:,1:]=(v[:,1:]-mu)/sd
        A=v[train];B=v[test];alpha=np.eye(A.shape[1])*10.;alpha[0,0]=0
        coef=np.linalg.solve(A.T@A+alpha,A.T@Y[train])
        out[test]=B@coef
    return float(np.mean((Y-out)**2))
mse={name:oos(source_X(name)) for name in result["original_predictor_covariates"]}
for name,n in mse.items():
    if not math.isclose(n,result["H7_disjoint_hash_5fold_ridge10_outcome_EQR_change_MSE"][name],rel_tol=1e-12,abs_tol=1e-12):
        raise RuntimeError("Independent source-pinned basin cross-validation differs for "+name)
C=collections.Counter(z["country"] for z in rows)
output={
 "schema":"H7_169EU_original_EQR_plant_S50_Lycaenidae_double-year_cochange_independent_qa_v1",
 "source_stage_original_sha":result["original_frozen_source_candidates_SHA256"],
 "H7_unique_basin_count":len(rows),
 "original_exact_stations":1023,
 "country_source_counts":dict(C),
 "majority_country_FR_fraction":C["FR"]/len(rows),
 "independent_exact_Spearman_recomputed":spears,
 "independent_all4_H7_hash_OOS_MSE_recomputed":mse,
 "all_coefficients_independent_script_reproduced":True,
 "not_independent_new_ecological_data":True,
 "true_source_2018_2022_method_same_station":True,
 "no_annual_actual_WWTP_discharge":True,
 "no_inferred_plant_pollination_ecosystem_service":True,
 "limitation":"This cannot support a universal serial wastewater→benthic→plant→pollinator mediation; 148/169 native H7 basins in France."
}
dst=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"H7_169EU_REPEAT_MACRO_DOUBLEYEAR_INDEPENDENT_QA.json"
dst.write_text(json.dumps(output,sort_keys=True,indent=2))
print("H7_169_EU_BASIN_MACRO_INDEPENDENT_RECHECK_START")
print(json.dumps(output,sort_keys=True))
print("H7_169_EU_BASIN_MACRO_INDEPENDENT_RECHECK_END")
