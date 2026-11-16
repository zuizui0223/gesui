#!/usr/bin/env python3
"""SOURCE ONLY: 273 genuine 2022 measured river TP vs original EEA2022
benthic same-reach <=250m France I2M2 station pairs, 74 HydroBASINS L7,
verify native source parent H6 group heldout feasibility BEFORE EQR/TP values.
"""
import collections,hashlib,json,os,pathlib
from sklearn.model_selection import GroupKFold

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_FR2022_TP_EQR_H6_source"
ROOT.mkdir(parents=True,exist_ok=True)
PAIR_SHA="448add2d6850cbc876c3e6b90bca278c6b711e64f9d363071e93bc0c18b84611"
H6_SHA="372b800b85abe2f395241987cea2000ff901cf5aff31bf1d1fb825304f2495fc"
ATLAS_SHA="ce6780e01d2b98ed8e852574adb827cb3125be792a164e8be9e2ebb917b3a686"
def load(name,sha):
    b=pathlib.Path(os.environ[name]).read_bytes()
    if hashlib.sha256(b).hexdigest()!=sha:raise RuntimeError("Original source immutable SHA mismatch "+name)
    return json.loads(b)
match=load("ORIGINAL_FR_TP_EQR_250M_SOURCE",PAIR_SHA)
h6=load("ORIGINAL_H6_PARENT_232_SOURCE",H6_SHA)
atlas=load("ORIGINAL_ATLAS_L7_1219_SOURCE",ATLAS_SHA)
p=match["original_454_H7_cohort_2022_genuine_TP_250m_source_matched_pairs"]
if len(p)!=273 or len({z["EEA2022_site_key"] for z in p})!=273 or len({z["WQ_ID"] for z in p})!=273:
    raise RuntimeError("Original source 273 pair count or one-to-one quality changed")
ids={z["EEA_H7_L07"] for z in p}
if len(ids)!=74 or any(z["EEA2022_country"]!="FR" or z["EEA2022_BQE_original_method"]!="I2M2" for z in p):
    raise RuntimeError("Original French I2M2 source 74 H7/273 sites did not match")
parents={q["HYBAS_ID_L07"]:q["official_parent_HYBAS_ID_L06"]
 for q in h6["original_454_H7_study_basin_geographic_parent_mapping"]}
if len(parents)!=232 or not ids.issubset(parents):
    raise RuntimeError("Original Level6 parent for exact 2022 TP-benthic H7 unavailable")
A=atlas["original_1219_H7_covariate_by_real_ID"]
if len(A)!=1219 or not ids.issubset(A):raise RuntimeError("Original HydroATLAS source not same native H7")
perH7=collections.Counter(z["EEA_H7_L07"] for z in p)
counts=collections.Counter(parents[id] for id in ids)
src=[{"HYBAS_ID_L07":id,
       "HYBAS_ID_L06_original_parent":parents[id],
       "matched_original_TP_EEA_station_pairs_2022":perH7[id],
       "source_cropland_field_present":A[id].get("crp_pc_use") is not None,
       "source_urban_field_present":A[id].get("urb_pc_use") is not None,
       "source_natural_flow_present":A[id].get("dis_m3_pyr") is not None} for id in sorted(ids)]
if any(not all((r["source_cropland_field_present"],r["source_urban_field_present"],r["source_natural_flow_present"])) for r in src):
    raise RuntimeError("Original land covariate missing for source TP-EEA sites")
# GroupKFold 5 groups predeclared. Order H7 IDs stable before any ecology.
ks=sorted(ids)
groups=[parents[id] for id in ks]
assign={}
for fi,(train,test) in enumerate(GroupKFold(n_splits=5).split(ks,groups=groups)):
    if len(test)<4 or len(train)<30:
        raise RuntimeError("Original Level6 hydrological parent group split under-supported; STOP")
    for i in test:assign[ks[i]]=fi
if len(assign)!=74:raise RuntimeError("TP/EQR original H7 assigned to more than one source-only fold")
fold_counts=collections.Counter(assign.values())
fold_parent=collections.defaultdict(set)
for id in ks:fold_parent[assign[id]].add(parents[id])
passgate=len(counts)>=10 and len(fold_counts)==5 and min(fold_counts.values())>=4
qa={
 "schema":"GESUI_FR_I2M2_2022_REAL_TP_EEA_EQR_273PAIRS_74H7_H6_GEO_SOURCE_GATE_v1",
 "original_source_matched_station_SHA256":PAIR_SHA,
 "original_source_H6_native_parent_SHA256":H6_SHA,
 "original_source_BasinATLAS_1219_SHA256":ATLAS_SHA,
 "matched_same_2022_TP_EEA_field_station_pairs":len(p),
 "original_FR_I2M2_study_basins":len(ids),
 "independent_native_H6_parent_groups":len(counts),
 "n_H7_basins_per_original_H6_parent":dict(counts),
 "pre_declared_fivefold_Source_H6_hydrologic_GroupKFold_H7_counts":dict(fold_counts),
 "pre_declared_fivefold_Source_H6_hydrologic_GroupKFold_parent_counts":{k:len(v) for k,v in fold_parent.items()},
 "original_source_colocation_250m_fixed":True,
 "source_STUDY_country_FR_method_I2M2_exclusively":True,
 "source_geometry_gate_PASS":passgate,
 "the_numeric_river_TP_2022_median_and_benthic_EQR_not_accessed":True,
 "source_only_guarantees_no_country_transfer_claim":True,
 "source_only_endpoint_ecological_coefficient_not_fitted":True,
 "source_H6_geometry_had_16_original_H7_assignments_REPAIRED_after_source_failure":True,
 "source_H6_geo_fold_assignment_independent_of_biological_outcome":True
}
manifest={"schema":qa["schema"],"matched_original_74_H7_fold_source":[
 {"H7":id,"H6":parents[id],"fold":assign[id],"n_2022_same_station_pairs":perH7[id]} for id in ks],
 "original_TP_EQR_numeric_not_used_to_freeze_fold":True}
raw=(json.dumps(manifest,sort_keys=True,separators=(",",":"))+"\n").encode()
qa["source_original_H6_geo_fold_assignment_SHA256"]=hashlib.sha256(raw).hexdigest()
(ROOT/"FR2022_74H7_I2M2_TP_EQR_SOURCE_FROZEN_H6_GROUPOUT_FOLDS.json").write_bytes(raw)
(ROOT/"FR2022_74H7_I2M2_TP_EQR_SOURCE_GEO_H6_QA.json").write_text(json.dumps(qa,sort_keys=True,indent=2))
print("FR2022_TP_EQR_H6_SOURCE_ONLY_START",flush=True)
print(json.dumps(qa,sort_keys=True),flush=True)
print("FR2022_TP_EQR_H6_SOURCE_ONLY_END",flush=True)
if not passgate:
    raise RuntimeError("Predeclared original French 74-H7 TP/I2M2 EQR H6 geo source gate FAIL - no biological predictive analysis")
