#!/usr/bin/env python3
"""Strict source-only 232-H7 original W/environmental covariate replay; then 29 new H7.

No observed EQR, NUTRI, HNS, butterfly response or fit opened. The old
454-row panel is read ONLY for identity and independent source covariate
reconstruction. Any mismatch is a hard stop.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from pathlib import Path
from collections import defaultdict

PANEL_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
ATLAS_SHA="ce6780e01d2b98ed8e852574adb827cb3125be792a164e8be9e2ebb917b3a686"
NATIVE_H7=1219
SOURCE_NEW_SCHEMA="GESUI_SOURCE_ONLY_EXTRA30_H7_BOTH_PERIOD_EIV_HNS_TRAIT_SUPPORT_V1"
SOURCE_GEOG_SCHEMA="GESUI_EEA_SINGLEERA_H7_OFFICIAL_H6_GEO_SOURCE_OPPORTUNITY_V1"
SOURCE_MASTER_NAME="HYBAS7_1219x2_GBIF_PLANTAE_LYCAENIDAE_ANIMALIA_HYDROWASTE_SOURCE_ONLY_MASTER.csv"
SOURCE_ATLAS_NAME="HYDROATLAS_V1_L07_1219_H7_BASIN_RAW_LANDUSE_ATTRIBUTES_SOURCE_ONLY.json"
STATIC_FLOW="HYDROWASTE_STATIC_UPSTREAM_TO_BASIN_OUTLET_MODELED_WASTE_DIS_M3_PER_DAY"
UP_AREA="UPSTREAM_TOTAL_AREA_KM2"
LOCAL_AREA="LOCAL_SUBBASIN_AREA_KM2"
EPS=1e-10

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def one(folder,name):
    matches=list(folder.rglob(name))
    if len(matches)!=1:raise RuntimeError("Expected one genuine original source "+name+": "+str(matches))
    return matches[0]
def rows(path):
    with path.open(encoding="utf-8-sig",newline="") as h:return list(csv.DictReader(h))
def put(path,data):
    if not data:return
    with path.open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(data[0]))
        writer.writeheader();writer.writerows(data)

def main(a):
    original=Path("results/20261115_sync/hns_macro/H7_454_HNS_panel.csv")
    if sha(original)!=PANEL_SHA:raise RuntimeError("Original 454 frozen basin-period panel SHA drift")
    panel=rows(original)
    old_h={r["h"] for r in panel}
    if len(panel)!=454 or len(old_h)!=232 or len({r["h6"] for r in panel})!=90:
        raise RuntimeError("Original reference units unexpected")
    atlasfile=one(a.atlas,SOURCE_ATLAS_NAME)
    if sha(atlasfile)!=ATLAS_SHA:raise RuntimeError("Previously published 1219 original HydroATLAS SHA changed")
    atlas=json.loads(atlasfile.read_text())
    meta=atlas.get("original_1219_H7_covariate_by_real_ID")
    if not isinstance(meta,dict) or len(meta)!=NATIVE_H7:
        raise RuntimeError("Original genuine 1219-H7 native HydroATLAS map count drift")
    masterfile=one(a.master,SOURCE_MASTER_NAME)
    master=rows(masterfile)
    if len(master)!=2*NATIVE_H7 or len({(r["HYBAS_ID_L07"],r["PERIOD"]) for r in master})!=2*NATIVE_H7:
        raise RuntimeError("Original 2438 native outlet-source master not unique")
    by={(r["HYBAS_ID_L07"],r["PERIOD"]):r for r in master}
    allids={q[0] for q in by}
    if len(allids)!=NATIVE_H7 or set(meta)!=allids:
        raise RuntimeError("Original W/HydroATLAS native H7 source mapping coverage changed")
    if not old_h.issubset(allids):
        raise RuntimeError("Original 232 source-model H7 not included in 1219 native master")
    records={}
    for hid in allids:
        v=meta[hid]
        try:
            crop=float(v["crp_pc_use"])
            urban=float(v["urb_pc_use"])
            natural=float(v["dis_m3_pyr"])
            oldw={}
            for time in ("2017-2019","2020-2022"):
                row=by[(hid,time)]
                wload=float(row[STATIC_FLOW])
                up=float(row[UP_AREA])
                area=float(row[LOCAL_AREA])
                if not (0<=crop<=100 and 0<=urban<=100 and natural>=0 and
                        wload>=0 and up>0 and area>0 and
                        all(math.isfinite(z) for z in (crop,urban,natural,wload,up,area))):
                    raise ValueError("Nonphysical genuine original source covariate "+hid)
                oldw[time]={
                  "W":math.log1p(max(0.,wload)/(up+1.)),
                  "area":math.log(area),
                  "crop":crop,"urban":urban,
                  "naturalout":math.log1p(natural)}
            for variable in ("W","area","crop","urban","naturalout"):
                if abs(oldw["2017-2019"][variable]-oldw["2020-2022"][variable])>EPS:
                    raise RuntimeError("Original native model source exposure changes across static periods: "+hid+":"+variable)
            records[hid]=oldw["2017-2019"]
        except Exception as e:
            raise RuntimeError("Original genuine source 1219 model input cannot be reconstructed H7="+hid+" "+str(e))
    mismatches=[]
    worst={k:0. for k in ("W","area","crop","urban","naturalout")}
    for row in panel:
        recreated=records[row["h"]]
        for k in worst:
            delta=abs(recreated[k]-float(row[k]))
            worst[k]=max(worst[k],delta)
            if delta>EPS:
                mismatches.append({"h":row["h"],"period":row["period"],"column":k,
                    "original":float(row[k]),"replayed":recreated[k],"absolute_error":delta})
    old_replay_success=not mismatches
    # Score quality membership and source-only GIS rank. Do not open outcomes.
    gate=json.loads(one(a.traits,"EXTRA30_H7_SOURCE_TRAIT_QC.json").read_text())
    if gate.get("schema")!=SOURCE_NEW_SCHEMA or not gate.get("source_full_strict_40_H6_opportunity_gate_with_old_traits"):
        raise RuntimeError("Actual new original EIV/HNS score support source gate not passed")
    numeric_set=set(gate["source_new_H7_selected_bothperiod_trait_supported_ids"])
    if len(numeric_set)!=29:
        raise RuntimeError("Previously frozen EIV/HNS two-interval 29 candidate H7 changed")
    geog=json.loads(one(a.geography,"EEA_SINGLEERA_H7_OFFICIAL_H6_SOURCE_GEO_GATE.json").read_text())
    if geog.get("schema")!=SOURCE_GEOG_SCHEMA or geog.get("target_existing_singleton_H6_with_geographic_additional_EEA_H7")!=12:
        raise RuntimeError("Official original 12 new H6 source geography QA missing")
    geo=rows(one(a.geography,"EEA_SINGLEERA_NEW_H7_IN_ORIGINAL_SINGLETON_H6_SOURCE_REGISTRY.csv"))
    if len(geo)!=40 or len({r["native_H7"] for r in geo})!=40:
        raise RuntimeError("Original 40 genuine EEA geographically possible H7 drift")
    perparent=defaultdict(list)
    checked=[]
    for r in geo:
        if r["native_H7"] not in numeric_set:continue
        hid=r["native_H7"]
        if hid in old_h:raise RuntimeError("New source-only H7 was already in old 232 H7 sample")
        if hid not in records:raise RuntimeError("New source EEA H7 not in genuine prior source master")
        rec=records[hid]
        if not all(math.isfinite(rec[x]) for x in ("W","area","crop","urban","naturalout")):
            raise RuntimeError("Additional H7 original static hydrology covariate not finite")
        record={"native_H7":hid,"new_official_H6":r["official_native_H6"],
                "source_country_codes":r["original_source_countries"],
                "source_EEA_monitoring_eras":r["EEA_years"],
                "source_predeclared_geo_rank":int(r["source_rank_within_original_H6"]),
                "old_NUTRI_HNS_scored_two_periods":True,
                **rec,"not_same_old_EEA_two_year_monitor_model":True}
        checked.append(record)
        perparent[r["official_native_H6"]].append(record)
    if len(checked)!=29 or len(perparent)!=10:
        raise RuntimeError("Genuine 29 EIV/HNS supported H7 /10 H6 missing from 1219 covariate source")
    chosen=[]
    for h6 in sorted(perparent,key=int):
        a=sorted(perparent[h6],key=lambda x:x["source_predeclared_geo_rank"])
        chosen.append(a[0])
    if len(chosen)!=10:
        raise RuntimeError("Exactly ten previously source-eligible singleton H6 need one new H7")
    output={
        "schema":"GESUI_20261115_H7_1219_NATIVE_SOURCE_W_ENVIRONMENT_REPLAY_V1",
        "original_454_panel_SHA256":PANEL_SHA,
        "original_native_1219_H7_HydroATLAS_SHA256":ATLAS_SHA,
        "original_native_source_master_SHA256":sha(masterfile),
        "original_source_master_rows":len(master),
        "original_H7_source_master_unique":len(allids),
        "original_old_232_H7_environment_exact_reproduction":old_replay_success,
        "original_454_H7period_environment_comparisons":len(panel),
        "original_W_area_crop_urban_naturalout_max_abs_replay_diff":worst,
        "original_old_H7_mismatch_n":len(mismatches),
        "original_old_H7_mismatch_examples":mismatches[:10],
        "new_EIV_HNS_source_supported_H7":len(checked),
        "new_H7_with_all_original_W_geography_environment_fields":len(checked),
        "new_original_singleton_H6_covered":len(perparent),
        "new_single_H7_per_H6_ranked_source_selected":len(chosen),
        "new_original_selected_H7_per_old_singleton_H6":{x["new_official_H6"]:x["native_H7"] for x in chosen},
        "strict_new_model_environment_source_gate_pass":bool(old_replay_success and len(chosen)==10),
        "true_H6_geography_30_previous_plus_10_new_source_eligible":True,
        "no_original_EEA_same_method_station_year_or_original_nEEA_riverflow_replay":True,
        "static_W_is_basin_outlet_modelled_WWTP_flow_not_measured_receptor_dose":True,
        "no_plant_NUTRI_or_butterfly_HNS_effect_or_coefficients_fitted":True,
        "new_study_cohort_must_be_separate_from_original_454_study":True}
    a.output.mkdir(parents=True,exist_ok=True)
    (a.output/"H7_1219_SOURCE_REPLAY_W_ENVIRONMENT_QA.json").write_text(
        json.dumps(output,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    put(a.output/"EXTRA_29_H7_NATIVE_W_ENVIRONMENT_SOURCE_ONLY.csv",checked)
    put(a.output/"NEW_TEN_H6_DETERMINISTIC_SOURCE_ONE_H7_EACH.csv",chosen)
    print("GESUI_H7_1219_COVARIATE_REPLAY="+json.dumps(output,sort_keys=True,ensure_ascii=False),flush=True)
    print("NEXT_MODEL_SOURCE_GATE="+("PASS_STATIC_ENV_ONLY_NEW_DESIGN_REQUIRED" if output["strict_new_model_environment_source_gate_pass"] else "STOP_W_ENV_SOURCE_INCONSISTENT"),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--master",type=Path,required=True)
    p.add_argument("--atlas",type=Path,required=True)
    p.add_argument("--traits",type=Path,required=True)
    p.add_argument("--geography",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
