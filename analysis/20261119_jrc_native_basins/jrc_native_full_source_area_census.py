#!/usr/bin/env python3
"""Full native JRC soilN + pointN source attribute census, NOT an ecological fit.

No numeric EQR/EIV/HNS. Never infer BasinID==HydroID solely from matching
integer strings: mark native-area-matching records as geometric candidates only.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,time
from pathlib import Path
from datetime import datetime,timezone
from jrc_soilN_waterpoint_native_inventory import CONFIG,fetch,validate_meta

PAGE=500
PAGE_LIMIT=200000
REQUIRE_SOIL_COMPLETE=1000
MIN_PILOT_AREA_CANDIDATE=100

def valid_number(raw, nonnegative=True):
    try:x=float(raw)
    except (ValueError,TypeError):return None
    if not math.isfinite(x) or (nonnegative and x<0):return None
    return x

def all_layer(name):
    base,layer,primary,value,areafield=CONFIG[name]
    endpoint=f"{base}/{layer}"
    meta,meta_sha=fetch(endpoint,{})
    validate_meta(meta,name,primary,value,areafield)
    obj,idhash=fetch(endpoint+"/query",{"where":"1=1","returnIdsOnly":"true"})
    counts,countsha=fetch(endpoint+"/query",{"where":"1=1","returnCountOnly":"true"})
    ids=obj.get("objectIds")
    n=counts.get("count")
    if not isinstance(ids,list) or type(n)is not int or len(ids)!=n or n>PAGE_LIMIT or len(set(ids))!=n:
        raise ValueError("Original full ID/count inventory failed for "+name)
    allids=sorted(ids)
    keys={}
    objects={}
    pagehash=[]
    nulls=0
    for i in range(0,len(allids),PAGE):
        batch=allids[i:i+PAGE]
        o,stamp=fetch(endpoint+"/query",{
          "objectIds":",".join(map(str,batch)),
          "outFields":",".join(["OBJECTID",primary,value,areafield]),
          "returnGeometry":"false"
        })
        data=[x.get("attributes",{}) for x in o.get("features",[])]
        if len(data)!=len(batch) or {z.get("OBJECTID") for z in data}!=set(batch) or o.get("exceededTransferLimit"):
            raise ValueError(f"Truncated page {name} {i}: {len(data)} vs {len(batch)}")
        for z in data:
            key=z.get(primary)
            if type(key) is not int or key in keys:
                raise ValueError(f"Ambiguous native key {name}: {key}")
            v=valid_number(z.get(value))
            area=valid_number(z.get(areafield))
            if area is None or area<=0:area=None
            if v is None:nulls+=1
            keys[key]=(v,area)
            objects[key]=z["OBJECTID"]
        pagehash.append(stamp)
        if (i//PAGE)%10==0 or i+PAGE>=n:
            print("JRC_NATIVE_FULL_ATTRIBUTE_PROGRESS",name,len(keys),n,flush=True)
        time.sleep(.04)
    if len(keys)!=n:raise ValueError("Incomplete original attribute inventory: "+name)
    manifest={
      "official_layer":endpoint,"original_features":n,
      "unique_native_keys":len(keys),"valid_model_value":n-nulls,
      "native_source_value_missing":nulls,
      "native_zero_count":sum(v==0 for v,_ in keys.values()),
      "area_nonmissing":sum(a is not None for _,a in keys.values()),
      "meta_sha256":meta_sha,"objectids_sha256":idhash,
      "count_sha256":countsha,
      "attribute_page_count":len(pagehash),
      "attribute_page_hashes_sha256":hashlib.sha256("\n".join(pagehash).encode()).hexdigest(),
      "native_ID_field":primary,"units_field":value,"area_field":areafield,
    }
    return keys,manifest

def compatible_areas(a,b):
    if a is None or b is None or a<=0 or b<=0:return False
    return abs(a-b)<=max(.05,.005*min(a,b))

def check_land_and_candidates(tables):
    names=["soil_atmosN_2015","soil_mineralN_2010_19",
           "soil_organicN_2010_19","soil_fixationN_2010_19"]
    primary=names[:3]
    ks=[set(tables[name]) for name in names]
    common4=set.intersection(*ks)
    common3=set.intersection(*ks[:3])
    point=tables["water_pointN_2014_18"]
    n_all=0;n_joint=0
    candidate=[]
    land_area_inconsistent=0
    for hid in sorted(common3):
        values=[tables[n][hid] for n in primary]
        measured=[x[0] for x in values]
        if any(v is None for v in measured):continue
        n_all+=1
        basearea=values[0][1]
        if not all(compatible_areas(basearea,x[1]) for x in values):
            land_area_inconsistent+=1
            continue
        # IMPORTANT: point BasinID==HydroID is only a candidate until polygons agree.
        if hid not in point:continue
        psn,parea=point[hid]
        if psn is None or not compatible_areas(basearea,parea):continue
        n_joint+=1
        values4=tables[names[3]].get(hid,(None,None))
        candidate.append({
          "candidate_ID_number_not_proof_of_same_basin":hid,
          "soil_hydroid":hid,"water_basinid":hid,
          "soil_area_sqkm":basearea,"point_area_sqkm":parea,
          "soil_atmos_N_mean_mg_m2":measured[0],
          "soil_mineral_N_mean_kg_ha":measured[1],
          "soil_organic_N_mean_kg_ha":measured[2],
          "soil_fixation_N_mean_kg_ha":values4[0],
          "water_point_N_modelled_t":psn,
          "N_source_2015_2010to19_2014to18_times_different":True
        })
    support={
      "soil_atmos_native_hydroID":len(ks[0]),
      "soil_mineral_native_hydroID":len(ks[1]),
      "soil_organic_native_hydroID":len(ks[2]),
      "soil_fixation_native_hydroID":len(ks[3]),
      "soil_3main_native_hydroID_intersection_before_missing":len(common3),
      "soil_4_native_hydroID_intersection_before_missing":len(common4),
      "soil_3main_values_nonmissing_hydroID":n_all,
      "soil_intersection_area_conflicts":land_area_inconsistent,
      "matching_integer_ID_and_native_area_only_geographic_candidates":n_joint,
      "point_and_soil_fully_geospatially_proved_same_HydroID":False,
      "source_census_source_only_go_for_geometry_audit":n_all>=REQUIRE_SOIL_COMPLETE and n_joint>=MIN_PILOT_AREA_CANDIDATE,
      "source_regional_correlation_not_computed_until_true_polygon_validation":True,
      "old_25NUTS2_6countries_biotic_joint_STOP_remains":True
    }
    return candidate,support

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output-dir",type=Path,required=True)
    a=p.parse_args()
    a.output_dir.mkdir(parents=True,exist_ok=True)
    data={};manifest={}
    for name in CONFIG:
        data[name],manifest[name]=all_layer(name)
        print("JRC_NATIVE_COMPLETE_SOURCE",name,json.dumps(manifest[name]),flush=True)
    candidates,support=check_land_and_candidates(data)
    out=a.output_dir/"JRC_NATIVE_POSSIBLE_SHARED_GEOMETRY_ONLY_SOURCE.csv"
    with out.open("w",newline="") as f:
        if candidates:
            wr=csv.DictWriter(f,fieldnames=list(candidates[0]))
            wr.writeheader()
            wr.writerows(candidates)
        else:f.write("no_source_shared_geometry_candidates\n")
    report={
      "schema":"gesui_JRC_native_riverbasin_HydroID_BasinID_FULL_SOURCE_ONLY_CENSUS_20261119_v1",
      "timestampUTC":datetime.now(timezone.utc).isoformat(),
      "layers":manifest,
      "source_geometry_match_precheck":support,
      "candidateCSVsha256":hashlib.sha256(out.read_bytes()).hexdigest(),
      "candidate_spatial_identity_unproven_until_original_GIS_polygon_overlay":True,
      "no_EEA_EQR_EIV_NUTRI_HNS_outcome_opened":True
    }
    (a.output_dir/"JRC_NATIVE_FIVE_LAYER_FULL_SOURCE_GATE.json").write_text(
      json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("JRC_NATIVE_FIVE_SOURCE_GEOMETRY_CANDIDATE_GATE",json.dumps(support,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
