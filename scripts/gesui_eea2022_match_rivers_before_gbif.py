#!/usr/bin/env python3
"""Pre-plant source-only matched EU rivers: actual reported EEA2022 P mass.

No GBIF, no measured river chemistry or plant outcomes inspected. Control
indicates zero REGISTERED EEA both-positive reported N/P load within 50km
hydrologically upstream (not sewage-free river).
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter
from pathlib import Path
import pandas as pd
SOURCE_SHA="79641d984b28fdbc11dbf02e32f32cda2818ec256f219439a1ccfae601da7dc3"
PLOAD="eea2022_upstream_P_measured_load_nativeunit"
EARTH_KM=6371.0088

def hav(lat1,lon1,lat2,lon2):
    a,b=math.radians(float(lat1)),math.radians(float(lat2))
    u=math.sin((b-a)/2)**2+math.cos(a)*math.cos(b)*math.sin(math.radians(float(lon2)-float(lon1))/2)**2
    return 2*EARTH_KM*math.asin(math.sqrt(max(0.,min(1.,u))))

def station(r):
    return {
        "station_id":str(r["WQ.ID"]),"latitude":float(r["WQ.lat"]),
        "longitude":float(r["WQ.lon"]),"river_id":int(r["HYRIV_ID"]),
        "block":str(r["spatial_3deg_block"]),
        "flow_m3s":float(r["DIS_AV_CMS"]),
        "upstream_area_km2":float(r["UPLAND_SKM"]),
        "cropland_fraction":float(r["crop2"]),
        "population_density":float(r["pop_dens"]),
        "eea2022_reported_measured_upstream_P_ton_yr":float(r[PLOAD]),
        "hydrowaste_prior_upstream_flow_volume":
            float(r["upstream_50km_ww_m3day"]) if "upstream_50km_ww_m3day" in r else None
    }

def run(args):
    digest=hashlib.sha256(args.source.read_bytes()).hexdigest()
    if digest!=SOURCE_SHA:raise ValueError("Pre-outcome EEA source table hash changed")
    full=pd.read_csv(args.source,dtype={"WQ.ID":str},low_memory=False)
    if len(full)!=1712 or full[["nutrient","WQ.ID"]].duplicated().any():
        raise ValueError("Source-station registry has unexpected size or duplicates")
    d=full[full.nutrient=="TP"].copy()
    assert len(d)==1418
    needed=["WQ.lat","WQ.lon","DIS_AV_CMS","UPLAND_SKM","crop2","pop_dens",PLOAD]
    for x in needed:d[x]=pd.to_numeric(d[x],errors="coerce")
    d=d[d[needed].notna().all(axis=1)&d.DIS_AV_CMS.gt(0)&d.UPLAND_SKM.gt(0)&
        d.pop_dens.ge(0)&d[PLOAD].ge(0)].copy()
    yes=d[d[PLOAD]>0].copy()
    no=d[d[PLOAD]==0].copy()
    if (len(d),len(yes))!=(510,46):
        raise ValueError("Source-only matching denominator changed")
    candidates=[]
    for _,a in yes.iterrows():
        controls=no[no.spatial_3deg_block==a.spatial_3deg_block]
        for _,b in controls.iterrows():
            if a.HYRIV_ID==b.HYRIV_ID:continue
            dist=hav(a["WQ.lat"],a["WQ.lon"],b["WQ.lat"],b["WQ.lon"])
            if not 5<=dist<=180:continue
            flow=abs(math.log(a.DIS_AV_CMS/b.DIS_AV_CMS))
            area=abs(math.log(a.UPLAND_SKM/b.UPLAND_SKM))
            crop=abs(float(a.crop2-b.crop2))
            pop=abs(math.log1p(a.pop_dens)-math.log1p(b.pop_dens))
            if flow>math.log(3) or area>math.log(3) or crop>.2 or pop>1.5:
                continue
            score=flow+area+2*crop+.5*pop+dist/200
            candidates.append((score,str(a["WQ.ID"]),str(b["WQ.ID"]),
                               dist,flow,area,crop,pop,a,b))
    candidates.sort(key=lambda z:(z[0],z[1],z[2]))
    used=set();locations=[];pairs=[]
    for score,aid,bid,dist,flow,area,crop,pop,a,b in candidates:
        if len(pairs)>=24:break
        if aid in used or bid in used:continue
        if any(hav(p["WQ.lat"],p["WQ.lon"],v["WQ.lat"],v["WQ.lon"])<5
               for p in (a,b) for v in locations):continue
        pair={
          "pair_id":"EUP22_"+str(len(pairs)+1).zfill(2),
          "source_exposed":station(a),
          "source_control":station(b),
          "spatial_block":str(a.spatial_3deg_block),
          "between_site_distance_km":float(dist),
          "source_matching_score":float(score),
          "covariate_distance":{"log_flow_ratio":float(flow),
               "log_upstream_area_ratio":float(area),
               "absolute_crop_fraction":float(crop),
               "log1p_population_density":float(pop)}
        }
        pairs.append(pair);used.update((aid,bid));locations.extend((a,b))
    if len(pairs)<12:raise ValueError("Fewer than 12 preregistered source-only pairs")
    summary={
      "schema":"gesui_EEA2022_source_positive_zero_plant_new_blind_match_v01",
      "source_only_registry_sha256":SOURCE_SHA,
      "source_positive_definition":"2022_EEA_measured_reported_P_tonnes_upstream_50km",
      "source_zero_definition":"NO_REPORTED_BOTH_POSITIVE_EEA_MEASURED_SOURCE_NOT_NATURAL_ZERO",
      "source_eligible_TP_stations":1418,"source_crop_population_complete":len(d),
      "source_positive_crop_population_complete":len(yes),
      "source_zero_crop_population_complete":len(no),
      "source_pair_candidates":len(candidates),
      "pairs_fixed":len(pairs),"plant_sites_fixed":len(pairs)*2,
      "geographic_3deg_blocks":dict(Counter(p["spatial_block"] for p in pairs)),
      "source_only_selection":True,
      "observed_Peters_2022_TP_concentrations_used_for_selection":False,
      "previous_GBIF_plant_community_used_for_selection":False,
      "new_GBIF_window":"2022,2024","new_GBIF_buffer_km":1,
      "pair_constraints":{"max_flow_ratio":3,"max_upstream_area_ratio":3,
         "max_absolute_crop_fraction_difference":.2,
         "max_log_population_density_difference":1.5,
         "pair_distance_range_km":[5,180],
         "all_site_center_minimum_distance_km":5},
      "primary_new_plant_record_depth":30,"sensitivity_depth":10,
      "not_a_causal_plant_or_2022_P_discharge_effect":True
    }
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/"SOURCE_MATCHED_PAIR_SELECTION.json").write_text(
       json.dumps({"source_design":summary,"pairs":pairs},sort_keys=True,indent=2)+"\n")
    (args.out/"MATCHING_QC.json").write_text(json.dumps(summary,sort_keys=True,indent=2)+"\n")
    print("GESUI_2022_EEA_P_SOURCE_MATCHED_BEFORE_GBIF",json.dumps(summary,sort_keys=True),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    run(p.parse_args())
