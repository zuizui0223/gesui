#!/usr/bin/env python3
"""First 16-outfall SOURCE-FROZEN directed river vs upstream Plantae GBIF study.

River geometry selected from pinned HydroWASTE/HydroRIVERS BEFORE inspecting GBIF.
This estimates differences in GBIF-recorded plant-species indices at river-reach
centers, not nutrient load, plant abundance, or a causal sewage effect.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime,timezone
import hashlib
import json
import math
from pathlib import Path
import statistics
import time
from urllib.parse import urlencode
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError

import numpy as np
from scipy.stats import spearmanr

SOURCE_RUN_ID=42364866919
SOURCE_OUTFALLS=16
SOURCE_PARQUET_SHA256="5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f"
SOURCE_HYDROWASTE_GIT_BLOB="0506f873f6541c5f324478894a9731290d722bce"
SOURCE_COUNTRIES=("AUT","CZE","DEU","DNK")
GBIF="https://api.gbif.org/v1/occurrence/search"
YEAR="2011,2024"
RADIUS_KM=1
PRIMARY_DEPTH=30
SENSITIVITY_DEPTH=10
MIN_SPECIES_PRIMARY=6
MIN_SPECIES_SENSITIVITY=4
MAX_SITE_TOTAL_RECORD_RATIO_PRIMARY=10.0
AGENT="gesui-directed-river-plantae-index-source-blind/0.1"

def sha256(path:Path)->str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def strict_source_manifest(design:dict,qc:dict):
    assert design["metadata"]["sample_selected_without_GBIF_response"] is True
    assert design["metadata"]["outfall_count"]==SOURCE_OUTFALLS
    assert qc["sources"]["hydrowaste_csv_git_blob"]==SOURCE_HYDROWASTE_GIT_BLOB
    assert qc["sources"]["hydrorivers_europe_shp_zip_sha256"]==SOURCE_PARQUET_SHA256
    assert qc["quality"]["selected_outfall_units"]==SOURCE_OUTFALLS
    assert qc["quality"]["hydrowaste_eligible_river_ids_found_in_europe_mirror"]==4713
    assert qc["sources"]["third_party_geoparquet_version_verified_against_official"] is False
    selected=design["selected"]
    assert len(selected)==SOURCE_OUTFALLS
    assert {v["country"] for v in selected}==set(SOURCE_COUNTRIES)
    assert {v["exposure_stratum"] for v in selected}=={"high","low"}
    assert len({v["wwtp_id"] for v in selected})==SOURCE_OUTFALLS
    assert len({v["river_id"] for v in selected})==SOURCE_OUTFALLS
    assert all(v["river_hydraulic_chain"]=="U2>U1>OUTFALL_REACH>D1" for v in selected)
    return selected

def site_registry(selected:list[dict]):
    sites=[]
    for w in selected:
        for role in ("u2","u1","d1"):
            sites.append({
                "site_id":str(w["wwtp_id"])+"_"+role,
                "wwtp_id":int(w["wwtp_id"]),
                "country":w["country"],
                "hydrobasin_l12":str(w["hydrobasin_l12"]),
                "river_role":role,
                "river_reach_id":int(w[role+"_river_id"]),
                "latitude":float(w[role+"_lat"]),
                "longitude":float(w[role+"_lon"]),
                "source_fraction":float(w["effluent_fraction_proxy"]),
                "source_fraction_stratum":w["exposure_stratum"],
                "distance_u2_u1_km":float(w["distance_u2_u1_km"]),
                "distance_u1_d1_km":float(w["distance_u1_d1_km"]),
                "river_ratio_u1d1_over_u2u1":float(w["distance_u1_d1_km"]/w["distance_u2_u1_km"]),
                "outfall_water_flow_m3day":float(w["estimated_ww_m3_day"]),
                "estimated_river_discharge_m3sec":float(w["estimated_river_m3_sec"]),
                "observed_N_mass_kg_year":None,
                "site_radius_km":RADIUS_KM,
                "gbif_year_range":YEAR,
                "source_run_id":SOURCE_RUN_ID,
            })
    assert len(sites)==SOURCE_OUTFALLS*3
    return sites

def gbif_query(params:dict):
    url=GBIF+"?"+urlencode(params)
    last=None
    for attempt in range(4):
        try:
            r=Request(url,headers={"User-Agent":AGENT,"Accept":"application/json"})
            with urlopen(r,timeout=65) as resp:
                raw=resp.read()
            obj=json.loads(raw)
            if "count" not in obj or "error" in obj:
                raise ValueError("GBIF response lacks count or has error")
            return obj
        except (HTTPError,URLError,TimeoutError,ValueError) as err:
            last=err
            if attempt!=3:time.sleep(3*(attempt+1))
    raise RuntimeError("GBIF source call failed: "+str(last))

def query_all_plantae(site:dict):
    params={
        "taxonKey":6, "year":YEAR,
        "geoDistance":f"{site['latitude']:.7f},{site['longitude']:.7f},{RADIUS_KM}km",
        "hasCoordinate":"true","hasGeospatialIssue":"false",
        "occurrenceStatus":"PRESENT","facet":"speciesKey","facetLimit":1000,"limit":0,
    }
    count=None
    keys={}
    pages=0
    while True:
        params["facetOffset"]=pages*1000
        data=gbif_query(params)
        new=int(data["count"])
        if count is None:count=new
        if new!=count:raise ValueError("GBIF record count changed between facet pages")
        facets=next((x.get("counts",[]) for x in data.get("facets",[])
                    if str(x.get("field","")).upper() in ("SPECIES_KEY","SPECIESKEY")),[])
        if count>0 and not facets and pages==0:
            raise ValueError("GBIF omitted species facet even though taxa are recorded")
        for item in facets:
            k=int(item["name"])
            if k in keys:raise ValueError("Species key duplicated between facet pages")
            keys[k]=int(item["count"])
        pages+=1
        if len(facets)<1000:break
        if pages>12:raise ValueError("Too many distinct species keys for fixed site radius")
        time.sleep(.12)
    n=sum(keys.values())
    if n>count:raise ValueError("GBIF species-keyed records exceed all Plantae records")
    return {**site,
        "total_plantae_occurrence_records":int(count),
        "indexed_species_record_count":int(n),
        "distinct_species_keys":len(keys),
        "species_key_record_counts":keys,
        "facet_pages":pages,
        "fraction_records_with_species_key":n/count if count>0 else None,
        "retrieved_utc":datetime.now(timezone.utc).isoformat()
    }

def expected_rarefaction(records,depth):
    c=[int(v) for v in records if int(v)>0]
    n=sum(c)
    if n<depth:return None
    def logcomb(n,k):
        return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)
    d=logcomb(n,depth)
    return float(sum(1-math.exp(logcomb(n-x,depth)-d) if n-x>=depth else 1 for x in c))

def jaccard(a:dict,b:dict):
    A=set(a);B=set(b)
    if not A and not B:return None
    return 1-len(A&B)/len(A|B)

def level(site:dict,depth:int,min_species:int):
    return site["indexed_species_record_count"]>=depth and site["distinct_species_keys"]>=min_species

def summarise(selected,observations):
    by={x["site_id"]:x for x in observations}
    output=[]
    for w in selected:
        wid=w["wwtp_id"]
        u2,u1,d1=[by[str(wid)+"_"+role] for role in ("u2","u1","d1")]
        rec=[z["indexed_species_record_count"] for z in (u2,u1,d1)]
        ratio=max(rec)/min(rec) if min(rec)>0 else float("inf")
        depth_results={}
        for depth,min_s in ((PRIMARY_DEPTH,MIN_SPECIES_PRIMARY),
                            (SENSITIVITY_DEPTH,MIN_SPECIES_SENSITIVITY)):
            eligible=all(level(z,depth,min_s) for z in (u2,u1,d1))
            vals=[expected_rarefaction(z["species_key_record_counts"].values(),depth)
                  if eligible else None for z in (u2,u1,d1)]
            if eligible:
                early=math.log1p(vals[1])-math.log1p(vals[0])
                late=math.log1p(vals[2])-math.log1p(vals[1])
                d2=late-early
                jl=jaccard(u1["species_key_record_counts"],d1["species_key_record_counts"])
                je=jaccard(u2["species_key_record_counts"],u1["species_key_record_counts"])
            else:early=late=d2=jl=je=None
            depth_results[str(depth)]={
                "eligible_by_records_and_species":eligible,
                "rarefied_expected_plant_species_u2_u1_d1":vals,
                "pre_outfall_upstream_log_richness_change":early,
                "across_outfall_downstream_log_richness_change":late,
                "extra_change_downstream_minus_prior_upstream":d2,
                "jaccard_dissimilarity_u2_u1":je,
                "jaccard_dissimilarity_u1_d1":jl,
                "turnover_increment_across_minus_upstream":jl-je if eligible else None,
            }
        primary=(depth_results[str(PRIMARY_DEPTH)]["eligible_by_records_and_species"]
                 and ratio<=MAX_SITE_TOTAL_RECORD_RATIO_PRIMARY)
        output.append({
            "wwtp_id":wid,"country":w["country"],
            "river_id":w["river_id"],
            "source_fraction":float(w["effluent_fraction_proxy"]),
            "source_fraction_stratum":w["exposure_stratum"],
            "distance_u2_u1_km":w["distance_u2_u1_km"],
            "distance_u1_d1_km":w["distance_u1_d1_km"],
            "site_species_key_record_counts":rec,
            "site_record_balance_ratio_maxmin":ratio if math.isfinite(ratio) else None,
            "strict_primary_eligible":bool(primary),
            "depth30":depth_results[str(PRIMARY_DEPTH)],
            "depth10_sensitivity":depth_results[str(SENSITIVITY_DEPTH)],
        })
    qualified=[x for x in output if x["strict_primary_eligible"]]
    permissive=[x for x in output if x["depth30"]["eligible_by_records_and_species"]]
    low=[x for x in output if x["depth10_sensitivity"]["eligible_by_records_and_species"]]
    def associations(xs,metric):
        if len(xs)<8 or len({x["country"] for x in xs})<3:
            return {"sufficient_support":False,"n":len(xs)}
        fr=[x["source_fraction"] for x in xs]
        y=[x[metric]["extra_change_downstream_minus_prior_upstream"] for x in xs]
        if len(set(y))<3:return {"sufficient_support":False,"n":len(xs)}
        r,p=spearmanr(fr,y)
        return {"sufficient_support":True,"n":len(xs),"spearman_rho":float(r),
                "two_sided_p_descriptive":float(p)}
    def high_low(xs,k):
        h=[x[k]["extra_change_downstream_minus_prior_upstream"]
           for x in xs if x["source_fraction_stratum"]=="high"]
        l=[x[k]["extra_change_downstream_minus_prior_upstream"]
           for x in xs if x["source_fraction_stratum"]=="low"]
        if not h or not l:return None
        return {"n_high":len(h),"n_low":len(l),
                "mean_high":float(np.mean(h)),"mean_low":float(np.mean(l)),
                "high_minus_low":float(np.mean(h)-np.mean(l))}
    return {"wwtp_outcomes":output,
            "site_count":len(observations),
            "source_selected_wwtps":len(selected),
            "all_sites_any_plantae_records":sum(s["total_plantae_occurrence_records"]>0 for s in observations),
            "primary_depth30_and_balance_n":len(qualified),
            "unbalanced_depth30_n":len(permissive),
            "sensitivity_depth10_n":len(low),
            "primary_source_fraction_association":associations(qualified,"depth30"),
            "depth30_all_balances_exploratory":associations(permissive,"depth30"),
            "depth10_association_sensitivity":associations(low,"depth10_sensitivity"),
            "primary_high_vs_low":high_low(qualified,"depth30"),
            "depth10_high_vs_low":high_low(low,"depth10_sensitivity"),
            "count_2011_2024_sites_by_country":dict(__import__("collections").Counter(x["country"] for x in observations)),
            "sample_window":YEAR,
            "limits":[
                "GBIF presence-only occurrence records, not local plant abundance or verified absence",
                "Kingdom Plantae includes non-vascular plants; no N trait or origin classification",
                "HydroRIVERS third-party GeoParquet version unverified against original official",
                "River reach centers differ in segment lengths, elevations, urbanization and recording effort",
                "Other WWTPs upstream, unrelated land use, irrigation or nitrogen may explain differences",
                "No measured wastewater N/P, delivered load, flow seasonality or causal attribution",
                "16 source-selected outfalls only, 4 countries and non-independent catchments possible",
                "Raw Jaccard turnover is sensitive to record effort and site distances",
            ]}

def run(args):
    design=json.loads(args.design.read_text())
    qc=json.loads(args.qc.read_text())
    selected=strict_source_manifest(design,qc)
    sites=site_registry(selected)
    args.out_dir.mkdir(parents=True,exist_ok=True)
    frozen={"source_geometry_original_run":SOURCE_RUN_ID,
            "geometry_sha256":sha256(args.design),
            "qc_sha256":sha256(args.qc),
            "date":"2026-11-07","taxon_kingdomKey":6,
            "year":YEAR,"radius_km":RADIUS_KM,
            "primary_depth":PRIMARY_DEPTH,
            "primary_min_species":MIN_SPECIES_PRIMARY,
            "primary_record_max_min_ratio":MAX_SITE_TOTAL_RECORD_RATIO_PRIMARY,
            "sensitivity_depth":SENSITIVITY_DEPTH,
            "site_units":sites}
    (args.out_dir/"preplant_sites_and_estimands.json").write_text(
        json.dumps(frozen,indent=2,sort_keys=True)+"\n")
    print("GESUI_DIRECTIONAL_FROZEN_SAMPLE",
          json.dumps({"outfalls":len(selected),"sites":len(sites),
                      "primary_rarefaction_records":PRIMARY_DEPTH,
                      "source_run":SOURCE_RUN_ID},sort_keys=True),flush=True)
    observations=[]
    for i,site in enumerate(sites,1):
        result=query_all_plantae(site)
        observations.append(result)
        print("GESUI_RIVER_PLANTAE_SITE",
              json.dumps({"index":i,"id":site["site_id"],"records":result["total_plantae_occurrence_records"],
                          "species":result["distinct_species_keys"],
                          "species_keyed_records":result["indexed_species_record_count"]},sort_keys=True),flush=True)
        time.sleep(.3)
    all_original={"schema":"gesui_source_fixed_directed_reach_gbif_plantae_facets_v01",
                  "source_run":SOURCE_RUN_ID,"year":YEAR,"radius_km":RADIUS_KM,
                  "facets":observations}
    (args.out_dir/"plantae_raw_species_facets.json").write_text(
        json.dumps(all_original,separators=(",",":"),ensure_ascii=False)+"\n")
    with (args.out_dir/"plantae_site_registry.csv").open("w",newline="") as out:
        cols=list(sites[0].keys())+["total_plantae_occurrence_records",
           "indexed_species_record_count","distinct_species_keys","fraction_records_with_species_key"]
        writer=csv.DictWriter(out,fieldnames=cols)
        writer.writeheader()
        for r in observations:
            writer.writerow({k:r.get(k) for k in cols})
    with (args.out_dir/"plantae_site_species_records.csv").open("w",newline="") as out:
        cols=["site_id","country","wwtp_id","river_role","river_reach_id",
              "gbif_speciesKey","speciesKey_record_count","gbif_year_range","source_run_id"]
        writer=csv.DictWriter(out,fieldnames=cols)
        writer.writeheader()
        for s in observations:
            for k,count in sorted(s["species_key_record_counts"].items()):
                writer.writerow({"site_id":s["site_id"],"country":s["country"],
                    "wwtp_id":s["wwtp_id"],"river_role":s["river_role"],
                    "river_reach_id":s["river_reach_id"],"gbif_speciesKey":k,
                    "speciesKey_record_count":count,"gbif_year_range":YEAR,"source_run_id":SOURCE_RUN_ID})
    summary=summarise(selected,observations)
    summary["source_run"]=SOURCE_RUN_ID
    summary["status"]="OBSERVATIONAL_GBIF_ALL_PLANTAE_DIRECTED_REACH_NOT_WASTEWATER_EFFECT"
    (args.out_dir/"directional_results.json").write_text(
        json.dumps(summary,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GESUI_RIVER_PLANTAE_DIRECTIONAL_SUMMARY",
          json.dumps({k:v for k,v in summary.items() if k not in ("wwtp_outcomes","limits")},
                     sort_keys=True),flush=True)
    return summary

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--design",type=Path,required=True)
    ap.add_argument("--qc",type=Path,required=True)
    ap.add_argument("--out-dir",type=Path,required=True)
    args=ap.parse_args()
    run(args)

if __name__=="__main__":
    main()
