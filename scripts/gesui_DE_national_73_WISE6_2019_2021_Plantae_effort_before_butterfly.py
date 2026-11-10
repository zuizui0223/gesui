#!/usr/bin/env python3
"""Freeze 2019-2021 plant photographic opportunity for original 73 DE source/zero WISE6 river pairs.

No butterfly observations (or river TP numeric concentrations) fetched here.
Reuses exact 32 prior sites, queries only the 68 new original station locations,
both preregistered 2km and 5km radii. Stage-0 source cohort must be on GitHub
before this script executes. Keep all unqualified source pairs visible.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import time
from collections import Counter, defaultdict
from datetime import datetime,timezone
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.parse import urlencode
from urllib.request import Request,urlopen

ROOT_FRAME = "results/lycaenidae_2022_DE_national_73pair_sourceonly/ORIGINAL_2022_73DE_WISE6_SOURCEZERO_PAIR_100SITE_PREPHOTO_FRAME.json"
PRIOR_PLANT = "results/lycaenidae_2022_wise6_directed_32_station_photo_sourceonly/2019_2021_ORIGINAL32_2KM_5KM_PLANT_ONLY_PREBUTTERFLY_GATE.json"
RADII=(2,5)
START="2019-01-01"
END="2021-12-31"
TAXON=47126
GERMANY=7207
SRC="registered_source_interposed"
ZERO="registered_qualifying_sourcezero"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def geo(a,b):
    r=math.pi/180
    dt=(a["lat"]-b["lat"])*r
    dl=(a["lon"]-b["lon"])*r
    z=math.sin(dt/2)**2+math.cos(a["lat"]*r)*math.cos(b["lat"]*r)*math.sin(dl/2)**2
    return 12742*math.asin(min(1,math.sqrt(z)))

def query(site,radius):
    q={
        "taxon_id": TAXON, "place_id":GERMANY,
        "lat":round(site["lat"],5),"lng":round(site["lon"],5),"radius":radius,
        "d1":START,"d2":END,
        "photos":"true","geo":"true",
        "geoprivacy":"open","taxon_geoprivacy":"open",
        "quality_grade":"research","captive":"false","acc_below":500,"per_page":1
    }
    url="https://api.inaturalist.org/v1/observations?"+urlencode(q)
    err=None
    for attempt in range(4):
        try:
            req=Request(url,headers={"User-Agent":"gesui-national-WISE6-2022-source-only-plantae/0.1 (zuizui0223/gesui)",
                                     "Accept":"application/json"})
            with urlopen(req,timeout=45) as f:
                raw=f.read()
                code=f.status
            parsed=json.loads(raw)
            n=parsed["total_results"]
            if code!=200 or type(n) is not int or n<0:
                raise ValueError("bad original API aggregate response")
            return {"status":200,"total_results":n,"original_JSON_SHA256":hashlib.sha256(raw).hexdigest(),
                    "original_public_API_url":url,"origin":"NEW_2026_11_10"}
        except (HTTPError,URLError,ValueError,TimeoutError) as e:
            err=str(e)
            if attempt!=3:time.sleep(5*(attempt+1))
    return {"status":"FAILED","total_results":None,"error":err,
            "original_public_API_url":url,"origin":"NEW_2026_11_10"}

def old_count(record,radius,site):
    """Use original 32-site identical time, taxon, circle filters and immutable URL."""
    p=record["independent_previous_year_plantae_photo_counts"][f"{radius}km"]
    url=p["original_API_url"]
    required=(f"taxon_id={TAXON}",f"place_id={GERMANY}",f"radius={radius}",
              f"d1={START}",f"d2={END}","acc_below=500",
              "photos=true","quality_grade=research")
    if not all(text in url for text in required):
        raise ValueError("Original 32-site historical source was not same query protocol")
    if abs(float(record["latitude"])-site["lat"])>1e-8:
        raise ValueError("Source-locked reused latitude mismatch")
    if abs(float(record["longitude"])-site["lon"])>1e-8:
        raise ValueError("Source-locked reused longitude mismatch")
    if p["status"]!=200 or type(p["count_public_2019_2021_research_grade_georeferenced_plantae_photos"]) is not int:
        raise ValueError("Prior 32-site source API did not succeed")
    return {
        "status":200,
        "total_results":int(p["count_public_2019_2021_research_grade_georeferenced_plantae_photos"]),
        "original_JSON_SHA256":p["original_JSON_response_SHA256"],
        "original_public_API_url":url,
        "origin":"REUSED_OLD32_2026_11_09"
    }

def qualifies_pair(pair,stations,obs,rad):
    u,d=pair["upstream_original_site_id_SHA256"],pair["downstream_original_site_id_SHA256"]
    uu,dd=obs[u][f"{rad}km"]["total_results"],obs[d][f"{rad}km"]["total_results"]
    distance=geo(stations[u],stations[d])
    okay=isinstance(uu,int) and isinstance(dd,int)
    ratio=min(uu,dd)/max(uu,dd) if okay and max(uu,dd)>0 else 0
    passed=bool(okay and distance>2*rad and min(uu,dd)>=5 and ratio>=.15)
    return {"eligible":passed,"upstream_plant_photos":uu,"downstream_plant_photos":dd,
            "within_pair_plant_record_ratio":ratio,
            "original_monitor_geodesic_km":distance,"nonoverlap_up_down":distance>2*rad}

def pair_ids(pair):
    return pair["upstream_original_site_id_SHA256"],pair["downstream_original_site_id_SHA256"]

def fit_nonoverlap(candidate,occupied,stations,rad):
    return all(geo(stations[s],stations[old])>2*rad for s in pair_ids(candidate) for old in occupied)

def choose_national(pairs,stations,rad):
    """Fully source/plant effort eligible, non-overlapping site-disjoint independent pairs."""
    ranked=sorted(pairs,key=lambda p:(p["network_full_reach_path_km"],p["pair_id_SHA256"]))
    selected=[]
    used=[]
    for p in ranked:
        if not fit_nonoverlap(p,used,stations,rad):continue
        selected.append(p)
        used.extend(pair_ids(p))
    return selected

def midpair(p,stations):
    u,d=pair_ids(p)
    return {"lat":(stations[u]["lat"]+stations[d]["lat"])/2,
            "lon":(stations[u]["lon"]+stations[d]["lon"])/2}

def choose_quartets(pairs,stations,rad):
    """Greedy fixed source/land/photo-effort feasible same-3° four-arm matches."""
    source=[x for x in pairs if x["exposure_class"]==SRC]
    zero=[x for x in pairs if x["exposure_class"]==ZERO]
    options=[]
    for p in source:
        for q in zero:
            if p["upstream_3deg_region"]!=q["upstream_3deg_region"]:continue
            if any(geo(stations[s],stations[t])<=2*rad for s in pair_ids(p) for t in pair_ids(q)):continue
            k=abs(math.log(p["network_full_reach_path_km"]/q["network_full_reach_path_km"]))
            md=geo(midpair(p,stations),midpair(q,stations))
            options.append((k,md,p["pair_id_SHA256"],q["pair_id_SHA256"],p,q))
    options.sort(key=lambda x:x[:4])
    chosen=[]
    used=[]
    for loss,mid,pid,qid,p,q in options:
        ids=list(pair_ids(p))+list(pair_ids(q))
        if not all(geo(stations[a],stations[b])>2*rad for a,b in itertools.combinations(ids,2)):
            continue
        if any(geo(stations[a],stations[b])<=2*rad for a in ids for b in used):continue
        chosen.append({"source_pair_SHA256":pid,"registeredzero_pair_SHA256":qid,
                       "original_3deg_block":p["upstream_3deg_region"],
                       "river_path_log_length_mismatch":loss,
                       "source_control_midpoint_distance_km":mid,
                       "four_station_hashes":ids,
                       "source_pair_reuses_prior32":p["reuses_any_prior_32WISE6_quartet_station"],
                       "zero_pair_reuses_prior32":q["reuses_any_prior_32WISE6_quartet_station"]})
        used.extend(ids)
    return chosen

def analyze_frame(frame,prior,root):
    assert frame["schema"]=="gesui_DE_national_2022_directedWISE6_73_sourcezero_pair_pre_plantae_pre_butterfly_frame_v01"
    assert frame["new_public_plant_photo_API_requests_performed"] is False
    assert frame["new_butterfly_observation_API_requests_performed"] is False
    assert len(frame["complete_original_fixed_candidate_pairs"])==73
    assert len(frame["original_unique_public_WISE6_monitor_site_metadata"])==100
    assert prior["successful_queries"]==prior["planned_queries"]==64
    assert prior["no_butterfly_species_or_photos_requested"]
    if prior["original_frozen_sourceSHA256"]["original_quartets_sha256"]!=frame["original_frozen_source_hashes_SHA256"]["prior_8quartets"]:
        raise ValueError("Old original source quartet hash changed")
    if prior["original_frozen_sourceSHA256"]["national_sites_sha256"]!=frame["original_frozen_source_hashes_SHA256"]["public_sites"]:
        raise ValueError("Old original national WISE6 station registry hash changed")
    stations={p["original_site_id_SHA256"]:p for p in frame["original_unique_public_WISE6_monitor_site_metadata"]}
    original={r["original_monitor_hash"]:r for r in prior["original_station_source_without_butterfly_outcomes"]}
    assert len(original)==32
    assert set(original).issubset(stations)
    new=[k for k in sorted(stations) if k not in original]
    assert len(new)==68
    obs={}
    for i,k in enumerate(sorted(stations),1):
        s=stations[k]
        obs[k]={}
        for r in RADII:
            if k in original:
                p=old_count(original[k],r,s)
            else:
                p=query(s,r)
                time.sleep(1.15)
            obs[k][f"{r}km"]=p
            print(f"[{i}/100, {r}km] origin={p['origin']} count={p['total_results']}",flush=True)
    out={
        "schema":"gesui_DE_national_2022_WISE6_all73_pair_2019_2021_plantae_prebutterfly_photo_effort_v01",
        "utc":datetime.now(timezone.utc).isoformat(),
        "original_73pair_100site_source_freeze_SHA256":sha(root/ROOT_FRAME),
        "original_32site_prior_plant_source_SHA256":sha(root/PRIOR_PLANT),
        "no_2022_butterfly_photo_counts_queried":True,
        "no_actual_river_TP_numeric_or_butterfly_substrate_exposure_measured":True,
        "original_sourcezero_not_absence_of_all_wastewater":True,
        "original_WISE6_station_count":len(stations),
        "original_reused_prior_sites_n":len(original),
        "new_station_2019_2021_plant_photo_queries_planned":len(new)*len(RADII),
        "independent_plant_photo_2022_butterfly_opportunity_years":[2019,2020,2021],
        "radii_km":list(RADII),
        "full_original_source_only_station_prior_year_plant_photo_receipts":{
            k:{"countryCode":"DE","lat":stations[k]["lat"],"lon":stations[k]["lon"],
               "original_WISE6_2022_monitored_dates":stations[k]["original_2022_distinct_TP_sample_dates_metadata"],
               "query_results":obs[k]} for k in sorted(stations)
        },
        "by_radius":{}
    }
    allpairs=frame["complete_original_fixed_candidate_pairs"]
    for rad in RADII:
        candidate=[]
        perpair={}
        for p in allpairs:
            evaluation=qualifies_pair(p,stations,obs,rad)
            perpair[p["pair_id_SHA256"]]=evaluation
            if evaluation["eligible"]:candidate.append(p)
        selected=choose_national(candidate,stations,rad)
        matched=choose_quartets(candidate,stations,rad)
        typ=Counter(p["exposure_class"] for p in selected)
        bothregion=defaultdict(set)
        classblocks=defaultdict(set)
        for p in selected:
            bothregion[p["upstream_3deg_region"]].add(p["exposure_class"])
            classblocks[p["exposure_class"]].add(p["same_2deg_upstream_region"])
        pcount=Counter(p["exposure_class"] for p in candidate)
        fourarm_GO=bool(len(matched)>=8 and len({m["original_3deg_block"] for m in matched})>=3)
        national_GO=bool(
            typ[SRC]>=10 and typ[ZERO]>=10
            and len(classblocks[SRC])>=4 and len(classblocks[ZERO])>=4
            and sum(len(v)==2 for v in bothregion.values())>=3
        )
        report={
            "all_original_candidate_pairs_n":73,
            "all_original_candidate_source_class_counts":dict(Counter(p["exposure_class"] for p in allpairs)),
            "plant_opportunity_valid_geometrically_nonoverlapping_pair_candidates_n":len(candidate),
            "photo_opportunity_pair_counts_by_class":dict(pcount),
            "new_pure_original_68site_only_pair_candidates":sum(not p["reuses_any_prior_32WISE6_quartet_station"] for p in candidate),
            "source_only_photo_opportunity_criteria":{"each_site_min_plantae_records":5,"up_down_min_max_ratio":.15,"distance_center_gt_2times_radius":True},
            "fully_site_disjoint_country_wide_directional_pairs":[
                {"pair_id_SHA256":p["pair_id_SHA256"],"class":p["exposure_class"],
                 "original_2deg_group":p["same_2deg_upstream_region"],
                 "original_3deg_group":p["upstream_3deg_region"],
                 "reuses_old32_site":p["reuses_any_prior_32WISE6_quartet_station"]} for p in selected],
            "national_pair_design_source_counts":dict(typ),
            "national_pair_design_geographic_source_two_degree_groups":{k:len(v) for k,v in classblocks.items()},
            "national_pair_design_three_degree_regions_with_both_statuses":sum(len(v)==2 for v in bothregion.values()),
            "national_pair_design_support_GATE":{
                "min_source_pairs":10,"min_registeredzero_pairs":10,
                "min_original_2deg_blocks_per_class":4,
                "min_3deg_regions_with_both_classes":3,
                "national_2022_butterfly_photograph_outcomes_permitted":national_GO
            },
            "original_3deg_source_vs_registeredzero_matched_fourarm_quartets":matched,
            "original_3deg_matched_quartets_n":len(matched),
            "original_3deg_matched_geographic_blocks":len({m["original_3deg_block"] for m in matched}),
            "three_degree_four_arm_design_support_GATE":{
                "min_four_station_disjoint_sourcezero_quartets":8,
                "min_original_3deg_geo_blocks":3,
                "matched_2022_butterfly_photo_outcomes_permitted":fourarm_GO
            },
            "all73_pair_plant_opportunity_and_geographic_nonoverlap_gates":perpair,
            "inference_note":"External iNat Plantae reports indicate public photographic opportunity, not unbiased plant biomass or butterfly counts. No butterfly/photo outcome accessed to select source geography or alter predeclared thresholds."
        }
        out["by_radius"][f"{rad}km"]=report
    out["actual_new_public_plant_photo_requests_success"]=sum(
        1 for k in new for r in RADII if obs[k][f"{r}km"]["status"]==200)
    out["actual_new_public_plant_photo_requests_failed"]=len(new)*len(RADII)-out["actual_new_public_plant_photo_requests_success"]
    out["failed_requests_are_missing_not_zero"]=True
    return out

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    frame=json.loads((root/ROOT_FRAME).read_text())
    prior=json.loads((root/PRIOR_PLANT).read_text())
    out=analyze_frame(frame,prior,root)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "new_public_2019_2021_plantae_source_success":out["actual_new_public_plant_photo_requests_success"],
        "new_public_request_failures":out["actual_new_public_plant_photo_requests_failed"],
        "radii_decisions":{k:{
            "photo_opportunity_pair_counts_by_class":v["photo_opportunity_pair_counts_by_class"],
            "matched_quartets":v["original_3deg_matched_quartets_n"],
            "matched_3deg_blocks":v["original_3deg_matched_geographic_blocks"],
            "national_pairs":v["national_pair_design_source_counts"],
            "fourarm_butterfly_gate":v["three_degree_four_arm_design_support_GATE"],
            "national_butterfly_gate":v["national_pair_design_support_GATE"]
        } for k,v in out["by_radius"].items()}
    },indent=2))

if __name__=="__main__":
    main()
