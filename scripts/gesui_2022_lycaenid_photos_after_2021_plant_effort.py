#!/usr/bin/env python3
"""2022 butterfly PHOTO REPORTING share in independently effort-selected 2022 river pairs.

This code must run ONLY after GitHub-committed 2021 Plantae effort screen,
source-only EEA WWTP pairs and explicit gates have been verified. These
are photo counts, not individual butterflies, mud-puddling, or sewage causality.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, time, urllib.request, urllib.parse, urllib.error
from collections import defaultdict,Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np

API="https://api.inaturalist.org/v1/observations"
TAXA={"Lycaenidae":47923,"Pieridae":48508,"Papilionidae":47223}
STATUS=("upstream_registered_WWTP","no_registered_qualifying_upstream_WWTP")
MIN_BUTTERFLY_TARGET_PHOTOS_EACH=5
MIN_PAIRS=12
MIN_COUNTRIES=3
RADIUS=5.0
MAX_GPS_ACCURACY_M=500
YEAR=2022
SEED=20261108
N_BOOT=1999

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def api_count(row,family,taxid):
    params={
      "lat":str(float(row["lat"])),"lng":str(float(row["lon"])),"radius":RADIUS,
      "taxon_id":taxid,"year":YEAR,"photos":"true","geo":"true",
      "geoprivacy":"open","taxon_geoprivacy":"open",
      "quality_grade":"research","captive":"false","acc_below":MAX_GPS_ACCURACY_M,
      "per_page":1}
    url=API+"?"+urllib.parse.urlencode(params)
    out={"url":url,"taxon_family":family,"http":None,"count":None}
    for attempt in range(4):
        try:
            req=urllib.request.Request(url,headers={
                "User-Agent":"gesui-2022-effort-locked-lycaenid-reporting-share/0.1",
                "Accept":"application/json"})
            with urllib.request.urlopen(req,timeout=35) as r:
                raw=r.read(500000)
                j=json.loads(raw)
                if "total_results" not in j:
                    raise ValueError("Missing source total_results")
                out.update({"http":r.status,"count":int(j["total_results"]),
                    "response_sha256":hashlib.sha256(raw).hexdigest()})
                return out
        except urllib.error.HTTPError as e:
            out.update({"http":e.code,"error":str(e)})
            if e.code not in (429,500,502,503) or attempt==3:break
        except Exception as e:
            out.update({"error":type(e).__name__+": "+str(e)})
            if attempt==3:break
        time.sleep(2.1*(attempt+1))
    return out

def main():
    ap=argparse.ArgumentParser()
    for p in ("source-only-qc","effort-qc","selected-pairs","out"):
        ap.add_argument("--"+p,type=Path,required=True)
    a=ap.parse_args()
    prior=json.loads(a.source_only_qc.read_text())
    effort=json.loads(a.effort_qc.read_text())
    if not prior.get("no_butterfly_or_other_ecological_outcomes_read"):
        raise ValueError("Non-source-only original river pair input")
    if not effort.get("NO_2022_BUTTERFLY_RESULT_QUERIED"):
        raise ValueError("Effort sample contaminated by same-year butterfly outcomes")
    if effort.get("source_registry_SHA256")!=prior["paired_sample_sha256"]:
        raise ValueError("Prior-year botanical screen not linked to original source-only registry")
    if effort.get("selected_pair_CSV_sha256")!=sha(a.selected_pairs):
        raise ValueError("2021 botanical effort-selected 2022 butterfly site set changed after freeze")
    source_pairs=list(csv.DictReader(a.selected_pairs.open()))
    if len(source_pairs)!=2*effort["selected_n"]:
        raise ValueError("2021 effort-selected pair registry size discrepancy")
    groups=defaultdict(dict)
    for r in source_pairs:groups[int(r["pair_id"])][r["status"]]=r
    if any(set(g)!=set(STATUS) for g in groups.values()):
        raise ValueError("Missing original source/control paired identities")
    pre_gate=effort["passes_independent_effort_gate_12_pairs_3_countries"]
    result={
        "schema":"gesui_2022_5km_butterfly_photo_reporting_share_after_2021_plant_effort_v01",
        "read_only_source_year":YEAR,
        "source_pair_original_SHA256":prior["paired_sample_sha256"],
        "prior_year_plant_effort_QC_SHA256":sha(a.effort_qc),
        "effort_selected_site_registry_SHA256":sha(a.selected_pairs),
        "prephoto_effort_gate":pre_gate,
        "both_photo_and_source_2022_year_matched":True,
        "photo_geographic_radius_km":RADIUS,
        "photo_max_reported_position_accuracy_m":MAX_GPS_ACCURACY_M,
        "taxa":TAXA,"no_actual_puddling_labels_or_images_read":True,
        "no_causal_WWTP_or_true_abundance_claim":True,
        "non_qualifying_source_zero_is_not_no_sewage":True,
        "date_utc":datetime.now(timezone.utc).isoformat()
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    if not pre_gate:
        result.update({"stage":"STOP_BEFORE_ANY_2022_BUTTERFLY_OUTCOME_QUERY",
                        "n_original_source_pairs":prior["country_3deg_matched_independent_pairs"],
                        "n_2021_plant_effort_selected_pairs":len(groups),
                        "needed_pairs":MIN_PAIRS,"needed_countries":MIN_COUNTRIES,
                        "per_site_counts":[],"eligible_pair_results":[],
                        "conclusion":"Independent previous-year nonbutterfly observation effort gate failed; no 2022 butterfly result."})
        a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
        print("GESUI_LYCAENIDAE_EFFORT_SCREEN_STOP_NO_BUTTERFLY_QUERIES",flush=True)
        return
    points=[]
    for i,r in enumerate(source_pairs):
        counts={name:api_count(r,name,taxid) for name,taxid in TAXA.items()}
        all_ok=all(q["count"] is not None for q in counts.values())
        n=sum(q["count"] for q in counts.values()) if all_ok else None
        points.append({"pair_id":int(r["pair_id"]),"status":r["status"],
          "country":r["countryCode"],"n_target_family_photo_observations":n,
          "family_photo_counts":counts,"all_API_count_queries_complete":all_ok})
        print("GESUI_LYCAENIDAE_PRIORYEAR_PLANT_EFFORT_LOCKED_2022_PHOTOS",
            json.dumps({"index":i+1,"of":len(source_pairs),"pair":r["pair_id"],
                         "status":r["status"],
                         "counts":{k:v["count"] for k,v in counts.items()}},sort_keys=True),flush=True)
        time.sleep(1.0)
    paired=defaultdict(dict)
    for x in points:paired[x["pair_id"]][x["status"]]=x
    valid=[]
    for pid,group in sorted(paired.items()):
        e=group[STATUS[0]];c=group[STATUS[1]]
        if not e["all_API_count_queries_complete"] or not c["all_API_count_queries_complete"]:
            continue
        if min(e["n_target_family_photo_observations"],c["n_target_family_photo_observations"])<MIN_BUTTERFLY_TARGET_PHOTOS_EACH:
            continue
        se=e["family_photo_counts"]["Lycaenidae"]["count"]/e["n_target_family_photo_observations"]
        sc=c["family_photo_counts"]["Lycaenidae"]["count"]/c["n_target_family_photo_observations"]
        valid.append({"pair_id":pid,"country":e["country"],
                      "source_registered_share":se,"no_registered_source_share":sc,
                      "difference":se-sc,
                      "source_registered_n":e["n_target_family_photo_observations"],
                      "no_registered_source_n":c["n_target_family_photo_observations"]})
    countries=sorted(set(x["country"] for x in valid))
    ready=len(valid)>=MIN_PAIRS and len(countries)>=MIN_COUNTRIES
    result.update({
       "n_source_pairs_2021_effort_selected":len(groups),
       "n_source_and_control_sites_queried":len(points),
       "n_2022_butterfly_API_count_queries":len(points)*len(TAXA),
       "all_complete_pair_n":sum(all(y["all_API_count_queries_complete"] for y in g.values()) for g in paired.values()),
       "eligible_targetgroup_photos_each_side_threshold":MIN_BUTTERFLY_TARGET_PHOTOS_EACH,
       "n_2022_butterfly_photo_support_eligible_pairs":len(valid),
       "n_independent_countries_eligible":len(countries),
       "eligible_country_counts":dict(Counter(x["country"] for x in valid)),
       "paired_reporting_share_gate_pass":ready,
       "per_site_counts":points,
       "eligible_pair_results":valid if ready else [],
       "data_support_without_effect":valid if not ready else None,
       "photos_are_observation_events_not_individuals":True
    })
    if ready:
        diffs=np.asarray([x["difference"] for x in valid],float)
        bycountry=defaultdict(list)
        for x in valid:bycountry[x["country"]].append(x["difference"])
        cmeans={c:float(np.mean(v)) for c,v in sorted(bycountry.items())}
        rng=np.random.default_rng(SEED)
        clusters=np.array(list(bycountry))
        boots=[]
        for _ in range(N_BOOT):
            sampled=rng.choice(clusters,size=len(clusters),replace=True)
            v=[p for country in sampled for p in bycountry[country]]
            boots.append(float(np.mean(v)))
        result["descriptive_reporting_share_only"]={
           "mean_station_pair_difference":float(np.mean(diffs)),
           "median_station_pair_difference":float(np.median(diffs)),
           "positive_pairs":int((diffs>0).sum()),"negative_pairs":int((diffs<0).sum()),
           "zero_pairs":int((diffs==0).sum()),
           "mean_equal_country_weight":float(np.mean(list(cmeans.values()))),
           "per_country_pair_difference_means":cmeans,
           "descriptive_country_resample_95pct_interval":[float(q) for q in np.quantile(boots,[0.025,0.975])],
           "bootstrap_is_posthoc_descriptive_not_caused_by_wastewater":True,
        }
    else:
        result["stage"]="STOP_BEFORE_INFERENTIAL_PHOTO_SHARE_MODEL"
        result["conclusion"]="2022 butterfly photo opportunity too sparse or concentrated; no wastewater association estimate is supported."
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_2022_LYCAENIDAE_EFFORT_RICH_PHOTO_COHORT_SOURCE_QC_COMPLETED",
       json.dumps({k:v for k,v in result.items() if k in (
           "n_source_pairs_2021_effort_selected","n_2022_butterfly_photo_support_eligible_pairs",
           "eligible_country_counts","paired_reporting_share_gate_pass",
           "descriptive_reporting_share_only")},sort_keys=True),flush=True)
if __name__=="__main__":main()
