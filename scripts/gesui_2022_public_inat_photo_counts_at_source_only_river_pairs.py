#!/usr/bin/env python3
"""Source-frozen 2022 EEA WWTP river sites vs iNaturalist butterfly photo opportunity.

The site/stratum registry was fixed and committed from EEA ONLY before
this script is run. Taxa are official iNaturalist IDs with descendants.
Counts are PUBLIC PHOTO OBSERVATION OPPORTUNITY, not butterfly density,
abundance, wet-soil use, puddling, or any wastewater mechanism.

No EQR data is opened. No photo contents are downloaded.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,time,urllib.parse,urllib.request,urllib.error
from collections import defaultdict,Counter
from datetime import datetime,timezone
from pathlib import Path

API="https://api.inaturalist.org/v1/observations"
TAXA={"Lycaenidae":47923,"Pieridae":48508,"Papilionidae":47223}
SOURCE_STATUS=("upstream_registered_WWTP","no_registered_qualifying_upstream_WWTP")
YEAR="2022"
RADIUS_KM=1.0
ACCURACY_CAP_M=500
MAX_PAIRS=24
MIN_PAIR_TARGET_GROUP_PHOTOS_EACH_SIDE=5
MIN_ELIGIBLE_PAIRS_FOR_DESCRIPTIVE_SHARE=12
PAUSE=1.12

def request(lat,lon,taxid):
    p={"lat":f"{float(lat):.7f}","lng":f"{float(lon):.7f}",
       "radius":str(RADIUS_KM),"taxon_id":str(taxid),
       "d1":"2022-01-01","d2":"2022-12-31",
       "photos":"true","geo":"true","geoprivacy":"open",
       "taxon_geoprivacy":"open",
       "quality_grade":"research","captive":"false",
       "acc_below":str(ACCURACY_CAP_M),"per_page":"1"}
    url=API+"?"+urllib.parse.urlencode(p)
    result={"public_api_url":url,"query":p,"http_status":None,
            "total_photo_observations":None,"error":None}
    for attempt in range(3):
        try:
            req=urllib.request.Request(url,headers={
                "User-Agent":"gesui-2022-sourcefrozen-butterfly-photo-count/0.1",
                "Accept":"application/json"})
            with urllib.request.urlopen(req,timeout=35) as r:
                raw=r.read(500000)
                data=json.loads(raw)
                if "total_results" not in data:
                    raise ValueError("iNat API did not expose observation count")
                result["http_status"]=r.status
                result["source_response_sha256"]=hashlib.sha256(raw).hexdigest()
                result["total_photo_observations"]=int(data["total_results"])
                # No observation IDs, positions, usernames or photo pixels stored.
                return result
        except urllib.error.HTTPError as e:
            result["error"]=f"HTTP {e.code}"
            result["http_status"]=e.code
            if e.code not in (429,500,502,503) or attempt==2:break
            time.sleep(3*(attempt+1))
        except Exception as e:
            result["error"]=type(e).__name__+": "+str(e)
            if attempt==2:break
            time.sleep(2)
    return result

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--paired-sites",required=True,type=Path)
    parser.add_argument("--source-qc",required=True,type=Path)
    parser.add_argument("--out",required=True,type=Path)
    args=parser.parse_args()
    prior=json.loads(args.source_qc.read_text())
    if prior["no_butterfly_or_other_ecological_outcomes_read"] is not True:
        raise ValueError("Potential data leakage from non-source-only stations")
    expected=prior["paired_sample_sha256"]
    if hashlib.sha256(args.paired_sites.read_bytes()).hexdigest()!=expected:
        raise ValueError("Prephoto site pairing changed")
    pairs=list(csv.DictReader(args.paired_sites.open()))
    if len(pairs)%2 or len(pairs)//2!=prior["country_3deg_matched_independent_pairs"] or len(pairs)>MAX_PAIRS*2:
        raise ValueError("Sample/paired source-only registry mutated")
    groups=defaultdict(list)
    for row in pairs:
        groups[row["pair_id"]].append(row)
    if not all(len(v)==2 and set(x["status"] for x in v)==set(SOURCE_STATUS) for v in groups.values()):
        raise ValueError("Unexpected exposure/control pairing")
    result=[]
    for row in pairs:
        counts={}
        for name, taxid in TAXA.items():
            counts[name]=request(row["lat"],row["lon"],taxid)
            print("GESUI_PUBLIC_2022_BUTTERFLY_SITE_PHOTO_COUNT",json.dumps({
                "pair":row["pair_id"],"source_status":row["status"],"taxon":name,
                "count":counts[name]["total_photo_observations"],
                "http":counts[name]["http_status"]},sort_keys=True),flush=True)
            time.sleep(PAUSE)
        v={k:row[k] for k in ("pair_id","site_key","status","countryCode","region_3deg_block",
                              "upstream_wwtp_count","matching_pair_distance_km")}
        v["taxa"]=counts
        available=[c["total_photo_observations"] for c in counts.values()]
        v["all_three_queries_success"]=all(n is not None for n in available)
        v["three_families_photorecord_total"]=(sum(available) if v["all_three_queries_success"] else None)
        result.append(v)
    by_pair=defaultdict(dict)
    for x in result:by_pair[x["pair_id"]][x["status"]]=x
    valid_pairs=[]
    all_ok=0
    for id,g in sorted(by_pair.items(),key=lambda q:int(q[0])):
        up=g[SOURCE_STATUS[0]]
        down=g[SOURCE_STATUS[1]]
        if up["all_three_queries_success"] and down["all_three_queries_success"]:
            all_ok+=1
            if min(up["three_families_photorecord_total"],down["three_families_photorecord_total"])>=MIN_PAIR_TARGET_GROUP_PHOTOS_EACH_SIDE:
                valid_pairs.append(id)
    pair_table=[]
    for id in valid_pairs:
        up=by_pair[id][SOURCE_STATUS[0]];down=by_pair[id][SOURCE_STATUS[1]]
        def share(t):
            x=t["taxa"]["Lycaenidae"]["total_photo_observations"]
            n=t["three_families_photorecord_total"]
            return x/n if n else None
        pair_table.append({"pair_id":int(id),"country":up["countryCode"],
           "registered_wwtp_lycaenidae_share":share(up),
           "unregistered_source_lycaenidae_share":share(down),
           "difference_upstream_minus_control":share(up)-share(down),
           "registered_wwtp_L_photo_count":up["taxa"]["Lycaenidae"]["total_photo_observations"],
           "registered_wwtp_target_group_photo_count":up["three_families_photorecord_total"],
           "no_registered_source_L_photo_count":down["taxa"]["Lycaenidae"]["total_photo_observations"],
           "no_registered_source_target_group_photo_count":down["three_families_photorecord_total"]})
    # A predeclared descriptive data-support gate, not post-hoc sign search.
    gate=len(valid_pairs)>=MIN_ELIGIBLE_PAIRS_FOR_DESCRIPTIVE_SHARE
    summary={
        "source_only_matched_2022_registered_WWTP_vs_zero_count":len(groups),
        "successful_public_API_complete_pairs":all_ok,
        "min_photos_per_side":MIN_PAIR_TARGET_GROUP_PHOTOS_EACH_SIDE,
        "eligible_photo_opportunity_pairs":len(valid_pairs),
        "eligible_pair_by_country":dict(Counter(row["country"] for row in pair_table)),
        "min_eligible_pairs":MIN_ELIGIBLE_PAIRS_FOR_DESCRIPTIVE_SHARE,
        "photo_opportunity_share_gate":gate,
        "not_behavior_or_puddling":True,
        "not_actual_butterfly_abundance_or_colonization":True,
        "not_wastewater_causality":True,
    }
    if gate:
        deltas=[z["difference_upstream_minus_control"] for z in pair_table]
        summary["source_vs_control_photo_share_descriptive_only"]={
            "mean_paired_difference":sum(deltas)/len(deltas),
            "median_pair_difference":sorted(deltas)[len(deltas)//2],
            "positive_difference_pairs":sum(v>0 for v in deltas),
            "negative_difference_pairs":sum(v<0 for v in deltas),
            "zero_difference_pairs":sum(v==0 for v in deltas),
            "not_hypothesis_test_and_no_puddling_claim":True,
        }
    output={
      "schema":"gesui_2022_prephoto_eea_matched_river_iNat_butterfly_photo_opportunity_v01",
      "source_year":YEAR,"iNat_photo_observation_year":YEAR,
      "site_matching_source_only_original_EEA_csv_sha256":expected,
      "site_match_qc_sha256":hashlib.sha256(args.source_qc.read_bytes()).hexdigest(),
      "search_radius_km":RADIUS_KM,"coordinate_accuracy_max_m":ACCURACY_CAP_M,
      "families":TAXA,"source":API,"retrieved_utc":datetime.now(timezone.utc).isoformat(),
      "summary":summary,"per_matched_site":result,"per_eligible_pair":pair_table,
      "no_images_downloaded":True,"no_image_behavior_classification":True,
      "no_observer_identification_archived":True,
      "caveats":[
        "Registered WWTP present/absent in source-limited 50km upstream EEA river graph, not actual effluent chemistry",
        "Site placement comes from public aquatic monitoring stations and could select peculiar habitats",
        "Controls are registered-source-zero not wastewater-free",
        "Within-pair urban land cover, host plants, season and observer effort are not matched",
        "Only research-grade georeferenced photographed iNaturalist observations counted, not butterflies",
        "No photos visually reviewed for puddling, inferred behavior forbidden",
        "Photo count near river monitor does not establish that a butterfly reached the connected river reach",
        "Same observer may supply several observations; photo count not independent individual butterflies",
        "Country coverage may be uneven and photo source sparsity may fail the gate"
      ],
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(output,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_LYCAENIDAE_2022_PUBLIC_PHOTO_OPPORTUNITY_PAIRS_COMPLETED",
      json.dumps(summary,ensure_ascii=False,sort_keys=True),flush=True)
if __name__=="__main__":main()
