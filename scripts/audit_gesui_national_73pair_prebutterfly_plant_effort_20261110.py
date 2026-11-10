#!/usr/bin/env python3
"""Independent offline source/effort/stop-gate QA, not an ecological-effect model."""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path

ROOT="results/lycaenidae_2022_DE_national_73pair_sourceonly/"
ORIGINAL=ROOT+"ORIGINAL_2022_73DE_WISE6_SOURCEZERO_PAIR_100SITE_PREPHOTO_FRAME.json"
EFFORT=ROOT+"ORIGINAL_100STATION_2KM_5KM_2019_2021_PLANTAE_EFFORT_BEFORE_2022_BUTTERFLIES.json"
SRC="registered_source_interposed"
ZERO="registered_qualifying_sourcezero"

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def geo(a,b):
    r=math.pi/180
    dlat=(a["lat"]-b["lat"])*r
    dlon=(a["lon"]-b["lon"])*r
    z=math.sin(dlat/2)**2+math.cos(a["lat"]*r)*math.cos(b["lat"]*r)*math.sin(dlon/2)**2
    return 12742*math.asin(min(1,math.sqrt(z)))

def mids(v):
    v=sorted(v)
    n=len(v)
    return (v[(n-1)//2]+v[n//2])/2

def audit(root):
    fp=root/ORIGINAL
    ep=root/EFFORT
    f=json.loads(fp.read_text())
    e=json.loads(ep.read_text())
    assert digest(fp)==e["original_73pair_100site_source_freeze_SHA256"]
    assert f["new_butterfly_observation_API_requests_performed"] is False
    assert f["new_public_plant_photo_API_requests_performed"] is False
    assert e["no_2022_butterfly_photo_counts_queried"] is True
    assert f["original_unique_monitoring_station_locations"]==100
    assert len(f["complete_original_fixed_candidate_pairs"])==73
    assert f["country_by_class_counts"]=={ZERO:40,SRC:33}
    assert f["spatial_positive_control_overlap_diagnostic"]["2deg"]["mathematical_sourcezero_matched_quartet_upper_bound_before_site_disjointness_and_plant_effort"]==7
    stations={x["original_site_id_SHA256"]:x for x in f["original_unique_public_WISE6_monitor_site_metadata"]}
    assert len(stations)==100
    pairs={p["pair_id_SHA256"]:p for p in f["complete_original_fixed_candidate_pairs"]}
    obs=e["full_original_source_only_station_prior_year_plant_photo_receipts"]
    assert set(obs)==set(stations)
    total=Counter()
    results={}
    for rad in (2,5):
        radius=f"{rad}km"
        expected={}
        count=[]
        for site_id,row in obs.items():
            assert row["lat"]==stations[site_id]["lat"]
            assert row["lon"]==stations[site_id]["lon"]
            query=row["query_results"][radius]
            assert query["status"]==200 and type(query["total_results"])==int
            url=query["original_public_API_url"]
            for part in ("taxon_id=47126","place_id=7207",f"radius={rad}","d1=2019-01-01",
                         "d2=2021-12-31","quality_grade=research","acc_below=500"):
                assert part in url
            assert len(query["original_JSON_SHA256"])==64
            count.append(query["total_results"])
            total[query["origin"]]+=1
        for pair_id,p in pairs.items():
            u,d=p["upstream_original_site_id_SHA256"],p["downstream_original_site_id_SHA256"]
            q=e["by_radius"][radius]["all73_pair_plant_opportunity_and_geographic_nonoverlap_gates"][pair_id]
            a=obs[u]["query_results"][radius]["total_results"]
            b=obs[d]["query_results"][radius]["total_results"]
            dst=geo(stations[u],stations[d])
            ratio=min(a,b)/max(a,b) if max(a,b)>0 else 0
            passgate=bool(a>=5 and b>=5 and ratio>=.15 and dst>2*rad)
            assert bool(q["eligible"])==passgate
            assert q["upstream_plant_photos"]==a and q["downstream_plant_photos"]==b
            expected[pair_id]=passgate
        qualified={p for p,v in expected.items() if v}
        rr=e["by_radius"][radius]
        assert rr["plant_opportunity_valid_geometrically_nonoverlapping_pair_candidates_n"]==len(qualified)
        assert Counter(pairs[k]["exposure_class"] for k in qualified)==Counter(rr["photo_opportunity_pair_counts_by_class"])
        quartets=rr["original_3deg_source_vs_registeredzero_matched_fourarm_quartets"]
        used=[]
        for quartet in quartets:
            sp=quartet["source_pair_SHA256"]
            zp=quartet["registeredzero_pair_SHA256"]
            assert sp in qualified and zp in qualified
            assert pairs[sp]["exposure_class"]==SRC and pairs[zp]["exposure_class"]==ZERO
            assert pairs[sp]["upstream_3deg_region"]==pairs[zp]["upstream_3deg_region"]==quartet["original_3deg_block"]
            ids=[pairs[k][side+"_original_site_id_SHA256"] for k in (sp,zp) for side in ("upstream","downstream")]
            assert sorted(ids)==sorted(quartet["four_station_hashes"]) and len(set(ids))==4
            assert all(geo(stations[x],stations[y])>2*rad for i,x in enumerate(ids) for y in ids[i+1:])
            assert all(geo(stations[x],stations[y])>2*rad for x in ids for y in used)
            used.extend(ids)
        selected=rr["fully_site_disjoint_country_wide_directional_pairs"]
        assert len(selected)==sum(rr["national_pair_design_source_counts"].values())
        national_used=[]
        byclass=Counter()
        shared=defaultdict(set)
        for r in selected:
            p=pairs[r["pair_id_SHA256"]]
            assert p["pair_id_SHA256"] in qualified
            assert r["class"]==p["exposure_class"]
            ids=[p["upstream_original_site_id_SHA256"],p["downstream_original_site_id_SHA256"]]
            assert all(geo(stations[x],stations[y])>2*rad for x in ids for y in national_used)
            national_used.extend(ids)
            byclass[p["exposure_class"]]+=1
            shared[p["upstream_3deg_region"]].add(p["exposure_class"])
        assert dict(byclass)==rr["national_pair_design_source_counts"]
        assert sum(len(v)==2 for v in shared.values())==rr["national_pair_design_three_degree_regions_with_both_statuses"]
        assert rr["national_pair_design_support_GATE"]["national_2022_butterfly_photograph_outcomes_permitted"] is False
        assert rr["three_degree_four_arm_design_support_GATE"]["matched_2022_butterfly_photo_outcomes_permitted"] is False
        results[radius]={
            "original_source_pairs_n":73,
            "photo_opportunity_source_zero_counts":dict(Counter(pairs[k]["exposure_class"] for k in qualified)),
            "photo_opportunity_eligible_pairs":len(qualified),
            "original_WISE6_independent_pair_geo_quartets":len(quartets),
            "original_3deg_quartet_blocks":len({q["original_3deg_block"] for q in quartets}),
            "original_countrywide_eligible_site_disjoint_pairs_by_class":dict(byclass),
            "nationwide_both_class_3deg_regions":sum(len(v)==2 for v in shared.values()),
            "plant_photo_station_min":min(count),"median":mids(count),"max":max(count),
            "station_plant_photos_zero_n":sum(v==0 for v in count),
            "station_plant_photos_below5_n":sum(v<5 for v in count),
            "butterfly_outcome_pre_source_gate":"STOP"
        }
    assert total["NEW_2026_11_10"]==136 and total["REUSED_OLD32_2026_11_09"]==64
    assert e["actual_new_public_plant_photo_requests_success"]==136
    assert e["actual_new_public_plant_photo_requests_failed"]==0
    return {"schema":"gesui_DE_national_73pair_2019_2021_plantae_prior_outcome_independent_offline_replay_v01",
            "source_frame_SHA256":digest(fp),"plant_200total_receipts_SHA256":digest(ep),
            "original_total_station_radius_requests":200,
            "new_2026_11_10_plant_requests_validated":136,
            "reused_2026_11_09_prior32_station_radius_plant_requests_validated":64,
            "no_butterfly_photo_numbers_opened_or_modeled":True,
            "all_predeclared_2km_and_5km_source_only_gates_fail":True,
            "no_test_of_lycaenidae_biological_negative_or_positive_effect":True,
            "by_radius":results}

if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args()
    root=Path(__file__).resolve().parents[1]
    result=audit(root)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print(json.dumps(result,indent=2))
