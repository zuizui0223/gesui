#!/usr/bin/env python3
"""Offline independent replay of pre-butterfly source/effort gates and asymmetry."""
from __future__ import annotations
import json
import math
import hashlib
import argparse
from collections import defaultdict, Counter
from itertools import product
from pathlib import Path


def h(lat1, lon1, lat2, lon2):
    r = math.pi/180
    a = math.sin((lat2-lat1)*r/2)**2 + math.cos(lat1*r)*math.cos(lat2*r)*math.sin((lon2-lon1)*r/2)**2
    return 12742 * math.asin(math.sqrt(min(1,a)))


def exact_signflip(vals):
    mean = sum(vals)/len(vals)
    possible = [sum(s*x for s,x in zip(signs,vals))/len(vals)
                for signs in product((-1,1),repeat=len(vals))]
    return {
        "n_original_quartets":len(vals), "observed_equal_quartet_mean":mean,
        "n_positive":sum(x>0 for x in vals),
        "exact_two_sided_signflip_reference_fraction":
            sum(abs(x)>=abs(mean)-1e-12 for x in possible)/len(possible),
        "not_a_randomized_source_effect_p_value":True,
    }


def audit(path):
    o=json.loads(path.read_text())
    rows=o["original_station_source_without_butterfly_outcomes"]
    assert o["successful_queries"]==o["planned_queries"]==64
    assert o["no_butterfly_species_or_photos_requested"]
    assert len(rows)==32 and len({r["original_monitor_hash"] for r in rows})==32
    assert Counter(r["original_quartet_number"] for r in rows)==Counter({i:4 for i in range(1,9)})
    out={"replay_matches_original_prebutterfly_gate":True,
         "original_plant_gate_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
         "record_total":len(rows),"original_quartets":8,"year_2019_2021_plantae_not_2022_butterflies":True,"by_radius":{}}
    for rad in (2,5):
        counts=[]
        candidates=[]
        contrasts=[]
        sources=[]
        for j in range(1,9):
            members=[r for r in rows if r["original_quartet_number"]==j]
            vals={r["original_2022_source_class"]+"_"+r["side_of_river_monitor_pair"]:
                  r["independent_previous_year_plantae_photo_counts"][f"{rad}km"]["count_public_2019_2021_research_grade_georeferenced_plantae_photos"]
                  for r in members}
            assert len(vals)==4 and all(isinstance(v,int) and v>=0 for v in vals.values())
            assert set(vals)=={
                "source_upstream","source_downstream",
                "registered_source_zero_upstream","registered_source_zero_downstream"
            }
            for r in members:
                v=r["independent_previous_year_plantae_photo_counts"][f"{rad}km"]
                assert v["status"]==200 and v["original_API_url"].startswith("https://api.inaturalist.org/v1/observations?")
                assert "&radius="+str(rad) in v["original_API_url"]
                assert "taxon_id=47126" in v["original_API_url"]
                assert "d1=2019-01-01" in v["original_API_url"]
                assert "d2=2021-12-31" in v["original_API_url"]
                counts.append(v["count_public_2019_2021_research_grade_georeferenced_plantae_photos"])
            distances=[
                h(x["latitude"],x["longitude"],y["latitude"],y["longitude"])
                for n,x in enumerate(members) for y in members[n+1:]
            ]
            ratio=min(vals.values())/max(vals.values()) if max(vals.values())>0 else 0.
            gate=bool(min(vals.values())>=5 and ratio>=.15 and min(distances)>2*rad)
            orig=o["prebutterfly_effort_gate"][f"{rad}km"]["disjoint_source_only_quartet_details"][j-1]["all_4_sampling_arm_source_photo_sufficient"]
            assert orig==gate, "Original source gate cannot be silently recalculated differently"
            if gate:candidates.append({"quartet":j,"block":members[0]["original_quartet_2deg_block"]})
            sf=math.log1p(vals["source_downstream"])-math.log1p(vals["source_upstream"])
            cf=math.log1p(vals["registered_source_zero_downstream"])-math.log1p(vals["registered_source_zero_upstream"])
            contrasts.append(sf-cf)
            sources.append((math.log1p(vals["source_downstream"])+math.log1p(vals["source_upstream"])
                           -math.log1p(vals["registered_source_zero_downstream"])
                           -math.log1p(vals["registered_source_zero_upstream"]))/2.)
        thresh=o["prebutterfly_effort_gate"][f"{rad}km"]["predeclared_gate"]
        expected=int(len(candidates)>=thresh["min_eligible_quartets"] and
                     len({q["block"] for q in candidates})>=thresh["min_geographic_2deg_blocks"])
        assert int(o["prebutterfly_effort_gate"][f"{rad}km"]["photo_outcome_query_allowed"])==expected
        counts.sort()
        out["by_radius"][f"{rad}km"]={
            "effort_source_only_measured_station_n":len(counts),
            "source_event_photo_count_min":counts[0],"median":(counts[15]+counts[16])/2,
            "max":counts[-1],"zero_sites":sum(v==0 for v in counts),
            "sites_fewer_than_5_plant_photos":sum(v<5 for v in counts),
            "eligible_quartets":len(candidates),
            "eligible_original_two_degree_blocks":len({q["block"] for q in candidates}),
            "eligible_quartet_indices":[q["quartet"] for q in candidates],
            "continue_to_2022_Lycaenidae_photo_response":bool(expected),
            "source_minus_registeredzero_log1p_plant_effort_equal_quartet":exact_signflip(sources),
            "source_vs_registeredzero_updown_difference_in_log_plant_effort":exact_signflip(contrasts),
            "these_2019_2021_plantae_records_are_photo_effort_not_control_for_all_nutrient_plant_response":True,
        }
    return out


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--out",type=Path,required=True)
    v=a.parse_args()
    root=Path(__file__).resolve().parents[1]
    path=root/"results/lycaenidae_2022_wise6_directed_32_station_photo_sourceonly/2019_2021_ORIGINAL32_2KM_5KM_PLANT_ONLY_PREBUTTERFLY_GATE.json"
    x=audit(path)
    assert x["by_radius"]["2km"]["eligible_quartets"]==1
    assert x["by_radius"]["5km"]["eligible_quartets"]==1
    assert not x["by_radius"]["2km"]["continue_to_2022_Lycaenidae_photo_response"]
    assert not x["by_radius"]["5km"]["continue_to_2022_Lycaenidae_photo_response"]
    v.out.parent.mkdir(parents=True,exist_ok=True)
    v.out.write_text(json.dumps(x,indent=2,sort_keys=True)+"\n")
    print(json.dumps(x,indent=2))


if __name__=="__main__":
    main()
