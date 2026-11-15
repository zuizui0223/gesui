#!/usr/bin/env python3
"""SOURCE ONLY: hard upper bound on extra same-era EEA H7 geographic replication.

Authentic original EEA18/EEA22/EPA19 H7 station source 6330 records, paired
with exact original 454 H7×period cohort /232 H7/90 genuine H6. No ecology
outcomes or coefficients are opened. This is NOT source gate relaxation.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from collections import Counter,defaultdict
from pathlib import Path

SOURCE_SHA="a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89"
SOURCE_NAME="EEA2018_EEA2022_EPA1819_NATIVE_HYBAS7_SOURCE_ONLY_GEOGRAPHY.json"
PANEL_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
H7PERIODS={"2017-2019","2020-2022"}
EXPECTED_EEA_BOTH_H7=234
EXPECTED_STUDY_H7=232
MIN_EXTRA_H6_NEEDED=10

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def original_only(folder,name):
    paths=list(folder.rglob(name))
    if len(paths)!=1:raise RuntimeError("Expected exactly one exact original source file: "+name+": "+repr([str(x) for x in paths]))
    return paths[0]

def main(args):
    path=original_only(args.station_source,SOURCE_NAME)
    if digest(path)!=SOURCE_SHA:raise RuntimeError("Original station-source original SHA256 mismatch")
    master=json.loads(path.read_text())
    rows=master.get("candidate_source_site_period_rows")
    if not isinstance(rows,list) or len(rows)!=6330:
        raise RuntimeError("Original six thousand three hundred thirty geographic station rows changed")
    if digest(args.study_panel)!=PANEL_SHA:
        raise RuntimeError("Original source-panel 454 H7-period SHA mismatch")
    with args.study_panel.open(encoding="utf-8-sig",newline="") as f:
        study=list(csv.DictReader(f))
    if len(study)!=454:raise RuntimeError("Not original 454 H7-period cohort")
    study_keys={(r["h"],r["period"]) for r in study}
    old_h7={r["h"] for r in study}
    if len(study_keys)!=454 or len(old_h7)!=232 or len({r["h6"] for r in study})!=90:
        raise RuntimeError("Study source unit identity/parent count changed")
    if {r["period"] for r in study}!=H7PERIODS:
        raise RuntimeError("Study source periods changed")
    # The source metadata has separate EEA18/EEA22 panels from the US EPA.
    seen_by_panel=defaultdict(set)
    original_countries=defaultdict(set)
    source_visit_rows=Counter()
    stationary=defaultdict(set)
    skipped=Counter()
    missing_id=0
    for r in rows:
        if r.get("HYBAS_ID_L07") is None:
            missing_id+=1
            continue
        h=str(r["HYBAS_ID_L07"])
        panel=r.get("panel")
        if panel not in ("EEA18","EEA22"):
            skipped[str(panel)]+=1
            continue
        srcper=str(r.get("period",""))
        if srcper not in H7PERIODS:
            raise RuntimeError("EEA original H7 station has unrecognized fixed period: "+srcper)
        if (panel=="EEA18" and srcper!="2017-2019") or (
                panel=="EEA22" and srcper!="2020-2022"):
            raise RuntimeError("Source panel, period mismatch; never reassign survey years")
        seen_by_panel[h].add(panel)
        source_visit_rows[(h,panel)]+=1
        original_countries[h].add(str(r.get("country","")))
        # Source original physical ids only; no biological values.
        if r.get("physical_id") is not None:
            stationary[(h,panel)].add(str(r["physical_id"]))
    any_EEA={h for h in seen_by_panel}
    both={h for h,p in seen_by_panel.items() if p=={"EEA18","EEA22"}}
    only18={h for h,p in seen_by_panel.items() if p=={"EEA18"}}
    only22={h for h,p in seen_by_panel.items() if p=={"EEA22"}}
    if len(both)!=EXPECTED_EEA_BOTH_H7:
        raise RuntimeError("Original source two-period 234 EEA H7 no longer exactly match")
    if not old_h7.issubset(both):
        raise RuntimeError("Fixed original 232 study H7 not all in original 234 both-period EEA frame")
    extra_both=sorted(both-old_h7)
    if len(extra_both)!=2:
        raise RuntimeError("Original two EEA H7 not included in 232 fixed science panel drift")
    # Independent H6 parent labels of extra basins are unknown unless native
    # full geometric Level7→Level6 rejoin is performed. The MAX possible
    # improvement is at most one singleton H6 promoted by each extra H7.
    original_eligible_2period_H6_multi=30
    upper=min(90,original_eligible_2period_H6_multi+len(extra_both))
    mathematically_impossible=upper<40
    if not mathematically_impossible:
        raise RuntimeError("Unexpected failure of strict hard upper bound")
    raw_2018=sorted(only18-old_h7)
    raw_2022=sorted(only22-old_h7)
    candidates=[]
    for k in sorted(any_EEA-old_h7):
        panels=seen_by_panel[k]
        candidates.append({"original_HYBAS_ID_L07":k,
              "source_EEA18_station_rows":source_visit_rows[(k,"EEA18")],
              "source_EEA22_station_rows":source_visit_rows[(k,"EEA22")],
              "source_EEA18_distinct_physical_stations":len(stationary[(k,"EEA18")]),
              "source_EEA22_distinct_physical_stations":len(stationary[(k,"EEA22")]),
              "source_EEA_countries":";".join(sorted(original_countries[k])),
              "source_both_EEA_eras":panels=={"EEA18","EEA22"},
              "source_one_EEA_era_only":len(panels)==1,
              "original_study_H7_excluded":True,
              "official_H6_parent_of_new_H7_not_yet_verified":True,
              "original_plantae_EIV_trait_and_butterfly_species_matrix_NOT_ACQUIRED":True})
    result={
      "schema":"GESUI_SOURCE_FIRST_EEA_2018_2022_H7_EXPANSION_HARD_CEILING_V1",
      "station_source_SHA256":SOURCE_SHA,"source_station_original_rows":len(rows),
      "original_EEA_H7_2018_only":len(only18),
      "original_EEA_H7_2022_only":len(only22),
      "original_EEA_H7_in_any_era":len(any_EEA),
      "original_EEA_H7_in_both_2018_and_2022":len(both),
      "original_454_basinperiod_study_H7":len(old_h7),
      "original_study_454_source_SHA256":PANEL_SHA,
      "original_BOTH_period_EEA_H7_not_in_study":len(extra_both),
      "original_BOTH_period_EEA_H7_not_in_study_ids":extra_both,
      "min_additional_multih7_H6_needed_for_40":MIN_EXTRA_H6_NEEDED,
      "original_previous_219_H7_supported_30_H6_multiH7":original_eligible_2period_H6_multi,
      "upper_bound_multiH7_H6_after_adding_ALL_remaining_two_EEA_period_H7":upper,
      "proof_strict_original_two_era_EEA_frame_CANNOT_reach_40_H6":mathematically_impossible,
      "single_era_EEA_H7_outside_original_study_count":len(raw_2018)+len(raw_2022),
      "single_era_EEA_H7_are_NOT_substitutes_for_two_era_EEA_observation":True,
      "EEA_both_period_does_not_guarantee_plants_or_butterfly_GBIF_support":True,
      "new_1219_or_all_new_H7_H6_geometry_mapping_not_done":True,
      "without_new_external_H7_cannot_reach_prespecified_local_H6_gate":True,
      "no_EQR_or_NUTRI_or_HNS_or_wastewater_effects_opened":True,
      "EPA_US_stations_skipped":dict(skipped),
      "not_an_EEA_plant_butterfly_observational_inference":True
    }
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/"EEA_NATIVE_SOURCE_H7_EXPANSION_HARD_CEILING.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    with (args.output/"EEA_NATIVE_SOURCE_ONEERA_AND_EXTRA_H7_REGISTRY.csv").open("w",newline="",encoding="utf-8") as f:
        if not candidates:raise RuntimeError("No source candidates: original geography changed")
        w=csv.DictWriter(f,fieldnames=list(candidates[0]))
        w.writeheader();w.writerows(candidates)
    print("GESUI_EEA_SOURCE_EXPANSION_CEILING="+json.dumps(result,sort_keys=True,ensure_ascii=False),flush=True)
    print("SOURCE_GEO_ADDITION_GATE=STOP_ORIGINAL_EEA_BOTH_ERA_HARD_UPPER_BOUND",flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--station-source",type=Path,required=True)
    p.add_argument("--study-panel",type=Path,default=Path("results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"))
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
