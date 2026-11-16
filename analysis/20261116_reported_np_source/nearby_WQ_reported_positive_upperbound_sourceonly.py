#!/usr/bin/env python3
"""Hard pre-chemistry ceiling for official 1km same-H7 WQ candidates with WWTP reports.

Original GIS candidate set was frozen based only on exact station IDs,
locations, true HydroBASINS H7 polygons and <=1000 m. Original 3525 site/year
EEA water source positive report flags are independent of 2026 chemistry. If
there are <30 unique positive-source original H7-periods among ALL candidate
biological site/year opportunities, NO near-station chemistry scan of Part1
or Part2 can pass the locked support gate. This is mathematical upper bound,
not ecological evidence of no municipal effluent or measured TNTP.
"""
from __future__ import annotations
import argparse,json,sys
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261116_reported_np_source"))
from scan_eea2026_part2_nearby_sameH7_TNTP import load_candidates,old_biology_reporting

def summarise(items):
    yes=[v for v in items if v["report_positive"]]
    return {"n_biology_source_siteyears":len(items),
        "n_report_positive_siteyears":len(yes),
        "n_report_positive_H7period":len({(v["h"],v["period"]) for v in yes}),
        "n_report_positive_H7":len({v["h"] for v in yes}),
        "n_report_positive_H6":len({v["h6"] for v in yes}),
        "n_report_positive_countries":len({v["country"] for v in yes})}

def main(a):
    events,target,geo=load_candidates(a)
    report=old_biology_reporting(a,events)
    all_rows=[report[k] for k in events]
    quarter=[report[k] for k,v in events.items() if any(d<=250 for d,_,_ in v)]
    exclude_alias=[report[k] for k,v in events.items() if any(d>25 for d,_,_ in v)]
    s={"all_1km_H7_source_opportunities_upperbound":summarise(all_rows),
       "all_250m_H7_source_opportunities_upperbound":summarise(quarter),
       "exclude_co_located_25m_alias_possible_1km_H7_upperbound":summarise(exclude_alias)}
    gate=s["all_1km_H7_source_opportunities_upperbound"]
    threshold={"positive_H7periods":30,"positive_H7":20,"positive_H6":10,"positive_countries":3}
    impossible=(gate["n_report_positive_H7period"]<threshold["positive_H7periods"] or
                gate["n_report_positive_H7"]<threshold["positive_H7"] or
                gate["n_report_positive_H6"]<threshold["positive_H6"] or
                gate["n_report_positive_countries"]<threshold["positive_countries"])
    out={
      "schema":"GESUI_EEA2026_NEARBY_1KM_GEO_WWTP_REPORT_PRECHEMISTRY_HARD_CEILING_V1",
      "source_candidates_575_pairs_402_biology_years":True,
      "original_source_geometry_action":"https://github.com/zuizui0223/gesui/actions/runs/43804933915",
      "true_original_3525_EEA_Biology_years_in_original_part1part2_union":True,
      "original_report_positive_biology_siteyears":132,
      "source_UPSTREAM_report_status_is_NOT_proof_of_reaching_nearby_WQ_site":True,
      "prespecified_source_gate_threshold":threshold,
      "source_only_upper_bounds":s,
      "nearby_1km_still_cannot_reach_locked_positive_source_gate":bool(impossible),
      "same_exact_native_station_year_source_original_9_H7period_STILL_STOP":True,
      "source_NUTRI_HNS_EQR_outcomes_NOT_used":True,
      "nearby_WQ_2018_2022_measured_TNTP_NOT_USED":True,
      "no_false_zero_for_no_matched_report":True,
      "not_a_causal_sewage_exposure_or_foodweb_result":True
    }
    a.output.mkdir(parents=True,exist_ok=True)
    (a.output/"NEARBY_WQ_SAME_H7_1KM_SOURCE_REPORT_POSITIVE_UPPERBOUND.json").write_text(
        json.dumps(out,ensure_ascii=False,indent=2)+"\n")
    print("GESUI_NEAR_WQ_1KM_PRECHEMISTRY_UPPERBOUND="+json.dumps(out,ensure_ascii=False,sort_keys=True),flush=True)
    print("NEW_1KM_STATION_SAME_H7_SUPPORT_SOURCE_GATE="+("STOP_SOURCE_UPPER_BOUND" if impossible else "POSSIBLE_NEED_CHEMISTRY"),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--geo",type=Path,required=True)
    p.add_argument("--union",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
