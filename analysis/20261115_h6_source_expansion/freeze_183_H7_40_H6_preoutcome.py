#!/usr/bin/env python3
"""Pre-outcome genuine H6 parent/H7 population gate: 183 unique H7 in 40 H6.

Source IDs only from the authentic original period-paired H7 EIV source cohort
and the original source-ranked ten new H7 with native 1219 environmental replay.
No NUTRI/HNS numeric or downstream ecological effect is opened here.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from collections import defaultdict
from pathlib import Path

PANEL_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
GATE_SCHEMA="GESUI_EU454_FULLSOURCE_PLANT_NUTRI_TRAIT_GATE_V2_HYBRID_QA"
ENV_SCHEMA="GESUI_20261115_H7_1219_NATIVE_SOURCE_W_ENVIRONMENT_REPLAY_V1"

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    with p.open("r",encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))
def write(p,rows):
    with p.open("w",encoding="utf-8",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
def one(d,name):
    all=list(d.rglob(name))
    if len(all)!=1:raise RuntimeError("Source file not unique "+name+": "+str(all))
    return all[0]

def main(a):
    base=Path("results/20261115_sync/hns_macro/H7_454_HNS_panel.csv")
    if sha(base)!=PANEL_SHA:raise RuntimeError("Original fixed H7 study panel SHA drift")
    main=read(base)
    if len(main)!=454 or len({x["h"] for x in main})!=232:
        raise RuntimeError("Original H7 population drift")
    gate=json.loads(one(a.prior,"EU454_FULL_EIV_SOURCE_GATE.json").read_text())
    if gate.get("schema")!=GATE_SCHEMA or not gate.get("hard_source_support_pass"):
        raise RuntimeError("Original EIV source gate did not complete")
    el=read(one(a.prior,"EU454_NUTRI_SOURCE_ELIGIBILITY_PER_H7.csv"))
    if len(el)!=454 or len({(x["h"],x["period"]) for x in el})!=454:
        raise RuntimeError("Original source EIV site-period qualification changed")
    full=defaultdict(list)
    for row in el:
        if str(row["full_source_trait_eligible"]).lower()=="true":
            full[row["h"]].append(row)
    paired={}
    for h,rows in full.items():
        if len(rows)!=2:continue
        if {z["period"] for z in rows}!={"2017-2019","2020-2022"}:
            raise RuntimeError("Original supported H7 has incorrect two periods")
        if len({z["h6"] for z in rows})!=1:raise RuntimeError("Original H7 changes H6 parent")
        paired[h]=rows[0]["h6"]
    if len(paired)!=219 or len(set(paired.values()))!=86:
        raise RuntimeError("Former 219 two-period H7 /86 H6 source status changed")
    byh6=defaultdict(list)
    for h,h6 in paired.items():byh6[h6].append(h)
    multi={h6 for h6,hs in byh6.items() if len(hs)>=2}
    singleton={h6 for h6,hs in byh6.items() if len(hs)==1}
    if len(multi)!=30 or len(singleton)!=56:
        raise RuntimeError("Former 30 multi-H7 /56 singleton H6 source status changed")
    env=json.loads(one(a.environment,"H7_1219_SOURCE_REPLAY_W_ENVIRONMENT_QA.json").read_text())
    if env.get("schema")!=ENV_SCHEMA or not env.get("strict_new_model_environment_source_gate_pass") or env.get("original_old_H7_mismatch_n")!=0:
        raise RuntimeError("Original 232 H7 W and HydroATLAS replay did not exactly pass")
    chosen=read(one(a.environment,"NEW_TEN_H6_DETERMINISTIC_SOURCE_ONE_H7_EACH.csv"))
    if len(chosen)!=10 or len({q["native_H7"] for q in chosen})!=10 or len({q["new_official_H6"] for q in chosen})!=10:
        raise RuntimeError("Original new H7 one-per-singleton-H6 source selection changed")
    addH6={z["new_official_H6"] for z in chosen}
    if not addH6.issubset(singleton) or addH6&multi:
        raise RuntimeError("Extra true H7 not inside original single-H7 H6")
    if any(z["native_H7"] in paired for z in chosen):
        raise RuntimeError("New H7 was already original study source")
    source=[]
    # H6 parents with two existing H7 remain exactly fixed.
    for h,h6 in sorted(paired.items(),key=lambda kv:int(kv[0])):
        if h6 in multi or h6 in addH6:
            source.append({"h":h,"h6":h6,"source":"original_old_H7","original_two_period_NUTRI_source_eligible":True})
    for z in chosen:
        source.append({"h":z["native_H7"],"h6":z["new_official_H6"],
            "source":"source_only_original_EEA_new_H7","original_two_period_NUTRI_source_eligible":True})
    groups=defaultdict(set)
    for z in source:groups[z["h6"]].add(z["h"])
    if len(source)!=183 or len({z["h"] for z in source})!=183 or len(groups)!=40 or min(map(len,groups.values()))<2:
        raise RuntimeError("True 183 distinct geographic H7 in 40 multi-H7 source groups not achieved")
    out={
      "schema":"GESUI_TERRESTRIAL_ONLY_NEW_WITHIN_H6_ORIGINAL_GEO_COHORT_BEFORE_OUTCOMES_V1",
      "original_H7_SHA256":PANEL_SHA,
      "source_original_two_period_H7":219,
      "source_original_two_period_true_H6":86,
      "source_old_eligible_multi_H6":30,
      "source_old_h7_inside_30_multi_H6":163,
      "source_old_singleton_H6":56,
      "source_selected_new_H7":len(chosen),
      "source_selected_new_true_H6":len(addH6),
      "source_final_true_H6":len(groups),
      "source_final_distinct_H7":len(source),
      "source_original_H7_periods_not_counted_as_independent_geography":True,
      "source_H6_min_two_different_H7":min(map(len,groups.values())),
      "source_new_H7_selection_preoutcome_original_station_rank":True,
      "no_new_original_EEA_both_year_EQR_inferred":True,
      "plant_N_and_butterfly_HNS_coefficients_not_yet_fitted":True,
      "source_support_gate_pass":True
    }
    a.output.mkdir(parents=True,exist_ok=True)
    write(a.output/"FIXED_183_H7_IN_40_NATIVE_H6_SOURCE_COHORT.csv",source)
    (a.output/"FIXED_NEW_40H6_SOURCE_COHORT_GATE.json").write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n")
    print("GESUI_TERRESTRIAL_ONLY_40H6_SOURCE_GATE="+json.dumps(out,ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--prior",required=True,type=Path)
    p.add_argument("--environment",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
