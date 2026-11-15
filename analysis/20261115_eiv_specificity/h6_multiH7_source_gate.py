#!/usr/bin/env python3
"""Source-only true H6 multi-H7 spatial opportunity gate; do NOT open outcomes."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd

PANEL_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
FACETS_SHA="251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4"
EIV_SHA="ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba"
PERIODS={"2017-2019","2020-2022"}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main(a):
    gate=json.loads((a.v2/"EU454_FULL_EIV_SOURCE_GATE.json").read_text())
    if (gate.get("schema")!="GESUI_EU454_FULLSOURCE_PLANT_NUTRI_TRAIT_GATE_V2_HYBRID_QA"
        or not gate.get("hard_source_support_pass")
        or gate.get("source_facets_SHA256")!=FACETS_SHA
        or gate.get("source_EIV",{}).get("original_source_SHA256")!=EIV_SHA
        or gate.get("official_confirmed_valid_unscorable_nonbinomial_hybrids")!=175
        or gate.get("source_hybrid_score_assigned")!=0):
        raise RuntimeError("Only the original source-v2 H7 site-based trait quality gate accepted")
    p=Path("results/20261115_sync/hns_macro/H7_454_HNS_panel.csv")
    if sha(p)!=PANEL_SHA:raise RuntimeError("Genuine source original 454 panel SHA mismatch")
    full=pd.read_csv(p,dtype={"h":str,"h6":str})
    if len(full)!=454 or full.h.nunique()!=232 or full.h6.nunique()!=90:
        raise RuntimeError("Original true H7/H6 source graph altered")
    el=pd.read_csv(a.v2/"EU454_NUTRI_SOURCE_ELIGIBILITY_PER_H7.csv",dtype={"h":str,"h6":str})
    if len(el)!=454 or el.duplicated(["h","period"]).any() or not (set(el.period)==PERIODS):
        raise RuntimeError("Original 454 original trait quality source periods mismatch")
    d=el.loc[el.full_source_trait_eligible.astype(str)=="True"].copy()
    if len(d)!=451 or d.h.nunique()!=232 or d.h6.nunique()!=90:
        raise RuntimeError("Original 451 source quality cohort changed")
    eligible=[]
    for h,g in d.groupby("h"):
        if len(g)!=2:continue
        if set(g.period)!=PERIODS:raise RuntimeError("Duplicate source periods")
        if g.h6.nunique()!=1 or g.country_holdout.nunique()!=1:
            raise RuntimeError("Same H7 changes H6 parent / country")
        eligible.append({"h":h,"h6":str(g.h6.iloc[0]),"country_holdout":str(g.country_holdout.iloc[0])})
    frame=pd.DataFrame(eligible)
    if len(frame)!=219 or frame.h6.nunique()!=86:
        raise RuntimeError("Expected original source-supported paired 219 H7/86 H6 panel changed")
    h6_counts=frame.groupby("h6").h.nunique()
    keep=set(h6_counts.index[h6_counts>=2])
    members=frame.loc[frame.h6.isin(keep)]
    support=(len(members)>=100 and len(keep)>=40 and members.country_holdout.nunique()>=5)
    out={"schema":"GESUI_EU454_TWO_PERIOD_TRUE_H6_MULTI_H7_SOURCE_OPPORTUNITY_V1",
         "original_panel_SHA256":PANEL_SHA,"source_original_facets_SHA256":FACETS_SHA,
         "source_EIV_SHA256":EIV_SHA,"original_supported_source_H7period":451,
         "exact_two_period_source_H7":len(frame),
         "exact_two_period_source_true_H6":int(frame.h6.nunique()),
         "H6_with_at_least_two_distinct_H7":len(keep),
         "distinct_H7_within_multiH7_H6":len(members),
         "source_country_groups_within_multiH7_H6":int(members.country_holdout.nunique()),
         "multiH7_H6_original_H7_count_distribution":{str(k):int(v) for k,v in h6_counts.value_counts().sort_index().items()},
         "source_gate_preoutcome_thresholds":{"minimum_distinct_H7":100,"minimum_true_H6":40,
             "minimum_majority_country_groups":5,"minimum_distinct_H7_within_parent_H6":2},
         "within_H6_spatial_association_sample_gate_pass":bool(support),
         "no_plantN_or_HNS_observation_opened":True,
         "false_independent_periods_not_counted":True,
         "not_independent_ecological_field_studies":True}
    a.output.mkdir(parents=True,exist_ok=True)
    members.to_csv(a.output/"WITHIN_H6_ELIGIBLE_DISTINCT_H7_KEYS.csv",index=False)
    (a.output/"WITHIN_H6_SOURCE_ONLY_GATE.json").write_text(
        json.dumps(out,ensure_ascii=False,indent=2)+"\n")
    print("GESUI_WITHIN_H6_SOURCE_ONLY_GATE="+json.dumps(out,ensure_ascii=False,sort_keys=True),flush=True)
    print("MODEL_GATE="+("PASS" if support else "STOP_INDEPENDENT_MULTIH7_H6_SUPPORT"),flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--v2",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    main(ap.parse_args())
