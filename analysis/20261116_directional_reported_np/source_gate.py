#!/usr/bin/env python3
"""Outcome-blind quality gate for genuine EU H7-period reported WWTP direction.

Only source geography/metadata and published plant EIV score-availability
flags. No plant NUTRI, HNS regression, actual receiving-river TN/TP or
true zero-sewage control is generated.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd

SOURCE_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
DIRECTION_BLOB_SHA="702131afc7d1a7bfb15e80357a506754cae7a70a"
SOURCE_SCHEMA="GESUI_EU454_FULLSOURCE_PLANT_NUTRI_TRAIT_GATE_V2_HYBRID_QA"
DIRECTION="results/20261115_np_direction/H7_454_DIRECTIONAL_SOURCE_AUDIT_PANEL.csv"
PANEL="results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
def gitblob(path):
    b=path.read_bytes()
    return hashlib.sha1(f"blob {len(b)}\0".encode()+b).hexdigest()

def main(a):
    prior=Path(PANEL);direction=Path(DIRECTION)
    if sha(prior)!=SOURCE_SHA:raise RuntimeError("Original 454 H7-period scientific source SHA drift")
    if gitblob(direction)!=DIRECTION_BLOB_SHA:raise RuntimeError("Frozen original PR92 454 up/down report source Git blob drift")
    df=pd.read_csv(prior,dtype={"h":str,"h6":str})
    exposure=pd.read_csv(direction,dtype={"h":str,"h6":str})
    if (len(df),df.h.nunique(),df.h6.nunique())!=(454,232,90) or len(exposure)!=454:
        raise RuntimeError("Original fixed 454 H7 periods 232 H7 90 genuine H6 changed")
    if df.duplicated(["h","period"]).any() or exposure.duplicated(["h","period"]).any():
        raise RuntimeError("Same original H7 period repeated")
    merged=df[["h","period","h6","country","W","primary_record_mean"]].merge(
        exposure,on=["h","period"],validate="one_to_one",how="outer",
        suffixes=("_orig","_direction"),indicator=True)
    if len(merged)!=454 or not (merged["_merge"]=="both").all():
        raise RuntimeError("Different original H7-year geography in 2018/22 WWTP source")
    if not (merged.h6_orig==merged.h6_direction).all():
        raise RuntimeError("Official H7 to true H6 source parent changed")
    if float((merged.W_orig-merged.W_direction).abs().max())>1e-10:
        raise RuntimeError("Fixed original modeled static WWTP source W differs")
    if float((merged.primary_record_mean_orig-merged.primary_record_mean_direction).abs().max())>1e-10:
        raise RuntimeError("Previously source-scored HNS record values disagree with PR92 report join")
    up=merged.upstream_report_positive.to_numpy(int)
    down=merged.downstream_report_positive.to_numpy(int)
    if set(up)!=set([0,1]) or set(down)!=set([0,1]):
        raise RuntimeError("Report-positive flags must stay 0/1 metadata")
    observed={(int(u),int(v)):int(((up==u)&(down==v)).sum()) for u in (0,1) for v in (0,1)}
    if observed!={(0,0):372,(1,0):21,(0,1):26,(1,1):35}:
        raise RuntimeError("Original 2x2 directional WWTP report status matrix changed: "+str(observed))
    v2=json.loads((a.eiv/"EU454_FULL_EIV_SOURCE_GATE.json").read_text())
    if (v2.get("schema")!=SOURCE_SCHEMA or not v2.get("hard_source_support_pass")
        or not v2.get("original_EIV_NUTRI_scored_keys_unchanged")
        or v2.get("official_confirmed_valid_unscorable_nonbinomial_hybrids")!=175
        or v2.get("source_hybrid_score_assigned")!=0):
        raise RuntimeError("Official frozen EIV species original hybrid source gate v2 is not authentic")
    eligible=pd.read_csv(a.eiv/"EU454_NUTRI_SOURCE_ELIGIBILITY_PER_H7.csv",dtype={"h":str,"h6":str})
    if len(eligible)!=454 or eligible.duplicated(["h","period"]).any():
        raise RuntimeError("Original EIV EEA/GBIF source eligibility rows changed")
    permitted=eligible.loc[eligible.full_source_trait_eligible.astype(str)=="True"].copy()
    if len(permitted)!=451 or permitted.h6.nunique()!=90:
        raise RuntimeError("Original 451 strict trait quality source count changed")
    s=merged.merge(permitted[["h","period","h6"]],left_on=["h","period","h6_orig"],
                   right_on=["h","period","h6"],how="inner",validate="one_to_one")
    if len(s)!=451:raise RuntimeError("The 451 original source trait cohort differs from directed WWTP source")
    s["report_fourcell"]=(s.upstream_report_positive.astype(str)+"_"+s.downstream_report_positive.astype(str))
    grid=s.report_fourcell.value_counts().to_dict()
    pos=s.loc[s.upstream_report_positive==1]
    dpos=s.loc[s.downstream_report_positive==1]
    for col in ("upstream_N_equiv_mgL","downstream_N_equiv_mgL",
                "upstream_P_equiv_mgL","downstream_P_equiv_mgL"):
        numeric=pd.to_numeric(s[col],errors="raise")
        if not numeric.notna().all() or (numeric<0).any():
            raise RuntimeError("Original reported annual potential delivered N/P source index invalid: "+col)
    support=(len(pos)>=30 and pos.h6_orig.nunique()>=20 and pos.country_orig.nunique()>=3
          and len(dpos)>=30 and all(grid.get(c,0)>=10 for c in ("1_0","0_1","1_1"))
          and s.h6_orig.nunique()==90)
    out={"schema":"GESUI_20261116_REPORT_POSITIVE_UPDOWN_EIV_PLANT_BUTTERFLY_SOURCE_GATE_V1",
         "original_H7panel_SHA256":SOURCE_SHA,
         "original_direction_panel_git_blob":DIRECTION_BLOB_SHA,
         "original_EIV_accepted_numeric_species_n":4011,
         "original_source_rows":454,"original_source_H7":232,"original_source_H6":90,
         "original_direction_fourcells_whole_original_454":{str(k):v for k,v in observed.items()},
         "original_direction_country_mismatch_rows":int((merged.country_orig!=merged.country_direction).sum()),
         "eligible_source_H7periods":len(s),"eligible_true_H6":s.h6_orig.nunique(),
         "eligible_source_fourcells":grid,
         "eligible_upstream_positive_H7periods":len(pos),
         "eligible_upstream_positive_H7":int(pos.h.nunique()),
         "eligible_upstream_positive_H6":int(pos.h6_orig.nunique()),
         "eligible_upstream_positive_country_labels":int(pos.country_orig.nunique()),
         "eligible_downstream_positive_H7periods":len(dpos),
         "eligible_downstream_positive_true_H6":int(dpos.h6_orig.nunique()),
         "source_gate_support_pass":bool(support),
         "hard_thresholds":{"up_positive_rows_min":30,"up_positive_H6_min":20,
                            "up_positive_countries_min":3,"down_positive_rows_min":30,
                            "each_non_neither_report_cell_min":10},
         "NO_MATCHED_POSITIVE_REPORT_is_NOT_zero_sewage":True,
         "original_reported_annual_NP_load_not_measured_receptor_water_chemistry":True,
         "no_biological_numeric_effect_fitted":True,
         "result_outcome_not_opened_before_source_gate":True}
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"UPDOWN_EIV_PLANT_BUTTERFLY_SOURCE_GATE.json").write_text(
        json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("GESUI_NEXT_SOURCE_UPDOWN_EIV_GATE="+json.dumps(out,ensure_ascii=False,sort_keys=True),flush=True)
    print("MODEL_GATE="+("PASS_REPORTED_SOURCE_ONLY" if support else "STOP_REPORTED_SOURCE_SUPPORT"),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--eiv",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    main(p.parse_args())
