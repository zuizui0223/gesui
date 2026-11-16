#!/usr/bin/env python3
"""Original-source-only 2018/2022 EEA reported N/P within-true-H6 opportunity.

DO NOT OPEN biological outputs. Three source indices are original byte-exact
Git blobs from PR90/PR92, full original H7 panel from integrated main,
frozen 183-H7/40-H6 population from exact original Actions source artifact.

Reported-positive versus no-positive-report is *reporting*, NOT wastewater
treatment control; N/P-equivalent values are mass / model flow, NOT measured
river concentration. New 10 H7 are MISSING historical reporting index, NOT 0.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from collections import defaultdict, Counter
from pathlib import Path
import pandas as pd

SOURCE_183_SHA="c0840b1c5c6e85342fea41d57dc40fa5ec295fabf6bcd6908a4072557db9a736"
SOURCE_454_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
REPORT_BLOB="45979613c1bc079f1c2346461eb7bdd31ac17e6a"
DIRECTION_BLOB="702131afc7d1a7bfb15e80357a506754cae7a70a"
QA_BLOB="dd283aef32eaa52a197bf70bd3bdf49c87a156aa"
PERIODS={"2017-2019","2020-2022"}
MIN_REPORT_POSITIVE_PERIODS=30
MIN_REPORT_POSITIVE_DISTINCT_H7=20
MIN_H6_WITH_TWO_POSITIVE_DISTINCT_H7=10
MIN_REPORT_POSITIVE_COUNTRY=5
MIN_BOTH_PERIOD_POSITIVE_H7=5

def sha(p:Path)->str:
    return hashlib.sha256(p.read_bytes()).hexdigest()

def blob(p:Path)->str:
    raw=p.read_bytes()
    return hashlib.sha1(("blob "+str(len(raw))+"\0").encode()+raw).hexdigest()

def rootfile(base:Path,name:str)->Path:
    found=list(base.rglob(name))
    if len(found)!=1:raise RuntimeError("Genuine original source file missing/ambiguous "+name+" "+repr([str(x) for x in found]))
    return found[0]

def read_csv_only(path,cols):
    return pd.read_csv(path,dtype={"h":str,"h6":str,"period":str,"country":str,
                       "original_source_cohort":str},usecols=cols)

def main(args):
    old_source=args.mainpanel
    if sha(old_source)!=SOURCE_454_SHA:
        raise RuntimeError("Original main 454 H7 science panel changed")
    main=read_csv_only(old_source,["h","h6","period","country"])
    if (len(main),main.h.nunique(),main.h6.nunique())!=(454,232,90) or main.duplicated(["h","period"]).any():
        raise RuntimeError("Original 454 source H7/H6 geography drift")
    if set(main.period)!=PERIODS:raise RuntimeError("Genuine original periods changed")
    frozen=rootfile(args.original183,"GENUINE_TERRESTRIAL_ONLY_183_H7_40_H6_PAIRED_H7_MEANS.csv")
    if sha(frozen)!=SOURCE_183_SHA:raise RuntimeError("Frozen 183 native H7 source population SHA drift")
    population=read_csv_only(frozen,["h","h6","original_source_cohort"])
    if (len(population),population.h.nunique(),population.h6.nunique())!=(183,183,40) or population.h.duplicated().any():
        raise RuntimeError("Fixed 183 H7/40 H6 population drift")
    nnew=population.original_source_cohort.ne("original_old_H7")
    if int(nnew.sum())!=10 or int((~nnew).sum())!=173:
        raise RuntimeError("Original immutable 173 old / 10 newly sourced H7 cohort changed")
    new=population[nnew];old=population[~nnew]
    source=pd.read_csv(args.reported,dtype={"h":str,"h6":str,"period":str,"country":str},
       usecols=["h","h6","period","country","reported_positive",
           "reporting_class","n_source","n_mg_l","p_mg_l"])
    direction=pd.read_csv(args.direction,dtype={"h":str,"h6":str,"period":str},
       usecols=["h","h6","period","upstream_report_positive","downstream_report_positive"])
    if blob(args.reported)!=REPORT_BLOB or blob(args.direction)!=DIRECTION_BLOB or blob(args.qa)!=QA_BLOB:
        raise RuntimeError("Copy of original source-only 454 reported-NP/updown Git blob changed")
    qa=json.loads(args.qa.read_text())
    if qa.get("source",{}).get("source_sha256")!="012ad79feabc8c11e0d74f803730b5b7913b96a64a6af26fb0c97263d9aa7c26":
        raise RuntimeError("Published original EEA 2018/22 reported mass source manifest SHA drift")
    for dat,label in [(source,"report"),(direction,"direction")]:
        if len(dat)!=454 or dat.h.nunique()!=232 or dat.duplicated(["h","period"]).any() or dat.h6.nunique()!=90:
            raise RuntimeError("Original 454 joint source cohort drift: "+label)
        if set(zip(dat.h,dat.period))!=set(zip(main.h,main.period)):
            raise RuntimeError("Original 454 basin-period key mismatch: "+label)
        test=main.merge(dat[["h","period","h6"]],on=["h","period"],validate="one_to_one")
        if not test.h6_x.eq(test.h6_y).all():raise RuntimeError("Original true H6 parent contradiction: "+label)
    if (len(source[source.reported_positive==1]),
        source.loc[source.reported_positive==1,"h"].nunique(),
        source.loc[source.reported_positive==1,"h6"].nunique())!=(56,38,31):
        raise RuntimeError("Original 56/38/31 report-positive source base changed")
    if sorted(source.groupby("period").reported_positive.sum().to_list())!=[28,28]:
        raise RuntimeError("Original balanced per-period positive report count changed")
    full=source.merge(direction,on=["h","h6","period"],validate="one_to_one")
    if len(full)!=454 or not full.reported_positive.eq(full.upstream_report_positive).all():
        raise RuntimeError("Source-positive status and independent upstream directional flag inconsistent")
    if not full.reported_positive.isin([0,1]).all() or not full.downstream_report_positive.isin([0,1]).all():
        raise RuntimeError("Source reporting state not boolean")
    if full["n_mg_l"].isna().any() or full["p_mg_l"].isna().any() or (
        (full[["n_mg_l","p_mg_l"]]<0).any().any()):
        raise RuntimeError("Reported source N/P theoretical mass/flow quantity invalid")
    if not ((full.reported_positive==1)==(full.n_mg_l.gt(0)&full.p_mg_l.gt(0))).all():
        raise RuntimeError("Reported source positive must agree with both N and P")
    sub=old.merge(full,on=["h","h6"],validate="one_to_many",how="left")
    if len(sub)!=346 or sub.h.nunique()!=173 or sub.h6.nunique()!=40:
        raise RuntimeError("Old 173 true H7 * two source periods not completely present")
    if sub[["period","reported_positive","upstream_report_positive","downstream_report_positive"]].isna().any().any():
        raise RuntimeError("Unreported source missing and incorrectly interpreted as a zero")
    if set(zip(sub.h,sub.period))!=set((h,t) for h in old.h for t in PERIODS):
        raise RuntimeError("Genuine selected 173 old H7 must have both report-source periods")
    if len(set(new.h)&set(source.h)):
        raise RuntimeError("New ten H7 unexpectedly have original 454 source-reporting data")
    if old.h6.nunique()!=40 or new.h6.nunique()!=10:
        raise RuntimeError("Source new ten true-H6 parents altered")
    positive=sub.loc[sub.reported_positive.eq(1)].copy()
    h7_flags=sub.groupby(["h6","h"],sort=False).agg(
        ever_up_positive=("upstream_report_positive","max"),
        ever_down_positive=("downstream_report_positive","max"),
        source_report_positive_periods=("reported_positive","sum")
    ).reset_index()
    byh6=h7_flags.groupby("h6").agg(
         distinct_old_h7_with_reporting_index=("h","nunique"),
         distinct_old_h7_upstream_ever_positive=("ever_up_positive","sum"),
         distinct_old_h7_downstream_ever_positive=("ever_down_positive","sum")
    ).reset_index()
    original_group_count=population.groupby("h6").h.nunique().rename("n_total_source_h7")
    byh6=byh6.merge(original_group_count,on="h6",validate="one_to_one")
    countries=(sub[["h","h6","country"]].drop_duplicates()
              .groupby("h6").country.agg(lambda v:";".join(sorted(set(v)))).rename("source_country_labels"))
    byh6=byh6.join(countries,on="h6")
    byh6["positive_upstream_H7_and_no_positive_report_H7"]=(
        (byh6.distinct_old_h7_upstream_ever_positive>=1) &
        (byh6.distinct_old_h7_with_reporting_index>byh6.distinct_old_h7_upstream_ever_positive))
    byh6["two_distinct_upstream_positive_H7"]=byh6.distinct_old_h7_upstream_ever_positive>=2
    byh6["new_H7_source_reporting_unknown"]=byh6.h6.isin(set(new.h6))
    byh6["source_info_is_not_unexposed_control"]=True
    pair_count=(sub.groupby("h").reported_positive.sum()==2).sum()
    positive_country=int(positive.country.nunique())
    pos_h6_two=int(byh6.two_distinct_upstream_positive_H7.sum())
    gate=bool(len(positive)>=MIN_REPORT_POSITIVE_PERIODS and
        positive.h.nunique()>=MIN_REPORT_POSITIVE_DISTINCT_H7 and
        pos_h6_two>=MIN_H6_WITH_TWO_POSITIVE_DISTINCT_H7 and
        positive_country>=MIN_REPORT_POSITIVE_COUNTRY and
        pair_count>=MIN_BOTH_PERIOD_POSITIVE_H7)
    mixed=int(byh6.positive_upstream_H7_and_no_positive_report_H7.sum())
    summary={
      "schema":"GESUI_20261116_SAME_183H7_40H6_EEA_REPORTED_NP_SOURCE_ONLY_GATING_V1",
      "original_454_H7_panel_SHA256":SOURCE_454_SHA,
      "original_183_H7_source_SHA256":SOURCE_183_SHA,
      "reported_index_original_git_blob_SHA1":REPORT_BLOB,
      "independent_direction_original_git_blob_SHA1":DIRECTION_BLOB,
      "source_original_EEA_upstream_NP_2018_2022_SHA256":"012ad79feabc8c11e0d74f803730b5b7913b96a64a6af26fb0c97263d9aa7c26",
      "original_population_distinct_H7":183,"original_population_true_H6":40,
      "source_old_H7_with_454_annual_source_indices":173,
      "source_new_H7_with_NO_corresponding_original_454_report_indices":10,
      "source_new_H7_index_missing_NOT_zero_exposure":True,
      "old_H7_original_two_period_report_rows":len(sub),
      "old_H7_upstream_report_positive_H7_periods":int(len(positive)),
      "old_H7_distinct_upstream_report_positive_H7":int(positive.h.nunique()),
      "old_H7_true_H6_with_any_upstream_report_positive":int(positive.h6.nunique()),
      "old_H7_country_labels_with_report_positive":positive_country,
      "old_H7_both_periods_positive_H7":int(pair_count),
      "true_H6_with_TWO_distinct_report_positive_old_H7":pos_h6_two,
      "true_H6_with_mixed_positive_and_no_matching_positive_source_REPORT":mixed,
      "true_H6_with_any_new_H7_report_index_unknown":int(byh6.new_H7_source_reporting_unknown.sum()),
      "old_H7_downstream_report_positive_periods":int(sub.downstream_report_positive.sum()),
      "original_all_454_report_positive_periods":56,
      "original_all_454_directional_upstream_positive_periods":56,
      "original_all_454_directional_downstream_positive_periods":int(direction.downstream_report_positive.sum()),
      "minimum_source_data_support_before_any_receptor_biology_model":{
           "upstream_positive_H7periods":MIN_REPORT_POSITIVE_PERIODS,
           "distinct_upstream_positive_H7":MIN_REPORT_POSITIVE_DISTINCT_H7,
           "distinct_true_H6_with_at_least_two_upstream_positive_H7":MIN_H6_WITH_TWO_POSITIVE_DISTINCT_H7,
           "report_positive_country_labels":MIN_REPORT_POSITIVE_COUNTRY,
           "both_periods_upstream_positive_H7":MIN_BOTH_PERIOD_POSITIVE_H7},
      "source_upstream_report_intensity_local_parent_replication_gate_pass":gate,
      "source_decision":"PASS_REPORTED_MASS_INDEX_ONLY_NOT_CAUSAL" if gate else "STOP_SAME_H6_REPORTED_INTENSITY_REPLICATION",
      "genuine_field_year_same_site_measured_TN_TP_available_from_this_source":False,
      "source_mass_to_modelled_flow_is_NOT_measured_ambient_nutrient_concentration":True,
      "GBIF_butterfly_plant_field_patches_not_exact_EEA_station":True,
      "reported_index_zero_means_NO_MATCHED_POSITIVE_REPORT_NOT_PRISTINE":True,
      "original_biological_HNS_EIV_EQR_W_COLUMNS_NOT_READ":True,
      "report_status_transitions_are_not_proven_facility_operating_changes":True,
      "no_effect_model_run":True,
      "not_a_causal_mediation_or_ecological_service_assessment":True}
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/"183H7_40H6_EEA_REPORTED_NP_SOURCE_ONLY_GATE.json").write_text(
        json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    byh6.sort_values("h6").to_csv(args.output/"183H7_40H6_REPORTING_SOURCE_H6_REPLICATION_AUDIT.csv",index=False)
    print("GESUI_NP_183H7_40H6_SOURCE_GATE="+json.dumps(summary,ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--original183",type=Path,required=True)
    parser.add_argument("--mainpanel",type=Path,default=Path("results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"))
    parser.add_argument("--reported",type=Path,default=Path("results/20261115_reported_np_hns/H7_454_REPORTED_SOURCE_INDEX_HNS_EQR_PLANT.csv"))
    parser.add_argument("--direction",type=Path,default=Path("results/20261115_np_direction/H7_454_DIRECTIONAL_SOURCE_AUDIT_PANEL.csv"))
    parser.add_argument("--qa",type=Path,default=Path("results/20261115_reported_np_hns/SOURCE_QA_AND_EXPOSURE.json"))
    parser.add_argument("--output",type=Path,required=True)
    main(parser.parse_args())
