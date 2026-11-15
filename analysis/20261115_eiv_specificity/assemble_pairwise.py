#!/usr/bin/env python3
"""Independent source-only matched NUTRI vs every one of five original EIV controls.

Posthoc *secondary* after prior six-indicator complete-case STOP.
Exact original source and the PR102 unchanged 4,011 accepted EIV-N matches.
At each site a pair uses precisely the SAME species for both indicators.
No W/HNS model is fitted here; source-only condition must pass per pair.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import re
from pathlib import Path
import numpy as np
import pandas as pd
import requests

EIV_SHA="ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba"
H7_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
FACET_SHA="251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4"
EIV_URL="https://www.sci.muni.cz/botany/juice/NET/sources/EIV_Europe.csv"
DIMENSIONS=["LIGHT","TEMPE","MOIST","REACT","NUTRI","SALIN"]
RIVALS=["LIGHT","TEMPE","MOIST","REACT","SALIN"]
RANGES={"LIGHT":(1,9),"TEMPE":(1,12),"MOIST":(1,12),"REACT":(1,9),"NUTRI":(1,9),"SALIN":(0,9)}
SOURCE_N_MATCHED=4011
SOURCE_PLANT_UNION=11511
MIN_N=10
MIN_FRAC=.40
MIN_H7PERIOD=100
MIN_H6=40
MIN_COUNTRY=5

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def one(folder,name):
    z=list(folder.rglob(name))
    if len(z)!=1:raise RuntimeError("Source file missing or duplicated "+name)
    return z[0]
def normal(s):
    return re.sub(r"\s+"," ",str(s).strip()).casefold() if pd.notna(s) else ""
def bino(s):
    p=s.split()
    return len(p)==2 and all(re.fullmatch(r"[a-z][a-z-]+",x) for x in p)

def traits():
    r=requests.get(EIV_URL,timeout=(15,120),headers={
        "User-Agent":"gesui-H7-original-EIV-NUTRI-pairwise-source-audit/20261115"})
    r.raise_for_status()
    if hashlib.sha256(r.content).hexdigest()!=EIV_SHA:
        raise RuntimeError("Published EIV source SHA drift")
    d=pd.read_csv(io.BytesIO(r.content),encoding="utf-8-sig")
    if len(d)!=8908 or list(d.columns)!=["Unnamed: 0"]+DIMENSIONS:
        raise RuntimeError("Published EIV six-column original structure changed")
    d["name"]=d["Unnamed: 0"].map(normal)
    d=d.loc[d["name"].map(bino)].copy()
    for col in DIMENSIONS:
        d[col]=pd.to_numeric(d[col],errors="coerce")
        low,high=RANGES[col]
        d.loc[~d[col].between(low,high),col]=np.nan
    ambiguous=set()
    # Freeze source name dropping for any duplicated conflicting trait,
    # without using any W/HNS ecological response.
    for col in DIMENSIONS:
        stat=d.groupby("name")[col].agg(["min","max"])
        ambiguous|=set(stat.index[stat["max"]-stat["min"]>1e-10])
    d=d.loc[~d.name.isin(ambiguous)].drop_duplicates("name")
    return {row["name"]:{col:float(row[col]) if pd.notna(row[col]) else None for col in DIMENSIONS}
            for _,row in d.iterrows()},len(ambiguous)

def main(args):
    v2=json.loads(one(args.v2,"EU454_FULL_EIV_SOURCE_GATE.json").read_text())
    if (v2.get("schema")!="GESUI_EU454_FULLSOURCE_PLANT_NUTRI_TRAIT_GATE_V2_HYBRID_QA"
       or not v2.get("hard_source_support_pass")
       or v2.get("official_confirmed_valid_unscorable_nonbinomial_hybrids")!=175
       or v2.get("source_hybrid_score_assigned")!=0
       or not v2.get("original_EIV_NUTRI_scored_keys_unchanged")):
        raise RuntimeError("Original source-v2 numerical-NUTRI quality hard gate required")
    cross=pd.read_csv(one(args.v2,"EU454_PLANT_11511_SOURCE_EIV_EXACT_CROSSWALK.csv"),
                      dtype={"speciesKey":str})
    if len(cross)!=SOURCE_PLANT_UNION or cross.speciesKey.duplicated().any():
        raise RuntimeError("Original 11511-source taxonomic keys changed")
    scored=cross.loc[cross.EIV_NUTRI_exact.notna()].copy()
    if len(scored)!=SOURCE_N_MATCHED:raise RuntimeError("Original unchanged 4011 N exact scores changed")
    published,nbad=traits()
    paired={control:{} for control in RIVALS}
    checked=0
    for row in scored.itertuples(index=False):
        name=normal(row.gbif_canonical_name)
        if not bino(name):raise RuntimeError("Unambiguous original GBIF accepted species key lost")
        trait=published.get(name)
        if trait is None or trait.get("NUTRI") is None:
            raise RuntimeError("Previously original exact score disappeared from EIV source "+str(row.speciesKey))
        n=float(row.EIV_NUTRI_exact)
        if abs(n-trait["NUTRI"])>1e-10:
            raise RuntimeError("Original source score disagreement for GBIF speciesKey "+str(row.speciesKey))
        checked+=1
        for c in RIVALS:
            if trait[c] is not None:
                paired[c][str(row.speciesKey)]=(n,trait[c])
    if checked!=SOURCE_N_MATCHED:raise RuntimeError("NUTRI species crosswalk check failed")
    panel_path=Path("results/20261115_sync/hns_macro/H7_454_HNS_panel.csv")
    if sha(panel_path)!=H7_SHA:raise RuntimeError("Original H7 cohort checksum changed")
    panel=pd.read_csv(panel_path,dtype={"h":str,"h6":str})
    if len(panel)!=454 or panel.h.nunique()!=232 or panel.h6.nunique()!=90:
        raise RuntimeError("Original H7/H6 source grain changed")
    elig=pd.read_csv(one(args.v2,"EU454_NUTRI_SOURCE_ELIGIBILITY_PER_H7.csv"),dtype={"h":str,"h6":str})
    if len(elig)!=454 or elig.duplicated(["h","period"]).any():
        raise RuntimeError("Original source-v2 prior trait cohort drift")
    elig_keys=set(zip(elig.loc[elig.full_source_trait_eligible,"h"],
                      elig.loc[elig.full_source_trait_eligible,"period"]))
    if len(elig_keys)!=451:raise RuntimeError("Original prior 451 approved H7-period cohort differs")
    path=one(args.facets,"GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json")
    if sha(path)!=FACET_SHA:raise RuntimeError("Original 2808 source facets SHA mismatch")
    src=json.loads(path.read_text())
    if len(src)!=2808 or any(x.get("source_status")!="ok" for x in src):
        raise RuntimeError("Original successful 2808 facet data changed")
    flora={(str(x["HYBAS_ID_L07"]),x["period"]):x
           for x in src if x["group"]=="Plantae"}
    if len(flora)!=468:raise RuntimeError("Original 468 plant facets absent")
    rows=[]
    for r in panel.itertuples(index=False):
        key=(r.h,r.period)
        f=flora[key]
        counts={str(k):int(v) for k,v in f["source_speciesKey_counts"].items() if int(v)>0}
        base={"h":r.h,"h6":r.h6,"period":r.period,"country_holdout":r.country_holdout,
              "all_plantae_scored_speciesKeys":len(counts),
              "all_original_plantae_records":int(f["source_query_occurrence_count"]),
              "original_v2_eligible_NUTRI":key in elig_keys}
        for c in RIVALS:
            keep=[(paired[c][k][0],paired[c][k][1],v)
                  for k,v in counts.items() if k in paired[c]]
            n=len(keep)
            cover=n/len(counts) if counts else 0.
            base["n_pair_"+c]=n
            base["cov_pair_"+c]=cover
            base["source_pair_eligible_"+c]=(key in elig_keys and n>=MIN_N
                 and cover>=MIN_FRAC and int(f["source_query_occurrence_count"])>=50)
            base["pair_N_"+c]=float(np.mean([z[0] for z in keep])) if keep else np.nan
            base["pair_"+c]=float(np.mean([z[1] for z in keep])) if keep else np.nan
        rows.append(base)
    frame=pd.DataFrame(rows)
    if len(frame)!=454 or frame.duplicated(["h","period"]).any():
        raise RuntimeError("Original 454 rows missing/duplicated")
    outputs={}
    for c in RIVALS:
        active=frame.loc[frame["source_pair_eligible_"+c]]
        outputs[c]={
            "scored_species_keys_original_11511_union":len(paired[c]),
            "paired_support_H7period":len(active),
            "paired_support_true_H6":int(active.h6.nunique()),
            "paired_support_country_groups":int(active.country_holdout.nunique()),
            "paired_coverage_q05_q50_q95":np.quantile(frame["cov_pair_"+c],[.05,.5,.95]).tolist(),
            "support_pass":bool(len(active)>=MIN_H7PERIOD and active.h6.nunique()>=MIN_H6
                                and active.country_holdout.nunique()>=MIN_COUNTRY)}
    summary={
        "schema":"GESUI_EIV_5_NUTRI_PAIRWISE_SAME_SCORED_SPECIES_SUPPORT_V1",
        "analysis_stage":"source_only_not_yet_ecological_regression",
        "stage_after_six_trait_stop_is_posthoc_secondary":True,
        "original_454_panel_SHA":H7_SHA,
        "original_2808_facets_SHA":FACET_SHA,
        "original_EIV_source_SHA":EIV_SHA,
        "original_NUTRI_scored_species_4011_verified":checked,
        "ambiguous_source_names_excluded":nbad,
        "pairwise_source":outputs,
        "original_six_alltrait_source_gate_STOP":True,
        "each_pair_on_exact_same_plants_for_both_N_and_other":True,
        "all_five_controls_prespecified":RIVALS,
        "joint_five_other_traits_not_tested":True,
        "not_actual_river_or_plant_nutrient_measurement":True}
    args.output.mkdir(parents=True,exist_ok=True)
    frame.to_csv(args.output/"EU454_EIV_NUTRI_FIVE_PAIRS_SOURCE.csv",index=False)
    (args.output/"EU454_EIV_NUTRI_FIVE_PAIR_SOURCE_GATE.json").write_text(
        json.dumps(summary,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_EIV_FIVE_PAIR_SOURCE_GATE="+json.dumps(summary,sort_keys=True,ensure_ascii=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--v2",required=True,type=Path)
    p.add_argument("--facets",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
