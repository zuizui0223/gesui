#!/usr/bin/env python3
"""US EPA NWCA2021 actual vegetation/water source availability BEFORE values.

Only reads original SITE_ID, visit and metadata from four official EPA CSVs.
Does not inspect H_ALL Shannon, FQAI, total NTL/total PTL measurements,
plant species or water concentrations until original source gate is frozen.
"""
from __future__ import annotations
import argparse,hashlib,json,math,re
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd

FILES={
 "SITEINFO":"nwca21_siteinfo-data.csv",
 "WATERCHEM":"nwca21_waterchem_wide-data.csv",
 "VEGMETRICS":"nwca21_vegmetrics-data.csv",
 "HYDROSOURCE":"nwca21_hydrology_sources_wide-data.csv"
}
SELECTED_SITE_FIELDS=("SITE_ID","UID","UNIQUE_ID","VISIT_NO","SITETYPE","XSTATUS","WGT_TP_CORE","YEAR",
 "STUDY","LAT_ANALYS","LON_ANALYS","PSTL_CODE","AG_ECO9","WETCLS_GRP","HYD_CLS")
SITE_GATE_MIN_N=300
STATE_GATE_MIN=15
ECO_GATE_MIN=4
def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()
def canon(x):return re.sub(r"[^A-Z0-9]","",str(x).upper())
def source(csv):
    # Original table header is read separately, then ONLY source IDs/metadata;
    # neither water chemistry nor biological response numerical columns.
    head=pd.read_csv(csv,nrows=0,low_memory=False)
    nms=list(head.columns)
    lookup={canon(x):x for x in nms}
    keys={z:lookup[canon(z)] for z in ("SITE_ID","VISIT_NO","UID") if canon(z) in lookup}
    if "SITE_ID" not in keys:raise ValueError("Official EPA table lacks SITE_ID key "+str(csv))
    allowed=SELECTED_SITE_FIELDS if "siteinfo" in str(csv).lower() else ("SITE_ID","VISIT_NO","UID")
    selected=[lookup[canon(z)] for z in allowed if canon(z) in lookup]
    x=pd.read_csv(csv,usecols=selected,low_memory=False,dtype={lookup["SITEID"]:"string"})
    x=x.rename(columns={col:canon(col) for col in x.columns})
    for k in ("SITEID","UID"):
        if k in x:x[k]=x[k].astype("string").str.strip()
    if "VISITNO" in x:
        x["VISITNO"]=pd.to_numeric(x["VISITNO"],errors="coerce")
    else:x["VISITNO"]=1
    return x,nms

def run(inputs):
    source_records={}
    matrices={}
    for key,p in inputs.items():
        x,headers=source(p)
        original_site_ids=set(x.SITEID.dropna().astype(str))
        primary=(x.SITEID.astype("string")+"|"+x.VISITNO.astype("Int64").astype("string"))
        if primary.duplicated().any():
            # Repeat source species/metric samples can be legitimate,
            # but join uses explicit SITE_ID + VISIT_NO after dedup.
            unique=primary.dropna().nunique()
        else:unique=len(primary)
        source_records[key]={
            "official_source_sha256":sha(p),
            "original_source_csv_bytes":p.stat().st_size,
            "original_source_file_rows":len(x),
            "original_distinct_SITE_IDs":len(original_site_ids),
            "distinct_SITE_ID_by_VISIT_NO_combinations":int(unique),
            "source_ID_fields":list(k for k in ("SITEID","VISITNO","UID") if k in x.columns),
            "source_response_FIELD_NAMES_ONLY_no_values_read":headers
        }
        matrices[key]=x
    site=matrices["SITEINFO"].copy()
    if "SITETYPE" not in site or "WGT TPCORE" in site:
        raise ValueError("Original US EPA site probability source schema invalid")
    def numeric(name):
        return pd.to_numeric(site[name],errors="coerce") if name in site else pd.Series([np.nan]*len(site),index=site.index)
    weights=numeric("WGT TPCORE") if "WGT TPCORE" in site else numeric("WGTTPCORE")
    lat=numeric("LATANALYS")
    lon=numeric("LONANALYS")
    year=numeric("YEAR")
    selected=(site.SITETYPE.astype("string").str.upper().eq("PROB").fillna(False)&
              weights.gt(0).fillna(False)&lat.between(24,50).fillna(False)&
              lon.between(-126,-66).fillna(False)&year.eq(2021).fillna(False)&
              site.VISITNO.eq(1).fillna(False))
    if "STUDY" in site:
        source_study=site.STUDY.astype("string")
        selected&=source_study.str.contains(r"NWCA",case=False,na=False)
    site["QUALIFIED_SOURCE_ONLY"]=selected
    eligible=site[site.QUALIFIED_SOURCE_ONLY].copy()
    keys=["SITEID","VISITNO"]
    if eligible.duplicated(keys).any():
        raise ValueError("Original EPA first-visit probability site rows repeated")
    other={}
    for k in ("WATERCHEM","VEGMETRICS","HYDROSOURCE"):
        other[k]=matrices[k][keys].dropna(subset=["SITEID"]).drop_duplicates()
    eligible_chem=eligible.merge(other["WATERCHEM"],on=keys,how="inner",validate="one_to_one")
    eligible_plants=eligible.merge(other["VEGMETRICS"],on=keys,how="inner",validate="one_to_one")
    eligible_both=eligible_chem.merge(other["VEGMETRICS"],on=keys,how="inner",validate="one_to_one")
    eligible_all=eligible_both.merge(other["HYDROSOURCE"],on=keys,how="inner",validate="one_to_one")
    states=Counter(eligible_both["PSTLCODE"].astype(str)) if "PSTLCODE" in eligible_both else {}
    eco=Counter(eligible_both["AGECO9"].astype(str)) if "AGECO9" in eligible_both else {}
    wetlands=Counter(eligible_both["WETCLSGRP"].astype(str)) if "WETCLSGRP" in eligible_both else {}
    hydro=Counter(eligible_both["HYDCLS"].astype(str)) if "HYDCLS" in eligible_both else {}
    gate=bool(len(eligible_both)>=SITE_GATE_MIN_N and len(states)>=STATE_GATE_MIN and len(eco)>=ECO_GATE_MIN)
    output={
        "schema":"gesui_US_EPA_NWCA2021_original_prob_wetland_plant_shannon_waterchem_source_only_join_GATE_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "official_EPA_file_urls":{k:"https://www.epa.gov/system/files/other-files/2024-05/"+v for k,v in FILES.items()},
        "official_EPA_original_source_csv_receipts":source_records,
        "source2021_site_qualifying_probability_weights_2021_visit1_latlon_n":len(eligible),
        "source_SITE_ID_and_VISIT_NO_matched_2021_chemistry_rows_n":len(eligible_chem),
        "source_SITE_ID_and_VISIT_NO_matched_2021_vegmetrics_rows_n":len(eligible_plants),
        "source_SITE_ID_and_VISIT_NO_matched_chemistry_vegmetrics_rows_n":len(eligible_both),
        "source_SITE_ID_and_VISIT_NO_matched_chemistry_vegmetrics_hydrology_rows_n":len(eligible_all),
        "source_matching_by_original_SITE_ID_VISIT_NO_not_by_response_values":True,
        "source_joint_plant_chemistry_states":dict(states),
        "source_joint_plant_chemistry_ecoregions":dict(eco),
        "source_joint_plant_chemistry_wetland_classes":dict(wetlands),
        "source_joint_plant_chemistry_hydrologic_classes":dict(hydro),
        "source_joined_original_SITE_ID_SHA256_not_plant_or_chemistry_values":[
            hashlib.sha256((str(r["SITEID"])+"|"+str(int(r["VISITNO"]))).encode()).hexdigest()
            for _,r in eligible_both.sort_values(keys).iterrows()],
        "precommitted_SOURCE_gate_before_measured_NTLP_and_plant_Shannon":{
            "min_original_joined_probability_NWCA2021_wetland_sites":SITE_GATE_MIN_N,
            "min_original_US_states":STATE_GATE_MIN,
            "min_original_NWCA_ecoregions":ECO_GATE_MIN,
            "may_open_2021_plant_diversity_and_measured_TNTP_results":gate},
        "all_original_2021_plant_H_ALL_FQAI_and_water_NTL_PTL_values_NOT_READ":True,
        "site2021_measured_TNTP_NOT_automatically_wastewater_derived":True,
        "NWCA_wetland_probability_site_sample_NOT_all_plant_ecosystems":True,
        "no_causal_sewage_biodiversity_claims":True
    }
    return output

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-dir",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    orig={key:a.source_dir/name for key,name in FILES.items()}
    out=run(orig)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({
        "source_probability_sites":out["source2021_site_qualifying_probability_weights_2021_visit1_latlon_n"],
        "source_both":out["source_SITE_ID_and_VISIT_NO_matched_chemistry_vegmetrics_rows_n"],
        "with_hydrology":out["source_SITE_ID_and_VISIT_NO_matched_chemistry_vegmetrics_hydrology_rows_n"],
        "distinct_states":len(out["source_joint_plant_chemistry_states"]),
        "distinct_ecoregions":len(out["source_joint_plant_chemistry_ecoregions"]),
        "gate":out["precommitted_SOURCE_gate_before_measured_NTLP_and_plant_Shannon"]
    },indent=2))
if __name__=="__main__":
    main()
