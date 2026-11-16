#!/usr/bin/env python3
"""Source-first, no ecological model: exact six-EIV species and EU454 common-score support.

Origin: same 2023 European EIV SHA, original H7 454-source GBIF source,
original 11,511 accepted-species-key crosswalk from PR102 v2.
Do not impute any missing plant score or source nothospecies/hybrid.
All six CWM indices use the identical site-by-species scored set.
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

FACET_SHA="251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4"
PANEL_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
EIV_SHA="ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba"
URL="https://www.sci.muni.cz/botany/juice/NET/sources/EIV_Europe.csv"
TRAITS=["LIGHT","TEMPE","MOIST","REACT","NUTRI","SALIN"]
BOUNDS={"LIGHT":(1,9),"TEMPE":(1,12),"MOIST":(1,12),"REACT":(1,9),"NUTRI":(1,9),"SALIN":(0,9)}
EXPECTED_NUTRI_MATCH=4011
MIN_SCORE=10
MIN_COVERAGE=.40
MIN_H7P=100
MIN_H6=40
MIN_COUNTRY=5

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def one(folder,name):
    hits=list(folder.rglob(name))
    if len(hits)!=1:raise RuntimeError("Exactly one original source file required: "+name)
    return hits[0]
def simplify(x):
    if not isinstance(x,str):return ""
    return re.sub(r"\s+"," ",x.strip()).casefold()
def exact_binomial(x):
    tokens=x.split()
    return len(tokens)==2 and all(re.fullmatch(r"[a-z][a-z-]+",t) for t in tokens)

def original_eiv():
    response=requests.get(URL,timeout=(15,120),headers={"User-Agent":"gesui-source-locked-six-trait-controls/20261115"})
    response.raise_for_status()
    raw=response.content
    digest=hashlib.sha256(raw).hexdigest()
    if digest!=EIV_SHA:raise RuntimeError("Source EIV checksum drift "+digest)
    d=pd.read_csv(io.BytesIO(raw),encoding="utf-8-sig")
    if list(d.columns)!=["Unnamed: 0"]+TRAITS or len(d)!=8908:
        raise RuntimeError("Source-only verified original EIV headers, shape drift")
    d["name_norm"]=d["Unnamed: 0"].map(simplify)
    d=d[d.name_norm.map(exact_binomial)].copy()
    for t in TRAITS:
        d[t]=pd.to_numeric(d[t],errors="coerce")
        low,high=BOUNDS[t]
        d.loc[~d[t].between(low,high),t]=np.nan
    # Numeric scores cannot be borrowed from a synonym / near-matching spelling.
    # When source duplicates disagree for any one of six, mark that exact
    # source name unusable in the primary complete-case six-score dictionary.
    ambiguous=set()
    for t in TRAITS:
        grp=d.groupby("name_norm")[t].agg(["min","max"])
        ambiguous.update(grp.index[(grp["max"]-grp["min"])>1e-10])
    d=d[~d.name_norm.isin(ambiguous)].drop_duplicates("name_norm")
    six=d.dropna(subset=TRAITS).copy()
    scores={r.name_norm:{t:float(getattr(r,t)) for t in TRAITS}
            for r in six.itertuples(index=False)}
    return scores,{
      "original_sha256":digest,
      "original_header":["Unnamed: 0"]+TRAITS,
      "original_rows":8908,
      "unambiguous_complete_six_trait_names":len(scores),
      "names_excluded_for_any_ambiguous_duplicate":len(ambiguous),
      "original_species_pairing_is_only_exact_accepted_vascular_binomial":True}

def main(args):
    v2=json.loads(one(args.v2,"EU454_FULL_EIV_SOURCE_GATE.json").read_text())
    if (v2.get("schema")!="GESUI_EU454_FULLSOURCE_PLANT_NUTRI_TRAIT_GATE_V2_HYBRID_QA"
        or not v2.get("hard_source_support_pass")
        or not v2.get("original_EIV_NUTRI_scored_keys_unchanged")
        or v2.get("official_confirmed_valid_unscorable_nonbinomial_hybrids")!=175
        or v2.get("source_hybrid_score_assigned")!=0):
        raise RuntimeError("Original true 11511 EIV source-v2 175-hybrid QA required")
    cross=pd.read_csv(one(args.v2,"EU454_PLANT_11511_SOURCE_EIV_EXACT_CROSSWALK.csv"),
                      dtype={"speciesKey":str})
    if len(cross)!=11511 or cross.speciesKey.duplicated().any():
        raise RuntimeError("Original 11,511 source speciesKey crosswalk drift")
    matched=cross.EIV_NUTRI_exact.notna()
    if int(matched.sum())!=EXPECTED_NUTRI_MATCH:
        raise RuntimeError("Original 4,011 exact NUTRI matches count drift")
    six_dict,eiv_info=original_eiv()
    scores={}
    for r in cross.itertuples(index=False):
        key=str(r.speciesKey)
        # Exact numeric NUTRI matches from original v2 are the only possible
        # recipients of a complete six-trait score.
        if not pd.notna(r.EIV_NUTRI_exact):
            continue
        name=simplify(r.gbif_canonical_name)
        if not exact_binomial(name):raise RuntimeError("Prior accepted exact N name invalid "+key)
        six=six_dict.get(name)
        if six is None:continue
        if abs(six["NUTRI"]-float(r.EIV_NUTRI_exact))>1e-10:
            raise RuntimeError("Source original 4,011 NUTRI changed when six traits matched "+key)
        scores[key]=six
    facets_path=one(args.facets,"GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json")
    if sha(facets_path)!=FACET_SHA:raise RuntimeError("Original 2808 source facets SHA mismatch")
    original=json.loads(facets_path.read_text())
    if len(original)!=2808 or any(x.get("source_status")!="ok" for x in original):
        raise RuntimeError("Not original successful 2808 source facet rows")
    flora={(str(r["HYBAS_ID_L07"]),r["period"]):r
           for r in original if r["group"]=="Plantae"}
    if len(flora)!=468:raise RuntimeError("Not original 468 Plantae source facets")
    p=Path("results/20261115_sync/hns_macro/H7_454_HNS_panel.csv")
    if sha(p)!=PANEL_SHA:raise RuntimeError("Original 454 H7 outcome and source panel SHA mismatch")
    panel=pd.read_csv(p,dtype={"h":str,"h6":str})
    if len(panel)!=454 or panel.h.nunique()!=232 or panel.h6.nunique()!=90:
        raise RuntimeError("Original geographic units changed")
    elig=pd.read_csv(one(args.v2,"EU454_NUTRI_SOURCE_ELIGIBILITY_PER_H7.csv"),dtype={"h":str,"h6":str})
    if len(elig)!=454:raise RuntimeError("Original source eligibility table changed")
    original_selected={(r.h,r.period) for r in panel.itertuples(index=False)}
    if len(original_selected)!=454:raise RuntimeError("Duplicated H7-period")
    records=[]
    for r in panel.itertuples(index=False):
        facet=flora[(r.h,r.period)]
        counts={str(k):int(v) for k,v in facet["source_speciesKey_counts"].items() if int(v)>0}
        n=len(counts)
        numeric=[(scores[k],v) for k,v in counts.items() if k in scores]
        cov=len(numeric)/n if n else 0.
        row={"h":r.h,"period":r.period,"h6":r.h6,"country_holdout":r.country_holdout,
             "all_original_Plantae_speciesKey_count":n,
             "source_raw_Plantae_occurrence_count":int(facet["source_query_occurrence_count"]),
             "six_EIV_common_scored_species":len(numeric),
             "six_EIV_common_species_coverage":cov,
             "six_trait_quality_eligible":bool(len(numeric)>=MIN_SCORE and cov>=MIN_COVERAGE
                and int(facet["source_query_occurrence_count"])>=50)}
        for t in TRAITS:
            row["CWM6_"+t]=float(np.mean([s[t] for s,_ in numeric])) if numeric else np.nan
            row["RECORD6_"+t]=float(np.average([s[t] for s,_ in numeric],weights=[v for _,v in numeric])) if numeric else np.nan
        records.append(row)
    f=pd.DataFrame(records)
    if len(f)!=454 or f.duplicated(["h","period"]).any():
        raise RuntimeError("Unchanged original 454 site grain not reproduced")
    support=f[f.six_trait_quality_eligible].copy()
    nH6=int(support.h6.nunique());nct=int(support.country_holdout.nunique())
    passing=len(support)>=MIN_H7P and nH6>=MIN_H6 and nct>=MIN_COUNTRY
    v2keys=set(cross.loc[matched,"speciesKey"])
    summary={
      "schema":"GESUI_EU454_EIV_SIXTRAIT_SAME_SPECIES_SOURCE_OPPORTUNITY_V1",
      "source_EIV":eiv_info,
      "original_H7panel_SHA256":PANEL_SHA,
      "original_2808_facets_SHA256":FACET_SHA,
      "original_EIV_trait_exact_NUTRI_species":len(v2keys),
      "original_six_complete_scored_speciesKeys":len(scores),
      "original_454_H7periods":454,
      "original_451_NUTRI_supported_H7periods":int(v2["eligible_H7_period"]),
      "all_six_source_supported_H7periods":len(support),
      "all_six_source_supported_H6":nH6,
      "all_six_source_supported_country_groups":nct,
      "all_six_coverage_q05_q50_q95":np.quantile(f.six_EIV_common_species_coverage,[.05,.5,.95]).tolist(),
      "original_same_univariate_scores_across_traits":True,
      "no_missing_trait_substitution":True,
      "all_six_preoutcome_gate_pass":bool(passing),
      "support_conditions":{"min_unique_species":MIN_SCORE,"min_species_fraction":MIN_COVERAGE,
                            "min_H7period":MIN_H7P,"min_H6":MIN_H6,"min_country":MIN_COUNTRY},
      "no_W_or_HNS_model_opened":True,
      "score_source_compatible_with_original_4011_NUTRI":True,
      "no_caused_sewage_nutrient_inference":True}
    args.output.mkdir(parents=True,exist_ok=True)
    f.to_csv(args.output/"EU454_EIV_SIXTRAIT_COMMON_SOURCE_SUPPORT.csv",index=False)
    (args.output/"EU454_EIV_SIXTRAIT_COMMON_SOURCE_GATE.json").write_text(
        json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("GESUI_EIV_SIXTRAIT_SOURCE_GATE="+json.dumps(summary,ensure_ascii=False,sort_keys=True),flush=True)
    print("MODEL_GATE="+("PASS_SAME_SCORED_SPECIES" if passing else "STOP_SIX_TRAIT_COMMON_SOURCE_SUPPORT"),flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--v2",type=Path,required=True)
    ap.add_argument("--facets",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    main(ap.parse_args())
