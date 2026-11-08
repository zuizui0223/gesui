#!/usr/bin/env python3
"""Exploratory global WWTP/outfall Plantae GBIF-index pilot, NOT a sewage-N effect.

Source-blind-to-outcome, fixed one high/low effluent-fraction pair per country.
Two nonoverlapping sites per WWTP: modelled outfall versus physical plant.
Only GBIF-indexed Plantae species and effort-standardized records are observed.
"""
from __future__ import annotations
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request,urlopen
from urllib.error import URLError,HTTPError
import argparse
import hashlib
import json
import math
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SOURCE_BLOB="0506f873f6541c5f324478894a9731290d722bce"
COUNTRIES=("AUT","BEL","CHE","CZE","DEU","DNK","FRA","GBR","NLD","POL","SWE")
EXCLUDED={"closed","decommissioned","non-operational","projected",
          "proposed","under construction"}
GBIF="https://api.gbif.org/v1/occurrence/search"
AGENT="gesui-plantae-wwtp-outfall-recording-pilot/0.1"
R=6371.0088
YEAR="2011,2024"
RADIUS_KM=2
MIN_CENTER_DISTANCE_KM=4.25
MIN_GBIF_SPECIES_RECORDS=30
RAREFACTION_RECORDS=30

def sha_blob(p):
    dat=Path(p).read_bytes()
    return hashlib.sha1(b"blob "+str(len(dat)).encode()+b"\0"+dat).hexdigest()

def haversine(la1,lo1,la2,lo2):
    a,b=math.radians(float(la1)),math.radians(float(la2))
    dlat=b-a
    dlon=math.radians(float(lo2)-float(lo1))
    q=math.sin(dlat/2)**2+math.cos(a)*math.cos(b)*math.sin(dlon/2)**2
    return 2*R*math.asin(math.sqrt(min(1,max(0,q))))

def water_fraction(ww_m3day,river_m3sec):
    if ww_m3day<=0 or river_m3sec<=0:return float("nan")
    return float(ww_m3day/(ww_m3day+86400*river_m3sec))

def select_before_gbif(df):
    """Pair at most one high and one low within each country, no GBIF used."""
    eligible=df[
        df.CNTRY_ISO.isin(COUNTRIES)
        & (~df.STATUS.astype(str).str.lower().isin(EXCLUDED))
        & (df.QUAL_LOC<=2)&(df.QUAL_WASTE<=3)
        & (df.WASTE_DIS>0)&(df.RIVER_DIS>=0.5)
        & (df.COAST_10KM==0)&df.HYRIV_ID.notna()
        & df.LAT_WWTP.between(-65,85)&df.LON_WWTP.between(-180,180)
        & df.LAT_OUT.between(-65,85)&df.LON_OUT.between(-180,180)
    ].copy()
    eligible["center_sep_km"]=eligible.apply(lambda x:haversine(
        x.LAT_OUT,x.LON_OUT,x.LAT_WWTP,x.LON_WWTP),axis=1)
    eligible=eligible[(eligible.center_sep_km>=MIN_CENTER_DISTANCE_KM)
                      &(eligible.center_sep_km<=35)].copy()
    eligible["water_fraction"]=eligible.apply(lambda x:water_fraction(
        float(x.WASTE_DIS),float(x.RIVER_DIS)),axis=1)
    eligible=eligible[np.isfinite(eligible.water_fraction)]
    pairs=[]
    for country in COUNTRIES:
        sub=eligible.loc[eligible.CNTRY_ISO==country].copy().sort_values(
            ["water_fraction","WASTE_ID"])
        if len(sub)<6:continue
        q=max(2,len(sub)//4)
        lo=sub.iloc[:q]
        hi=sub.iloc[-q:]
        opts=[]
        for h in hi.itertuples(index=False):
            for l in lo.itertuples(index=False):
                dist=haversine(h.LAT_OUT,h.LON_OUT,l.LAT_OUT,l.LON_OUT)
                if not (30<=dist<=350):continue
                ratio=max(h.RIVER_DIS,l.RIVER_DIS)/min(h.RIVER_DIS,l.RIVER_DIS)
                if ratio>6:continue
                score=abs(math.log(h.RIVER_DIS/l.RIVER_DIS)) + dist/700
                opts.append((score,int(h.WASTE_ID),int(l.WASTE_ID),h,l,dist))
        if not opts:continue
        _,_,_,h,l,d=sorted(opts,key=lambda z:z[:3])[0]
        pairs.append({"pair_id":country+"-01","country":country,
                      "pair_distance_km":d,"high":h,"low":l})
    summary={"eligible_outfall_wwtp":len(eligible),
             "countries_with_pairs":len(pairs),"matched_pairs":len(pairs),
             "selected_wwtps":len(pairs)*2,
             "buffer_km":RADIUS_KM,
             "source_selection_used_no_gbif":True}
    return pairs,summary

def gbif_query(params, retries=4):
    url=GBIF+"?"+urlencode(params)
    error=None
    for retry in range(retries):
        try:
            req=Request(url,headers={"User-Agent":AGENT,"Accept":"application/json"})
            with urlopen(req,timeout=65) as resp:
                body=resp.read()
            parsed=json.loads(body)
            if "count" not in parsed:raise ValueError("GBIF response lacks count")
            return parsed
        except (HTTPError,URLError,TimeoutError,ValueError) as e:
            error=e
            if retry+1<retries:time.sleep(4*(retry+1))
    raise RuntimeError(f"GBIF request failed after {retries} attempts: {error}")

def facet_site(lat,lon,role,wwtp_id):
    params={
       "taxonKey":6,"year":YEAR,
       "geoDistance":f"{lat:.7f},{lon:.7f},{RADIUS_KM}km",
       "hasCoordinate":"true","hasGeospatialIssue":"false",
       "occurrenceStatus":"PRESENT",
       "facet":"speciesKey","facetLimit":1000,"limit":0
    }
    keys={}
    total=None
    pages=0
    while True:
        params["facetOffset"]=pages*1000
        data=gbif_query(params)
        if total is None:total=int(data["count"])
        if int(data["count"])!=total:
            raise ValueError("GBIF count drift during facet pagination")
        facets=next((x.get("counts",[]) for x in data.get("facets",[])
                    if str(x.get("field","")).upper() in ("SPECIES_KEY","SPECIESKEY")),[])
        if total and not facets:
            raise ValueError("GBIF species facet absent despite records")
        for obj in facets:
            key=int(obj["name"])
            if key in keys:raise ValueError("Duplicate species facet page")
            keys[key]=int(obj["count"])
        pages+=1
        if len(facets)<1000:break
        if pages>=10:raise ValueError("Facet pagination above 10,000 species per 2km")
        time.sleep(.25)
    n=sum(keys.values())
    if n>total:raise ValueError("Faceted species counts exceed overall plant records")
    return {"role":role,"wwtp_id":wwtp_id,"lat":float(lat),"lon":float(lon),
            "total_plant_records":total,"species_record_count":n,
            "distinct_indexed_species":len(keys),"species_facets":keys,
            "facet_pages":pages}

def rarefy_expected(counts, depth=RAREFACTION_RECORDS):
    arr=[int(x) for x in counts if int(x)>0]
    n=sum(arr)
    if n<depth:return None
    def log_comb(total,k):
        return math.lgamma(total+1)-math.lgamma(k+1)-math.lgamma(total-k+1)
    den=log_comb(n,depth)
    return sum(1-math.exp(log_comb(n-k,depth)-den)
               if n-k>=depth else 1.0 for k in arr)

def summarize_facets(rows):
    by_site={}
    for row in rows:
        copy=dict(row)
        copy["expected_richness_rarefied_30"]=rarefy_expected(
            row["species_facets"].values())
        by_site[(row["wwtp_id"],row["role"])]=copy
    return by_site

def model(pairs,sites):
    detail=[]
    pair_deltas=[]
    for p in pairs:
        cells={}
        for stratum in ("high","low"):
            w=p[stratum]
            w_id=int(w.WASTE_ID)
            o=sites[(w_id,"outfall")]
            t=sites[(w_id,"wwtp")]
            qualified=(o["species_record_count"]>=MIN_GBIF_SPECIES_RECORDS and
                       t["species_record_count"]>=MIN_GBIF_SPECIES_RECORDS and
                       o["distinct_indexed_species"]>=6 and t["distinct_indexed_species"]>=6)
            delta=(math.log1p(o["expected_richness_rarefied_30"])
                   -math.log1p(t["expected_richness_rarefied_30"])) if qualified else None
            cells[stratum]=delta
            detail.append({
                "pair_id":p["pair_id"],"country":p["country"],
                "exposure_stratum":stratum,"WASTE_ID":w_id,
                "water_fraction_proxy":float(w.water_fraction),
                "wwtp_outfall_center_km":float(w.center_sep_km),
                "eligible_both_site_records":qualified,
                "plant_records_outfall":o["total_plant_records"],
                "plant_records_wwtp":t["total_plant_records"],
                "index_species_outfall":o["distinct_indexed_species"],
                "index_species_wwtp":t["distinct_indexed_species"],
                "expected_richness_30_outfall":o["expected_richness_rarefied_30"],
                "expected_richness_30_wwtp":t["expected_richness_rarefied_30"],
                "delta_log_rarefied_outfall_minus_wwtp":delta
            })
        if all(v is not None for v in cells.values()):
            pair_deltas.append({"pair_id":p["pair_id"],"country":p["country"],
                 "high_minus_low_outfall_wwtp_delta":cells["high"]-cells["low"]})
    qual=[x for x in detail if x["eligible_both_site_records"]]
    cor=None
    if len(qual)>=8 and len(set(x["water_fraction_proxy"] for x in qual))>=8:
        r,pval=spearmanr(
            [x["water_fraction_proxy"] for x in qual],
            [x["delta_log_rarefied_outfall_minus_wwtp"] for x in qual])
        cor={"n":len(qual),"spearman_rho":float(r),"two_sided_p_descriptive":float(pval)}
    return {"facility_level":detail,"paired_country_differences":pair_deltas,
            "facility_count_quality":len(qual),"matched_country_pairs_quality":len(pair_deltas),
            "spearman_of_fraction_vs_outfall_facility_delta":cor,
            "mean_high_minus_low_if_any":float(np.mean([x["high_minus_low_outfall_wwtp_delta"]
              for x in pair_deltas])) if pair_deltas else None}

def run(source,output):
    if sha_blob(source)!=SOURCE_BLOB:raise ValueError("Original HydroWASTE source mismatch")
    df=pd.read_csv(source,encoding="latin1",low_memory=False)
    if len(df)!=58502 or df.WASTE_ID.nunique()!=58502:
        raise ValueError("Original source facility count changed")
    for k in ("QUAL_LOC","QUAL_WASTE","WASTE_DIS","RIVER_DIS","COAST_10KM",
              "HYRIV_ID","LAT_OUT","LON_OUT","LAT_WWTP","LON_WWTP"):
        df[k]=pd.to_numeric(df[k],errors="coerce")
    pairs,design=select_before_gbif(df)
    output.mkdir(parents=True,exist_ok=True)
    (output/"preoutcome_source_design.json").write_text(json.dumps({
        "source_blob":SOURCE_BLOB,"sample_design":design,
        "facilities":[{"country":p["country"],"pair_id":p["pair_id"],
                        "pair_distance_km":p["pair_distance_km"],
                        "strata":{q:{"WWTP_ID":int(p[q].WASTE_ID),
                          "nominal_water_fraction":float(p[q].water_fraction)}
                                  for q in ("high","low")}}
                      for p in pairs]},indent=2)+"\n")
    print("GESUI_OUTFALL_PLANT_PREOUTCOME_SELECTION",
          json.dumps(design,sort_keys=True),flush=True)
    if len(pairs)<4:
        raise ValueError("Fewer than 4 source-matched country pairs; no GBIF queries")
    rows=[]
    for k,p in enumerate(pairs):
        for group in ("high","low"):
            w=p[group]
            for role,lat,lon in [
                ("outfall",w.LAT_OUT,w.LON_OUT),
                ("wwtp",w.LAT_WWTP,w.LON_WWTP)]:
                result=facet_site(float(lat),float(lon),role,int(w.WASTE_ID))
                rows.append(result)
                print("GESUI_GBIF_PLANT_INDEX_SITE",json.dumps({
                    "pair":p["pair_id"],"group":group,"role":role,
                    "records":result["total_plant_records"],
                    "species":result["distinct_indexed_species"]}),flush=True)
                time.sleep(.6)
    output_full={"source":"GBIF occurrence search live 2011-2024 Plantae speciesKey facets",
                 "retrieved_utc":datetime.now(timezone.utc).isoformat(),
                 "record_window":YEAR,"facet_radius_km":RADIUS_KM,
                 "source_design":design,
                 "sites":rows,"data_not_host_subset":True}
    (output/"raw_gbif_facets.json").write_text(
        json.dumps(output_full,ensure_ascii=False,separators=(",",":"))+"\n")
    result=model(pairs,summarize_facets(rows))
    result.update({"status":"EXPLORATORY_LIVE_GBIF_PLANTAE_RECORD_INDEX_NOT_WASTEWATER_N_EFFECT",
                   "hydrowaste_git_blob":SOURCE_BLOB,
                   "gbif_window":YEAR,"gbif_sites":len(rows),
                   "source_design":design,
                   "cautions":["outfall is modeled, not sampled N exposure",
                      "WWTP coordinate is urban facility, NOT a hydrologically matched upstream control",
                      "Plantae facets include nonvascular plants; not a proven all-vascular-plant sample",
                      "richness is GBIF recorded species and rarefied record diversity, not plant abundance",
                      "not directionally downstream and not a causal effluent effect",
                      "source-only sampling may yield low GBIF coverage",
                      "country sample intentionally selected one high-low pair, not globally representative"]})
    (output/"plant_observational_pilot.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_FIRST_LIVE_PLANT_OBSERVATIONAL_PILOT",
          json.dumps({k:v for k,v in result.items() if k not in
                     ("facility_level","paired_country_differences")},sort_keys=True),flush=True)
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--hydrowaste",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    run(a.hydrowaste,a.out)

if __name__=="__main__":
    main()
