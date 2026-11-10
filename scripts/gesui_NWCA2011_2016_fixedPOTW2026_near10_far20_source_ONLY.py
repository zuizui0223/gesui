#!/usr/bin/env python3
"""EPA NWCA2011/2016 exact source-only original plots × fixed 2026 POTW registry.

Uses already SHA-frozen pre-H_ALL original 2011 n967, 2016 n966 source keys,
the *same* Nov2026 EPA ECHO national POTW discharge-point ZIP SHA, registered
sewage code and 2021 plausible permit-date proxy as previous NWCA2021.
Only source site lat/lon and IDs are read; NO vegetation metric values.
"""
from __future__ import annotations
import argparse,hashlib,json
from datetime import datetime,timezone
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd

from gesui_NWCA2021_568source_EPA_NPDES_registered_sewage_outfall_geoONLY import (
    ndpes_outfall_source,assess_sites
)
SRC=Path("results/us_NWCA2021_plant_TNTP_20261110/NWCA2011_2016_ORIGINAL_SITE_VEG_SOURCE_BEFORE_PLANT_OUTCOMES.json")
ECHO_SHA="f0d1310b1c07085c42527e2ede5cf97ab83d517d50548b7dadb7621a44aa6160"
MIN_NEAR=50;MIN_FAR=60;MIN_STATES=10;MIN_CELLS=20;MIN_PERMITS=15
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(2**20),b""):h.update(b)
    return h.hexdigest()
def site_hash(site,visit):
    return hashlib.sha256((str(site)+"|"+str(int(visit))).encode()).hexdigest()
def parse_year(s):
    return pd.to_datetime(s,format="mixed",errors="coerce").dt.year

def load(year,root,source_dir,manifest):
    origin=manifest["source_years"][year]
    if origin["source_gate_STATUS"]!="SOURCE_SITE_VEG_PASS":
        return None,"SOURCE_VEGETATION_SITE_GATE_FAILED"
    p=source_dir/(year+"_siteinfo.csv")
    v=source_dir/(year+"_vegmetrics.csv")
    if sha(p)!=origin["original_official_siteinfo_SHA256"] or sha(v)!=origin["original_official_vegetation_metrics_SHA256"]:
        raise ValueError("Original independently pinned historical EPA NWCA source file bytes changed")
    col=origin["source_siteinfo_ID_columns_resolved"]
    required=("SITE_ID","VISIT_NO","LAT","LON","SITETYPE","WEIGHT","YEAR","STATE","ECO","WETCLS","HYDCLS")
    if any(col[x] is None for x in required):
        raise ValueError("Previous historical site source schema incomplete")
    d=pd.read_csv(p,usecols=list({col[x] for x in required}),low_memory=False)
    d=d.rename(columns={value:key for key,value in col.items() if value})
    d["SITE_ID"]=d.SITE_ID.astype("string").str.strip()
    d["VISIT_NO"]=pd.to_numeric(d.VISIT_NO,errors="coerce")
    for x in ("LAT","LON","WEIGHT"):
        d[x]=pd.to_numeric(d[x],errors="coerce")
    d["YEAR"]=parse_year(d.YEAR)
    source=(d.SITETYPE.astype("string").str.strip().str.upper().eq("PROB")&
        d.WEIGHT.gt(0)&d.LAT.between(24,50)&d.LON.between(-126,-66)&
        d.YEAR.eq(int(year))&d.VISIT_NO.eq(1)&
        d.STATE.notna()&d.ECO.notna()&d.WETCLS.notna()&d.HYDCLS.notna())
    d=d.loc[source].copy()
    d["site_hash"]=[site_hash(site,visit) for site,visit in zip(d.SITE_ID,d.VISIT_NO)]
    selected=set(origin["source_site_hashes_only"])
    d=d.loc[d.site_hash.isin(selected)].copy()
    if len(d)!=origin["source_vegmetric_SITE_ID_VISIT_NO_join_matched_n"] or set(d.site_hash)!=selected:
        raise ValueError("Original fixed historical source-only 967/966 vegetation-site keys changed")
    d=d.rename(columns={"LAT":"LAT_ANALYS","LON":"LON_ANALYS","STATE":"PSTL_CODE",
             "ECO":"AG_ECO9","WETCLS":"WETCLS_GRP","HYDCLS":"HYD_CLS",
             "site_hash":"source_site_visit_SHA256"})
    return d,None

def audit(year,d,sewage):
    original_receipts,stats=assess_sites(d,sewage,"strict")
    near=[x for x in original_receipts if x["outfall_rows_within_10km"]>0]
    far=[x for x in original_receipts if x["nearest_outfall_geodesic_km_within25km_else_null"] is None
       or x["nearest_outfall_geodesic_km_within25km_else_null"]>=20.]
    near_group={(x["original_state"],x["original_ecoregion"]) for x in near}
    far_group={(x["original_state"],x["original_ecoregion"]) for x in far}
    shared=near_group&far_group
    selected_near=[x for x in near if (x["original_state"],x["original_ecoregion"]) in shared]
    selected_far=[x for x in far if (x["original_state"],x["original_ecoregion"]) in shared]
    states={r["original_state"] for r in selected_near}
    permits={p for x in selected_near for p in x["independent_10km_NPDES_permits_SHA256"]}
    # 2026 registry might not contain historic 2011/2016 permit periods,
    # so report factual source opportunity, NEVER treat missing as true absence.
    pass_=(len(selected_near)>=MIN_NEAR and len(selected_far)>=MIN_FAR and
           len(states)>=MIN_STATES and len(shared)>=MIN_CELLS and len(permits)>=MIN_PERMITS)
    return {
       "original_historical_year":int(year),
       "original_PROB_firstvisit_fullveg_source_n":len(d),
       "strict_2026_ECHO_2021plausible_2021permit_feature_source_geography":stats,
       "original_shared_state_and_ecoregion_near10_far20_source_n":{
           "near10":len(selected_near),"far20":len(selected_far),
           "shared_original_states":len(states),"shared_original_state_ecoregion_strata":len(shared),
           "unique_administrative_strict_NPDES_permits_in_matched_near":len(permits)
       },
       "source_precommitted_independent_year_gate":{
           "min_near10_sites":MIN_NEAR,"min_far20_sites":MIN_FAR,
           "min_states":MIN_STATES,"min_state_eco_strata":MIN_CELLS,
           "min_distinct_registered_outfall_NPDES_permit_ids":MIN_PERMITS,
           "may_open_original_NWCA_H_ALL_and_fit_2011_2016_comparison":bool(pass_)
       },
       "original_historical_source_matched_site_date_weight_and_outfall_receipts":[
         {"source_year":int(year),"site_SHA256":x["NWCA2021_source_site_visit_SHA256"],
          "original_state":x["original_state"],"original_ecoregion":x["original_ecoregion"],
          "original_wetland_group":x["original_wetland_group"],
          "original_hydro_class":x["original_hydro_class"],
          "source_group":("near10" if x in selected_near else "far20"),
          "strict_outfall_nearest_km_within25km":x["nearest_outfall_geodesic_km_within25km_else_null"],
          "number_original_NPDES_registered_permit_ids_within10km":x["unique_NPDES_permit_IDs_10km"],
          "source_permit_IDs_SHA256_not_actual_discharge":x["independent_10km_NPDES_permits_SHA256"]}
         for x in sorted(selected_near+selected_far,key=lambda z:z["NWCA2021_source_site_visit_SHA256"])
       ],
       "all_original_source_site_plant_metrics_NOT_ACCESSED_before_2026_NPDES_geo_classification":True,
       "year2011_2016_discharge_NOT_VERIFIED_from_current_2026_NPDES_register":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-dir",type=Path,required=True)
    p.add_argument("--npdes-zip",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    manifest=json.loads((root/SRC).read_text())
    if sha(a.npdes_zip)!=ECHO_SHA:
        raise ValueError("Original original 2026-11-08 50MB official EPA ECHO source SHA changed; cannot reuse 2021 provenance")
    sewage,meta=ndpes_outfall_source(a.npdes_zip)
    if (len(sewage)!=58231 or int(meta["sewage_strict_both_original_feature_geometry_and_dates_n"])!=10865):
        raise ValueError("Original NPDES source normalized 2021 status changed since 2021 cohort")
    result={"schema":"gesui_US_NWCA_2011_2016_before_HALL_2026_registered_POTW_GEO_2cycle_v01",
        "utc":datetime.now(timezone.utc).isoformat(),
        "original_2011_2016_vegetation_source_only_manifest_SHA256":sha(root/SRC),
        "original_2026_NPDES_official_ECHO_2021permit_PROXY_ZIP_SHA256":ECHO_SHA,
        "same_original_2021_strict_sewage_source_filter_reused_without_tuning":True,
        "2026_NPDES_can_miss_actual_2011_or_2016_historic_outfalls":True,
        "source_2011_2016_years":{}}
    for year in ("2011","2016"):
        sites,stop=load(year,root,a.source_dir,manifest)
        if stop:result["source_2011_2016_years"][year]={"stop":stop}
        else:result["source_2011_2016_years"][year]=audit(year,sites,sewage)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({year:{"source_n":v.get("original_PROB_firstvisit_fullveg_source_n"),
                "groups":v.get("original_shared_state_and_ecoregion_near10_far20_source_n"),
                "gate":v.get("source_precommitted_independent_year_gate"),
                "stop":v.get("stop")} for year,v in result["source_2011_2016_years"].items()},indent=2))
if __name__=="__main__":main()
