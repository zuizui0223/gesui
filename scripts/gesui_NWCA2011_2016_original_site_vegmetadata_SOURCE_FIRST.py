#!/usr/bin/env python3
"""EPA NWCA 2011/2016 original site+vegetation metadata ONLY.

No H_ALL/H_NAT/FQAI individual values, no N/P. Source opportunity is frozen
before any new historical wetland-plant relation. Report schema barriers
rather than guessing 2021 column identities or cherry-picking sites.
"""
from __future__ import annotations
import argparse,hashlib,json,re
from pathlib import Path
from datetime import datetime,timezone
from collections import Counter
import pandas as pd
import numpy as np

URLS={
  "2011":{"siteinfo":"https://www.epa.gov/system/files/other-files/2023-04/nwca11_siteinfo.csv",
          "vegmetrics":"https://www.epa.gov/sites/default/files/2016-10/nwca2011_vegmetrics.csv"},
  "2016":{"siteinfo":"https://www.epa.gov/system/files/other-files/2022-04/nwca-2016-site-information-data_0.csv",
          "vegmetrics":"https://www.epa.gov/system/files/other-files/2023-03/nwca_2016_veg_metrics.csv"},
}
CANON_ALIASES={
  "SITE_ID":["SITE_ID","SITEID","SITE-ID"],
  "VISIT_NO":["VISIT_NO","VISITNO","VISIT"],
  "LAT":["LAT_ANALYS","LATITUDE","LAT","LAT_DD","LATITUDE83","LATDD"],
  "LON":["LON_ANALYS","LONGITUDE","LON","LONG","LON_DD","LONGITUDE83","LONDD"],
  "STATE":["PSTL_CODE","STATE","STATE_CODE","PSTLCODE","STATE_ABBR"],
  "ECO":["AG_ECO9","AG_ECO","ECO9","ECORE_GRP","AG_ECOREGION"],
  "SITETYPE":["SITETYPE","SITE_TYPE","SAMPLE_TYPE","PANEL_TYPE"],
  "WEIGHT":["WGT_TP_CORE","WGT_TP","WGT","WGT_CORE","WGT_TP_ALL"],
  "YEAR":["YEAR","SAMPLE_YEAR","SAMP_YEAR","SURVEY_YEAR"],
  "STUDY":["STUDY"],
  "WETCLS":["WETCLS_GRP","WETCLS","NWCA_WETCLS_GRP"],
  "HYDCLS":["HYD_CLS","HYDCLS"]
}
METRIC_COLS=("H_ALL","HALL","H_NAT","HNAT","H_ALIEN","HALIEN","FQAI_ALL","FQAIALL","FQAI_NAT","FQAINAT")
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(2**20),b""):h.update(b)
    return h.hexdigest()
def norm(s):return re.sub("[^A-Z0-9]","",str(s).upper())
def choose(headers,k):
    lookup={norm(c):c for c in headers}
    for a in CANON_ALIASES[k]:
        if norm(a) in lookup:return lookup[norm(a)]
    return None
def header(p):
    return pd.read_csv(p,nrows=0,low_memory=False).columns.tolist()
def nrows(p):
    with p.open("rb") as f:return max(0,sum(1 for _ in f)-1)

def inspect(year,sitepath,metricpath):
    years={"2011":2011,"2016":2016}
    a=header(sitepath);b=header(metricpath)
    map_source={k:choose(a,k) for k in CANON_ALIASES}
    map_metric={k:choose(b,k) for k in ("SITE_ID","VISIT_NO")}
    result={
      "original_EPA_siteinfo_source_URL":URLS[year]["siteinfo"],
      "original_EPA_vegetation_metrics_source_URL":URLS[year]["vegmetrics"],
      "original_official_siteinfo_SHA256":sha(sitepath),
      "original_official_vegetation_metrics_SHA256":sha(metricpath),
      "source_siteinfo_columns":a,
      "source_vegetation_metric_column_NAMES_ONLY":b,
      "source_siteinfo_ID_columns_resolved":map_source,
      "source_metric_ID_columns_resolved":map_metric,
      "source_year":years[year],
      "original_siteinfo_source_row_count":nrows(sitepath),
      "original_vegetation_metrics_source_row_count":nrows(metricpath),
      "H_ALL_value_NOT_READ":True,
      "FQAI_and_biodiversity_metric_values_NOT_READ":True
    }
    req=("SITE_ID","LAT","LON","STATE","ECO","SITETYPE","WEIGHT","YEAR","WETCLS","HYDCLS")
    missing=[k for k in req if map_source[k] is None]
    if map_metric["SITE_ID"] is None:missing.append("METRIC.SITE_ID")
    if missing:
        result["source_gate_STATUS"]="SCHEMA_INSUFFICIENT_STOP"
        result["source_gate_missing_fields"]=missing
        return result
    # Read only originally fixed source site ID + design metadata, not biology.
    names=[v for v in map_source.values() if v]
    raw=pd.read_csv(sitepath,usecols=names,low_memory=False)
    c={k:raw[v] for k,v in map_source.items() if v}
    vtab=pd.read_csv(metricpath,usecols=[v for v in map_metric.values() if v],low_memory=False)
    ids=c["SITE_ID"].astype("string").str.strip()
    metids=vtab[map_metric["SITE_ID"]].astype("string").str.strip()
    if map_source["VISIT_NO"]:
        visit=pd.to_numeric(c["VISIT_NO"],errors="coerce")
    else:
        visit=pd.Series(1,index=raw.index,dtype=int)
    if map_metric["VISIT_NO"]:
        mvisit=pd.to_numeric(vtab[map_metric["VISIT_NO"]],errors="coerce")
    else:mvisit=pd.Series(1,index=vtab.index,dtype=int)
    kind=c["SITETYPE"].astype("string").str.strip().str.upper()
    weight=pd.to_numeric(c["WEIGHT"],errors="coerce")
    lat=pd.to_numeric(c["LAT"],errors="coerce")
    lon=pd.to_numeric(c["LON"],errors="coerce")
    year_num=pd.to_numeric(c["YEAR"],errors="coerce")
    original=pd.DataFrame({
        "id":ids,"visit":visit,"type":kind,"weight":weight,"lat":lat,"lon":lon,
        "year":year_num,"state":c["STATE"].astype("string").str.strip(),
        "eco":c["ECO"].astype("string").str.strip(),
        "wetcls":c["WETCLS"].astype("string").str.strip(),
        "hydcls":c["HYDCLS"].astype("string").str.strip()
    })
    result["source_site_type_counts"]=dict(Counter(kind.dropna()))
    result["source_year_counts"]=dict(Counter(year_num.dropna()))
    result["source_positive_weight_count"]=int(weight.gt(0).sum())
    result["source_valid_geographic_count"]=int((lat.between(24,50)&lon.between(-126,-66)).sum())
    result["metric_has_2011or2016_H_ALL_header"] = any(norm(s)=="HALL" for s in b)
    result["source_metric_visit_nonnull"] = int(mvisit.notna().sum())
    good=(original.type.eq("PROB")&
          original.weight.gt(0)&
          original.lat.between(24,50)&original.lon.between(-126,-66)&
          original.year.eq(years[year])&original.visit.eq(1)&
          original.state.notna()&original.eco.notna()&
          original.wetcls.notna()&original.hydcls.notna())
    original=original.loc[good].copy()
    if original.duplicated(["id","visit"]).any():
        result["source_gate_STATUS"]="DUPLICATE_SITE_VISIT_STOP"
        return result
    vkey=pd.DataFrame({"id":metids,"visit":mvisit}).dropna(subset=["id"]).drop_duplicates()
    common=original.merge(vkey,on=["id","visit"],how="inner",validate="one_to_one")
    source_countries=set(common.state)
    source_eco=set(common.eco)
    n=len(common)
    gate=n>=500 and len(source_countries)>=25 and len(source_eco)>=5
    result.update({
        "source_original_PROB_firstvisit_positiveweight_2021_style_sites_n":len(original),
        "source_vegmetric_SITE_ID_VISIT_NO_join_matched_n":n,
        "source_vegmetric_join_unique_states":len(source_countries),
        "source_vegmetric_join_ecoregions":len(source_eco),
        "source_site_NEAR_FAR_NOT_YET_COMPUTED":True,
        "source_precommitted_gate_A":{"at_least_500_source_weight_eligible_plant_visits":n>=500,
          "at_least_25_states":len(source_countries)>=25,
          "at_least_5_ecoregions":len(source_eco)>=5,
          "SOURCE_2011_2016_SITE_AND_VEGMETRICS_POPULATION_PASS":bool(gate)},
        "source_site_hashes_only":[hashlib.sha256((str(r.id)+"|"+str(int(r.visit))).encode()).hexdigest()
                                   for r in common.itertuples()],
        "source_group_counts_by_state":dict(Counter(common.state)),
        "source_group_counts_by_ecoregion":dict(Counter(common.eco)),
        "source_gate_STATUS":"SOURCE_SITE_VEG_PASS" if gate else "SOURCE_SITE_VEG_INSUFFICIENT_STOP"
    })
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-dir",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    result={
      "schema":"gesui_US_NWCA2011_2016_PROB_sites_vegmetrics_BEFORE_HALL_OUTCOME_v01",
      "utc":datetime.now(timezone.utc).isoformat(),
      "previous_2021_weighted_HALL_near10_positive_known_POSTHOC":True,
      "original_2021_registered_WWTP_proximity_not_actual_2011_2016_water_delivery":True,
      "source_years":{},
    }
    for year in ("2011","2016"):
        paths=[a.source_dir/(f"{year}_{kind}.csv") for kind in ("siteinfo","vegmetrics")]
        result["source_years"][year]=inspect(year,*paths)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({y:{k:v for k,v in x.items() if k in
         ("source_gate_STATUS","source_gate_missing_fields","source_site_type_counts",
          "source_year_counts","source_original_PROB_firstvisit_positiveweight_2021_style_sites_n",
          "source_vegmetric_SITE_ID_VISIT_NO_join_matched_n",
          "source_vegmetric_join_unique_states","source_vegmetric_join_ecoregions",
          "source_siteinfo_ID_columns_resolved","source_metric_ID_columns_resolved")}
         for y,x in result["source_years"].items()},indent=2))
if __name__=="__main__":main()
