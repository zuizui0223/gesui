#!/usr/bin/env python3
"""Descriptive SOURCE-ONLY JRC river-basin N geographic associations.

Uses 9739 *candidate* HydroID/BasinID with exact native area and a prior
independent 30/30 polygon geometry pilot. Does NOT assume native source basins
are independent hydrological replicates. No new biology outcomes or p-values.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import pandas as pd

LAND={
  "atmosN":"soil_atmos_N_mean_mg_m2",
  "mineralN":"soil_mineral_N_mean_kg_ha",
  "organicN":"soil_organic_N_mean_kg_ha",
  "fixationN":"soil_fixation_N_mean_kg_ha",
}
ID="candidate_ID_number_not_proof_of_same_basin"

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def pearson(x,y):
    x=np.asarray(x,float);y=np.asarray(y,float)
    if len(x)<3 or np.std(x)<1e-12 or np.std(y)<1e-12:return None
    return float(np.corrcoef(x,y)[0,1])

def rankcorr(x,y):
    return pearson(pd.Series(x).rank(method="average").to_numpy(),
                   pd.Series(y).rank(method="average").to_numpy())

def source_input(full_csv,full_json,geom_json):
    original=json.loads(Path(full_json).read_text())
    verified=json.loads(Path(geom_json).read_text())
    info=original["source_geometry_match_precheck"]
    polygon=verified["source_polygon_geometry_support"]
    if info["matching_integer_ID_and_native_area_only_geographic_candidates"]!=9739:
        raise ValueError("Original native 9739 sources changed")
    if not polygon.get("source_polygon_pilot_geospatial_match_PASS") or polygon.get("original_polygon_sample_tested")!=30:
        raise ValueError("Fixed 30 native polygons did not pass preread source geometry audit")
    if original["candidateCSVsha256"]!=sha(full_csv) or verified["full_JRC_candidateCSV_SHA256"]!=sha(full_csv):
        raise ValueError("Full JRC original candidate CSV digest mismatch between source and geometry")
    if not original.get("no_EEA_EQR_EIV_NUTRI_HNS_outcome_opened") or not verified.get("no_biological_EQR_EIV_butterfly_values_opened"):
        raise ValueError("Source-only contract violated")
    x=pd.read_csv(full_csv)
    required=[ID,"soil_area_sqkm","point_area_sqkm","water_point_N_modelled_t",*LAND.values()]
    if not set(required).issubset(x.columns) or len(x)!=9739 or x[ID].duplicated().any():
        raise ValueError("Original native source census changed")
    for field in required:
        x[field]=pd.to_numeric(x[field],errors="coerce")
    if x[required].isna().any().any() and any(x[k].isna().any() for k in required if k!="soil_fixation_N_mean_kg_ha"):
        raise ValueError("Three primary source fields missing")
    for col in ["soil_area_sqkm","point_area_sqkm","water_point_N_modelled_t",*LAND.values()]:
        if (x[col].dropna()<0).any():raise ValueError("Negative JRC modeled native source value")
    if not (x.soil_area_sqkm>0).all():raise ValueError("Invalid basin native area")
    return x

def summarize(x):
    x=x.copy()
    x["pointN_t_per_km2"]=x.water_point_N_modelled_t/x.soil_area_sqkm
    x["pointN_density_log"]=np.log1p(x.pointN_t_per_km2)
    x["pointN_total_log"]=np.log1p(x.water_point_N_modelled_t)
    quart=pd.qcut(x.soil_area_sqkm,4,duplicates="raise",labels=False)
    x["area_quartile"]=quart+1
    overall={}
    for name,col in LAND.items():
        y=np.log1p(x[col].astype(float))
        overall[name]={
           "complete_pairs":int(y.notna().sum()),
           "Pearson_log_pointN_density_vs_log_soilN_mean":pearson(x.pointN_density_log,y),
           "Spearman_pointN_density_vs_soilN_mean":rankcorr(x.pointN_t_per_km2,x[col]),
           "Pearson_log_pointN_total_vs_log_soilN_mean":pearson(x.pointN_total_log,y),
           "within_area_quartile_Pearson_log_density":{}
        }
        for q,g in x.groupby("area_quartile"):
            v=np.log1p(g[col].to_numpy(float))
            ok=np.isfinite(v)
            overall[name]["within_area_quartile_Pearson_log_density"][str(int(q))]={
              "n":int(sum(ok)),"r":pearson(g.pointN_density_log.to_numpy()[ok],v[ok])
            }
        positive=x.water_point_N_modelled_t>0
        overall[name]["point_positive_only_descriptive_r"]=pearson(
          x.loc[positive & y.notna(),"pointN_density_log"],
          y[positive&y.notna()])
    return {
      "native_candidate_count":len(x),
      "native_model_pointN_positive_n":int((x.water_point_N_modelled_t>0).sum()),
      "native_model_pointN_true_zero_n":int((x.water_point_N_modelled_t==0).sum()),
      "native_soil_fixation_missing_n":int(x[LAND["fixationN"]].isna().sum()),
      "native_basin_area_sqkm_q05_q50_q95":[float(v) for v in np.quantile(x.soil_area_sqkm,[.05,.50,.95])],
      "native_all_candidate_area_quartiles":{str(int(k)):{
        "native_basin_n":len(v),
        "native_model_pointN_positive_n":int((v.water_point_N_modelled_t>0).sum())}
        for k,v in x.groupby("area_quartile")},
      "descriptive_source_correlations":overall,
      "correlations_not_independent_samples_or_source_causes":True,
      "native_basin_poly_geo_30of30_PASS_but_not_all_geometries_verified":True
    }

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--native-source-csv",type=Path,required=True)
    a.add_argument("--native-source-report",type=Path,required=True)
    a.add_argument("--geom30-report",type=Path,required=True)
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    x=source_input(args.native_source_csv,args.native_source_report,args.geom30_report)
    results=summarize(x)
    o={
      "schema":"gesui_JRC_9739_native_candidate_basin_sourceN_after_polygon30_descriptive_v1",
      "UTC":datetime.now(timezone.utc).isoformat(),
      "native_candidate_sourceCSV_SHA256":sha(args.native_source_csv),
      "native_polygon30_original_report_SHA256":sha(args.geom30_report),
      "point_source_year":"GREEN average2014-2018",
      "land_atmosN_year":"EMEP 2015",
      "land_fertilizer_year":"DAYCENT 2010-2019",
      "source_values_not_ecological_response":True,
      "source_only_descriptive":results,
      "no_EQR_EIV_NUTRI_HNS_opened":True,
      "original_2018_water_2017to19_land_joint_2x2_source_STOP":True
    }
    (args.out/"JRC_9739_NATIVE_SOURCE_POINTN_SOILN_DESCRIPTIVE_20261119.json").write_text(
      json.dumps(o,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("JRC_NATIVE_9739_SOURCE_DESCRIPTIVE",json.dumps(results,ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":main()
