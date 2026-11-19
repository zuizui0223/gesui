#!/usr/bin/env python3
"""Adversarial SOURCE-ONLY check of water point-N versus soil N-input geography.

Analyses no EEA EQR, plant EIV or butterfly HNS. The original water-by-land
2x2 biotic source gate stays STOP and cannot be reopened by this script.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from datetime import datetime, timezone

import geopandas as gpd
import numpy as np
import pandas as pd

SOURCE_KEYS = [
    "water_point_N_t",
    "water_point_P_t",
    "land_atmos_N_mg_m2",
    "land_mineral_N_kg_ha",
    "land_organic_N_kg_ha",
    "land_fixation_N_kg_ha",
]
LAND = SOURCE_KEYS[2:]
PRIMARY_LAND = LAND[:3]
EXPECTED_SOURCE_ROW_N = 269
EXPECTED_COUNTRY_N = 28
BOOT_N = 3999
SEED = 20261119
GPKG_SHA256 = "7a233d3d683dd48241b4a9e5a9ef6c5f3e294008934c5bc7f1a465eb5cbf9c0b"

def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def clean(v):
    if isinstance(v, dict):
        return {str(k): clean(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [clean(x) for x in v]
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating, float)):
        return float(v) if math.isfinite(float(v)) else None
    return v

def pearson(x, y, weight=None):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if len(x) != len(y) or len(x) < 3 or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("Invalid common source sample")
    w = np.ones(len(x), float) if weight is None else np.asarray(weight, float)
    if len(w) != len(x) or not np.isfinite(w).all() or (w < 0).any() or w.sum() <= 0:
        raise ValueError("Invalid weighting")
    w = w / w.sum()
    xx = x - np.sum(x * w)
    yy = y - np.sum(y * w)
    den = np.sqrt(np.sum(w * xx * xx) * np.sum(w * yy * yy))
    return float(np.sum(w * xx * yy) / den) if den > 1e-12 else None

def group_demean(x, country):
    x = np.asarray(x, float)
    g = np.asarray(country)
    z = x.copy()
    for c in sorted(set(g)):
        mask = g == c
        z[mask] -= x[mask].mean()
    return z

def country_equal_weights(country):
    c = np.asarray(country)
    out = np.zeros(len(c), float)
    for key in sorted(set(c)):
        ix = c == key
        out[ix] = 1.0 / ix.sum()
    return out

def within_country_r(x, y, country, equal_country=True):
    x = group_demean(x, country)
    y = group_demean(y, country)
    w = country_equal_weights(country) if equal_country else None
    return pearson(x, y, w)

def moments_by_country(x, y, country):
    x=np.asarray(x,float)
    y=np.asarray(y,float)
    g=np.asarray(country)
    out=[]
    for c in sorted(set(g)):
        a=x[g==c]
        b=y[g==c]
        da=a-a.mean()
        db=b-b.mean()
        out.append({"country":str(c),"n":len(a),
                    "xx":float(np.mean(da*da)),
                    "yy":float(np.mean(db*db)),
                    "xy":float(np.mean(da*db)),
                    "r_within_country":pearson(a,b) if len(a)>=3 else None})
    return out

def bootstrap_country_r(moment, seed=SEED, nboot=BOOT_N):
    x=np.array([z["xx"] for z in moment])
    y=np.array([z["yy"] for z in moment])
    cross=np.array([z["xy"] for z in moment])
    if len(x)<2 or min(x.sum(),y.sum())<=0:
        return [None,None]
    rng=np.random.default_rng(seed)
    draws=rng.integers(0,len(moment),size=(nboot,len(moment)))
    numerator=cross[draws].sum(axis=1)
    den=np.sqrt(x[draws].sum(axis=1)*y[draws].sum(axis=1))
    ok=den>1e-12
    if ok.sum()<nboot*.9:
        return [None,None]
    return np.quantile(numerator[ok]/den[ok],[.025,.975]).tolist()

def country_loo(moment):
    output=[]
    for i,z in enumerate(moment):
        without=[v for j,v in enumerate(moment) if j!=i]
        numerator=sum(q["xy"] for q in without)
        den=np.sqrt(sum(q["xx"] for q in without)*sum(q["yy"] for q in without))
        output.append({"excluded_country":z["country"],"r_within_others":float(numerator/den) if den>1e-12 else None})
    return output

def residual(y, x):
    y=np.asarray(y,float)
    x=np.asarray(x,float)
    if not np.isfinite(x).all() or not np.isfinite(y).all():raise ValueError("Incomplete design")
    coef=np.linalg.lstsq(x,y,rcond=None)[0]
    return y-x@coef

def fit_source_multivariate(point, predictors, baseline):
    rp=residual(point,baseline)
    rd=np.column_stack([residual(p,baseline) for p in predictors])
    if (rp@rp)<1e-10:raise ValueError("PointN residual source variation zero")
    reduced=rp@rp
    errors=residual(rp,rd)
    full=errors@errors
    r2=1-float(full/reduced)
    if r2 < -1e-8 or r2 > 1.0000001:raise ValueError("Unexpected source partial R2")
    r2=float(np.clip(r2,0,1))
    vif=float(1/(1-r2)) if r2<1 else None
    vifs=[]
    for i in range(rd.shape[1]):
        tar=rd[:,i]
        other=rd[:,[j for j in range(rd.shape[1]) if j!=i]]
        err=residual(tar,other)
        if tar@tar <1e-10:
            vifs.append(None)
        else:
            r=1-float((err@err)/(tar@tar))
            vifs.append(float(1/(1-r)) if r<1 else None)
    return {"partial_R2_pointN_given_landN":r2,
            "VIF_pointN_given_landN":vif,
            "VIF_each_landN_given_other_landN":vifs}

def load_panel(source_csv, source_report, gisco_gpkg):
    report=json.loads(Path(source_report).read_text())
    if report.get("NO_EEA_EQR_EIV_GBIF_BIOLOGICAL_VALUES_READ") is not True:
        raise ValueError("Original JRC source-only report not attested")
    if sha(source_csv)!=report["csv_sha256"]:
        raise ValueError("Original JRC NUTS2 model values CSV digest changed")
    if sha(gisco_gpkg)!=GPKG_SHA256:
        raise ValueError("Official GISCO 2016 01M reference digest changed")
    df=pd.read_csv(source_csv,dtype={"NUTS2_country_id":str})
    if not set(["NUTS2_country_id",*SOURCE_KEYS]).issubset(df.columns):
        raise ValueError("Missing source-only NUTS2 fields")
    if df.NUTS2_country_id.duplicated().any():raise ValueError("Duplicate GISCO NUTS2 sources")
    for k in SOURCE_KEYS:
        df[k]=pd.to_numeric(df[k],errors="coerce")
    df=df.dropna(subset=SOURCE_KEYS).copy()
    if len(df)!=EXPECTED_SOURCE_ROW_N or df.NUTS2_country_id.str.split("|",regex=False).str[0].nunique()!=EXPECTED_COUNTRY_N:
        raise ValueError("Fixed 269/28 source cohort changed")
    if (df[SOURCE_KEYS]<0).any().any():raise ValueError("Negative N/P model source")
    polygons=gpd.read_file(gisco_gpkg)
    if not set(["NUTS_ID","CNTR_CODE","geometry"]).issubset(polygons.columns):
        raise ValueError("Official 2016 GISCO NUTS2 schema absent")
    if "LEVL_CODE" in polygons.columns:
        polygons=polygons.loc[polygons.LEVL_CODE.astype(str).eq("2")].copy()
    polygons["NUTS2_country_id"]=polygons.CNTR_CODE.astype(str).str.upper()+"|"+polygons.NUTS_ID.astype(str).str.upper()
    if polygons.NUTS2_country_id.duplicated().any():raise ValueError("Duplicate GISCO polygon ID")
    pol=polygons.to_crs(3035)
    polygons["area_km2"]=(pol.geometry.area/1e6).to_numpy()
    points=gpd.GeoSeries(pol.geometry.centroid,crs=3035).to_crs(4326)
    polygons["lon"]=points.x.to_numpy()
    polygons["lat"]=points.y.to_numpy()
    df=df.merge(polygons[["NUTS2_country_id","area_km2","lat","lon"]],on="NUTS2_country_id",how="left",validate="one_to_one")
    if df[["area_km2","lat","lon"]].isna().any().any():
        raise ValueError("Official GISCO 2016 01M 269-area coverage missing; no improvised replacement")
    if not (df.area_km2>0).all():raise ValueError("Invalid 2016 regional GIS area")
    df["country"]=df.NUTS2_country_id.str.split("|",regex=False).str[0]
    df["point_raw_log"]=np.log1p(df["water_point_N_t"])
    df["point_density_log"]=np.log1p(df["water_point_N_t"]/df.area_km2)
    for k in LAND:df[k+"_log"]=np.log1p(df[k])
    df["logarea"]=np.log1p(df.area_km2)
    return df, report, len(polygons)

def analyze(df):
    countries=df.country.to_numpy()
    xraw=df.point_raw_log.to_numpy()
    x=df.point_density_log.to_numpy()
    dum=pd.get_dummies(df.country,dtype=float).to_numpy()
    geographic=np.column_stack([dum,df[["logarea","lat","lon"]].to_numpy(float)])
    records={}
    for i,k in enumerate(LAND):
        y=df[k+"_log"].to_numpy()
        moments=moments_by_country(x,y,countries)
        within_equal=within_country_r(x,y,countries,True)
        ci=bootstrap_country_r(moments,seed=SEED+i*100)
        loo=country_loo(moments)
        validloo=[v["r_within_others"] for v in loo if v["r_within_others"] is not None]
        direction=np.sign(within_equal) if within_equal is not None else 0
        maintained=sum(v*direction>0 and abs(v)>=.30 for v in validloo)/len(validloo) if validloo and direction else None
        records[k]={
            "raw_269_Pearson_log_mass_vs_log_landinput":pearson(xraw,y),
            "density_269_Pearson_log_t_per_km2_vs_log_landinput":pearson(x,y),
            "density_269_Spearman_rank":pearson(pd.Series(x).rank().to_numpy(),pd.Series(y).rank().to_numpy()),
            "density_within_country_Pearson_region_equal":within_country_r(x,y,countries,False),
            "density_within_country_Pearson_country_equal":within_equal,
            "density_country_equal_3999_cluster_bootstrap95_descriptive":ci,
            "country_area_lat_lon_adjusted_partial_Pearson":pearson(residual(x,geographic),residual(y,geographic)),
            "countries_with_distinct_observed_values":sum(z["xx"]>0 and z["yy"]>0 for z in moments),
            "individual_country_univariate_within_r":moments,
            "leave_one_country_out":loo,
            "LOCO_same_sign_and_abs_r_ge_0p30_fraction":maintained,
        }
    logland=[df[k+"_log"].to_numpy() for k in LAND]
    three=fit_source_multivariate(x,logland[:3],geographic)
    four=fit_source_multivariate(x,logland,geographic)
    dcm=np.column_stack(logland)
    cross={a:{b:pearson(dcm[:,i],dcm[:,j]) for j,b in enumerate(LAND) if j>i} for i,a in enumerate(LAND)}
    strong=any((r["density_within_country_Pearson_country_equal"] is not None
                 and abs(r["density_within_country_Pearson_country_equal"])>=.50
                 and r["LOCO_same_sign_and_abs_r_ge_0p30_fraction"] is not None
                 and r["LOCO_same_sign_and_abs_r_ge_0p30_fraction"]>=.70) for r in records.values())
    weak=(all((records[k]["density_within_country_Pearson_country_equal"] is not None
               and abs(records[k]["density_within_country_Pearson_country_equal"])<.30
               and records[k]["density_country_equal_3999_cluster_bootstrap95_descriptive"][0] is not None
               and records[k]["density_country_equal_3999_cluster_bootstrap95_descriptive"][0]>-.30
               and records[k]["density_country_equal_3999_cluster_bootstrap95_descriptive"][1]<.30)
              for k in PRIMARY_LAND)
          and three["VIF_pointN_given_landN"] is not None
          and three["VIF_pointN_given_landN"]<5)
    classification="STRONG_WITHIN_COUNTRY_COVARIATION_COUNTEREVIDENCE" if strong else (
        "WEAK_LINEAR_REDUNDANCY_SUPPORTED_WITH_SOURCE_ONLY_CAVEAT" if weak else
        "INDETERMINATE_SOURCE_GEOGRAPHIC_COUPLING")
    return records, {"source_partial_collinearity_three_landN":three,
                     "source_partial_collinearity_four_landN":four,
                     "inter_land_source_Pearson_loginput":cross,
                     "source_geographic_pattern_classification":classification,
                     "biotic_2x2_pr125_source_STOP_remains":True}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-csv",type=Path,required=True)
    ap.add_argument("--source-report",type=Path,required=True)
    ap.add_argument("--gisco-gpkg",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    df,source_receipt,polygon_n=load_panel(a.source_csv,a.source_report,a.gisco_gpkg)
    records,drivers=analyze(df)
    report={
        "schema":"gesui_20161119_NUTS2_water_pointN_vs_terrestrial_N_source_only_adversarial_v01",
        "utc":datetime.now(timezone.utc).isoformat(),
        "n_NUTS2":len(df),"n_countries":df.country.nunique(),
        "official_GISCO2016_polygon_count":polygon_n,
        "original_SOURCE_CSV_sha256":sha(a.source_csv),
        "original_GISCO2016_GPKG_sha256":sha(a.gisco_gpkg),
        "source_variables_not_measured_ecological_dose":True,
        "source_country_number_with_two_or_more_regions":int((df.country.value_counts()>=2).sum()),
        "source_country_number_with_three_or_more_regions":int((df.country.value_counts()>=3).sum()),
        "diagnostics":records,"joint_diagnostics":drivers,
        "interpretation":"Posthoc source-only geographic comparison; not wastewater-specific effect, no aquatic/terrestrial biological outcomes read; original 2x2 source gate remains STOP.",
    }
    reportpath=a.out/"POINT_LAND_N_WITHIN_COUNTRY_AREA_ADVERSARIAL_20261119.json"
    reportpath.write_text(json.dumps(clean(report),indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    cols=["NUTS2_country_id","country","area_km2","lat","lon","point_raw_log","point_density_log",*SOURCE_KEYS]
    df[cols].to_csv(a.out/"POINT_LAND_N_NUTS2_SOURCE_GIS_AREA_COMPLETE_269.csv",index=False,float_format="%.12g")
    print("GESUI_POINT_LAND_N_ADVERSARIAL_SOURCE_ONLY",json.dumps(clean({
        "n_regions":len(df),"n_countries":df.country.nunique(),
        "source_correlations":{k:{
          key:records[k][key] for key in [
            "raw_269_Pearson_log_mass_vs_log_landinput",
            "density_269_Pearson_log_t_per_km2_vs_log_landinput",
            "density_within_country_Pearson_country_equal",
            "density_country_equal_3999_cluster_bootstrap95_descriptive",
            "country_area_lat_lon_adjusted_partial_Pearson",
            "LOCO_same_sign_and_abs_r_ge_0p30_fraction"
          ]} for k in LAND},
        "joint":drivers
    }),ensure_ascii=False),flush=True)

if __name__=="__main__":main()
