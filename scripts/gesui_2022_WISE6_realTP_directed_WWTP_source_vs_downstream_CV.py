#!/usr/bin/env python3
"""POSTHOC new WISE6 river station-level directed WWTP→ACTUAL TP source test.

Source-only strictly directed network eligibility was Git-frozen before this
script can read the numeric site-year TP medians. Same original 2022 WISE6
numeric quality QC, 2021 WorldCover country-land classes, 2022 registered
WWTP annual measured N/P, third-party converted HydroRIVERS v1 EU network.
Original 2degree spatial heldout CV with downstream source placebo.
No causality or Lycaenidae population/feeding/host result can be inferred.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np

SOURCE_FIELDS={
  "strict_upstream_reported_P_mass_t_year":
      "strict_upstream_qualified_WWTP_reported_P_native_t_year",
  "strict_downstream_reported_P_mass_placebo":
      "strict_downstream_50km_WWTP_placebo_P_native_t_year",
  "strict_upstream_reported_N_mass_t_year":
      "strict_upstream_qualified_WWTP_reported_N_native_t_year",
  "strict_upstream_qualified_WWTP_count":
      "strict_upstream_qualified_WWTP_count_50km",
  "source_reported_P_over_modeled_receptor_flow_mgL":
      "source_apparent_P_load_over_modeled_avg_flow_mg_L",
}
MIN_EXPOSED=20
MIN_ZERO=20
MIN_POSITIVE_EACH_COUNTRY=5
MIN_POSITIVE_BLOCKS=5
MIN_POSITIVE_REACH=15
SOURCE_UPSTREAM_PRIMARY="strict_upstream_reported_P_mass_t_year"
SOURCE_DOWNSTREAM_PLACEBO="strict_downstream_reported_P_mass_placebo"
BOOTSTRAP_N=999
SEED=20261109

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for q in iter(lambda:f.read(1048576),b""):h.update(q)
    return h.hexdigest()

def keys(x):
    return (x["countryCode"],int(float(x["country_1deg_lat_floor"])),int(float(x["country_1deg_lon_floor"])))

def geo_block(x):
    return x["countryCode"]+"|"+str(x["original_country2degree_block"])

def response(rows):
    return np.asarray([math.log1p(float(x["TP_original_2022_annual_site_median_mg_P_per_L"])) for x in rows],float)

def values(rows,fn):
    return np.asarray([fn(x) for x in rows],float)

def design(train,test,predictor):
    funcs={
      "country_AT":lambda x:float(x["countryCode"]=="AT"),
      "latitude":lambda x:float(x["public_monitored_WGS84_lat"]),
      "longitude":lambda x:float(x["public_monitored_WGS84_lon"]),
      "modeled_discharge":lambda x:math.log1p(float(x["HydroRIVERS_original_modeled_DIS_AV_CMS"])),
      "modeled_basin_area":lambda x:math.log1p(float(x["HydroRIVERS_original_upstream_basin_UPLAND_SKM"])),
      "worldcover2021_built_country_fraction":lambda x:float(x["WC2021_512_built_land_fraction"]),
      "worldcover2021_crop_country_fraction":lambda x:float(x["WC2021_512_crop_land_fraction"]),
      "worldcover2021_landarea_km2":lambda x:math.log1p(float(x["WC2021_512_country_terrestrial_area_km2"]))
    }
    cols_train=[np.ones(len(train),float)]
    cols_test=[np.ones(len(test),float)]
    for name,fn in funcs.items():
        aa=values(train,fn);bb=values(test,fn)
        if name=="country_AT":
            cols_train.append(aa);cols_test.append(bb);continue
        m=float(aa.mean());std=float(aa.std())
        if std<1e-12:raise ValueError("No training variance in fixed original hydrological / landcover predictor "+name)
        cols_train.append((aa-m)/std);cols_test.append((bb-m)/std)
    X0=np.column_stack(cols_train);T0=np.column_stack(cols_test)
    name=SOURCE_FIELDS[predictor]
    a=values(train,lambda x:math.log1p(float(x[name])))
    b=values(test,lambda x:math.log1p(float(x[name])))
    za=np.zeros(len(train));zb=np.zeros(len(test))
    for country in ("DE","AT"):
        it=np.asarray([x["countryCode"]==country for x in train])
        ie=np.asarray([x["countryCode"]==country for x in test])
        if it.sum()<3:raise ValueError("Insufficient 2022 TP source station country train records")
        sub=a[it];m=float(sub.mean());std=float(sub.std())
        if std<1e-12:raise ValueError("No source variation in original country train region "+country+" "+predictor)
        za[it]=(sub-m)/std;zb[ie]=(b[ie]-m)/std
    X1=np.column_stack([X0,za]);T1=np.column_stack([T0,zb])
    return X0,T0,X1,T1

def OLS(X,T,Y):
    if np.linalg.matrix_rank(X)<X.shape[1]:
        raise ValueError("Numerical rank deficiency in original training-only blocked OLS design")
    return T@np.linalg.lstsq(X,Y,rcond=None)[0]

def blocked_cv(rows,source,include_predictions=False):
    groups=sorted(set(geo_block(x) for x in rows))
    predictions=[]
    for k in groups:
        train=[x for x in rows if geo_block(x)!=k]
        test=[x for x in rows if geo_block(x)==k]
        X0,T0,X1,T1=design(train,test,source)
        yt=response(train);y=response(test)
        p0=OLS(X0,T0,yt);p1=OLS(X1,T1,yt)
        for r,actual,a,b in zip(test,y,p0,p1):
            predictions.append({
              "WISE6_original_station_country_scheme_site_SHA256":r["original_WISE6_country_scheme_site_SHA256"],
              "country":r["countryCode"],
              "original_2degree_geographic_block":k,
              "original_HydroRIVERS_reach_ID":r["HYRIV_ID"],
              "original_WISE2022_TP_mgP_per_L":r["TP_original_2022_annual_site_median_mg_P_per_L"],
              "actual_TP_log1p":float(actual),
              "baseline_predicted_log1p_TP":float(a),
              "plus_source_predicted_log1p_TP":float(b),
              "baseline_squared_error":float((actual-a)**2),
              "source_squared_error":float((actual-b)**2),
              "baseline_absolute_error":float(abs(actual-a)),
              "source_absolute_error":float(abs(actual-b)),
              "original_strict_upstream_WWTP_P_mass":r["strict_upstream_qualified_WWTP_reported_P_native_t_year"],
              "original_strict_downstream_placebo_P_mass":r["strict_downstream_50km_WWTP_placebo_P_native_t_year"],
            })
    m0=float(np.mean([x["baseline_squared_error"] for x in predictions]))
    m1=float(np.mean([x["source_squared_error"] for x in predictions]))
    out={
      "source":source,"original_2022_TP_station_n":len(rows),
      "original_country_2degree_geographic_block_n":len(groups),
      "baseline_spatial_holdout_MSE":m0,
      "plus_source_spatial_holdout_MSE":m1,
      "relative_source_added_MSE_improvement_pct":100*(m0-m1)/m0,
      "baseline_spatial_holdout_MAE":float(np.mean([x["baseline_absolute_error"] for x in predictions])),
      "plus_source_spatial_holdout_MAE":float(np.mean([x["source_absolute_error"] for x in predictions])),
      "source_is_2022_reported_annual_WWTP_mass_not_independent_butterfly_nutrition":True,
      "source_importance_not_causal_estimate":True,
    }
    out["country_specific"]={}
    for country in ("DE","AT"):
        countryres=[x for x in predictions if x["country"]==country]
        a=float(np.mean([x["baseline_squared_error"] for x in countryres]))
        b=float(np.mean([x["source_squared_error"] for x in countryres]))
        out["country_specific"][country]={
         "n":len(countryres),"baseline_MSE":a,"plus_source_MSE":b,
         "source_added_MSE_improvement_pct":100*(a-b)/a
        }
    if include_predictions:
        out["original_2022_TP_country2degree_block_holdout_site_predictions"]=predictions
        grouped=defaultdict(list)
        for x in predictions:grouped[x["original_2degree_geographic_block"]].append(x)
        keys=list(grouped)
        rng=np.random.default_rng(SEED)
        deltas=[]
        for _ in range(BOOTSTRAP_N):
            sample=rng.choice(keys,size=len(keys),replace=True)
            z=[x for key in sample for x in grouped[key]]
            aa=float(np.mean([x["baseline_squared_error"] for x in z]))
            bb=float(np.mean([x["source_squared_error"] for x in z]))
            if aa>0:deltas.append(100*(aa-bb)/aa)
        out["source_2degree_error_block_bootstrap_descriptive"]={
            "seed":SEED,"n":len(deltas),
            "CI_95_percentile":[float(x) for x in np.quantile(deltas,[0.025,0.975])],
            "not_a_randomized_treatment_or_independent_preregistered_CI":True,
        }
    return out

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--source-only-directed-WWTP",type=Path,required=True)
    parser.add_argument("--original-WISE6-numeric-TP",type=Path,required=True)
    parser.add_argument("--original-WISE6-PREnumeric-TP-source",type=Path,required=True)
    parser.add_argument("--WorldCover-2021-photo-original",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    a=parser.parse_args()
    strict=json.loads(a.source_only_directed_WWTP.read_text())
    nums=json.loads(a.original_WISE6_numeric_TP.read_text())
    pre=json.loads(a.original_WISE6_PREnumeric_TP_source.read_text())
    wc=json.loads(a.WorldCover_2021_photo_original.read_text())
    if not strict["no_numeric_TP_opened_source_gate_pass"]:
        raise ValueError("Original real 2022 TP directed source-only sample gate DID NOT PASS; stop")
    if strict["original_2022_WISE6_metadata_QC_SHA256"]!=sha(a.original_WISE6_PREnumeric_TP_source):
        raise ValueError("Original directed WWTP source not linked to exact original metadata-only WISE TP station frame")
    if nums["original_PREnumeric_2022_TP_site_metadata_source_SHA256"]!=sha(a.original_WISE6_PREnumeric_TP_source):
        raise ValueError("Original numerical TP station values not linked to original nonnumeric TP source frame")
    if not nums["numeric_P_water_TP_to_Lycaenidae_photo_geographic_gate_pass"]:
        raise ValueError("Original TP numerical 2022 LOQ gate failed")
    if wc["status"]!="POSTHOC_LANDCOVER_OUTOFREGION_PREDICTION_AND_PLACEBO_COMPLETED":
        raise ValueError("Original WorldCover 2021 land source not fixed")
    seasource={x["original_station_siteKey_scheme_country_SHA256"]:x
         for x in nums["original_source_187_sites_all_numeric_and_censor_QC"]}
    world={(x["countryCode"],int(float(x["lat_floor"])),int(float(x["lon_floor"]))):x
          for x in wc["all_original_same_sample_cells_2022_photo_counts_and_WorldCover_categorical_2021_covariates"]}
    if len(seasource)!=187 or len(strict["original_187_WISE6_TP_stations_with_public_source_reach_exposures"])!=187:
        raise ValueError("Original same 187 TP stations/source-only directed exposure registry altered")
    rows=[];exclusions=Counter();bad=[]
    for e in strict["original_187_WISE6_TP_stations_with_public_source_reach_exposures"]:
        idhash=e["original_WISE6_country_scheme_site_SHA256"]
        v=seasource.get(idhash)
        if not v:raise ValueError("Unknown original WISE6 station SHA in directed source")
        if not v["numeric_site_support_ge4days_ge3months_LOQfraction_lehalf"]:
            exclusions["original_2022_numeric_TP_LOQ_or_repeated_days_fail"]+=1
            continue
        if e["HYRIV_ID"] is None:
            exclusions["original_WISE6_2022_no_unique_snap_to_directed_HydroRIVERS_500m"]+=1
            continue
        c=e["countryCode"]
        key=(c,int(float(e["country_1deg_lat_floor"])),int(float(e["country_1deg_lon_floor"])))
        if key not in world:
            raise ValueError("TP station cell not original 2021 WorldCover land-cover source")
        w=world[key]
        flow=e["HydroRIVERS_original_modeled_DIS_AV_CMS"]
        basin=e["HydroRIVERS_original_upstream_basin_UPLAND_SKM"]
        if not (isinstance(flow,(int,float)) and isinstance(basin,(int,float)) and
                 math.isfinite(flow) and flow>0 and math.isfinite(basin) and basin>0):
            exclusions["missing_or_nonpositive_original_HydroRIVERS_modeled_flow_or_basin"]+=1
            continue
        originalTP=v["original_EEA2022_TP_station_median_of_daily_medians_mg_P_per_L"]
        if originalTP is None or not math.isfinite(originalTP) or originalTP<0:
            raise ValueError("Previously QA-passed original 2022 measured TP numeric invalid")
        row=dict(e)
        row["TP_original_2022_annual_site_median_mg_P_per_L"]=float(originalTP)
        row["WC2021_512_built_land_fraction"]=float(w["WC2021_512_built_land_fraction"])
        row["WC2021_512_crop_land_fraction"]=float(w["WC2021_512_crop_land_fraction"])
        row["WC2021_512_country_terrestrial_area_km2"]=float(w["WC2021_512_country_terrestrial_area_km2"])
        rows.append(row)
    positive=[r for r in rows if r["strict_upstream_qualified_WWTP_count_50km"]>0]
    zero=[r for r in rows if r["strict_upstream_qualified_WWTP_count_50km"]==0]
    cs=Counter(r["countryCode"] for r in positive)
    regions=len(set(geo_block(r) for r in positive))
    reaches=len(set(r["HYRIV_ID"] for r in positive))
    ready=(len(positive)>=MIN_EXPOSED and len(zero)>=MIN_ZERO and
           cs["DE"]>=MIN_POSITIVE_EACH_COUNTRY and cs["AT"]>=MIN_POSITIVE_EACH_COUNTRY and
           regions>=MIN_POSITIVE_BLOCKS and reaches>=MIN_POSITIVE_REACH)
    report={
       "schema":"gesui_2022_actual_WISE6_riverTP_strict_upstream_2022_EEA_WWTP_HydroRIVERS_site_CV_v01",
       "UTC":datetime.now(timezone.utc).isoformat(),
       "posthoc_relative_to_2022_source_photo_geography":True,
       "source_only_directed_WWTP_source_gate_SHA256":sha(a.source_only_directed_WWTP),
       "2022_WISE6_original_PREnumeric_site_source_SHA256":sha(a.original_WISE6_PREnumeric_TP_source),
       "2022_WISE6_original_numeric_TP_site_year_source_SHA256":sha(a.original_WISE6_numeric_TP),
       "WorldCover2021_original_frozen_land_SHA256":sha(a.WorldCover_2021_photo_original),
       "source_mirror_HydroRIVERS_converted_topology_independently_official_verified":False,
       "source_hydrological_directional_distance_max_km":50,
       "source_station_snap_max_m":500,
       "original_source_187_station_count":187,
       "original_2022_TP_numeric_and_directed_reach_eligible_station_count":len(rows),
       "source_baseline_missingness_or_TP_quality_exclusion":dict(exclusions),
       "2022_strict_upstream_qualified_WWTP_P_positive_station_n":len(positive),
       "2022_registered_qualifying_upstream_source_zero_station_n":len(zero),
       "positive_station_counts_by_country":dict(cs),
       "positive_independent_river_reaches":reaches,
       "positive_original_country2degree_blocks":regions,
       "source_and_numeric_support_regression_gate_pass":ready,
       "source_zero_is_not_no_WWTP_or_no_river_P":True,
       "no_butterfly_2022_photos_opened_here":True,
       "source_annual_WWTP_reported_mass_not_measured_chemical_dose_at_monitoring_receptor":True,
       "original_2022_river_2021_land_2022_EEA_site_source_model_results":{},
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    if not ready:
        report["STOP"]="NO_DIRECTED_SOURCE_EXPOSURE_COUNTRY_AND_SPATIAL_SUPPORT_AFTER_NUMERIC_TP_QUALITY; NO EFFECT MODEL"
        a.out.write_text(json.dumps(report,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
        print("GESUI_WISE6_2022_DIRECTED_WWTP_TP_SOURCE_CV_STOP",json.dumps({
           "n":len(rows),"positive":len(positive),"negative":len(zero),
           "positive_country":dict(cs),"blocks":regions,"reaches":reaches},sort_keys=True),flush=True)
        return
    report["directed_source_station_model_input_records"]=[{
       "original_site_SHA256":r["original_WISE6_country_scheme_site_SHA256"],
       "country":r["countryCode"],"country2degree_block":r["original_country2degree_block"],
       "HYRIV_ID":r["HYRIV_ID"],"original_2022_WISE6_TP_mgP_per_L":r["TP_original_2022_annual_site_median_mg_P_per_L"],
       "original_2022_EEA_strict_upstream_P_t_year":r["strict_upstream_qualified_WWTP_reported_P_native_t_year"],
       "original_2022_EEA_strict_downstream_P_t_year":r["strict_downstream_50km_WWTP_placebo_P_native_t_year"],
       "modeled_original_river_flow_cms":r["HydroRIVERS_original_modeled_DIS_AV_CMS"],
       "original_WC2021_built_fraction":r["WC2021_512_built_land_fraction"],
       "original_WC2021_crop_fraction":r["WC2021_512_crop_land_fraction"],
       "no_causal_claim":True,
    } for r in rows]
    models={}
    for source in SOURCE_FIELDS:
        models[source]=blocked_cv(rows,source,include_predictions=True)
        v=models[source]
        print("GESUI_2022_ACTUAL_WISE6_TP_DIRECTED_WWTP_RECEPTOR_GEOBLOCK_CV",
          json.dumps({"source":source,
            "TP_site_years":len(rows),"regions":v["original_country_2degree_geographic_block_n"],
            "MSE_baseline":v["baseline_spatial_holdout_MSE"],
            "MSE_with_source":v["plus_source_spatial_holdout_MSE"],
            "improvement_percent":v["relative_source_added_MSE_improvement_pct"],
            "DE_AT":v["country_specific"]},sort_keys=True),flush=True)
    report["original_2022_river_2021_land_2022_EEA_site_source_model_results"]=models
    report["upstream_P_minus_downstream_placebo_incremental_CV_MSE_gain_pct"]=(models[SOURCE_UPSTREAM_PRIMARY]["relative_source_added_MSE_improvement_pct"]-
                     models[SOURCE_DOWNSTREAM_PLACEBO]["relative_source_added_MSE_improvement_pct"])
    report["not_yet_a_randomized_WWTP_water_quality_treatment_effect"]=True
    report["not_Lycaenidae_organicN_or_puddling_or_population_response"]=True
    # Original 500m station/receptor and strict directed 50km WWTP matching,
    # but only source-reported mass, not actual annual WWTP-to-station dose.
    a.out.write_text(json.dumps(report,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GESUI_WISE6_2022_DIRECTED_WWTP_VS_REAL_TP_SOURCE_PRIMARY_AND_DOWNSTREAM_PLACEBO_COMPLETED",
      json.dumps({"n":len(rows),"upstream_P_gain_pct":models[SOURCE_UPSTREAM_PRIMARY]["relative_source_added_MSE_improvement_pct"],
       "downstream_P_placebo_gain_pct":models[SOURCE_DOWNSTREAM_PLACEBO]["relative_source_added_MSE_improvement_pct"],
       "flow_normalized_P_gain_pct":models["source_reported_P_over_modeled_receptor_flow_mgL"]["relative_source_added_MSE_improvement_pct"],
       "n_reported_EEA2022_sources_only":True},sort_keys=True),flush=True)

if __name__=="__main__":main()
