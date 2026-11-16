#!/usr/bin/env python3
"""Source-only 2018 EEA measured WWTP N/P mass → historical 2018 river BQE sites.

Original independent 2017–2018 EEA UWWTD (NOT 2022 source); 2018 BQE sites
and five geographic fold assignments from separate pre-outcome source-only
run. All numerical EQR response fields deliberately absent from input.
2018 WWTP source and 2018 biological-source sample are compared BEFORE
opening any 2018 biological scores.
"""
from __future__ import annotations
import argparse,hashlib,io,json,zipfile
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
import geopandas as gpd
import gesui_eea_nutrient_load_river_crosswalk as eea

SROOT="eea_t_uwwtd-waterbase_p_2017-2018_v07_r01/Waterbase_UWWTD_v8_csv/"
SOURCE_NAMES=("UWWTPS.csv","DischargePoints.csv","ReportPeriod.csv")
YEAR="2018"
BQE_CODE={
 "EEA_13-01-4":"benthic_invertebrates",
 "EEA_124-04-9":"phytobenthos",
 "EEA_14-02-8":"fish"}
BBOX=(-12.,44.,26.,64.)
EXPECTED_HISTORICAL_MEASURED_BOTH=2507
EXPECTED_HISTORICAL_FW_BBOX=1889
RIVER_SHA="5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f"
SOURCE_2018_QC_RUN=42441852933
SITE_2018_SOURCE_RUN=42441464939
EXACT_SNAP_METERS=500

def sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()

def load_original_source(archive,receipt):
    q=json.loads(receipt.read_text())
    saved=q["source_QC"]
    if saved["source_2018_NP_both_measured_positive_plant_rows"]!=EXPECTED_HISTORICAL_MEASURED_BOTH:
        raise ValueError("Original 2018 source both-positive measured quality gate changed")
    if saved["source_2018_positive_both_measured_unique_outlet_direct_FW_in_BBOX"]!=EXPECTED_HISTORICAL_FW_BBOX:
        raise ValueError("Original fixed geographic official FW source candidate count changed")
    out={}
    with zipfile.ZipFile(archive) as z:
        for name in SOURCE_NAMES:
            full=SROOT+name
            raw=z.read(full)
            expected=saved["source_table_original_SHA256"][name]["source_uncompressed_SHA256"]
            digest=hashlib.sha256(raw).hexdigest()
            if expected!=digest:raise ValueError("EEA 2018 source table bytes changed "+name)
            out[name]=pd.read_csv(io.BytesIO(raw),dtype="string",low_memory=False)
    return out,saved

def source_outfalls(tables):
    plants=tables["UWWTPS.csv"].copy()
    outlets=tables["DischargePoints.csv"].copy()
    periods=tables["ReportPeriod.csv"].copy()
    period_id=periods[["repCode","repReportedPeriod","repSituationAt"]].copy()
    for col in ("repReportedPeriod","repSituationAt","repCode"):
        period_id[col]=period_id[col].astype("string").str.strip()
    if period_id.repCode.duplicated().any():
        raise ValueError("Historical EEA reporter source has multiple incompatible year reference rows")
    valid_period=period_id[
       period_id.repReportedPeriod.eq(YEAR)&
       period_id.repSituationAt.str.contains("2018",regex=False,na=False)]
    if len(valid_period)!=30:
        raise ValueError("Original 30 countries must explicitly report 2018 reference year")
    allowed=set(valid_period.repCode)
    plants=plants[plants.repCode.isin(allowed)].copy()
    outlets=outlets[outlets.repCode.isin(allowed)].copy()
    for tbl in (plants,outlets):
        for c in ("repCode","uwwCode"):
            tbl[c]=tbl[c].astype("string").str.strip()
        tbl["source_key"]=tbl.repCode.fillna("")+"|"+tbl.uwwCode.fillna("")
        tbl.loc[tbl.repCode.isna()|tbl.uwwCode.isna()|
                tbl.repCode.eq("")|tbl.uwwCode.eq(""),"source_key"]=pd.NA
    for col in ("uwwNDischargeMeasured","uwwPDischargeMeasured"):
        plants[col]=pd.to_numeric(plants[col],errors="coerce")
    # No pandas nullable .all() false positives:
    vals=plants[["uwwNDischargeMeasured","uwwPDischargeMeasured"]].to_numpy(
        dtype=float,na_value=np.nan)
    good=np.isfinite(vals)&(vals>0)
    plants=plants.loc[good.all(axis=1)&plants.source_key.notna()].copy()
    if len(plants)!=EXPECTED_HISTORICAL_MEASURED_BOTH:
        raise ValueError("Original 2018 both-positive measured plant selection changed")
    if plants.source_key.duplicated().any():
        raise ValueError("2018 unique measured plant source IDs not actually unique")
    outlets=outlets[outlets.source_key.notna()].copy()
    counts=outlets.source_key.value_counts()
    plants["count_matching_outlets"]=plants.source_key.map(counts).fillna(0).astype(int)
    eligible=plants[plants.count_matching_outlets.eq(1)].merge(
      outlets.drop_duplicates("source_key"),on="source_key",
      how="inner",validate="one_to_one",suffixes=("_plant","_outlet"))
    for c in ("dcpLatitude","dcpLongitude"):
        eligible[c]=pd.to_numeric(eligible[c],errors="coerce")
    eligible=eligible[
      eligible.dcpLatitude.between(BBOX[1],BBOX[3])&
      eligible.dcpLongitude.between(BBOX[0],BBOX[2])&
      eligible.dcpWaterBodyType.astype("string").str.strip().str.upper().eq("FW")
    ].copy()
    if len(eligible)!=EXPECTED_HISTORICAL_FW_BBOX:
        raise ValueError("Original FW freshwater unique outlet number shifted")
    # Primary historical source has 2018 measured year and direct freshwater,
    # but operational-status and surfaceWaters are reported for audit rather
    # than excluded by inspecting any EQR outcomes.
    quality={
     "2018_original_both_positive_NP_measured_plant_n":len(plants),
     "2018_FW_direct_freshwater_unique_outlets_in_bbox":len(eligible),
     "2018_plant_outlet_status_counts":{
       str(k):int(v) for k,v in eligible.uwwState.fillna("MISSING").value_counts().items()},
     "2018_discharge_surfacewater_flag_counts":{
       str(k):int(v) for k,v in eligible.dcpSurfaceWaters.fillna("MISSING").value_counts().items()},
     "2018_original_reporter_code_counts":{
       str(k):int(v) for k,v in eligible["repCode_plant"].fillna("MISSING").value_counts().items()},
     "2018_all_30_report_years_explicitly_verified":True,
     "2018_source_measured_mass_not_calculated_or_estimated":True,
    }
    return eligible,quality

def river_source(path):
    if sha(path)!=RIVER_SHA:raise ValueError("HydroRIVERS original previously pinned third-party mirror changed")
    g=gpd.read_parquet(path).to_crs("EPSG:4326")
    g=g.rename(columns={x:x.upper() for x in g.columns if x!=g.geometry.name})
    g=g.cx[BBOX[0]:BBOX[2],BBOX[1]:BBOX[3]].copy()
    for c in ("HYRIV_ID","NEXT_DOWN","LENGTH_KM","DIS_AV_CMS","UPLAND_SKM"):
        g[c]=pd.to_numeric(g[c],errors="coerce")
    g=g[g.HYRIV_ID.notna()].drop_duplicates("HYRIV_ID")
    rivers={int(x.HYRIV_ID):(int(x.NEXT_DOWN) if pd.notna(x.NEXT_DOWN) else 0,
                             float(x.LENGTH_KM))
            for x in g[["HYRIV_ID","NEXT_DOWN","LENGTH_KM"]].itertuples(index=False)}
    return g,rivers

def compute(args):
    tables,prior=load_original_source(args.historical_zip,args.historical_qc)
    fixed=json.loads(args.original_sites_qc.read_text())
    if fixed["independent_year_fixed"]!="2018" or fixed["pre_2018_Biology_EQR_values_read"]:
        raise ValueError("Original independent 2018 biological sites not outcome-blind")
    sites=pd.read_csv(args.original_sites,dtype="string",low_memory=False)
    if len(sites)!=7034 or sites.duplicated(["site_key","observedPropertyDeterminandBiologyEQRCode"]).any():
        raise ValueError("Original 2018 BQE source-only geography changed")
    ww,sourceqc=source_outfalls(tables)
    rivers,graph=river_source(args.river_parquet)
    ww=ww.reset_index(drop=True)
    ww["source_row_id"]=np.arange(len(ww),dtype=int)
    pts=eea.gdf_points(ww,"dcpLatitude","dcpLongitude")
    snapped=eea.snap(pts,rivers,EXACT_SNAP_METERS)
    snapped=snapped[snapped.HYRIV_ID.notna()].copy()
    snapped["HYRIV_ID"]=snapped.HYRIV_ID.astype("int64")
    snapped=snapped[snapped.HYRIV_ID.isin(graph)].copy()
    if len(snapped)>len(ww):raise ValueError("Duplicate source rows after nearest HydroRIVERS join")
    totals,byreach,links=eea.accum_upstream(snapped,graph)
    # The original 2018 biological data are source metadata without any EQR.
    site=sites.drop_duplicates("site_key").copy()
    site["lat"]=pd.to_numeric(site.lat,errors="raise")
    site["lon"]=pd.to_numeric(site.lon,errors="raise")
    site=site.reset_index(drop=True)
    site["source_row_id"]=np.arange(len(site),dtype=int)
    geom=eea.snap(eea.gdf_points(site,"lat","lon"),rivers,EXACT_SNAP_METERS)
    geom=geom[geom.HYRIV_ID.notna()].copy()
    geom["HYRIV_ID"]=geom.HYRIV_ID.astype("int64")
    geom=geom[geom.HYRIV_ID.isin(graph)].copy()
    river_meta=rivers[["HYRIV_ID","DIS_AV_CMS","UPLAND_SKM"]].drop_duplicates("HYRIV_ID")
    geom=geom.merge(river_meta,on="HYRIV_ID",how="left",validate="many_to_one")
    if geom.site_key.duplicated().any():raise ValueError("Repeated 2018 biological monitoring site")
    for nutrient,idx in (("N",1),("P",2)):
        geom["2018_measured_upstream_"+nutrient+"_t_per_year"]=geom.HYRIV_ID.map(
          lambda k:totals.get(int(k),[0,0.,0,0])[idx])
    geom["2018_upstream_qualifying_WWTP_count"]=geom.HYRIV_ID.map(
      lambda k:totals.get(int(k),[0,0.,0,0])[0])
    placebo=[eea.downstream_placebo(int(k),graph,byreach) for k in geom.HYRIV_ID]
    geom["2018_downstream_WWTP_count_placebo"]=[v[0] for v in placebo]
    geom["2018_measured_downstream_N_t_per_year_placebo"]=[v[1] for v in placebo]
    geom["2018_measured_downstream_P_t_per_year_placebo"]=[v[2] for v in placebo]
    join=sites.merge(geom[["site_key","HYRIV_ID","snap_distance_m",
                 "DIS_AV_CMS","UPLAND_SKM",
                 "2018_measured_upstream_N_t_per_year","2018_measured_upstream_P_t_per_year",
                 "2018_upstream_qualifying_WWTP_count",
                 "2018_downstream_WWTP_count_placebo",
                 "2018_measured_downstream_N_t_per_year_placebo",
                 "2018_measured_downstream_P_t_per_year_placebo"]],
                 on="site_key",how="inner",validate="many_to_one")
    rows={}
    for code,name in BQE_CODE.items():
        f=join[join.observedPropertyDeterminandBiologyEQRCode.eq(code)]
        covered=f[f["2018_upstream_qualifying_WWTP_count"]>0]
        rows[name]={
            "official_2018_BQE_pre_outcome_geolocated_sites":int(
                (sites.observedPropertyDeterminandBiologyEQRCode.eq(code)).sum()),
            "2018_BQE_monitoring_sites_snap_to_river_500m":len(f),
            "upstream_2018_true_measured_both_NP_source_sites":len(covered),
            "source_positive_countries":int(covered.countryCode.nunique()),
            "all_monitoring_countries":int(f.countryCode.nunique()),
            "3degree_geographic_blocks":int(f.region_3deg_block.nunique()),
            "source_positive_countries_counts":{
                str(k):int(v) for k,v in covered.countryCode.value_counts().items()},
            "source_zero_unregistered_50km_sites":int(f["2018_upstream_qualifying_WWTP_count"].eq(0).sum()),
            "source_sufficient_before_opening_2018_EQR":bool(
                len(f)>=100 and len(covered)>=20 and
                f.countryCode.nunique()>=5 and covered.countryCode.nunique()>=3
                and f.region_3deg_block.nunique()>=10)
        }
        print("GESUI_2018_INDEPENDENT_WWTP_BQE_PREOUTCOME_GEOGRAPHY",
              json.dumps({"BQE":name,**rows[name]},sort_keys=True),flush=True)
    args.out.mkdir(parents=True,exist_ok=True)
    keep=["site_key","countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme",
          "observedPropertyDeterminandBiologyEQRCode","procedureClassificationSystem",
          "phenomenonTimeReferenceYear","resultObservationStatus",
          "waterBodyIdentifier","surfaceWaterBodyTypeCode","lat","lon",
          "region_3deg_block","pre_2018_source_only_spatial_fold",
          "HYRIV_ID","DIS_AV_CMS","UPLAND_SKM","snap_distance_m",
          "2018_measured_upstream_N_t_per_year","2018_measured_upstream_P_t_per_year",
          "2018_upstream_qualifying_WWTP_count",
          "2018_downstream_WWTP_count_placebo",
          "2018_measured_downstream_N_t_per_year_placebo",
          "2018_measured_downstream_P_t_per_year_placebo"]
    join[keep].to_csv(args.out/"2018_SOURCE_ONLY_BQE_WWTP_NP_UPSTREAM_AND_NEGATIVE_CONTROL.csv",
                      index=False,float_format="%.12g")
    report={
      "schema":"gesui_2018_independent_year_real_measured_UWWTD_WWTP_Biology_River_Source_v01",
      "source_state":"SOURCE_FROZEN_BEFORE_2018_NUMERIC_EQR_OUTCOMES",
      "timestamp_utc":datetime.now(timezone.utc).isoformat(),
      "original_2018_UWWTD_report_qc_run":SOURCE_2018_QC_RUN,
      "original_2018_BQE_source_only_run":SITE_2018_SOURCE_RUN,
      "source_original_archive_data_sha256":{
         x:prior["source_table_original_SHA256"][x]["source_uncompressed_SHA256"]
         for x in SOURCE_NAMES},
      "original_fixed_2018_BQE_source_csv_sha256":sha(args.original_sites),
      "original_fixed_2018_biology_source_metadata_sha256":sha(args.original_sites_qc),
      "2018_annual_WWTP_measured_both_NP_source":sourceqc,
      "source_EEA_Fw_2018_outfalls_snapped_500m":len(snapped),
      "source_2018_unique_biological_site_keys_snapped_500m":int(geom.site_key.nunique()),
      "source_2018_BQE_site_element_rows_after_river_graph":len(join),
      "source_2018_WWTP_reported_load_directional_graph_links":links,
      "source_2018_BQE":rows,
      "2018_EEA_Biology_normalized_EQR_values_read":False,
      "2018_original_WWTP_measured_mass_units":"tons per year (verify original EEA Art15 measured field documentation)",
      "source_Europe_HydroRIVERS_SHA256":RIVER_SHA,
      "official_source_HydroRIVERS_independent_validation":False,
      "source_zero_only_registered_2018_qualifying_WWTPs":True,
      "caveats":[
       "WWTP measured annual source mass is not an independently measured delivered pollutant dose at BQE monitoring station",
       "Old archived UWWTD full 2018 annual reporting differs by nation and source selection from 2022",
       "Only one FW outlet per reporting plant; other reported source load kinds and multiple-outlet plants omitted",
       "2018 biological sites differ from 2022 in country coverage, data quality and river-reach geometry",
       "Third-party HydroRIVERS Europe v1 derived GeoParquet retains original 2022 study quality limits",
       "No 2018 measured EQR numeric data accessed or used in spatial site/exposure selection",
      ]}
    (args.out/"2018_MEASURED_SOURCE_BIOLOGY_COHORT_QC.json").write_text(
      json.dumps(report,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GESUI_2018_NP_MEASURED_TEMPORAL_SOURCE_READY",
      json.dumps({"verified_source_year":2018,
        "qualifying_both_measured_2018_plant_n":EXPECTED_HISTORICAL_MEASURED_BOTH,
        "direct_FW_source_outlets_in_bbox":len(ww),
        "snap_source_outlets":len(snapped),
        "biology_group_support":{n:{"total":v["2018_BQE_monitoring_sites_snap_to_river_500m"],
             "upstream_real_measured":v["upstream_2018_true_measured_both_NP_source_sites"],
             "gate":v["source_sufficient_before_opening_2018_EQR"]}
             for n,v in rows.items()},
        "numeric_2018_EQR_opened":False},sort_keys=True),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--historical-zip",type=Path,required=True)
    p.add_argument("--historical-qc",type=Path,required=True)
    p.add_argument("--original-sites",type=Path,required=True)
    p.add_argument("--original-sites-qc",type=Path,required=True)
    p.add_argument("--river-parquet",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    compute(p.parse_args())
