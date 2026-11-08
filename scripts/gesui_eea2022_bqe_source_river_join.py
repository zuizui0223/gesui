#!/usr/bin/env python3
"""Source-only crosswalk: official EEA 2022 WWTP measured N/P mass -> 2022 river BQE sites.

EEA Waterbase Biology 2026 source tables are read WITHOUT the numerical
normalized EQR or biological status outcome columns. Only monitoring-site IDs,
2022 observation existence, BQE taxonomy/method and coordinates are read.
Official EEA UWWTD 2022 measured annual source loads are joined to unique FW
outlets, projected nearest HydroRIVERS and directed upstream <=50 km.

No biological EQR effect or sign is evaluated by this source gate.
"""
from __future__ import annotations
import argparse,hashlib,io,json,math,sys,zipfile
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import geopandas as gpd
import numpy as np,pandas as pd
import gesui_eea_nutrient_load_river_crosswalk as eea

BIO_SHA={
 "EQR":"f11001a57172b618bee95f58fd0ec51c8e71be7b9c01a70927dde21b7ca6e2a9",
 "SPATIAL":"f6e2448ca2d99f528c21cdb9496673711b1d349e04576cf9024c6165c9a52008",
}
PREFIX="eea_t_waterbase-biology-2026_p_1990-2024_v01_r00/"
EQR=PREFIX+"Waterbase_v2025_1_T_WISE2_BiologyEQRData.csv"
SPATIAL=PREFIX+"Waterbase_v2025_1_S_WISE2_SpatialObject_DerivedData.csv"
PRIMARY_YEAR="2022"
BQE_CODES={
 "EEA_13-01-4":"benthic_invertebrates",
 "EEA_124-04-9":"phytobenthos",
 "EEA_14-02-8":"fish",
}
RIVER_SHA="5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f"
MAX_SNAP_M=500
LIMIT_KM=50
BBOX=(-12,44,26,64)

def uncompressed_digest(arc,name):
    h=hashlib.sha256()
    with arc.open(name) as fd:
        for block in iter(lambda:fd.read(1024*1024),b""):
            h.update(block)
    return h.hexdigest()

def source_keys(df):
    cols=("countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme")
    for c in cols:
        if c not in df.columns:raise ValueError("Official biology monitoring ID field missing")
        df[c]=df[c].astype("string").str.strip()
    df["site_key"]=df[list(cols)].fillna("").agg("|".join,axis=1)
    df.loc[df.monitoringSiteIdentifier.isna()|df.monitoringSiteIdentifier.eq(""),"site_key"]=pd.NA
    return df

def read_fixed_biological_sites(path):
    with zipfile.ZipFile(path) as z:
        for name,hash_expected in ((EQR,BIO_SHA["EQR"]),(SPATIAL,BIO_SHA["SPATIAL"])):
            if uncompressed_digest(z,name)!=hash_expected:
                raise ValueError("EEA Biology 2026 uncompressed scientific source SHA mismatch")
        # IMPORTANT: no resultEQRValue, resultNormalisedEQRValue,
        # resultEcologicalStatusClassValue or other biological response columns
        # are accessed in the source-only match.
        needed=["countryCode","monitoringSiteIdentifier",
            "monitoringSiteIdentifierScheme","parameterWaterBodyCategory",
            "observedPropertyDeterminandBiologyEQRCode",
            "observedPropertyDeterminandLabel","procedureClassificationSystem",
            "phenomenonTimeReferenceYear","resultObservationStatus"]
        e=pd.read_csv(z.open(EQR),usecols=needed,dtype="string",low_memory=False)
        p=pd.read_csv(z.open(SPATIAL),dtype="string",low_memory=False)
    e=source_keys(e);p=source_keys(p)
    e=e[e.phenomenonTimeReferenceYear.eq(PRIMARY_YEAR)&
        e.parameterWaterBodyCategory.eq("RW")&
        e.observedPropertyDeterminandBiologyEQRCode.isin(BQE_CODES)].copy()
    original_source_rows=len(e)
    for field in ("lat","lon"):p[field]=pd.to_numeric(p[field],errors="coerce")
    p["has_coords"]=p.lat.between(-90,90)&p.lon.between(-180,180)
    p["coord_signature"]=p[["lat","lon"]].astype("string").fillna("MISSING").agg("|".join,axis=1)
    p=p.dropna(subset=["site_key"])
    groups=p.groupby("site_key").agg(location_count=("coord_signature","nunique"),
                                     geometry_present=("has_coords","any"),
                                     waterbody_count=("waterBodyIdentifier","nunique"))
    unique=set(groups.index[groups.location_count.eq(1)&groups.geometry_present])
    sites=p[p.site_key.isin(unique)].drop_duplicates("site_key").copy()
    sites=sites[sites.lat.between(BBOX[1],BBOX[3])&sites.lon.between(BBOX[0],BBOX[2])]
    # Geometry and identifiers only; not even the biological score exists
    # in this table. No site selection from sign/magnitude of EQR.
    meta=e.merge(sites[["site_key","countryCode","lat","lon","waterBodyIdentifier",
                        "waterBodyIdentifierScheme","surfaceWaterBodyTypeCode"]],
                 how="inner",on=["site_key","countryCode"],validate="many_to_one")
    quality={"original_2022_rw_fixed_codes_rows":original_source_rows,
             "2022_rw_source_rows_with_complete_official_spatial_reference":len(meta),
             "distinct_RW_site_keys_with_coords":int(meta.site_key.nunique()),
             "distinct_RW_waterbodies_metadata":int(meta.waterBodyIdentifier.nunique()),
             "source_EQR_columns_loaded_only":needed,
             "original_BQE_codes_selected":BQE_CODES}
    return meta,quality

def source_eea_to_reaches(river_gdf,river):
    w,pages_w=eea.download_layer(eea.PLANT,eea.BOTH_WHERE,eea.PLANT_FIELDS,5492)
    d,pages_d=eea.download_layer(eea.OUTFALL,"1=1",eea.OUT_FIELDS,24416)
    w=eea.normalize_key(w)
    d=eea.normalize_key(d)
    w=w.dropna(subset=["compound_join_key"])
    d=d.dropna(subset=["compound_join_key"])
    duplications=w.compound_join_key.value_counts()
    w=w[w.compound_join_key.map(duplications).eq(1)].copy()
    counts=d.compound_join_key.value_counts()
    w["outlets_per_plant"]=w.compound_join_key.map(counts).fillna(0).astype(int)
    unique=w[w.outlets_per_plant.eq(1)].merge(
        d.drop_duplicates("compound_join_key"),
        on="compound_join_key",how="inner",suffixes=("_treatment","_dcp"),validate="one_to_one")
    for field in ("dcpLatitude","dcpLongitude","uwwPDischargeMeasured","uwwNDischargeMeasured"):
        unique[field]=pd.to_numeric(unique[field],errors="coerce")
    unique=unique[unique.dcpWaterBodyType.astype("string").str.strip().str.upper().eq("FW")]
    unique=unique[unique.dcpLatitude.between(BBOX[1],BBOX[3])&
                  unique.dcpLongitude.between(BBOX[0],BBOX[2])]
    unique=unique.reset_index(drop=True)
    unique["source_row_id"]=np.arange(len(unique),dtype=int)
    matched=eea.snap(eea.gdf_points(unique,"dcpLatitude","dcpLongitude"),
                      river_gdf,MAX_SNAP_M)
    matched=matched[matched.HYRIV_ID.notna()].copy()
    matched["HYRIV_ID"]=matched.HYRIV_ID.astype("int64")
    matched=matched[matched.HYRIV_ID.isin(river)].copy()
    upstream,sourcebyreach,edges=eea.accum_upstream(matched,river)
    qc={"source_original_2022_both_measured_plants":5492,
        "source_2022_dcp_all":24416,
        "unique_fw_dcp_in_bbox":len(unique),
        "unique_fw_source_snapped_to_river_500m":len(matched),
        "source_traversed_graph_reach_links":edges,
        "original_source_plant_pages":len(pages_w),
        "original_dcp_pages":len(pages_d),
        "EEA_UPSTREAM_2022_mass_tonnes_per_year":True,
        "outlet_ambiguous_multiple_points_excluded":True}
    return upstream,sourcebyreach,qc

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--biology-zip",type=Path,required=True)
    ap.add_argument("--river-parquet",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    meta,bioqc=read_fixed_biological_sites(a.biology_zip)
    if eea.sha256(a.river_parquet)!=RIVER_SHA:raise ValueError("Fixed 2022 river source sha changed")
    g=gpd.read_parquet(a.river_parquet).to_crs("EPSG:4326")
    g=g.rename(columns={x:x.upper() for x in g.columns if x!=g.geometry.name})
    g=g.cx[BBOX[0]:BBOX[2],BBOX[1]:BBOX[3]].copy()
    for c in ("HYRIV_ID","NEXT_DOWN","LENGTH_KM"):
        g[c]=pd.to_numeric(g[c],errors="coerce")
    g=g[g.HYRIV_ID.notna()].drop_duplicates("HYRIV_ID")
    river={int(r.HYRIV_ID):(int(r.NEXT_DOWN) if pd.notna(r.NEXT_DOWN) else 0,
                             float(r.LENGTH_KM))
           for r in g[["HYRIV_ID","NEXT_DOWN","LENGTH_KM"]].itertuples(index=False)}
    upstream,byreach,eeaqc=source_eea_to_reaches(g,river)
    station=meta.drop_duplicates("site_key").copy().reset_index(drop=True)
    station["source_row_id"]=np.arange(len(station),dtype=int)
    spatial=eea.snap(eea.gdf_points(station,"lat","lon"),g,MAX_SNAP_M)
    spatial=spatial[spatial.HYRIV_ID.notna()].copy()
    spatial["HYRIV_ID"]=spatial.HYRIV_ID.astype("int64")
    spatial=spatial[spatial.HYRIV_ID.isin(river)].copy()
    # Only source-published hydrologic predictors, before EQR outcomes.
    for col in ("DIS_AV_CMS","UPLAND_SKM"):
        g[col]=pd.to_numeric(g[col],errors="coerce")
    river_meta=g[["HYRIV_ID","DIS_AV_CMS","UPLAND_SKM"]].drop_duplicates("HYRIV_ID")
    spatial=spatial.merge(river_meta,on="HYRIV_ID",how="left",validate="many_to_one")
    for nutrient,idx in (("N",1),("P",2)):
        spatial[f"upstream_2022_EEA_{nutrient}_measured_t_per_yr"]=spatial.HYRIV_ID.map(
          lambda x:upstream.get(int(x),[0,0.,0.,0])[idx])
    spatial["upstream_2022_unique_EEA_both_load_WWTP_count"]=spatial.HYRIV_ID.map(
      lambda x:upstream.get(int(x),[0,0.,0.,0])[0])
    downstream=[eea.downstream_placebo(int(x),river,byreach) for x in spatial.HYRIV_ID]
    spatial["downstream_2022_WWTP_count_placebo"]=[v[0] for v in downstream]
    spatial["downstream_2022_EEA_N_t_per_yr_placebo"]=[v[1] for v in downstream]
    spatial["downstream_2022_EEA_P_t_per_yr_placebo"]=[v[2] for v in downstream]
    spatial["region_3deg_block"]=(np.floor(spatial.lat.astype(float)/3).astype("int64").astype(str)+"_"+
                                  np.floor(spatial.lon.astype(float)/3).astype("int64").astype(str))
    matched=meta.merge(spatial[["site_key","HYRIV_ID","snap_distance_m","region_3deg_block",
       "DIS_AV_CMS","UPLAND_SKM",
       "upstream_2022_EEA_N_measured_t_per_yr","upstream_2022_EEA_P_measured_t_per_yr",
       "upstream_2022_unique_EEA_both_load_WWTP_count",
       "downstream_2022_WWTP_count_placebo",
       "downstream_2022_EEA_N_t_per_yr_placebo","downstream_2022_EEA_P_t_per_yr_placebo"]],
       on="site_key",how="inner",validate="many_to_one")
    # Fixed outcome-blind spatial folds: assign whole 3degree regions by
    # descending sample count then region ID; choose currently least-loaded
    # fold and resolve ties by fold number. Independent of EQR values and
    # of scikit-learn GroupKFold version-specific tie-breaking.
    matched["fixed_spatial_cv_fold"]=0
    folds_manifest={}
    for code,label in BQE_CODES.items():
        mask=matched.observedPropertyDeterminandBiologyEQRCode.eq(code)
        size=matched.loc[mask].groupby("region_3deg_block").size().to_dict()
        bins=[0,0,0,0,0]
        assignments={}
        for region,n in sorted(size.items(),key=lambda x:(-x[1],str(x[0]))):
            fold=min(range(5),key=lambda idx:(bins[idx],idx))
            assignments[str(region)]=fold+1
            bins[fold]+=int(n)
        matched.loc[mask,"fixed_spatial_cv_fold"]=matched.loc[mask,"region_3deg_block"].map(assignments)
        folds_manifest[label]={"region_to_fold":assignments,"test_sites_per_fold":bins,
           "from_biological_response_values":False}
    stats={}
    for code,label in BQE_CODES.items():
        rows=matched[matched.observedPropertyDeterminandBiologyEQRCode.eq(code)].copy()
        dedup=rows.drop_duplicates("site_key")
        pos=dedup[dedup.upstream_2022_unique_EEA_both_load_WWTP_count.gt(0)]
        src_countries=int(pos.countryCode.nunique())
        ctrys=int(dedup.countryCode.nunique())
        regions=int(dedup.region_3deg_block.nunique())
        bqe={
            "biological_quality_element":label,
            "2022_original_source_2022_RW_rows_with_valid_station_coordinates":
                int((meta.observedPropertyDeterminandBiologyEQRCode==code).sum()),
            "2022_BQE_sites_snapped_to_RIVER_within_500m":len(dedup),
            "2022_BQE_rows_snapped":len(rows),
            "2022_distinct_countries":ctrys,
            "2022_upstream_EEA2022_both_measured_FW_source_sites":len(pos),
            "2022_positive_EEA_source_countries":src_countries,
            "2022_region_3degree_blocks":regions,
            "source_2022_downstream_placebo_sites":int(dedup.downstream_2022_WWTP_count_placebo.gt(0).sum()),
            "source_2022_zero_registered_upstream_sites":int((dedup.upstream_2022_unique_EEA_both_load_WWTP_count==0).sum()),
            "source_dated_2022_spatial_association_gate_pass":bool(
              len(dedup)>=100 and len(pos)>=20 and ctrys>=5 and src_countries>=3 and regions>=10),
            "primary_thresholds_predeclared":{
                "min_monitoring_sites":100,"min_upstream_exposed":20,
                "min_reporting_countries":5,"min_upstream_source_countries":3,
                "min_3degree_blocks":10},
            "source_positive_country_counts":{str(k):int(v) for k,v in pos.countryCode.value_counts().items()},
            "source_all_country_counts":{str(k):int(v) for k,v in dedup.countryCode.value_counts().items()},
            "classification_procedures_with_source":int(pos.procedureClassificationSystem.nunique()),
        }
        stats[label]=bqe
        print("GESUI_NEW_EEA_BIOLOGICAL_SOURCE_READY",json.dumps(bqe,sort_keys=True),flush=True)
    report={
        "schema":"gesui_EEA_Biology_2026_EEA2022_FW_loads_hydrorivers_2022_BQE_source_only_v01",
        "status":"SOURCE_ONLY_NO_BIOLOGICAL_SCORE_VALUES_ACCESSED",
        "timestamp_utc":datetime.now(timezone.utc).isoformat(),
        "official_Biology_2026_EQR_SHA256":BIO_SHA["EQR"],
        "official_Biology_2026_Spatial_SHA256":BIO_SHA["SPATIAL"],
        "EEA2022_UWWTD_measured_NP_Fw_source":True,
        "HydroRIVERS_derived_Europe_SHA256":RIVER_SHA,
        "river_original_official_topology_independently_verified":False,
        "source_2022_BQE_and_spatial_coverage":bioqc,
        "EEA2022_official_FW_source_join":eeaqc,
        "2022_biological_element_monitoring_source_support":stats,
       "fixed_2022_geographic_fold_manifest":folds_manifest,
        "no_normalized_EQR_values_loaded_or_used":True,
        "no_GBIF_plant_biodiversity_values_loaded":True,
        "future_EQR_primary_year":"2022",
        "source_zero_is_no_reported_both_positive_direct_FW_WWTP_NOT_pollution_free":True,
        "caveats":[
          "EEA biological site code is official but nearest HydroRIVERS reach is still a modeled link requiring catchment validation",
          "Actual 2022 EQR biological class values were not read and cannot be inferred from source station counts",
          "EEA WWTP measured load data only covers both-positive N/P sources with unique FW discharge point",
          "Biological classification procedures differ by country and method; do not fit pooled raw EQR without standardization",
          "No measured agricultural nutrient flux, hydrological flow seasonality or sewage-specific organism dispersal in this gate"
        ]}
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"BQE_SOURCE_EXPOSURE_COVERAGE.json").write_text(
       json.dumps(report,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    keep=["countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme",
      "site_key","parameterWaterBodyCategory","observedPropertyDeterminandBiologyEQRCode",
      "observedPropertyDeterminandLabel","procedureClassificationSystem",
      "phenomenonTimeReferenceYear","resultObservationStatus","lat","lon",
      "waterBodyIdentifier","surfaceWaterBodyTypeCode","HYRIV_ID","snap_distance_m",
      "region_3deg_block","fixed_spatial_cv_fold","DIS_AV_CMS","UPLAND_SKM",
      "upstream_2022_EEA_N_measured_t_per_yr",
      "upstream_2022_EEA_P_measured_t_per_yr","upstream_2022_unique_EEA_both_load_WWTP_count",
      "downstream_2022_WWTP_count_placebo",
      "downstream_2022_EEA_N_t_per_yr_placebo","downstream_2022_EEA_P_t_per_yr_placebo"]
    matched[keep].to_csv(a.out/"EQR_2022_river_source_only_station_registry.csv",
                         index=False,float_format="%.9g")
    print("GESUI_EEA_BQE_SOURCE_ONLY_FINAL",json.dumps({
       "2022_biological_sites_before_river_snap":bioqc["distinct_RW_site_keys_with_coords"],
       "2022_biological_rows_after_graph_join":len(matched),
       "2022_biological_stations_after_graph_join":int(matched.site_key.nunique()),
       "bqe_gate":{"name":{k:v for k,v in q.items()
           if k.startswith("2022_") or k=="source_2022_zero_registered_upstream_sites"}
           for name,q in stats.items()},
       "no_BQE_outcomes_examined":True},sort_keys=True),flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--biology-zip",type=Path,required=True)
    ap.add_argument("--river-parquet",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    main()
