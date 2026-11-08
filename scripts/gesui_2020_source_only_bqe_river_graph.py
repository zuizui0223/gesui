#!/usr/bin/env python3
"""Outcome-blind 2020 UWWTD reporter-year N/P loads to directed EU river BQE sites.

Uses ONLY the independently downloaded 2019–2020 public original GPKGs'
exact 2020 reporter-code + unique FW outlet source table, the previously
frozen 2020 Biology *site metadata without scores*, historical 2018/2022
source-only method catalogs, and SHA-locked HydroRIVERS EU graph.
No 2020 numeric biological EQR is opened. This is a PRE-OUTCOME source gate.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
from collections import Counter
import geopandas as gpd
import numpy as np,pandas as pd
import gesui_eea_nutrient_load_river_crosswalk as eea

YEAR="2020"
RIVER_SHA="5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f"
BBOX=(-12.,44.,26.,64.)
SNAP_M=500
BQE={"EEA_13-01-4":"benthic_invertebrates",
     "EEA_124-04-9":"phytobenthos",
     "EEA_14-02-8":"fish"}

def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for ch in iter(lambda:f.read(1024*1024),b""):h.update(ch)
    return h.hexdigest()

def data_only(args):
    source_qc=json.loads(args.source_qc.read_text())
    expected=source_qc["source_candidate_file_sha256"]
    if sha(args.source_candidates)!=expected:
        raise ValueError("Original 2020 public WWTP/outfall source candidate SHA changed")
    if not source_qc.get("official_2020_reporter_suffix_enforced_in_candidate_CSV"):
        raise ValueError("2022 reporting-code contamination was not excluded")
    ww=pd.read_csv(args.source_candidates,dtype={"source_key":str,"repCode":str,"uwwCode":str},low_memory=False)
    if len(ww)!=source_qc["official_both_measured_positive_NP_unique_FW_outlets_in_study_bbox"]:
        raise ValueError("Original public 2020 candidate inventory count changed")
    if ww.source_key.duplicated().any():raise ValueError("Duplicate official full plant ID")
    if not ww.repCode.astype(str).str.strip().str.endswith("_2020").all():
        raise ValueError("Non-2020 reporting code in historical source")
    if not ww.dcpWaterBodyType.eq("FW").all():raise ValueError("Nonfreshwater source")
    ww["uwwNDischargeMeasured"]=pd.to_numeric(
      ww["uwwNDischargeMeasured_t_per_year"],errors="raise")
    ww["uwwPDischargeMeasured"]=pd.to_numeric(
      ww["uwwPDischargeMeasured_t_per_year"],errors="raise")
    if not (ww[["uwwNDischargeMeasured","uwwPDischargeMeasured"]].to_numpy(float)>0).all():
        raise ValueError("Invalid reported measured annual source mass")
    for k in ("dcpLatitude","dcpLongitude"):
        ww[k]=pd.to_numeric(ww[k],errors="raise")
    if not (ww.dcpLatitude.between(BBOX[1],BBOX[3]) &
            ww.dcpLongitude.between(BBOX[0],BBOX[2])).all():
        raise ValueError("2020 source outside fixed EU region")

    site_qc=json.loads(args.site_qc.read_text())
    if site_qc["independent_year_fixed"]!="2020" or site_qc["pre_2020_Biology_EQR_values_read"]:
        raise ValueError("2020 station cohort not frozen before outcome")
    sites=pd.read_csv(args.frozen_sites,dtype="string",low_memory=False)
    if len(sites)!=site_qc["source_only_geolocated_2020_site_element_rows"]:
        raise ValueError("Original 2020 site metadata count changed")
    if sites.duplicated(["site_key","observedPropertyDeterminandBiologyEQRCode"]).any():
        raise ValueError("2020 original station has multiple element observations")

    # Only SOURCE-ONLY archives for past years; no EQR outcome ever loaded.
    method_cols=["countryCode","procedureClassificationSystem","observedPropertyDeterminandBiologyEQRCode"]
    past18=pd.read_csv(args.methods2018,usecols=method_cols,dtype="string",low_memory=False)
    past22=pd.read_csv(args.methods2022,usecols=method_cols,dtype="string",low_memory=False)
    keycols=["countryCode","procedureClassificationSystem"]
    catalog={}
    for year,past in [("2018",past18),("2022",past22)]:
        catalog[year]={}
        for code in BQE:
            x=past.loc[past.observedPropertyDeterminandBiologyEQRCode.eq(code),keycols].fillna("UNKNOWN").astype(str)
            catalog[year][code]=set(x.countryCode+"|"+x.procedureClassificationSystem)
    return ww,sites,catalog,site_qc,source_qc

def river_graph(path):
    if sha(path)!=RIVER_SHA:raise ValueError("Changed HydroRIVERS EU parquet original mirror")
    rivers=gpd.read_parquet(path).to_crs("EPSG:4326")
    rivers=rivers.rename(columns={x:x.upper() for x in rivers.columns if x!=rivers.geometry.name})
    rivers=rivers.cx[BBOX[0]:BBOX[2],BBOX[1]:BBOX[3]].copy()
    for x in ("HYRIV_ID","NEXT_DOWN","LENGTH_KM","DIS_AV_CMS","UPLAND_SKM"):
        rivers[x]=pd.to_numeric(rivers[x],errors="coerce")
    rivers=rivers[rivers.HYRIV_ID.notna()].drop_duplicates("HYRIV_ID")
    graph={int(x.HYRIV_ID):(int(x.NEXT_DOWN) if pd.notna(x.NEXT_DOWN) else 0,float(x.LENGTH_KM))
           for x in rivers[["HYRIV_ID","NEXT_DOWN","LENGTH_KM"]].itertuples(index=False)}
    return rivers,graph

def main():
    parser=argparse.ArgumentParser()
    for name in ("source-candidates","source-qc","frozen-sites","site-qc",
                 "methods2018","methods2022","river-parquet","out"):
        parser.add_argument("--"+name,required=True,type=Path)
    a=parser.parse_args()
    ww,sites,catalog,siteqc,sourceqc=data_only(a)
    rivers,graph=river_graph(a.river_parquet)

    ww=ww.reset_index(drop=True)
    ww["source_row_id"]=np.arange(len(ww),dtype=int)
    source_nearest=eea.snap(eea.gdf_points(ww,"dcpLatitude","dcpLongitude"),rivers,SNAP_M)
    source_nearest=source_nearest[source_nearest.HYRIV_ID.notna()].copy()
    source_nearest["HYRIV_ID"]=source_nearest.HYRIV_ID.astype("int64")
    source_nearest=source_nearest[source_nearest.HYRIV_ID.isin(graph)].copy()
    if source_nearest.source_key.duplicated().any():raise ValueError("Duplicate source snapped")
    # Graph helper excludes source and receptor ordering within the same reach,
    # and uses directed NEXT_DOWN within <=50km.
    upstream,byreach,edges=eea.accum_upstream(source_nearest,graph)

    unique=sites.drop_duplicates("site_key").copy()
    for x in ("lat","lon"):unique[x]=pd.to_numeric(unique[x],errors="raise")
    unique=unique.reset_index(drop=True)
    unique["source_row_id"]=np.arange(len(unique))
    geom=eea.snap(eea.gdf_points(unique,"lat","lon"),rivers,SNAP_M)
    geom=geom[geom.HYRIV_ID.notna()].copy()
    geom["HYRIV_ID"]=geom.HYRIV_ID.astype("int64")
    geom=geom[geom.HYRIV_ID.isin(graph)].copy()
    meta=rivers[["HYRIV_ID","DIS_AV_CMS","UPLAND_SKM"]].drop_duplicates("HYRIV_ID")
    geom=geom.merge(meta,on="HYRIV_ID",how="left",validate="many_to_one")
    if geom.site_key.duplicated().any():raise ValueError("Duplicate biological site snapped")
    for nutrient,index in (("N",1),("P",2)):
        geom[f"2020_measured_upstream_{nutrient}_t_per_year"]=geom.HYRIV_ID.map(
            lambda k:upstream.get(int(k),[0,0.,0.,0])[index])
    geom["2020_upstream_qualifying_WWTP_count"]=geom.HYRIV_ID.map(
            lambda k:upstream.get(int(k),[0,0.,0.,0])[0])
    placebo=[eea.downstream_placebo(int(x),graph,byreach) for x in geom.HYRIV_ID]
    geom["2020_downstream_WWTP_count_placebo"]=[x[0] for x in placebo]
    geom["2020_measured_downstream_N_t_per_year_placebo"]=[x[1] for x in placebo]
    geom["2020_measured_downstream_P_t_per_year_placebo"]=[x[2] for x in placebo]

    values=["site_key","HYRIV_ID","snap_distance_m","DIS_AV_CMS","UPLAND_SKM",
            "2020_upstream_qualifying_WWTP_count","2020_downstream_WWTP_count_placebo",
            "2020_measured_upstream_N_t_per_year",
            "2020_measured_upstream_P_t_per_year",
            "2020_measured_downstream_N_t_per_year_placebo",
            "2020_measured_downstream_P_t_per_year_placebo"]
    joined=sites.merge(geom[values],on="site_key",how="inner",validate="many_to_one")
    counts={}
    for code,label in BQE.items():
        s=joined[joined.observedPropertyDeterminandBiologyEQRCode.eq(code)].copy()
        keys=s.countryCode.fillna("UNKNOWN").astype(str)+"|"+s.procedureClassificationSystem.fillna("UNKNOWN").astype(str)
        src=s["2020_upstream_qualifying_WWTP_count"].gt(0)
        src18=keys.isin(catalog["2018"][code])
        src22=keys.isin(catalog["2022"][code])
        both=src18 & src22
        either=src18 | src22
        exposed=s.loc[src]
        eligible=s.loc[src & either]
        counts[label]={
          "2020_original_site_elements_preoutcome":int(
            (sites.observedPropertyDeterminandBiologyEQRCode==code).sum()),
          "2020_sites_500m_nearest_river":len(s),
          "reported_2020_both_NP_upstream_source_sites":int(src.sum()),
          "registered_source_positive_countries":int(exposed.countryCode.nunique()),
          "regional_3deg_blocks":int(s.region_3deg_block.nunique()),
          "countries_total":int(s.countryCode.nunique()),
          "source_positive_by_country":{str(k):int(v) for k,v in
                     exposed.countryCode.value_counts().items()},
          "source_positive_exact_country_method_in_2018":int((src & src18).sum()),
          "source_positive_exact_country_method_in_2022":int((src & src22).sum()),
          "source_positive_exact_country_method_in_both":int((src & both).sum()),
          "source_positive_exact_country_method_in_either":int((src & either).sum()),
          "source_positive_compatible_countries_either":int(eligible.countryCode.nunique()),
          "2020_compatible_site_total":int(either.sum()),
          "predeclared_generalization_gate_50_source_positive_and_5_countries":bool(
              int((src & either).sum())>=50 and eligible.countryCode.nunique()>=5),
          "predeclared_2020_basic_EQR_model_coverage_gate":bool(
              len(s)>=100 and int(src.sum())>=20 and s.countryCode.nunique()>=5
              and exposed.countryCode.nunique()>=3 and s.region_3deg_block.nunique()>=10),
        }
        print("GESUI_2020_SOURCE_ONLY_RIVER_RECEPTOR_GATE",
               json.dumps({"BQE":label,**counts[label]},sort_keys=True),flush=True)
    a.out.mkdir(parents=True,exist_ok=True)
    keep=list(sites.columns)+[x for x in values if x!="site_key"]
    joined[keep].to_csv(a.out/"2020_SOURCE_ONLY_BQE_WWTP_NP_UPSTREAM_AND_NEGATIVE_CONTROL.csv",
              index=False,float_format="%.12g")
    csvpath=a.out/"2020_SOURCE_ONLY_BQE_WWTP_NP_UPSTREAM_AND_NEGATIVE_CONTROL.csv"
    receipt={
      "schema":"gesui_2020_reporter_year_source_only_BQE_graph_v01",
      "year":"2020","source_2020_reporter_code_verified":True,
      "individual_WWTP_measurement_dates_not_available":True,
      "source_and_biology_site_cohorts_frozen_without_2020_numeric_EQR":True,
      "source_candidate_SHA256":sha(a.source_candidates),
      "source_QC_SHA256":sha(a.source_qc),
      "frozen_2020_biology_site_metadata_SHA256":sha(a.frozen_sites),
      "frozen_2020_biology_site_QC_SHA256":sha(a.site_qc),
      "2018_source_only_method_table_SHA256":sha(a.methods2018),
      "2022_source_only_method_table_SHA256":sha(a.methods2022),
      "HydroRIVERS_EU_SHA256":RIVER_SHA,
      "source_original_public_GPKG_SHA256":{
          k:sourceqc["datasets"][k]["source"]["source_SHA256"] for k in ["WWTP","DCP"]},
      "source_only_2020_FW_both_NP_eligible_GPKG_plant_n":len(ww),
      "source_snapped_to_HydroRIVERS_500m":len(source_nearest),
      "source_snapping_median_m":float(source_nearest.snap_distance_m.median()) if len(source_nearest) else None,
      "source_2020_outlet_within_50km_downstream_graph_link_count":edges,
      "source_only_2020_distinct_biological_sites_snapped_to_river":int(geom.site_key.nunique()),
      "source_only_biological_site_element_rows":len(joined),
      "bqe_by_group":counts,
      "2020_numeric_EQR_values_read":False,
      "source_only_receptor_CSV_SHA256":sha(csvpath),
      "causal_not_identified":True,
      "registered_no_source_does_not_mean_no_sewage":True,
      "caveats":[
        "WWTP annual discharge masses are reported measured fields, individual collection dates not provided",
        "Different report-source inventories by country and year; 2020 means official reporter code suffix, not proven identical date of laboratory sample",
        "HydroRIVERS EU public third-party mirror not original official byte-equivalent proven",
        "Upstream river graph assumes connectivity and 50km cutoff, no within-reach ordering",
        "Missing registered source != confirmed sewage-free control",
        "No 2020 numeric biological EQR opened before source-site selection and fold freeze"
      ]
    }
    (a.out/"2020_SOURCE_ONLY_BQE_WWTP_NP_RIVER_QC.json").write_text(
       json.dumps(receipt,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_2020_BIOLOGICAL_EQR_SOURCE_ONLY_GRAPH_FREEZE_COMPLETED",json.dumps({
        "source_2020_n":len(ww),"source_snapped_n":len(source_nearest),
        "BQE_method_gate":{k:v["predeclared_generalization_gate_50_source_positive_and_5_countries"]
                           for k,v in counts.items()},
        "no_2020_EQR_scores":True},sort_keys=True),flush=True)

if __name__=="__main__":main()
