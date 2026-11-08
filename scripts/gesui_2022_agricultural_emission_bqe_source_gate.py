#!/usr/bin/env python3
"""EEA 2022 agricultural TN/TP emissions availability at frozen river BQE stations.

Independent source-only crosswalk: source EEA Emissions category NP1 is
Diffuse - Agricultural emissions (WISE SoE emissions code list). Output
does NOT fit ecological EQR, use any biological scores or rank sources by
their predictive sign. Annual masses are reported by RBD/subunit, not a local
monitoring-reach nutrient delivery measurement.
"""
from __future__ import annotations
import argparse,hashlib,json,zipfile
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import pandas as pd
import numpy as np

ARCH_ROOT="eea_t_waterbase-emissions_p_1977-2025_v01_r00/"
TABLE=ARCH_ROOT+"Waterbase_v2025_1_T_WISE1_Emissions.csv"
SPATIAL=ARCH_ROOT+"Waterbase_v2025_1_S_WISE6_SpatialObject_DerivedData.csv"
SOURCE_BQE_REGISTRY_SHA="e12df511d64c4c580ef2422dccdbbd11a1309c84d29d467fd4527206e41fa636"
CATEGORIES={
 "NP1":"agricultural_diffuse",
 "NP2":"atmospheric_deposition",
 "NP4":"urban_runoff"
}
PARAMETERS={"Total nitrogen":"N","Total phosphorus":"P"}
UNITS_TO_TONNES={"t/a":1.,"kg/a":.001}
SITE_CODES={"EEA_13-01-4":"benthic_invertebrates",
            "EEA_124-04-9":"phytobenthos",
            "EEA_14-02-8":"fish"}
VALID_SCHEMES={"euSubUnitCode":"subUnitIdentifier",
               "euRBDCode":"rbdIdentifier"}

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for x in iter(lambda:f.read(1024*1024),b""):h.update(x)
    return h.hexdigest()

def site_waterbody_metadata(original):
    keep=("countryCode","waterBodyIdentifier","subUnitIdentifier","rbdIdentifier")
    df=original[list(keep)].copy()
    for col in keep:df[col]=df[col].astype("string").str.strip()
    df=df.dropna(subset=["countryCode","waterBodyIdentifier"])
    bad={}
    for field in ("subUnitIdentifier","rbdIdentifier"):
        group=df.groupby(["countryCode","waterBodyIdentifier"])[field].nunique(dropna=True)
        bad[field]=int((group>1).sum())
    # Do not assign location metadata to waterbody IDs that map to >1 region.
    valid_keys=df.groupby(["countryCode","waterBodyIdentifier"]).filter(
        lambda group:all(group[x].dropna().nunique()<=1
                         for x in ("subUnitIdentifier","rbdIdentifier")))
    valid=valid_keys.drop_duplicates(["countryCode","waterBodyIdentifier"])
    if valid.duplicated(["countryCode","waterBodyIdentifier"]).any():
        raise ValueError("Nonunique waterbody source identity")
    return valid,bad

def normalized_emission_source(df):
    rows=df.loc[
       df.phenomenonTimeReferencePeriod.astype("string").str.strip().eq("2022")&
       df.observedPropertyDeterminandLabel.isin(PARAMETERS)&
       df.parameterEmissionsSourceCategory.isin(CATEGORIES)&
       df.spatialUnitIdentifierScheme.isin(VALID_SCHEMES)
    ].copy()
    for col in ("countryCode","spatialUnitIdentifier","spatialUnitIdentifierScheme",
                "parameterEmissionsSourceCategory","observedPropertyDeterminandLabel",
                "resultEmissionsUom","procedureEmissionsMethod"):
        rows[col]=rows[col].astype("string").str.strip()
    rows["mass_tonnes_2022"]=pd.to_numeric(rows.resultEmissionsValue,errors="coerce")
    factors=rows.resultEmissionsUom.map(UNITS_TO_TONNES)
    rows["original_annual_mass_tonnes"]=rows.mass_tonnes_2022*factors
    rows["numeric_mass_usable"]=np.isfinite(rows.original_annual_mass_tonnes.astype(float))&rows.original_annual_mass_tonnes.ge(0)
    # No pooling across reporting methods or hierarchy, and no choosing max.
    groupcols=("countryCode","spatialUnitIdentifierScheme","spatialUnitIdentifier",
               "parameterEmissionsSourceCategory","observedPropertyDeterminandLabel")
    repeats=rows.groupby(list(groupcols),dropna=False).size()
    ambiguous={tuple(map(str,k)) for k,v in repeats.items() if v>1}
    def ambiguous_row(x):
        return tuple(map(str,(x.countryCode,x.spatialUnitIdentifierScheme,
             x.spatialUnitIdentifier,x.parameterEmissionsSourceCategory,
             x.observedPropertyDeterminandLabel))) in ambiguous
    rows["duplicate_source_grain"]=rows.apply(ambiguous_row,axis=1)
    valid=rows[rows.numeric_mass_usable&~rows.duplicate_source_grain].copy()
    valid["country_source"]=valid.countryCode.astype(str)
    valid["grained_annual_source_tons"]=valid.original_annual_mass_tonnes.astype(float)
    return valid,rows,{
      "original_2022_NP1_NP2_NP4_totalN_totalP_rows_at_region_grain":len(rows),
      "distinct_original_multi_report_duplicate_keys":len(ambiguous),
      "unusable_numeric_or_missing_unit":int((~rows.numeric_mass_usable).sum()),
      "original_valid_unambiguous_annual_mass_rows":len(valid),
      "original_source_category_rows":{
         str(k):int(v) for k,v in rows.parameterEmissionsSourceCategory.value_counts().items()},
      "source_category_determinand_counts":{
         f"{cat}:{param}":int(n)
         for (cat,param),n in rows.groupby(["parameterEmissionsSourceCategory",
                                            "observedPropertyDeterminandLabel"]).size().items()},
      "source_category_method_counts":{
         f"{cat}:{method}":int(n)
         for (cat,method),n in rows.groupby(["parameterEmissionsSourceCategory",
                                            "procedureEmissionsMethod"],dropna=False).size().items()},
      "original_unit_counts":{str(k):int(v) for k,v in rows.resultEmissionsUom.value_counts().items()},
      "valid_annual_source_countries":{
         str(k):int(v) for k,v in valid.countryCode.value_counts().items()},
    }

def one_source_unit(site,source,category,param):
    """Prefer a reported subunit over its reporting parent RBD; NEVER add both."""
    country=str(site.countryCode)
    for grain,field in VALID_SCHEMES.items():
        # Order explicitly defined below instead of relying on dictionary.
        pass
    for scheme,field in (("euSubUnitCode","subUnitIdentifier"),("euRBDCode","rbdIdentifier")):
        v=getattr(site,field)
        if pd.isna(v) or str(v).strip() in ("","<NA>"):continue
        key=(country,scheme,str(v),category,param)
        rec=source.get(key)
        if rec is not None:return rec
    return None

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--original-archive",type=Path,required=True)
    parser.add_argument("--registry",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    a=parser.parse_args()
    if sha(a.registry)!=SOURCE_BQE_REGISTRY_SHA:
        raise ValueError("Original BQE source-only registry changed")
    bqe=pd.read_csv(a.registry,dtype="string",low_memory=False)
    if len(bqe)!=8245:raise ValueError("Original 2022 BQE source site denominator changed")
    with zipfile.ZipFile(a.original_archive) as z:
        emissions=pd.read_csv(z.open(TABLE),dtype="string",low_memory=False)
        spatial=pd.read_csv(z.open(SPATIAL),dtype="string",low_memory=False)
    source_rows,selected_source,qc=normalized_emission_source(emissions)
    region,bad=site_waterbody_metadata(spatial)
    metadata=bqe.merge(region,how="left",on=["countryCode","waterBodyIdentifier"],
                       validate="many_to_one")
    metadata["waterbody_region_mapped"]=metadata.subUnitIdentifier.notna()|metadata.rbdIdentifier.notna()
    src={}
    for r in source_rows.itertuples(index=False):
        key=(str(r.countryCode),str(r.spatialUnitIdentifierScheme),
             str(r.spatialUnitIdentifier),
             str(r.parameterEmissionsSourceCategory),
             str(r.observedPropertyDeterminandLabel))
        if key in src:raise ValueError("Non-unique annual nutrient-source input")
        src[key]={
          "reporting_method":str(r.procedureEmissionsMethod),
          "annual_tonnes":float(r.grained_annual_source_tons),
          "spatial_grain":str(r.spatialUnitIdentifierScheme),
        }
    per_site=[]
    for s in metadata.itertuples(index=False):
        record={
          "site_key":str(s.site_key),"biological_code":str(s.observedPropertyDeterminandBiologyEQRCode),
          "country":str(s.countryCode),
          "original_EEA_WFD_waterbody_id":str(s.waterBodyIdentifier),
          "upstream_EEA_measured_2022_P_WWTP_t_per_yr":
               float(s.upstream_2022_EEA_P_measured_t_per_yr),
          "upstream_EEA_measured_2022_N_WWTP_t_per_yr":
               float(s.upstream_2022_EEA_N_measured_t_per_yr),
          "source_2022_EEA_WWTP_positive":float(s.upstream_2022_unique_EEA_both_load_WWTP_count)>0,
          "source_reporting_region_resolved":bool(s.waterbody_region_mapped),
        }
        for cat,name in CATEGORIES.items():
            for param,sym in PARAMETERS.items():
                entry=one_source_unit(s,src,cat,param)
                key=name+"_"+sym
                record[key+"_2022_t_per_yr"]=entry["annual_tonnes"] if entry else None
                record[key+"_method"]=entry["reporting_method"] if entry else None
                record[key+"_spatial_grain"]=entry["spatial_grain"] if entry else None
        per_site.append(record)
    matched=pd.DataFrame(per_site)
    groupqc={}
    for code,name in SITE_CODES.items():
        part=matched[matched.biological_code==code].copy()
        p=part["agricultural_diffuse_P_2022_t_per_yr"].notna()
        n=part["agricultural_diffuse_N_2022_t_per_yr"].notna()
        positives=part.source_2022_EEA_WWTP_positive
        both=p&n
        groupqc[name]={
          "source_only_original_monitoring_site_n":len(part),
          "source_upstream_WWTP_reported_NP_positive_sites":int(positives.sum()),
          "source_2022_agricultural_P_reported_sites":int(p.sum()),
          "source_2022_agricultural_N_reported_sites":int(n.sum()),
          "source_2022_agricultural_PN_both_reported_sites":int(both.sum()),
          "source_2022_agricultural_P_plus_positive_WWTP_sites":int((p&positives).sum()),
          "source_2022_agricultural_PN_plus_positive_WWTP_sites":int((both&positives).sum()),
          "agricultural_P_and_WWTP_positive_countries":
             int(part.loc[p&positives,"country"].nunique()),
          "agricultural_P_reporting_countries":int(part.loc[p,"country"].nunique()),
          "agricultural_N_reporting_countries":int(part.loc[n,"country"].nunique()),
          "agricultural_P_reporting_methods":{
              str(k):int(v) for k,v in part.loc[p,"agricultural_diffuse_P_method"].value_counts().items()},
          "agricultural_P_reporting_grain":{
              str(k):int(v) for k,v in part.loc[p,"agricultural_diffuse_P_spatial_grain"].value_counts().items()},
          "country_agricultural_P_and_WWTP_positive_source":{
              str(k):int(v) for k,v in part.loc[p&positives,"country"].value_counts().items()},
        }
        print("GESUI_INDEPENDENT_AGRICULTURE_P_SOURCE_COVERAGE",
              json.dumps({"BQE":name,**groupqc[name]},sort_keys=True),flush=True)
    result={
      "schema":"gesui_EEA2026_emissions_2022_NP1_diffuse_PN_source_only_BQE_coverage_v01",
      "original_emissions_source":"https://www.eea.europa.eu/en/datahub/datahubitem-view/5bd950b1-52b9-42f3-b0d7-bf51ff3056b7",
      "timestamp_utc":datetime.now(timezone.utc).isoformat(),
      "original_source_only_BQE_registry_SHA256":SOURCE_BQE_REGISTRY_SHA,
      "spatial_2022_source_grain":"RBD or subunit, not local monitoring station reach",
      "official_code_list":"https://dd.eionet.europa.eu/datasets/3351",
      "agricultural_source_category_code":"NP1",
      "original_source_quality":qc,
      "ambiguous_source_to_waterbody_area_mapping":bad,
      "unique_source_WISE_BQE_site_waterbody_records":int(metadata.site_key.nunique()),
      "source_region_mapped_EEA_monitor_rows":int(metadata.waterbody_region_mapped.sum()),
      "source_2022_diffuse_NP1_by_BQE":groupqc,
      "no_2022_EQR_biological_outcome_read":True,
      "no_source_tonnes_added_across_hierarchical_RBD_and_subunit":True,
      "no_missing_agriculture_data_imputed_as_zero":True,
      "duplicate_or_unverifiable_emissions_excluded":True,
      "further_limits":[
         "RBD/subunit source area grain is different from EEA WWTP receiving-water reach grain",
         "Annual reported diffuse agricultural emissions may be modelled/estimated; method preserved",
         "Source coverage selection can be nonrandom with respect to reporting country and 2022 aquatic EQR",
         "A shared regional agricultural source value for many river sites does not give 1:1 local exposure",
         "No biological source-to-EQR sign/causation claimed; first report source validity and geographic denominator",
         "Spatial join uses source metadata exact country/waterBody/subunit/RBD, not nearest-neighbor pollution assumption"
      ]}
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"SOURCE_2022_AGRICULTURAL_TNTP_COVERAGE.json").write_text(
       json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    matched.to_csv(a.out/"SOURCE_ONLY_BQE_MONITOR_AGRICULTURAL_NP1_2022.csv",
                   index=False,float_format="%.10g")
    print("GESUI_NP1_2022_SOURCE_GATE_FINAL",
          json.dumps({"matched_source_site_n":len(matched),
            "source_biology_numeric_EQR_not_accessed":True,
            "agriculture_P_site_by_BQE":{x:v["source_2022_agricultural_P_reported_sites"]
                 for x,v in groupqc.items()},
            "agriculture_P_and_WWTP_positive":{x:v["source_2022_agricultural_P_plus_positive_WWTP_sites"]
                 for x,v in groupqc.items()}},sort_keys=True),flush=True)

if __name__=="__main__":main()
