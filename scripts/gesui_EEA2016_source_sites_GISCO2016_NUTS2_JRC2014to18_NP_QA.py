#!/usr/bin/env python3
"""Source-first geo overlap: 2016 EEA BQE station IDs × official GISCO NUTS2016 × JRC source.

The actual 2016 numeric EQR outcomes are NOT opened. Reproduce fixed 5823
2016 river biological source-only site-element rows from EEA Biology 2026,
spatially match to official Eurostat GISCO 1:1M NUTS 2016 level2 polygons,
then *exact* NUTS2 IDs and area to JRC 2014–18 modeled N+P point-source
tonnes. Checks independent region counts and country support before EQR.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import geopandas as gpd
import numpy as np,pandas as pd
from shapely.geometry import Point

ROOT=Path("results/eea2016_thirdyear_source_20261111")
JRC=ROOT/"JRC_GREEN_OFFICIAL_2014_18_MODELLED_POINT_NP_COUNTRY_NUTS2_BASIN_SOURCE.json"
BQE=ROOT/"2016_BQE_PREOUTCOME_SOURCE_ONLY/2016_BQE_PREOUTCOME_SOURCE_QC.json"
SITES=ROOT/"2016_BQE_PREOUTCOME_SOURCE_ONLY/2016_FIXED_BQE_MONITORING_SITE_METADATA_NO_EQR.csv"
GISCO="https://gisco-services.ec.europa.eu/distribution/v2/nuts/geojson/NUTS_RG_01M_2016_4326_LEVL_2.geojson"
BQE_CODE={"EEA_13-01-4":"benthic_invertebrates",
          "EEA_124-04-9":"phytobenthos",
          "EEA_14-02-8":"fish"}
MIN={"benthic_sites":900,"countries":8,"NUTS2_regions":60,"spatial_blocks":15,"countries_with_NUTS2_between_region_modelNP_density_variation":5}
def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(2*1024*1024),b""):
            h.update(b)
    return h.hexdigest()
def source(root,gisco):
    j=json.loads((root/JRC).read_text())
    b=json.loads((root/BQE).read_text())
    if j["schema"]!="gesui_OFFICIAL_JRC_GREEN_MAPSERVER_original_2014to2018_POINT_SOURCE_NP_3GEO_SCALES_SOURCE_v01":
        raise ValueError("Original public JRC 2014–18 source not independently frozen")
    if not j["source_gate_to_consider_EEA2016_NUTS2_spatial_join"]:
        raise ValueError("2014–18 official NUTS2 JRC N/P source quality support failed")
    if b["pre_2016_Biology_EQR_values_read"] is not False:
        raise ValueError("2016 EEA normalized biological EQR values leaked into source-only cohort")
    records=j["original_JRC_GREEN_2014_2018_source_data"]["NUTS2"]
    if records["original_modelled_NP_both_nonmissing_regions_n"]!=337:
        raise ValueError("Official JRC full NUTS2 N/P source count changed")
    rr={r["original_geography_ID"].upper():r for r in records["source_all_original_2014to2018_JRC_modelled_point_NP_t_region_aggregates"]}
    if len(rr)!=337:raise ValueError("Duplicated JRC NUTS2 model units")
    site=pd.read_csv(root/SITES,dtype={"site_key":"string","countryCode":"string",
             "procedureClassificationSystem":"string","region_3deg_block":"string"},
             low_memory=False)
    if len(site)!=b["source_only_geolocated_2016_site_element_rows"] or len(site)!=5823:
        raise ValueError("Original EEA2016 no-outcome site source changed")
    if site.duplicated(["site_key","observedPropertyDeterminandBiologyEQRCode"]).any():
        raise ValueError("Original EEA source has duplicated biology site/element pair")
    if not site.phenomenonTimeReferenceYear.eq(2016).all():
        raise ValueError("2016 biological survey source year changed")
    if sha(gisco) in ("",None):
        raise ValueError("Original official GISCO polygon SHA missing")
    polys=gpd.read_file(gisco)
    columns=set(polys.columns)
    needed={"NUTS_ID","CNTR_CODE","geometry"}
    if not needed.issubset(columns):
        raise ValueError("Original Eurostat GISCO NUTS2016 L2 polygon identifiers unavailable: "+str(columns))
    if polys.crs is None:raise ValueError("Official GISCO polygon lacks CRS")
    if polys.crs.to_epsg()!=4326:
        raise ValueError("Original Eurostat GISCO official NUTS source not WGS84")
    if "LEVL_CODE" in columns:
        polys=polys[polys.LEVL_CODE.astype(str).eq("2")].copy()
    if polys.NUTS_ID.duplicated().any() or len(polys)<250:
        raise ValueError("Official 2016 polygon NUTS2 census incomplete or ambiguous")
    polys["NUTS_ID"]=polys.NUTS_ID.astype(str).str.upper()
    polys["CNTR_CODE"]=polys.CNTR_CODE.astype(str).str.upper()
    area=polys.to_crs(3035).geometry.area/1_000_000
    if (area<=0).any() or area.isna().any():
        raise ValueError("Official 2016 GISCO polygon area invalid")
    polys["official_GISCO_NUTS2_area_km2"]=area.to_numpy(float)
    site["lat"]=pd.to_numeric(site.lat,errors="coerce")
    site["lon"]=pd.to_numeric(site.lon,errors="coerce")
    if site.lat.isna().any() or site.lon.isna().any():
        raise ValueError("Original source-only EEA2016 coordinates invalid")
    pts=gpd.GeoDataFrame(site,geometry=gpd.points_from_xy(site.lon,site.lat),crs=4326)
    idx=gpd.sjoin(pts,polys[["CNTR_CODE","NUTS_ID","official_GISCO_NUTS2_area_km2","geometry"]],
                  how="left",predicate="within",lsuffix="EEA",rsuffix="GISCO")
    # Explicitly distinguish missing and multi-polygon ambiguous station IDs.
    collision=idx.groupby(["site_key","observedPropertyDeterminandBiologyEQRCode"]).size()
    ambiguous={k for k,v in collision.items() if v>1}
    idx["ambiguous_polygon"]=[
       (a,b) in ambiguous for a,b in zip(idx.site_key,idx.observedPropertyDeterminandBiologyEQRCode)]
    idx["country_identifier_agrees"]=idx.countryCode.astype(str).str.upper().eq(idx.CNTR_CODE.astype(str).str.upper())
    idx["official_JRC_join_id"]=idx.CNTR_CODE.astype(str).str.upper()+"|"+idx.NUTS_ID.astype(str).str.upper()
    idx["exposure"]=idx.official_JRC_join_id.map(rr)
    keep=idx.loc[
       ~idx.ambiguous_polygon&
       idx.country_identifier_agrees&
       idx.NUTS_ID.notna()&
       idx.exposure.notna()
    ].copy()
    keep["PsN"]=keep.exposure.map(lambda x:x.get("original_GREEN2014_18_pointN_modelled_t"))
    keep["PsP"]=keep.exposure.map(lambda x:x.get("original_GREEN2014_18_pointP_modelled_t"))
    keep=keep.loc[keep.PsN.notna()&keep.PsP.notna()].copy()
    if keep.duplicated(["site_key","observedPropertyDeterminandBiologyEQRCode"]).any():
        raise ValueError("True 2016 river source site lacks unambiguous NUTS2 region")
    if len(keep)>len(site):raise ValueError("SJOIN generated more source sites than original")
    keep["log_modelled_PsN_t_per_km2"]=np.log1p(keep.PsN.astype(float)/keep.official_GISCO_NUTS2_area_km2)
    keep["log_modelled_PsP_t_per_km2"]=np.log1p(keep.PsP.astype(float)/keep.official_GISCO_NUTS2_area_km2)
    keep["site_source_original_sha256"]=keep.site_key.map(lambda x:hashlib.sha256(x.encode()).hexdigest())
    samples=[]
    report={}
    for code,name in BQE_CODE.items():
        gg=keep.loc[keep.observedPropertyDeterminandBiologyEQRCode.eq(code)].copy()
        siteKeys={x.site_key for x in gg.itertuples(index=False)}
        cells=gg.official_JRC_join_id.nunique()
        country_unique=gg.countryCode.nunique()
        uniqB=gg.region_3deg_block.nunique()
        by_country=gg.groupby("countryCode")
        varied=0
        for country,z in by_country:
            if z.official_JRC_join_id.nunique()<2:continue
            byregion=z.drop_duplicates("official_JRC_join_id")
            if byregion.log_modelled_PsN_t_per_km2.nunique()>1 and byregion.log_modelled_PsP_t_per_km2.nunique()>1:
                varied+=1
        report[name]={
           "original_2016_geolocated_EEA_source_sites_pre_NUTS2":int((site.observedPropertyDeterminandBiologyEQRCode==code).sum()),
           "2016_source_sites_with_exact_OFFICIAL_GISCO_2016_NUTS2_and_JRC_NP":len(gg),
           "2016_original_JRC_numeric_NP_independent_NUTS2_regions_n":cells,
           "2016_original_EEA_BQE_countries_source_covered":country_unique,
           "2016_original_EEA_3degree_source_blocks":uniqB,
           "2016_EEA_countries_with_interregional_N_and_P_modelled_source_density_variation_n":varied,
           "original_NUTS2_source_country_counts":{str(k):int(v) for k,v in Counter(gg.countryCode).items()},
           "original_NUTS2_source_unique_region_country_counts":{str(k):int(v) for k,v in gg.groupby("countryCode").official_JRC_join_id.nunique().items()},
           "original_NUTS2_SOURCE_3deg_fixed_blocks_by_country":{str(k):int(v) for k,v in gg.groupby("countryCode").region_3deg_block.nunique().items()},
           "2016_real_EEA_EQR_normalized_value_NOT_READ":True
        }
        for x in gg.itertuples(index=False):
            samples.append({
              "original_2016_BQE_site_key_SHA256":x.site_source_original_sha256,
              "original_2016_BQE_group":name,
              "original_2016_EEA_country":str(x.countryCode),
              "original_procedureClassificationSystem":str(x.procedureClassificationSystem),
              "original_fixed_2016_3deg_geo_block":str(x.region_3deg_block),
              "original_fixed_2016_3deg_BQE_fold":str(x.pre_2016_source_only_spatial_fold),
              "original_GISCO2016_NUTS2":str(x.NUTS_ID),
              "original_GISCO2016_NUTS2_km2":float(x.official_GISCO_NUTS2_area_km2),
              "original_JRC2014to18_modelled_PsN_t":float(x.PsN),
              "original_JRC2014to18_modelled_PsP_t":float(x.PsP),
              "log1p_modelled_POINT_N_tonnes_per_km2":float(x.log_modelled_PsN_t_per_km2),
              "log1p_modelled_POINT_P_tonnes_per_km2":float(x.log_modelled_PsP_t_per_km2)
            })
    bent=report["benthic_invertebrates"]
    gate=bool(
        bent["2016_source_sites_with_exact_OFFICIAL_GISCO_2016_NUTS2_and_JRC_NP"]>=MIN["benthic_sites"] and
        bent["2016_original_EEA_BQE_countries_source_covered"]>=MIN["countries"] and
        bent["2016_original_JRC_numeric_NP_independent_NUTS2_regions_n"]>=MIN["NUTS2_regions"] and
        bent["2016_original_EEA_3degree_source_blocks"]>=MIN["spatial_blocks"] and
        bent["2016_EEA_countries_with_interregional_N_and_P_modelled_source_density_variation_n"]>=MIN["countries_with_NUTS2_between_region_modelNP_density_variation"])
    return {
      "schema":"gesui_official_2016_EEA_benthic_source_NUTS2016_geopolygons_2014to2018_JRC_GREEN_modelled_point_NP_SOURCE_v01",
      "UTC":datetime.now(timezone.utc).isoformat(),
      "original_EEA2016_BQE_site_metadata_SHA256":sha(root/SITES),
      "original_EEA2016_BQE_source_gate_SHA256":sha(root/BQE),
      "official_JRC_337NUTS2_sourcevalues_SHA256":sha(root/JRC),
      "official_Eurostat_GISCO_NUTS2016_01M_L2_source_URL":GISCO,
      "original_official_Eurostat_GISCO_source_SHA256":sha(gisco),
      "GISCO_original_NUTS2016_polygon_records_n":len(polys),
      "original_2016_BQE_site_elements_before_GISCO_join":len(site),
      "original_GISCO_geographic_source_missing":int(idx.NUTS_ID.isna().sum()),
      "original_GISCO_multi_polygon_ambiguous_source_rows":int(idx.ambiguous_polygon.sum()),
      "original_original_EEA_country_and_GISCO_country_conflicting_rows":int((~idx.country_identifier_agrees).sum()),
      "original_source_sites_that_MATCH_2016_GISCO2016_and_JRC_GREEN_NP":report,
      "original_2016_source_site_and_modelled_NUTS2_PN_values_HASHED_IDs_NO_EQR":samples,
      "precommitted_2016_benthic_NP_modelled_NUTS2_SUPPORT_GATE":{"thresholds":MIN,
              "can_open_2016_EEA_numeric_normalised_EQR_and_test_POSTHOC_correlation":gate},
      "modelled_JRC_2014to2018_green_pointsource_includes_DOMESTIC_AND_INDUSTRIAL":True,
      "not_2016_actual_measured_WWTP_discharge_to_EEA_station":True,
      "all_original_2016_EEA_BQE_EQR_results_NOT_READ":True,
      "original_2016_waterbody_site_coordinates_NOT_REPUBLISHED":True
    }
def main():
    a=argparse.ArgumentParser();a.add_argument("--gisco",type=Path,required=True);a.add_argument("--out",type=Path,required=True);v=a.parse_args()
    root=Path(__file__).resolve().parents[1]
    output=source(root,v.gisco)
    v.out.parent.mkdir(parents=True,exist_ok=True)
    v.out.write_text(json.dumps(output,sort_keys=True,indent=2)+"\n")
    print("JRC2016_NUTS2_STATION_SOURCE",json.dumps({"bybiology":output["original_source_sites_that_MATCH_2016_GISCO2016_and_JRC_GREEN_NP"],
          "gate":output["precommitted_2016_benthic_NP_modelled_NUTS2_SUPPORT_GATE"]},indent=2),flush=True)
if __name__=="__main__":main()
