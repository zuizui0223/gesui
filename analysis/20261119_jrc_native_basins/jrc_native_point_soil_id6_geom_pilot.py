#!/usr/bin/env python3
"""One-record JRC HydroID=6 / BasinID=6 geospatial diagnostic only.

Chosen from prior published source schema sample BEFORE any polygon results.
Never claim all JRC IDs equivalent based on this one positive example.
"""
import argparse, json
from pathlib import Path
from datetime import datetime,timezone
from jrc_native_point_soil_polygon30_gate import original_geom

ID=6

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    soil,sattr,soilha=original_geom("soil_atmosN_2015",ID)
    point,pattr,pointsha=original_geom("water_pointN_2014_18",ID)
    overlap=soil.intersection(point).area
    fracsoil=overlap/soil.area
    fracpoint=overlap/point.area
    report={
      "schema":"gesui_JRC_20261119_native_HydroID6_BasinID6_SINGLE_illustrative_geom_source_v1",
      "UTC":datetime.now(timezone.utc).isoformat(),
      "sample_source_prechosen_ID":ID,
      "soil_source_native_field":"Basin_3035_HydroID",
      "point_source_native_field":"BasinID",
      "soil_original_area_km2":sattr["Basin_3035_AreaSqKm"],
      "point_original_area_km2":pattr["AreaSqKm"],
      "soil_polygon_area_km2":soil.area/1e6,
      "point_polygon_area_km2":point.area/1e6,
      "mutual_polygon_overlap_fraction_of_soil":fracsoil,
      "mutual_polygon_overlap_fraction_of_point":fracpoint,
      "single_case_geometry_agreement_at_least95":min(fracsoil,fracpoint)>=.95,
      "original_soil_geometry_sha256":soilha,
      "original_point_geometry_sha256":pointsha,
      "one_ID_does_NOT_validate_others":True,
      "biotic_data_NOT_read":True,
    }
    (a.out/"JRC_SOIL_HYDROID6_POINT_BASINID6_SINGLE_GEOMETRY.json").write_text(
      json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("JRC_SINGLE_ID6_POLYGON_SOURCE_PILOT",json.dumps(report,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
