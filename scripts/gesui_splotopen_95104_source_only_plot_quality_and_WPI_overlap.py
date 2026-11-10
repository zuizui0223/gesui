#!/usr/bin/env python3
"""Genuine global sPlotOpen 95,104 plot SITE SOURCE/geometry gate, NO species.

Exact 2024-version mirror LFS header SHA; observational years 2000-2015,
full vascular flora, plot size 1-1000m2, coordinate uncertainty <=1000m.
Simple nearest WPI 1055 port distances are mere geography/support and never
assigned as WWp exposure. Biological species/cover source DT never opened.
"""
from __future__ import annotations
import argparse,hashlib,json,math,re
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

SRC=Path("results/global_all_Plantae_WWp2015_20261110/ORIGINAL_ALL1055_42COUNTRY_WWp2015_BEFORE_NEW808_GBIF_SOURCE.json")
HEADER_SHA="b20e21c70bf1c7fc2d79db403797e018bacbc38f74bb851dcf28825bc131977a"
HEADER_N=95104
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def clean_name(s):return re.sub(r"[^a-z0-9]","",str(s).lower())
def find_col(columns,*options,required=True):
    by={clean_name(x):x for x in columns}
    for opt in options:
        if clean_name(opt) in by:return by[clean_name(opt)]
    if required:raise ValueError("Required original sPlotOpen source column absent: "+repr(options)+" in "+str(list(columns)))
    return None

def sphere(lat,lon):
    a=np.deg2rad(lat.astype(float))
    b=np.deg2rad(lon.astype(float))
    return np.column_stack((np.cos(a)*np.cos(b),np.cos(a)*np.sin(b),np.sin(a)))
def full_frame(path):
    if sha(path)!=HEADER_SHA or path.stat().st_size!=24975505:
        raise ValueError("Secondary public mirror LFS 2024 sPlotOpen original header checksum changed")
    df=pd.read_csv(path,sep="\t",encoding="utf-8",na_values=["NA"],low_memory=False)
    if len(df)!=HEADER_N:raise ValueError("Original sPlotOpen mirror did not contain 95,104 plot headers")
    lat=find_col(df.columns,"Latitude","lat")
    lon=find_col(df.columns,"Longitude","lon","long")
    uncert=find_col(df.columns,"Location_uncertainty","Location uncertainty","coordinateUncertaintyInMeters")
    area=find_col(df.columns,"Releve_area","Relevé area","Plot_area","Plot_size")
    date=find_col(df.columns,"Date_of_recording","Date","Recording_date")
    scope=find_col(df.columns,"Plant_recorded","Plants_recorded","Plant_recording")
    country=find_col(df.columns,"Country","Country_code")
    biome=find_col(df.columns,"Biome",required=False)
    dataset=find_col(df.columns,"Dataset",required=False)
    plot_id=find_col(df.columns,"PlotID","Plot_ID","Plot_ID_sPlotOpen","ID","PlotID_sPlot","sPlotOpen_ID",required=False)
    coords=pd.DataFrame({
      "lat":pd.to_numeric(df[lat],errors="coerce"),
      "lon":pd.to_numeric(df[lon],errors="coerce"),
      "uncert_m":pd.to_numeric(df[uncert],errors="coerce"),
      "area_m2":pd.to_numeric(df[area],errors="coerce"),
      "date":pd.to_datetime(df[date],format="mixed",dayfirst=True,errors="coerce"),
      "plant_scope":df[scope].astype("string").str.strip(),
      "country":df[country].astype("string").str.strip(),
      "biome":df[biome].astype("string").str.strip() if biome else pd.Series(["UNKNOWN"]*len(df)),
      "dataset":df[dataset].astype("string").str.strip() if dataset else pd.Series(["UNKNOWN"]*len(df)),
    })
    coords["year"]=coords.date.dt.year
    valid_coord=(coords.lat.between(-90,90,inclusive="both")&
                 coords.lon.between(-180,180,inclusive="both"))
    geo=(valid_coord&coords.uncert_m.between(0,1000,inclusive="both"))
    tempor=(coords.year.between(2000,2015,inclusive="both"))
    full=coords.plant_scope.str.casefold().eq("all vascular plants").fillna(False)
    standardized=coords.area_m2.between(1,1000,inclusive="both")
    eligible=(geo&tempor&full&standardized&coords.country.notna())
    coords["SOURCE_GATE_PASS"]=eligible
    info={
        "column_names_original":list(df.columns),
        "source_identified_columns":{"latitude":lat,"longitude":lon,"location_uncertainty":uncert,
            "plot_area":area,"date":date,"plant_scope":scope,"country":country,
            "biome":biome,"source_dataset":dataset,"plot_id":plot_id},
        "source_header_plot_rows":len(df),
        "source_valid_wgs84_coordinates_n":int(valid_coord.sum()),
        "source_location_uncertainty_le1000_with_valid_coordinates_n":int(geo.sum()),
        "source_year_2000_2015_n":int(tempor.sum()),
        "source_all_vascular_plants_protocol_n":int(full.sum()),
        "source_original_plot_area_1_1000m2_n":int(standardized.sum()),
        "source_full_gate_2000to2015_full_flora_area_location_n":int(eligible.sum()),
        "source_full_gate_countries_n":int(coords.loc[eligible,"country"].nunique()),
        "source_full_gate_biomes":dict(Counter(coords.loc[eligible,"biome"].dropna())),
        "source_full_gate_countries":dict(Counter(coords.loc[eligible,"country"].dropna())),
        "source_full_gate_source_subdatasets":dict(Counter(coords.loc[eligible,"dataset"].dropna())),
        "source_full_gate_years":{str(int(y)):int(n) for y,n in Counter(coords.loc[eligible,"year"].dropna()).items()},
        "source_area_median_m2":float(coords.loc[eligible,"area_m2"].median()) if eligible.any() else None,
        "source_uncertainty_median_m":float(coords.loc[eligible,"uncert_m"].median()) if eligible.any() else None,
        "original_sPlotOpen_header_plot_ID_field_available":plot_id is not None,
        "all_95104_sPlotOpen_species_cover_OR_community_richness_NOT_ACCESSED":True
    }
    return df,coords,info

def analyze(root,header):
    old=json.loads((root/SRC).read_text(encoding="utf-8"))
    if old["all_original_country_WWp_source_sites_n"]!=1055 or old["original_country_count"]!=42:
        raise ValueError("Original all-1055 WPI geo source cohort missing")
    if not old["GBIF_new808_species_records_NOT_OBSERVED_in_stage0"]:
        raise ValueError("Initial WPI source was not independent of new plant records")
    _,plot,q=full_frame(header)
    wpi=old["source_all_1055_public_country_port_site_metadata_no_new_GBIF_response"]
    x=np.array([[v["original_port_lat"],v["original_port_lon"]] for v in wpi])
    valid=np.asarray(plot.SOURCE_GATE_PASS)
    p=np.column_stack((plot.loc[valid,"lat"].to_numpy(dtype=float),
                       plot.loc[valid,"lon"].to_numpy(dtype=float)))
    radius={}
    distances=np.array([])
    if len(p):
        distance,cidx=cKDTree(sphere(x[:,0],x[:,1])).query(sphere(p[:,0],p[:,1]),k=1)
        km=2*6371.0088*np.arcsin(np.clip(distance/2,0,1))
        distances=km
        for r in (25,50,100):
            mask=km<=r
            region_country=plot.loc[valid,"country"].reset_index(drop=True)
            radius[str(r)+"km"]={
               "original_source_quality_vegetation_plot_centroids":int(mask.sum()),
               "original_countries_with_any_plot_near_a_WPI_source_port":int(region_country[mask].nunique()),
               "nearest_WPI_port_not_actual_2015_wastewater_plot_dose":True
            }
    gate={
        "min_complete_source_plots_year2000_2015_all_vascular_area1_1000_uncertainty1000":2000,
        "min_source_countries":15,
        "plot_centered_PANGAEA_2015_WWp_source_geometry_allowed":bool(
            q["source_full_gate_2000to2015_full_flora_area_location_n"]>=2000 and
            q["source_full_gate_countries_n"]>=15)
    }
    return {
        "schema":"gesui_global_sPlotOpen2024_mirror_95104_header_2000to2015_original_fullflora_2015WWp_SOURCE_ONLY_v01",
        "utc":datetime.now(timezone.utc).isoformat(),
        "original_public_2024_mirror_sPlotOpen_header_LFS_SHA256":sha(header),
        "original_sPlotOpen_header_2024_mirror_file_size_bytes":header.stat().st_size,
        "original_global_WPI1055_positiveWWp_port_source_SHA256":hashlib.sha256((root/SRC).read_bytes()).hexdigest(),
        "sPlotOpen_2024_mirror_data_not_verified_bit_identical_to_2021_primary_idiv_release":True,
        "source_vegetation_plot_eligibility":q,
        "original_nearest_WPI_1055port_geographical_overlap":radius,
        "original_nearest_WPI_distance_median_km":float(np.median(distances)) if len(distances) else None,
        "source_gate_before_any_2015_PANGAEA_raster_or_sPlot_species":gate,
        "year2015_modelled_WWp_production_not_vegetation_nutrient_dose":True,
        "sPlotOpen_sample_excluded_wetlands_and_anthropogenic_vegetation":True,
        "no_original_species_cover_abundance_or_richness_values_read_in_this_stage":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--header",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    out=analyze(root,args.header)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "source":{k:v for k,v in out["source_vegetation_plot_eligibility"].items()
           if k not in ("source_full_gate_countries","source_full_gate_biomes","source_full_gate_source_subdatasets","column_names_original")},
        "geometry":out["original_nearest_WPI_1055port_geographical_overlap"],
        "next_source_gate":out["source_gate_before_any_2015_PANGAEA_raster_or_sPlot_species"]
    },indent=2))

if __name__=="__main__":
    main()
