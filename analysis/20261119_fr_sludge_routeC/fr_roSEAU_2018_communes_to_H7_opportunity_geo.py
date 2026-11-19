#!/usr/bin/env python3
"""Source-only 2018 French ROSEAU candidate communes x original H7 WKT.

Historical Etalab 2018 municipality boundaries are used only to screen sites
for possible field follow-up; NOT a parcel-level sludge exposure, negative
control, nitrogen mass, or plant/butterfly ecological effect.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urlsplit,parse_qs

import geopandas as gpd
import pandas as pd
from shapely import wkt,make_valid
from shapely.ops import unary_union

FACETS_SHA="251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4"
PANEL_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
CANDIDATES_SHA="971e20215a0788a229c5f0a18663e1eee9ddf70b6161fc4f70e78c8da9621333"
ORIGINAL_EIV_TOTAL=451
ORIGINAL_FIRST_PERIOD_H7=226
ORIGINAL_2018_CANDIDATE_CODES=72
CRS="EPSG:3035"
M2_PER_HA=10000.0
MIN_MEANINGFUL_HA=10.

def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):h.update(chunk)
    return h.hexdigest()

def candidates(path):
    if digest(path)!=CANDIDATES_SHA:raise ValueError("Original 130 validated candidate INSEE CSV SHA changed")
    with Path(path).open(newline="",encoding="utf-8") as f:
        records=list(csv.DictReader(f))
    if len(records)!=130 or len({r["INSEE_code"] for r in records})!=130:
        raise ValueError("Original six-year unique municipality candidate source changed")
    selected={r["INSEE_code"] for r in records if "2018" in r["candidate_source_years"].split(";")}
    if len(selected)!=ORIGINAL_2018_CANDIDATE_CODES:
        raise ValueError("Original validated 2018 candidate INSEE source not 72")
    if any(len(k)!=5 or not k.isdigit() for k in selected):
        raise ValueError("Original candidate code not valid 5-digit INSEE (mainland only)")
    return selected

def flora_h7(panel_path,eiv_path,facets_path):
    if digest(panel_path)!=PANEL_SHA or digest(facets_path)!=FACETS_SHA:
        raise ValueError("Original H7 or GBIF polygon source SHA changed")
    cohort=pd.read_csv(panel_path,usecols=["h","h6","period","country"],dtype=str)
    eiv=pd.read_csv(eiv_path,usecols=["h","period","full_source_trait_eligible"],dtype=str)
    if len(cohort)!=454 or len(eiv)!=454 or cohort.duplicated(["h","period"]).any() or eiv.duplicated(["h","period"]).any():
        raise ValueError("Original H7×period roster invalid or EIV duplication")
    eiv["eligible"]=eiv.full_source_trait_eligible.str.strip().str.lower().map({"true":True,"false":False})
    if eiv.eligible.isna().any() or int(eiv.eligible.sum())!=ORIGINAL_EIV_TOTAL:
        raise ValueError("Original 451 quality-passed NUTRI metadata changed")
    valid=cohort[cohort.period.eq("2017-2019")].merge(
        eiv[["h","period","eligible"]],on=["h","period"],validate="one_to_one")
    if len(valid)!=227 or int(valid.eligible.sum())!=ORIGINAL_FIRST_PERIOD_H7:
        raise ValueError("Original first-period EIV quality supports no longer 226")
    fr=valid.loc[valid.eligible&valid.country.eq("FR"),["h","h6","country"]].copy()
    if fr.h.duplicated().any() or fr.empty:raise ValueError("No original French flora-quality H7 source records")
    selected=set(fr.h)
    raw=json.loads(Path(facets_path).read_text())
    if len(raw)!=2808:raise ValueError("Original 2808 GBIF source facet roster changed")
    polygons={}
    for r in raw:
        if r.get("group")!="Plantae" or r.get("period")!="2017-2019":continue
        h=str(r.get("HYBAS_ID_L07"))
        if h not in selected:continue
        if h in polygons:raise ValueError("Duplicate original GBIF H7 geometry receipt")
        parts=[]
        for request in r["source_request_receipts"]:
            geomstr=parse_qs(urlsplit(request["url"]).query).get("geometry",[])
            if len(geomstr)!=1 or hashlib.sha256(geomstr[0].encode()).hexdigest()!=request["component_SHA"]:
                raise ValueError("Original H7 WKT query component SHA mismatch")
            parts.append(wkt.loads(geomstr[0]))
        if not parts:raise ValueError("Missing genuine H7 polygon geometry")
        g=unary_union(parts)
        if g.is_empty or not g.is_valid or g.area<=0:
            raise ValueError("Original GBIF H7 polygon invalid")
        polygons[h]=g
    if set(polygons)!=selected:
        raise ValueError("Original French H7 GBIF geometry source incomplete")
    gdf=gpd.GeoDataFrame(fr,geometry=gpd.GeoSeries(
        [polygons[h] for h in fr.h],index=fr.index,crs="EPSG:4326"),crs="EPSG:4326")
    return gdf.to_crs(CRS)

def historical_communes(file,codes):
    with gzip.open(file,"rt",encoding="utf-8") as f:
        obj=json.load(f)
    if obj.get("type")!="FeatureCollection":
        raise ValueError("Historical Etalab municipal shapefile not GeoJSON FeatureCollection")
    gdf=gpd.GeoDataFrame.from_features(obj["features"],crs="EPSG:4326")
    props=set(gdf.columns)
    idcol=next((k for k in ("code","code_insee","INSEE_COM","insee") if k in props),None)
    if idcol is None:raise ValueError("Original historical commune GeoJSON lacks identifiable code field: "+str(sorted(props)[:50]))
    ids=gdf[idcol].astype(str).str.zfill(5)
    gdf["insee_code"]=ids
    if gdf.insee_code.duplicated().any():
        raise ValueError("Original historic 2018 municipality polygons share same code")
    matched=gdf.loc[gdf.insee_code.isin(codes)].copy()
    missing=sorted(codes-set(matched.insee_code))
    if len(matched)+len(missing)!=len(codes):
        raise ValueError("Historic candidate municipality matches ambiguous")
    if len(matched)<10:
        raise ValueError("Historical 2018 municipality source has fewer than 10 matches, likely wrong schema")
    gdf=matched.to_crs(CRS)
    if any(g.is_empty or not g.is_valid or g.area<=0 for g in gdf.geometry):
        raise ValueError("Selected municipality polygon invalid in original Etalab 2018 source")
    return gdf,{"original_etatlab_geojson_feature_count":len(obj["features"]),
                "INSEE_original_code_field":idcol,
                "2018_candidate_codes":len(codes),
                "codes_found_in_2018_historical_geometries":len(matched),
                "candidate_codes_without_2018_polygon":missing}

def overlap(communes,h7):
    hits=[]
    index=h7.sindex
    for c in communes.itertuples():
        for idx in index.query(c.geometry,predicate="intersects"):
            h=h7.iloc[idx]
            intersect=c.geometry.intersection(h.geometry)
            if intersect.is_empty or intersect.area<=0:continue
            area=float(intersect.area/M2_PER_HA)
            hits.append({
               "year_of_commune_boundary":"2018",
               "original_candidate_INSEE_code":str(c.insee_code),
               "original_H7":str(h.h),"original_true_H6":str(h.h6),
               "intersection_ha":area,
               "pct_of_commune_area_inside_H7":100*intersect.area/c.geometry.area,
               "pct_of_H7_area_intersected_by_commune":100*intersect.area/h.geometry.area,
               "minimum_10ha_opportunity_screen":area>=MIN_MEANINGFUL_HA,
               "source_interpretation":"municipality_listed_in_positive_total_STEU_spreading_CANDIDATE; NOT confirmed field application or H7 nutrient exposure"
            })
    if len({(h["original_candidate_INSEE_code"],h["original_H7"]) for h in hits})!=len(hits):
        raise ValueError("Duplicate H7 municipality overlay pairing")
    strong=[z for z in hits if z["minimum_10ha_opportunity_screen"]]
    def counts(rows):
        return {
           "candidate_commune_codes_intersecting_FR_H7":len({x["original_candidate_INSEE_code"] for x in rows}),
           "unique_H7_with_candidate_commune_overlap":len({x["original_H7"] for x in rows}),
           "unique_true_H6_with_candidate_commune_overlap":len({x["original_true_H6"] for x in rows}),
           "source_overlap_pairs":len(rows),
           "sum_intersection_ha_NOT_sludge_treated_area":sum(x["intersection_ha"] for x in rows)
        }
    return hits,{"positive_area_intersection":counts(hits),"intersection_ge_10ha_screen":counts(strong)}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--candidates",type=Path,required=True)
    parser.add_argument("--municipal-2018-gz",type=Path,required=True)
    parser.add_argument("--h7-panel",type=Path,required=True)
    parser.add_argument("--eiv-eligibility",type=Path,required=True)
    parser.add_argument("--original-h7-facets",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    codes=candidates(args.candidates)
    h7=flora_h7(args.h7_panel,args.eiv_eligibility,args.original_h7_facets)
    communes,cm=historical_communes(args.municipal_2018_gz,codes)
    pairs,screen=overlap(communes,h7)
    result={
       "schema":"gesui_FR_ROSEAU_2018_COG_2018_50m_GBIF_H7_CANDIDATE_GEOGRAPHIC_ONLY_v1",
       "utc":datetime.now(timezone.utc).isoformat(),
       "2018_historic_commune_dataset":"https://etalab-datasets.geo.data.gouv.fr/contours-administratifs/2018/geojson/communes-50m.geojson.gz",
       "source_sha256":{"original_130_COG_candidate_codes":digest(args.candidates),
                        "original_2018_historical_municipality_geometry_gz":digest(args.municipal_2018_gz),
                        "original_H7_454_panel":digest(args.h7_panel),
                        "original_451_quality_pass_metadata":digest(args.eiv_eligibility),
                        "original_2808_facets_true_H7_wkt":digest(args.original_h7_facets)},
       "source_only_FR_plants_2017to19_H7_cohort_count":len(h7),
       "2018_historical_commune_source":cm,
       "SOURCE_ONLY_geographical_opportunity":screen,
       "not_original_FR139_shared_H7_72_cohort_unless_independently_verified":True,
       "no_plant_NUTRI_butterfly_HNS_river_EQR_outcomes_read":True,
       "not_confirmed_sludge_field_application_or_exposure":True,
       "do_NOT_use_zero_overlaps_as_unexposed_controls":True,
       "not_allowed_to_override_original_biotic_source_STOP":True,
    }
    (args.out/"FR_ROSEAU_2018_CANDIDATES_X_H7_OPPORTUNITY_SOURCE_GIS.json").write_text(
        json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    with (args.out/"FR_ROSEAU_2018_CANDIDATE_COMMUNES_X_ORIGINAL_FR_H7_ONLY.csv").open("w",newline="") as f:
        fields=["year_of_commune_boundary","original_candidate_INSEE_code","original_H7","original_true_H6","intersection_ha",
                "pct_of_commune_area_inside_H7","pct_of_H7_area_intersected_by_commune",
                "minimum_10ha_opportunity_screen","source_interpretation"]
        writer=csv.DictWriter(f,fieldnames=fields)
        writer.writeheader();writer.writerows(pairs)
    print("GESUI_FR_SLUDGE_2018_COMMUNES_H7_OPPORTUNITY_ONLY",json.dumps({
       "candidate_codes":len(codes),"matched_historical_codes":cm["codes_found_in_2018_historical_geometries"],
       "FR_plant_h7":len(h7),"screen":screen,
       "SLUDGE_BIOLOGICAL_EXPOSURE_CONFIRMED":False
    },sort_keys=True),flush=True)

if __name__=="__main__":
    main()
