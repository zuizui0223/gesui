#!/usr/bin/env python3
"""Only check terrestrial-plant H7 source-cohort overlap with JRC NUTS2 P/D inputs.

No EQR, plant NUTRI or butterfly HNS outcomes are read. Source geometry
functions are verbatim adapted from earlier independently tested source gate.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from pathlib import Path
from datetime import datetime,timezone
from collections import defaultdict,Counter
from urllib.parse import urlparse,parse_qs
import geopandas as gpd
import pandas as pd
import numpy as np
from shapely import make_valid,wkt
from shapely.ops import unary_union

FACETS_SHA="251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4"
PANEL_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
GPKG_SHA="7a233d3d683dd48241b4a9e5a9ef6c5f3e294008934c5bc7f1a465eb5cbf9c0b"
AREA_CUTOFF=.90
PANEL=Path("results/20261115_sync/hns_macro/H7_454_HNS_panel.csv")

def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""):h.update(block)
    return h.hexdigest()

def read_land(eligibility_path):
    if digest(PANEL)!=PANEL_SHA:raise ValueError("Original 454 H7 panel SHA mismatch")
    d=pd.read_csv(PANEL,usecols=["h","h6","period","country"],dtype=str)
    if len(d)!=454 or d.h.nunique()!=232 or d.h6.nunique()!=90:
        raise ValueError("Original H7 study cohort or H6 parent mapping changed")
    e=pd.read_csv(eligibility_path,usecols=["h","period","full_source_trait_eligible"],dtype=str)
    if len(e)!=454 or e.duplicated(["h","period"]).any():
        raise ValueError("Original EIV source quality rows changed")
    ok=e.full_source_trait_eligible.str.lower().map({"true":True,"false":False})
    if ok.isna().any() or int(ok.sum())!=451:raise ValueError("Prior 451 EIV trait-support gate changed")
    e["original_eiv_source_eligible"]=ok
    d=d[d.period=="2017-2019"].merge(
        e[["h","period","original_eiv_source_eligible"]],on=["h","period"],validate="one_to_one")
    if len(d)!=227:raise ValueError("Original 2017-19 H7 cohort changed")
    return d[d.original_eiv_source_eligible].copy()

def read_source(path,receipt):
    j=json.loads(Path(receipt).read_text())
    if not j.get("NO_EEA_EQR_EIV_GBIF_BIOLOGICAL_VALUES_READ",False):
        raise ValueError("Original JRC source-only quality receipt failed")
    if digest(path)!=j["csv_sha256"]:raise ValueError("Original JRC269 source SHA drift")
    df=pd.read_csv(path,dtype={"NUTS2_country_id":str})
    cols=["NUTS2_country_id","water_point_N_t","land_atmos_N_mg_m2",
          "land_mineral_N_kg_ha","land_organic_N_kg_ha"]
    df=df[cols].copy()
    for k in cols[1:]:df[k]=pd.to_numeric(df[k],errors="coerce")
    df=df.dropna(subset=cols[1:])
    if len(df)<269:raise ValueError("JRC fixed 269 complete-source NUTS2 unavailable")
    if df.NUTS2_country_id.duplicated().any():raise ValueError("Duplicate original JRC NUTS2")
    return df.set_index("NUTS2_country_id")

def geometries(gisco_path,hybas_path,land):
    nuts=gpd.read_file(gisco_path)
    if "NUTS_ID" not in nuts.columns and "id" in nuts.columns:
        nuts["NUTS_ID"]=nuts["id"].astype(str)
    if "CNTR_CODE" not in nuts.columns and "NUTS_ID" in nuts.columns:
        nuts["CNTR_CODE"]=nuts["NUTS_ID"].astype(str).str[:2]
    if not {"CNTR_CODE","NUTS_ID","geometry"}.issubset(nuts.columns):
        raise ValueError("2016 NUTS source columns missing: "+str(list(nuts.columns)))
    nuts=nuts.to_crs(3035)
    if "LEVL_CODE" in nuts:
        nuts=nuts.loc[nuts.LEVL_CODE.astype(str).eq("2")].copy()
    nuts["geo"]=nuts.CNTR_CODE.astype(str).str.upper()+"|"+nuts.NUTS_ID.astype(str).str.upper()
    if nuts.geo.duplicated().any():raise ValueError("Multiple polygons per native NUTS2")
    nuts["geometry"]=nuts.geometry.map(lambda g: make_valid(g) if not g.is_valid else g)
    full_geom_count=len(nuts)
    # Official HydroBASINS L7 zip OR the original GBIF H7 WKT query receipts.
    # The latter retains genuine historical geometry instead of centroid proxies.
    if str(hybas_path).endswith(".zip"):
        bas=gpd.read_file("zip://"+str(hybas_path))
        if "HYBAS_ID" not in bas.columns or bas.crs is None:
            raise ValueError("Official HydroBASINS L7 shape source or CRS invalid")
        bas["h"]=bas.HYBAS_ID.astype(str)
        bas=bas.loc[bas.h.isin(set(land.h))].copy()
    else:
        if digest(hybas_path)!=FACETS_SHA:
            raise ValueError("Original 2808-facet WKT source SHA mismatch")
        facets=json.loads(Path(hybas_path).read_text())
        if len(facets)!=2808:raise ValueError("Original source facet count differs")
        verts={}
        for record in facets:
            if record.get("group")!="Plantae" or record.get("period")!="2017-2019":continue
            h=str(record["HYBAS_ID_L07"])
            if h not in set(land.h):continue
            if h in verts:raise ValueError("Duplicated true H7 original WKT records")
            pieces=[]
            for rec in record["source_request_receipts"]:
                raw=parse_qs(urlparse(rec["url"]).query).get("geometry",[])
                if len(raw)!=1 or hashlib.sha256(raw[0].encode()).hexdigest()!=rec["component_SHA"]:
                    raise ValueError("Original GBIF queried H7 polygon receipt changed")
                pieces.append(wkt.loads(raw[0]))
            if not pieces:raise ValueError("H7 source polygon has zero geometry parts")
            verts[h]=unary_union(pieces)
        bas=gpd.GeoDataFrame([{"h":h,"geometry":g} for h,g in verts.items()],
                             geometry="geometry",crs=4326)
    
    if bas.h.duplicated().any() or len(bas)!=len(land):
        raise ValueError(f"All original 227 H7 genuine polygons required; found {len(bas)}")
    bas=bas.to_crs(3035)
    bas["geometry"]=bas.geometry.map(lambda g: make_valid(g) if not g.is_valid else g)
    if any(g.is_empty or not math.isfinite(g.area) or g.area<=0 for g in bas.geometry):
        raise ValueError("True original H7 basin polygon geometry invalid")
    return nuts,bas,full_geom_count

def map_water(water,nuts):
    pts=gpd.GeoDataFrame(water.copy(), geometry=gpd.points_from_xy(water.lon,water.lat),crs=4326).to_crs(3035)
    join=gpd.sjoin(pts,nuts[["geo","CNTR_CODE","geometry"]],how="left",predicate="within")
    counts=join.groupby("site_key",dropna=False).size()
    dup=set(counts[counts>1].index)
    join["ambiguous"]=join.site_key.isin(dup)
    ok=join.geo.notna()&~join.ambiguous&join.countryCode.eq(join.CNTR_CODE)
    return join[ok].copy(),{"water_pre_sites":len(water),"mapped_sites":int(ok.sum()),
                             "ambiguous_sites":len(dup),
                             "outside_source_polygons":int(join.geo.isna().sum()),
                             "country_conflict_rows":int((join.geo.notna()&~join.countryCode.eq(join.CNTR_CODE)).sum())}

def map_land(land,bas,nuts):
    idx=nuts.sindex
    rows=[]
    for h,g in bas.set_index("h").geometry.items():
        gids=list(idx.query(g,predicate="intersects"))
        area=g.area
        portions=[]
        for j in gids:
            rec=nuts.iloc[j]
            overlap=g.intersection(rec.geometry).area/area
            if overlap>0:portions.append((str(rec.geo),float(overlap)))
        portions.sort(key=lambda p:(-p[1],p[0]))
        maxgeo,frac=portions[0] if portions else (None,0.0)
        assigned=maxgeo if frac>=AREA_CUTOFF and sum(p[1]>=AREA_CUTOFF for p in portions)==1 else None
        rows.append({"h":h,"assigned_NUTS2":assigned,"max_overlap":frac,"crossed_nuts2_n":len(portions)})
    hmap=pd.DataFrame(rows)
    out=land.merge(hmap,on="h",validate="one_to_one")
    out["country_geom"]=out.assigned_NUTS2.str.split("|",regex=False).str[0]
    return out,{"land_original_H7_period_primary":len(land),
                 "H7_strict_assigned_over90pct":int(out.assigned_NUTS2.notna().sum()),
                 "H7_cross_nuts2":int((out.crossed_nuts2_n>1).sum()),
                 "H7_not_strict_assigned":int(out.assigned_NUTS2.isna().sum()),
                 "H7_overlap_ratio_minmedianmax":[float(out.max_overlap.min()),float(out.max_overlap.median()),float(out.max_overlap.max())]}

def components(frame):
    """Graph nodes are true H6 and NUTS2; H7 links them. Prevent fold leakage."""
    parents={}
    def find(x):
        parents.setdefault(x,x)
        while parents[x]!=x:
            parents[x]=parents[parents[x]]
            x=parents[x]
        return x
    def union(a,b):
        ra,rb=find(a),find(b)
        if ra!=rb:parents[max(ra,rb)]=min(ra,rb)
    for r in frame.itertuples():
        union("h6:"+str(r.h6),"nuts2:"+str(r.assigned_NUTS2))
    group=defaultdict(list)
    for r in frame.itertuples():
        group[find("h6:"+str(r.h6))].append(r)
    fold_loads=[0]*5
    folds={}
    for key,z in sorted(group.items(),key=lambda x:(-len(x[1]),x[0])):
        f=min(range(5),key=lambda k:(fold_loads[k],k))
        folds[key]=f; fold_loads[f]+=len(z)
    return len(group),[sum(f==i for f in folds.values()) for i in range(5)],fold_loads

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--gisco",type=Path,required=True)
    p.add_argument("--facets",type=Path,required=True)
    p.add_argument("--eiv-eligibility",type=Path,required=True)
    p.add_argument("--jrc-source",type=Path,required=True)
    p.add_argument("--jrc-report",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    if digest(args.gisco)!=GPKG_SHA:raise ValueError("GISCO2016 official original1M source SHA changed")
    land=read_land(args.eiv_eligibility)
    source=read_source(args.jrc_source,args.jrc_report)
    nuts,bas,ngeo=geometries(args.gisco,args.facets,land)
    mapped,mqc=map_land(land,bas,nuts)
    eligible=mapped[mapped.assigned_NUTS2.isin(source.index)].copy()
    keyset=set(eligible.assigned_NUTS2)
    countries=Counter(k.split("|",1)[0] for k in keyset)
    area={z.geo:float(z.geometry.area/1e6) for z in nuts.itertuples()}
    varied=[]
    for c,n in sorted(countries.items()):
        if n<3:continue
        keys=sorted(k for k in keyset if k.startswith(c+"|"))
        pvals={round(math.log1p(float(source.at[k,"water_point_N_t"])/area[k]),10) for k in keys}
        dvals={round(math.log1p(float(source.at[k,"land_atmos_N_mg_m2"])),10) for k in keys}
        if len(pvals)>=2 and len(dvals)>=2:varied.append(c)
    ncomp,fold_component,fold_H7=components(eligible)
    decisions={
      "NUTS2_ge_55":len(keyset)>=55,
      "countries_ge_8":len(countries)>=8,
      "true_H6_ge_40":eligible.h6.nunique()>=40,
      "countries_with_three_NUTS2_P_D_variation_ge_5":len(varied)>=5,
      "connected_components_ge_20":ncomp>=20,
      "each_5fold_ge_4_components":all(x>=4 for x in fold_component)
    }
    results={
      "schema":"gesui_2017to19_terrestrial_EIV_quality_H7_JRC269_P_D_SOURCEONLY_v1",
      "UTC":datetime.now(timezone.utc).isoformat(),
      "NUTS2_country_IDs":sorted(keyset),
      "NUTS2_n":len(keyset),"country_n":len(countries),
      "countries":sorted(countries),
      "H7_eiv_source_quality_eligible_n":len(land),
      "H7_with_JRC_source_complete_and_overlap90_n":len(eligible),
      "true_H6_n":eligible.h6.nunique(),
      "countries_three_regions_and_P_D_variation":varied,
      "countries_three_regions_and_P_D_variation_n":len(varied),
      "independent_connected_geo_components_n":ncomp,
      "geo_components_per_fixed5fold":fold_component,
      "H7_per_fixed5fold":fold_H7,
      "JRC_269_source_original_SHA":digest(args.jrc_source),
      "original_H7_panel_SHA":digest(PANEL),
      "original_H7_facets_SHA":digest(args.facets),
      "original_2016_GISCO01M_sha":digest(args.gisco),
      "original_EIV_451_quality_only_SHA":digest(args.eiv_eligibility),
      "gis_source_original_polygon_count":ngeo,
      "geographic_quality":mqc,
      "source_gate_conditions":decisions,
      "source_only_GO_for_new_land_only_biotic_model":all(decisions.values()),
      "original_water_land_joint_2x2_SOURCE_STOP_remains":True,
      "no_plant_EIV_NUTRI_numeric_butterfly_HNS_or_EQR_read":True,
    }
    (args.out/"TERRESTRIAL_ONLY_EIV_SOURCE_P_D_NUTS2_ELIGIBILITY_20261119.json").write_text(
        json.dumps(results,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    with (args.out/"TERRESTRIAL_ONLY_NUTS2_SOURCE_SUPPORT.csv").open("w",newline="") as f:
        writer=csv.writer(f)
        writer.writerow(["NUTS2_country_ID","H7_count","true_H6_count","pointN_model_t","atmosN_model_mg_m2"])
        for k in sorted(keyset):
            gg=eligible[eligible.assigned_NUTS2==k]
            writer.writerow([k,len(gg),gg.h6.nunique(),source.at[k,"water_point_N_t"],source.at[k,"land_atmos_N_mg_m2"]])
    print("GESUI_TERRESTRIAL_SOURCE_ONLY_POINT_LANDN",json.dumps({
        k:results[k] for k in ("NUTS2_n","country_n","countries","H7_eiv_source_quality_eligible_n",
            "H7_with_JRC_source_complete_and_overlap90_n","true_H6_n",
            "countries_three_regions_and_P_D_variation_n",
            "independent_connected_geo_components_n","geo_components_per_fixed5fold",
            "source_gate_conditions","source_only_GO_for_new_land_only_biotic_model")
    },ensure_ascii=False),flush=True)

if __name__=="__main__":main()
