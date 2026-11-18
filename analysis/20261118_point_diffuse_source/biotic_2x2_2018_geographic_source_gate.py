#!/usr/bin/env python3
"""Pre-outcome geographic eligibility for the frozen point-vs-landN 2x2.

Uses immutable source-only aquatic registry, only geography/effort columns
of existing H7 panel, GISCO 2016 original NUTS2 geometries, HydroBASINS L7
polygons, and archived JRC source-only input values. Explicitly never
reads EQR, EIV NUTRI or butterfly HNS numerical outcome columns.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math
from collections import defaultdict, Counter
from pathlib import Path
from datetime import datetime, timezone
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import make_valid
from shapely.geometry import Point

ROOT=Path(__file__).resolve().parents[2]
FROZEN=ROOT/"analysis/20261118_point_diffuse_source/BIOTIC_2X2_FROZEN_CONFIG_20261118.json"
AQUATIC=ROOT/"results/eea2018_independent_source_frozen_20261108/2018_SOURCE_ONLY_BQE_WWTP_NP_UPSTREAM_AND_NEGATIVE_CONTROL.csv"
AQUATIC_QC=ROOT/"results/eea2018_independent_source_frozen_20261108/2018_MEASURED_SOURCE_BIOLOGY_COHORT_QC.json"
H7_META=ROOT/"results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"
EXPECTED_H7_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
HYBAS_SHA="b9b810f4db271f0fa6b39fecf01784041b948944b3c3f296a393766388e86b44"
SITE_COLS=("site_key","countryCode","observedPropertyDeterminandBiologyEQRCode","procedureClassificationSystem","phenomenonTimeReferenceYear","lat","lon")
H7_COLS=("h","h6","period","country","nPlant","nButterfly")
SOURCE_FIELDS=("NUTS2_country_id","water_point_N_t","water_point_P_t","land_atmos_N_mg_m2","land_mineral_N_kg_ha","land_organic_N_kg_ha","land_fixation_N_kg_ha")
AREA_CUTOFF=0.90
BENTHIC="EEA_13-01-4"

def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as fp:
        for buf in iter(lambda:fp.read(1<<20),b""):h.update(buf)
    return h.hexdigest()

def strict_main_source(path, expected, cols):
    header=pd.read_csv(path,nrows=0).columns.tolist()
    if not set(cols).issubset(header):
        raise ValueError("Original source-only fields unavailable: "+repr(set(cols)-set(header)))
    if expected and digest(path)!=expected:
        raise ValueError(f"Hash-locked cohort bytes have changed: {path}")
    return pd.read_csv(path,usecols=list(cols),dtype=str,low_memory=False)

def read_sources(path, report, config):
    raw=pd.read_csv(path,dtype=str)
    j=json.loads(Path(report).read_text())
    if not j.get("NO_EEA_EQR_EIV_GBIF_BIOLOGICAL_VALUES_READ",False):
        raise ValueError("Original source-only receipt lacks outcome-exclusion attestation")
    if digest(path)!=j["csv_sha256"]:
        raise ValueError("Archived exact JRC source-only CSV SHA changed")
    if not set(SOURCE_FIELDS).issubset(raw.columns):
        raise ValueError("JRC source-only 6 variables not all present")
    raw=raw[list(SOURCE_FIELDS)].copy()
    if raw.NUTS2_country_id.duplicated().any(): raise ValueError("Duplicated official country|NUTS2")
    for x in SOURCE_FIELDS[1:]:raw[x]=pd.to_numeric(raw[x],errors="coerce")
    complete=raw.dropna(subset=list(SOURCE_FIELDS)).copy()
    if len(complete)!=269 or len(complete.NUTS2_country_id.str.split("|",regex=False).str[0].unique())!=28:
        raise ValueError("269 NUTS2/28 country original source cohort no longer matches frozen receipt")
    if (complete[list(SOURCE_FIELDS[1:])]<0).any().any():
        raise ValueError("Negative original model source value")
    return complete.set_index("NUTS2_country_id")

def read_water():
    if not AQUATIC.exists(): raise ValueError("Original 2018 water SOURCE-only archive not present")
    orig=strict_main_source(AQUATIC,None,SITE_COLS)
    qc=json.loads(AQUATIC_QC.read_text())
    if len(orig)!=qc["source_2018_BQE_site_element_rows_after_river_graph"]:
        raise ValueError("Original 2018 aquatic source row total changed")
    sub=orig[(orig.observedPropertyDeterminandBiologyEQRCode==BENTHIC)&(orig.phenomenonTimeReferenceYear=="2018")].copy()
    if len(sub)!=2639 or sub.site_key.duplicated().any():raise ValueError("Original 2018 benthic source site roster drift")
    sub["lat"]=pd.to_numeric(sub.lat,errors="coerce")
    sub["lon"]=pd.to_numeric(sub.lon,errors="coerce")
    if sub[["lat","lon"]].isna().any().any(): raise ValueError("2018 original site lat/lon invalid")
    return sub

def read_land(eiv_eligibility):
    if digest(H7_META)!=EXPECTED_H7_SHA:
        raise ValueError("Original H7 original outcome-containing panel SHA changed (will not read outcomes)")
    d=strict_main_source(H7_META,EXPECTED_H7_SHA,H7_COLS)
    if len(d)!=454 or d.h.nunique()!=232 or d.h6.nunique()!=90:
        raise ValueError("Original full 454 H7 cohort or H6 parent IDs changed")
    if d.duplicated(["h","period"]).any():raise ValueError("Original H7period duplicate")
    one=d[d.period.eq("2017-2019")].copy()
    if len(one)!=227 or one.h.nunique()!=227:
        raise ValueError("Primary 2017-2019 source cohort drift")
    # Retrieve original v2 source-only trait eligibility exactly (no NUTRI outcomes).
    elig=strict_main_source(eiv_eligibility,None,("h","period","full_source_trait_eligible"))
    if len(elig)!=454 or elig.duplicated(["h","period"]).any():
        raise ValueError("Original EIV 454 source eligibility table missing or duplicated")
    ok=elig.full_source_trait_eligible.str.strip().str.lower().map({"true":True,"false":False})
    if ok.isna().any() or int(ok.sum())!=451:
        raise ValueError("Prior source-only EIV 451 pass count did not reproduce")
    elig["quality_eligible"]=ok
    one=one.merge(elig[["h","period","quality_eligible"]],on=["h","period"],validate="one_to_one",how="left")
    if one.quality_eligible.isna().any():raise ValueError("2017-2019 EIV eligibility source gaps")
    return one[one.quality_eligible].copy()

def geometries(gisco_path,hybas_path,land):
    nuts=gpd.read_file(gisco_path)
    if not {"CNTR_CODE","NUTS_ID","geometry"}.issubset(nuts.columns):
        raise ValueError("GISCO 2016 source columns missing")
    nuts=nuts.to_crs(3035)
    if "LEVL_CODE" in nuts:
        nuts=nuts.loc[nuts.LEVL_CODE.astype(str).eq("2")].copy()
    nuts["geo"]=nuts.CNTR_CODE.astype(str).str.upper()+"|"+nuts.NUTS_ID.astype(str).str.upper()
    if nuts.geo.duplicated().any():raise ValueError("Multiple polygons per native NUTS2")
    nuts["geometry"]=nuts.geometry.map(lambda g: make_valid(g) if not g.is_valid else g)
    full_geom_count=len(nuts)
    # Read all original EU HydroBASINS L7. This is not the WKT centroid shortcut.
    bas=gpd.read_file("zip://"+str(hybas_path))
    if "HYBAS_ID" not in bas.columns or bas.crs is None:
        raise ValueError("Official HydroBASINS L7 shape source or CRS invalid")
    bas["h"]=bas.HYBAS_ID.astype(str)
    bas=bas.loc[bas.h.isin(set(land.h))].copy()
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

def gate(source,water,land):
    a=set(water.geo)&set(source.index)
    b=set(land.assigned_NUTS2.dropna())&set(source.index)
    common=sorted(a & b)
    w=water[water.geo.isin(common)]
    l=land[land.assigned_NUTS2.isin(common)]
    # EIV quality must be exact row-level; currently only 451/454 global.
    country=Counter(s.split("|")[0] for s in common)
    country_var={}
    for c,n in country.items():
        if n<3:continue
        data=source.loc[[s for s in common if s.startswith(c+"|")]]
        p=data.water_point_N_t.nunique()>1
        dep=data.land_atmos_N_mg_m2.nunique()>1
        country_var[c]={"shared_NUTS2":n,"point_N_variation":bool(p),"atmosN_variation":bool(dep)}
    countries_ok=sum(z["point_N_variation"] and z["atmosN_variation"] for z in country_var.values())
    ncomp,fold_components,fold_h7=components(l)
    core={
        "common_water_and_land_strict_NUTS2_source_upper_bound":len(common),
        "common_country_n":len(country),
        "common_countries":sorted(country),
        "common_true_H6":l.h6.nunique(),
        "common_distinct_H7":l.h.nunique(),
        "common_2018_EEA_benthic_sites":w.site_key.nunique(),
        "countries_with_three_shared_NUTS2_and_both_P_D_variation":countries_ok,
        "linked_H6_NUTS2_geographic_components":ncomp,
        "components_per_fixed_5fold":fold_components,
        "H7_per_fixed_5fold":fold_h7,
        "status_per_frozen_thresholds":{
            "NUTS2>=40":len(common)>=40,
            "countries>=8":len(country)>=8,
            "H6>=30":l.h6.nunique()>=30,
            ">=5_countries_3_regions_P_D_change":countries_ok>=5,
            "geo_components>=20":ncomp>=20,
            "each_fold>=4_components":all(n>=4 for n in fold_components)
        },
        "biotic_trait_451_H7period_exact_ids_recovered":True,
        "background_covariate_complete_cases_NOT_ASSESSED":True
    }
    thresholds_pass=all(core["status_per_frozen_thresholds"].values())
    core["geographic_upper_bound_thresholds_pass"]=thresholds_pass
    core["final_biotic_source_gate_pass"]=False
    core["decision"]="SOURCE_STOP_GEOGRAPHIC_UPPER_BOUND" if not thresholds_pass else "GEOGRAPHIC_EIV_QUALITY_PASS_BACKGROUND_COVARIATES_PENDING"
    return core

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--gisco",type=Path,required=True)
    p.add_argument("--hybas",type=Path,required=True)
    p.add_argument("--jrc-source",type=Path,required=True)
    p.add_argument("--jrc-report",type=Path,required=True)
    p.add_argument("--eiv-eligibility",type=Path,required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    args=p.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=True)
    if digest(args.hybas)!=HYBAS_SHA:raise ValueError("HydroBASINS official L7 archive SHA differs from prior source hash")
    freeze=json.loads(FROZEN.read_text())
    if freeze.get("schema")!="gesui_20261118_point_vs_landN_biotic_2x2_preoutcome_v1":
        raise ValueError("2x2 preread frozen analysis contract missing")
    source=read_sources(args.jrc_source,args.jrc_report,freeze)
    water=read_water()
    land=read_land(args.eiv_eligibility)
    nuts,bas,n_all=geometries(args.gisco,args.hybas,land)
    water_map,wdiag=map_water(water,nuts)
    land_map,ldiag=map_land(land,bas,nuts)
    core=gate(source,water_map,land_map)
    # no site coordinates, EQR, NUTRI or butterfly outcome values written out
    derived={
      "schema":"gesui_2x2_biotic_source_geography_upperbound_v1",
      "utc":datetime.now(timezone.utc).isoformat(),
      "only_original_source_metadata_opened":True,
      "outcome_numeric_fields_not_read":["resultNormalisedEQRValue","Plant_N_incidence","primary_incidence_mean","primary_record_mean"],
      "source_SHA256":{"H7_original_panel":digest(H7_META),"2018_EEA_source_noEQR":digest(AQUATIC),"2016_GISCO_L2_polygon":digest(args.gisco),"HydroBASINS_L7_official":digest(args.hybas),"archived_JRC_NUTS2_CSV":digest(args.jrc_source),"EIV_SOURCE_ONLY_451_IDS":digest(args.eiv_eligibility)},
      "orig_gisco_polygon_count":n_all,
      "water_mapping":wdiag,
      "land_mapping":ldiag,
      "gates":core,
      "critical_limitations":[
         "All 451 source-eligible EIV H7-period IDs are restored; only 2017-2019 strict source records used.",
         "H7 polygon 90% overlap is not observed riparian plant location.",
         "EEA site numeric EQR availability and complete climate/land-use/effort covariates not yet checked.",
         "Modelled GREEN pointN combines municipal/industrial and does not measure wastewater fraction reaching vegetation.",
         "Country blocked and H6-NUTS2 geography needed before evaluation; absence of new effects is not a null result."
      ]
    }
    report=args.output_dir/"BIOTIC_2X2_2018_2017TO19_GEOGRAPHIC_SOURCE_ONLY_GATE.json"
    report.write_text(json.dumps(derived,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    summary=args.output_dir/"common_nuts2_anonymized_source_support.csv"
    eligible=sorted(set(water_map.geo)&set(land_map.assigned_NUTS2.dropna())&set(source.index))
    with summary.open("w",newline="") as f:
        wr=csv.writer(f);wr.writerow(["country_NUTS2","EEA2018_benthic_site_count","H7_2017_2019_count","genuine_H6_count","modelled_point_N","modelled_land_atmos_N"])
        for k in eligible:
            land_sel=land_map[land_map.assigned_NUTS2==k]
            wr.writerow([k,int(water_map[water_map.geo==k].site_key.nunique()),len(land_sel),int(land_sel.h6.nunique()),
                         float(source.at[k,"water_point_N_t"]),float(source.at[k,"land_atmos_N_mg_m2"])])
    print("GESUI_POINT_LAND_SOURCE_ONLY_BIOTIC_GEOGRAPHIC_GATE",json.dumps(core,ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__": main()
