#!/usr/bin/env python3
"""Source-only actual official Level7→Level6 polygon mapping of EEA single-year H7.

Parent polygons matched by >=95% of native H7 polygon area overlap and an
unambiguous 0.5 gap over second best, matching prior genuine 232 H7 source.
NO EQR, butterfly, plant NUTRI numeric, treatment or slope opened. Only H7,
H6, source monitoring station counts, and original eligibility booleans.
"""
from __future__ import annotations
import argparse
import collections
import csv
import hashlib
import io
import json
from pathlib import Path
import time
import urllib.request
import zipfile

import shapefile
from shapely import make_valid
from shapely.geometry import shape as shaped
from shapely.ops import unary_union
from shapely.strtree import STRtree

GEO_SOURCE_SHA="a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89"
GEO_SOURCE_NAME="EEA2018_EEA2022_EPA1819_NATIVE_HYBAS7_SOURCE_ONLY_GEOGRAPHY.json"
PANEL_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
EIV_SHA="ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba"
L7_ZIP_SHA="b9b810f4db271f0fa6b39fecf01784041b948944b3c3f296a393766388e86b44"
PERIODS={"2017-2019","2020-2022"}
NEED_NEW_SINGLETON_H6=10

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def one(directory,name):
    files=list(directory.rglob(name))
    if len(files)!=1:raise RuntimeError("Exactly one original file required: "+name+" "+repr([str(x) for x in files]))
    return files[0]
def fetch(lev):
    url=f"https://data.hydrosheds.org/file/hydrobasins/standard/hybas_eu_lev{lev:02d}_v1c.zip"
    err=None
    for k in range(3):
        try:
            request=urllib.request.Request(url,headers={
                "User-Agent":"gesui-source-only-EU501-H7-to-H6-official-geometry-20261115/1",
                "Accept":"application/zip,*/*"})
            with urllib.request.urlopen(request,timeout=90) as response:
                raw=response.read(75000000)
            if len(raw)>=75000000 or not zipfile.is_zipfile(io.BytesIO(raw)):
                raise RuntimeError("Original official HydroBASINS ZIP not valid")
            sha=hashlib.sha256(raw).hexdigest()
            if lev==7 and sha!=L7_ZIP_SHA:
                raise RuntimeError("Published official Lev7 HydroBASINS original source SHA mismatch")
            return raw,{"url":url,"sha256":sha,"size_bytes":len(raw)}
        except Exception as ex:
            err=ex;time.sleep(2+k*2)
    raise RuntimeError("Could not verify genuine official native HYBAS geometry Level"+str(lev)+": "+str(err))
def geom_records(data,lev,selected):
    z=zipfile.ZipFile(io.BytesIO(data))
    groups={}
    for name in z.namelist():
        if f"lev{lev:02d}" not in name.lower():continue
        extension=name.rsplit(".",1)[-1].lower()
        if extension not in ("shp","shx","dbf"):continue
        base=name.rsplit(".",1)[0]
        groups.setdefault(base,{})[extension]=name
    grp=[v for v in groups.values() if all(k in v for k in ("shp","shx","dbf"))]
    if len(grp)!=1:raise RuntimeError("Genuine source shapefile not unique at level "+str(lev))
    fieldsgrp=grp[0]
    sh=shapefile.Reader(shp=io.BytesIO(z.read(fieldsgrp["shp"])),
        shx=io.BytesIO(z.read(fieldsgrp["shx"])),
        dbf=io.BytesIO(z.read(fieldsgrp["dbf"])))
    fields=[x[0] for x in sh.fields[1:]]
    if "HYBAS_ID" not in fields:raise RuntimeError("Official source does not have HYBAS_ID")
    ids=[];geoms=[];invalid=0;repaired=0;failed=0
    for rec in sh.iterShapeRecords():
        vals=dict(zip(fields,rec.record))
        h=str(vals["HYBAS_ID"])
        if lev==7 and h not in selected:continue
        try:
            g=shaped(rec.shape.__geo_interface__)
            if not g.is_valid or g.is_empty or g.area<=0:
                invalid+=1
                if lev==6:
                    good=make_valid(g)
                    if good.geom_type=="GeometryCollection":
                        pool=[x for x in good.geoms if x.geom_type in ("Polygon","MultiPolygon") and not x.is_empty]
                        good=unary_union(pool) if pool else good
                    if not good.is_empty and good.is_valid and good.area>0:
                        g=good;repaired+=1
                    else:
                        failed+=1;continue
                else:
                    failed+=1;continue
            ids.append(h);geoms.append(g)
        except Exception:
            failed+=1
    qa={"n_original_source_geometries_after_repair":len(ids),
        "n_invalid":invalid,"n_repaired_L6_only":repaired,
        "n_unusable":failed,"HYBAS_ID_column":True}
    if len(ids)!=len(set(ids)):raise RuntimeError("Official source geometry IDs duplicated")
    if lev==7 and len(ids)+failed!=len(selected):
        raise RuntimeError("Official level7 source missing requested EEA H7")
    if lev==6 and len(ids)<200:
        raise RuntimeError("Official H6 source has too few polygons")
    return ids,geoms,qa

def main(args):
    station=one(args.station_dir,GEO_SOURCE_NAME)
    if digest(station)!=GEO_SOURCE_SHA:raise RuntimeError("Original EEA station source SHA mismatch")
    stationmeta=json.loads(station.read_text())
    rows=stationmeta.get("candidate_source_site_period_rows")
    if not isinstance(rows,list) or len(rows)!=6330:
        raise RuntimeError("Authentic EEA/US station dataset rowcount drift")
    source=collections.defaultdict(set)
    stations=collections.defaultdict(set)
    ct=collections.defaultdict(set)
    for r in rows:
        if r.get("HYBAS_ID_L07") is None or r.get("panel") not in ("EEA18","EEA22"):
            continue
        h=str(r["HYBAS_ID_L07"])
        panel=r["panel"]
        period=str(r.get("period"))
        if (panel=="EEA18" and period!="2017-2019") or (panel=="EEA22" and period!="2020-2022"):
            raise RuntimeError("Source EEA original panel period assignment mismatch")
        source[h].add(panel)
        if r.get("physical_id") is not None:stations[(h,panel)].add(str(r["physical_id"]))
        if r.get("country") is not None:ct[h].add(str(r["country"]))
    if len(source)!=501:raise RuntimeError("Original 501 EEA H7 source opportunity drift")
    seen_both={h for h,p in source.items() if p=={"EEA18","EEA22"}}
    seen_single={h for h,p in source.items() if len(p)==1}
    if len(seen_both)!=234 or len(seen_single)!=267:
        raise RuntimeError("Original EEA 234 both and 267 single period source numbers changed")
    panelpath=args.panel
    if digest(panelpath)!=PANEL_SHA:raise RuntimeError("Original H7 454-panel SHA drift")
    with panelpath.open(encoding="utf-8-sig",newline="") as f:
        panel=list(csv.DictReader(f))
    if len(panel)!=454 or len({(x["h"],x["period"]) for x in panel})!=454:
        raise RuntimeError("Original 454 H7-period observations changed")
    original_h7=set(x["h"] for x in panel)
    original_h6={h:x["h6"] for x in panel for h in [x["h"]]}
    if len(original_h7)!=232 or len(set(original_h6.values()))!=90:
        raise RuntimeError("Original study H7/H6 parent geography changed")
    if not original_h7.issubset(seen_both):
        raise RuntimeError("Original 232 study H7 not subset original EEA both-period H7")
    eligibility=one(args.trait_dir,"EU454_NUTRI_SOURCE_ELIGIBILITY_PER_H7.csv")
    with eligibility.open(encoding="utf-8-sig",newline="") as f:
        el=list(csv.DictReader(f))
    if len(el)!=454 or len({(x["h"],x["period"]) for x in el})!=454:
        raise RuntimeError("EIV original source rows duplicate or missing")
    active=collections.defaultdict(set)
    parent=collections.defaultdict(set)
    country=collections.defaultdict(set)
    for r in el:
        if r["full_source_trait_eligible"].lower()=="true":
            active[r["h"]].add(r["period"])
            parent[r["h"]].add(r["h6"])
            country[r["h"]].add(r["country_holdout"])
    bothscore={h for h,v in active.items() if v==PERIODS}
    if len(bothscore)!=219 or len({next(iter(parent[h])) for h in bothscore})!=86:
        raise RuntimeError("Previously frozen EIV two-period source cohort 219 H7/86 H6 changed")
    inH6=collections.defaultdict(set)
    for h in bothscore:inH6[next(iter(parent[h]))].add(h)
    singleton={h6 for h6,s in inH6.items() if len(s)==1}
    multi={h6 for h6,s in inH6.items() if len(s)>=2}
    if len(singleton)!=56 or len(multi)!=30:
        raise RuntimeError("Original 56 singleton and 30 multi-H7 H6 groups changed")
    # Authentic overlap-based parent mapping. No source responses opened.
    polygons={};report={}
    for lev in (6,7):
        blob,meta=fetch(lev)
        ids,geoms,qa=geom_records(blob,lev,set(source))
        polygons[lev]=(ids,geoms)
        report[str(lev)]={**meta,**qa}
        print("OFFICIAL_HYDROBASINS_GEOMETRY="+json.dumps({"level":lev,**report[str(lev)]},
            ensure_ascii=False,sort_keys=True),flush=True)
    h6,geo6=polygons[6]
    h7,geo7=polygons[7]
    tree=STRtree(geo6)
    mapping={}
    unassigned=[]
    for hid,g in zip(h7,geo7):
        pairs=[]
        for i in tree.query(g,predicate="intersects"):
            j=int(i);part=g.intersection(geo6[j]).area
            if part>0:pairs.append((h6[j],float(part/g.area)))
        pairs.sort(key=lambda q:(-q[1],q[0]))
        top=pairs[0] if pairs else (None,0.)
        second=pairs[1][1] if len(pairs)>1 else 0.
        ok=bool(top[1]>=.95 and top[1]-second>.5)
        mapping[hid]={"parent_h6":top[0] if ok else None,
                     "max_L7_area_inside_L6_fraction":top[1],
                     "second_best_fraction":second,
                     "source_geometry_parent_95pct":ok}
        if not ok:unassigned.append(hid)
    if len(mapping)+report["7"]["n_unusable"]!=501:
        raise RuntimeError("Did not exhaust original 501 EEA H7 geometries")
    # EXACT original repaired 232 H7-to-H6 original ID agreement required.
    mismatch=[{"H7":h,"original_H6":orig,"geometry_H6":mapping.get(h,{}).get("parent_h6")}
         for h,orig in sorted(original_h6.items()) if mapping.get(h,{}).get("parent_h6")!=orig]
    if mismatch:
        raise RuntimeError("Official H6 source mismatch with 232 original repaired parents: "+
            json.dumps(mismatch[:10]))
    verified=collections.defaultdict(list)
    for hid in sorted(seen_single| (seen_both-original_h7),key=int):
        rec=mapping.get(hid)
        if rec is None or not rec["source_geometry_parent_95pct"]:continue
        h6id=rec["parent_h6"]
        if h6id in singleton:
            verified[h6id].append({
                "native_H7":hid,"official_native_H6":h6id,
                "EEA_years":"both" if hid in seen_both else
                     ("2018_only" if source[hid]=={"EEA18"} else "2022_only"),
                "n_physical_EEA18":len(stations[(hid,"EEA18")]),
                "n_physical_EEA22":len(stations[(hid,"EEA22")]),
                "original_source_countries":";".join(sorted(ct[hid])),
                "source_max_area_overlap_fraction":rec["max_L7_area_inside_L6_fraction"],
                "source_second_best_fraction":rec["second_best_fraction"],
                "not_new_EIV_or_HNS_eligible":True})
    number_targets=len(verified)
    both_targets=len([h for h,rs in verified.items() if any(x["EEA_years"]=="both" for x in rs)])
    single_targets=len([h for h,rs in verified.items() if any(x["EEA_years"]!="both" for x in rs)])
    # Candidate selection is determined entirely by source presence, NOT ecology.
    chosen=[]
    allrows=[]
    for target in sorted(verified,key=int):
        candidates=sorted(verified[target],key=lambda x:(
            0 if x["EEA_years"]=="both" else 1,
            -(x["n_physical_EEA18"]+x["n_physical_EEA22"]),int(x["native_H7"])))
        for idx,rec in enumerate(candidates):
            allrows.append({**rec,"source_rank_within_original_H6":idx+1,
                            "source_selected_1_per_parent":idx==0})
        chosen.append(candidates[0])
    if number_targets>56 or any(h not in singleton for h in verified):
        raise RuntimeError("Unrecognized target H6 included")
    reportsummary={
       "schema":"GESUI_EEA_SINGLEERA_H7_OFFICIAL_H6_GEO_SOURCE_OPPORTUNITY_V1",
       "original_6330_source_SHA256":GEO_SOURCE_SHA,"original_454_panel_SHA256":PANEL_SHA,
       "official_geometry_sources":report,
       "original_501_eea_H7_all":len(source),
       "original_234_eea_twoera_H7":len(seen_both),
       "original_267_eea_singleera_H7":len(seen_single),
       "original_232_study_H7_source_parent_matches_geo":True,
       "original_219_twoperiod_plant_eiv_supported_H7":len(bothscore),
       "original_86_twoperiod_plant_eiv_supported_H6":len(inH6),
       "original_single_H7_parent_H6":len(singleton),
       "original_multi_H7_parent_H6":len(multi),
       "original_EEA_H7_geometric_unassigned":len(unassigned),
       "unassigned_source_H7_examples":unassigned[:12],
       "target_existing_singleton_H6_with_geographic_additional_EEA_H7":number_targets,
       "target_singleton_H6_with_at_least_one_additional_EEA_bothperiod":both_targets,
       "target_singleton_H6_with_at_least_one_additional_EEA_singleperiod":single_targets,
       "unique_extra_EEA_H7_candidates_in_target_singleton_H6":len(allrows),
       "source_one_preadvance_candidate_per_target_H6":len(chosen),
       "source_target_geo_gate_min_10_distinct_singleton_H6":number_targets>=10,
       "possible_geographically_multiH7_H6_upper_bound_30_plus_targets":30+number_targets,
       "source_singleera_EEA_is_not_twoera_comparable_2018_and_2022":True,
       "source_H7_GBIF_plant_butterfly_pair_opportunity_not_verified":True,
       "source_geography_result_does_not_pass_trait_or_observed_dose_gate":True,
       "no_EQR_NUTRI_HNS_W_or_butterfly_response_opened":True
    }
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/"EEA_SINGLEERA_H7_OFFICIAL_H6_SOURCE_GEO_GATE.json").write_text(
       json.dumps(reportsummary,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    if not allrows:
        allrows=[{"native_H7":"","official_native_H6":"","EEA_years":"NO_CANDIDATE",
                  "n_physical_EEA18":0,"n_physical_EEA22":0,"original_source_countries":"",
                  "source_max_area_overlap_fraction":None,
                  "source_second_best_fraction":None,
                  "not_new_EIV_or_HNS_eligible":True,
                  "source_rank_within_original_H6":None,
                  "source_selected_1_per_parent":False}]
    with (args.output/"EEA_SINGLEERA_NEW_H7_IN_ORIGINAL_SINGLETON_H6_SOURCE_REGISTRY.csv").open(
           "w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(allrows[0]));w.writeheader();w.writerows(allrows)
    print("GESUI_EEA_SINGLEERA_H6_GEOGRAPHY="+json.dumps(reportsummary,ensure_ascii=False,sort_keys=True),flush=True)
    print("GEO_OPPORTUNITY_GATE="+("PASS_ONLY_GEO_NEED_SEPARATE_GBIF_SOURCE" if number_targets>=10
         else "STOP_NOT_ENOUGH_SINGLETON_H6_SOURCE"),flush=True)

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--station-dir",required=True,type=Path)
    parser.add_argument("--trait-dir",required=True,type=Path)
    parser.add_argument("--panel",type=Path,
        default=Path("results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"))
    parser.add_argument("--output",required=True,type=Path)
    main(parser.parse_args())
