#!/usr/bin/env python3
"""SOURCE ONLY: original full HydroBASINS v1c Europe Level6 vs Level7 H7
native polygon intersection parent mapping. Exactly the same original 232
native catchments from the fixed 454xperiod original field study.

Do not derive source H6 ID from PFAF string prefixes, H3 bins or centroid
neighbors. Do not read biological response coefficients or adjust source
geography based on any ecological signs.
"""
import collections,csv,hashlib,io,json,os,pathlib,statistics,time,urllib.request,zipfile
import shapefile
from shapely import make_valid
from shapely.ops import unary_union
from shapely.geometry import shape as shaped
from shapely.strtree import STRtree

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_HYBAS6_GEOGRAPHY_SOURCE"
ROOT.mkdir(parents=True,exist_ok=True)
EXPECTED_H7_SHA="b9b810f4db271f0fa6b39fecf01784041b948944b3c3f296a393766388e86b44"
GEO_STATION_SHA="a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89"
GEO_SITES=pathlib.Path(os.environ["SOURCE_2018_22_ORIGINAL_NATIVE_H7"])
if hashlib.sha256(GEO_SITES.read_bytes()).hexdigest()!=GEO_STATION_SHA:raise RuntimeError("Original EEA/USA 6330 native source SHA changed")
site_manifest=json.loads(GEO_SITES.read_text())
site_ids={str(q["HYBAS_ID_L07"]) for q in site_manifest["candidate_source_site_period_rows"]
          if q["HYBAS_ID_L07"] is not None and q["panel"] in ("EEA18","EEA22")}
path=pathlib.Path(os.environ["ORIGINAL_454_EEA_GBIF_SOURCE_COHORT"])
with path.open(encoding="utf-8-sig",newline="") as f:
    study=list(csv.DictReader(f))
studyids={z["HYBAS_ID_L07"] for z in study}
if len(study)!=454 or len(studyids)!=232 or \
   len({(z["HYBAS_ID_L07"],z["PERIOD"]) for z in study})!=454 or \
   not studyids.issubset(site_ids):
    raise RuntimeError("Original fixed 454 H7×year river ecology source IDs or population changed")
if set(z["PERIOD"] for z in study)!={"2017-2019","2020-2022"}:
    raise RuntimeError("Original source period changed")
def fetch(lev):
    url=f"https://data.hydrosheds.org/file/hydrobasins/standard/hybas_eu_lev{lev:02d}_v1c.zip"
    err=None
    for k in range(3):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"gesui-europe-HYBAS6-source-geographic-transfer-20261113/1","Accept":"application/zip,*/*"})
            with urllib.request.urlopen(req,timeout=80) as r:
                raw=r.read(75000000)
            if len(raw)>=75000000 or not zipfile.is_zipfile(io.BytesIO(raw)):raise RuntimeError("Not original HydroBASINS source ZIP")
            if lev==7 and hashlib.sha256(raw).hexdigest()!=EXPECTED_H7_SHA:
                raise RuntimeError("Previously source-frozen HydroBASINS7 original ZIP changed")
            return raw,url
        except Exception as e:
            err=e;time.sleep(1+k)
    raise RuntimeError("Authentic HydroBASINS original level"+str(lev)+" source inaccessible "+str(err))
def parse(raw,lev):
    archive=zipfile.ZipFile(io.BytesIO(raw));n=archive.namelist()
    groups={}
    for path in n:
        if f"lev{lev:02d}" in path.lower() and path.rsplit(".",1)[-1].lower() in ("shp","shx","dbf"):
            base,ext=path.rsplit(".",1)
            groups.setdefault(base,{})[ext.lower()]=path
    options=[g for g in groups.values() if set(g)=={"shp","shx","dbf"}]
    if len(options)!=1:raise RuntimeError("Original HydroBASINS original shapefile not exactly one in Level"+str(lev))
    d=options[0]
    sh=shapefile.Reader(shp=io.BytesIO(archive.read(d["shp"])),
      shx=io.BytesIO(archive.read(d["shx"])),
      dbf=io.BytesIO(archive.read(d["dbf"])))
    fields=[x[0] for x in sh.fields[1:]]
    if "HYBAS_ID" not in fields:raise RuntimeError("Original Level "+str(lev)+" source lacks native HYBAS_ID")
    ids=[];geoms=[];invalid=0;repaired=0;repair_failed=0;repaired_ids=set()
    for rec in sh.iterShapeRecords():
        vals=dict(zip(fields,rec.record))
        hid=str(vals["HYBAS_ID"])
        if lev==7 and hid not in studyids:continue
        g=shaped(rec.shape.__geo_interface__)
        if g.is_empty or not g.is_valid or g.area<=0:
            invalid+=1
            # Explicit post-source failure sensitivity. Apply GEOS make_valid
            # ONLY to originally invalid source H6 polygons, never to H7 or
            # already-valid H6 geometry. Preserve original ZIP/CRC and IDs.
            if lev==6:
                try:
                    gg=make_valid(g)
                    if gg.geom_type=="GeometryCollection":
                        polygons=[q for q in gg.geoms if q.geom_type in ("Polygon","MultiPolygon") and not q.is_empty]
                        gg=unary_union(polygons) if polygons else gg
                    if not gg.is_empty and gg.is_valid and gg.area>0:
                        g=gg
                        repaired+=1
                        repaired_ids.add(hid)
                    else:
                        repair_failed+=1
                        continue
                except Exception:
                    repair_failed+=1
                    continue
            else:
                continue
        ids.append(hid);geoms.append(g)
    if lev==7 and (len(ids)!=232 or len(set(ids))!=232):
        raise RuntimeError("All original H7 study targets not found in original valid native source")
    if lev==6 and len(ids)<200:
        raise RuntimeError("Official source Level6 unexpectedly few basin polygons")
    return ids,geoms,{"original_source_polygons_after_repair":len(ids),"original_invalid_source_polygons":invalid,"source_only_repaired_polygons":repaired,"source_only_repair_failed_polygons":repair_failed,"repair_applied_ONLY_invalid_L6_polygons":lev==6,"original_columns":fields,"_repaired_IDs":list(sorted(repaired_ids))}
data={};source_report={}
for lev in [6,7]:
    raw,url=fetch(lev)
    ids,geoms,qa=parse(raw,lev)
    repairedIDs=set(qa.pop("_repaired_IDs"))
    if lev==6:originalRepairedH6=repairedIDs
    source_report[str(lev)]={"original_url":url,"original_zip_SHA256":hashlib.sha256(raw).hexdigest(),
                              "source_bytes":len(raw),**qa}
    data[lev]=(ids,geoms)
    print("original source H",lev,"valid polygons",len(ids),"source_SHA256",source_report[str(lev)]["original_zip_SHA256"],flush=True)
ids6,geom6=data[6];ids7,geom7=data[7]
tree=STRtree(geom6)
record=[]
amb=[]
for h7,g in zip(ids7,geom7):
    ix=tree.query(g,predicate="intersects")
    intersects=[]
    for i in ix:
        k=int(i)
        area=g.intersection(geom6[k]).area
        if area>0:
            intersects.append((ids6[k],float(area/g.area),float(area)))
    intersects.sort(key=lambda p:(-p[1],p[0]))
    if intersects:
        top=intersects[0]
        second=intersects[1][1] if len(intersects)>1 else 0
        ok=bool(top[1]>=.95 and (top[1]-second)>.5)
    else:ok=False;top=(None,0,0);second=0
    rec={
     "HYBAS_ID_L07":h7,
     "official_parent_HYBAS_ID_L06":top[0] if ok else None,
     "original_L7_geometry_to_L6_parent_area_fraction":top[1],
     "original_next_best_L6_overlap_fraction":second,
     "original_H6_candidate_geometries_intersected":len(intersects),
     "source_parent_unambiguously_geometric_95pct":ok,
     "parent_source_was_invalid_and_required_geometry_repair":bool(ok and top[0] in originalRepairedH6),
     "explicit_geometric_sensitivity_post_original_source_gate_failure":True,
     "group_geo_not_outcome_driven":True
    }
    record.append(rec)
    if not ok:amb.append(rec)
if len(record)!=232 or len({x["HYBAS_ID_L07"] for x in record})!=232:
    raise RuntimeError("Original H7 source geometry one-to-one mapping failed")
valid=[r for r in record if r["source_parent_unambiguously_geometric_95pct"]]
groups=collections.Counter(r["official_parent_HYBAS_ID_L06"] for r in valid)
folds=collections.defaultdict(set)
for parent in groups:
    fold=int(hashlib.sha256(parent.encode()).hexdigest()[:8],16)%5
    folds[str(fold)].add(parent)
sourcepass=len(valid)==232 and len(groups)>=30 and all(len(folds[str(i)])>=5 for i in range(5))
out={
 "schema":"GESUI_232_EU_H7_TO_AUTHENTIC_HYBAS6_EXPLICIT_INVALID_GEOMETRY_REPAIR_SOURCE_SENSITIVITY_v2",
 "original_native_H7_source_station_geo_SHA256":GEO_STATION_SHA,
 "fixed_prior_original_454_H7_ecology_rows":454,
 "fixed_prior_original_232_unique_H7":232,
 "source_original_ZIP_metadata":source_report,
 "n_H7_unique_parent_geometrically_assigned_95pct":len(valid),
 "n_H7_unassignable_after_explicit_original_H6_geometry_repair":len(amb),
 "H7_parent_only_possible_through_repaired_H6_n":sum(x["parent_source_was_invalid_and_required_geometry_repair"] for x in valid),
 "original_NO_REPAIR_H6_geography_GATE_FAILED_216_of_232":True,
 "unique_native_Level6_parent_geo_groups":len(groups),
 "Level6_parents_H7_member_range":[min(groups.values()),max(groups.values())] if groups else None,
 "Level6_parent_H7_member_median":statistics.median(groups.values()) if groups else None,
 "Level6_parent_source_folds_n_parents":{str(k):len(v) for k,v in folds.items()},
 "H6_spatial_model_gate_pass_original_full_454_same_sample":sourcepass,
 "source_geometry_unassigned_examples":amb[:15],
 "fixed_model_geography_no_ecological_values_used":True,
 "post_failure_REPAIR_source_only_not_a_preregistered_original_GEOMETRY_test":True
}
payload={
 "schema":out["schema"],"source_original_level6_archive_SHA256":source_report["6"]["original_zip_SHA256"],
 "source_original_level7_archive_SHA256":source_report["7"]["original_zip_SHA256"],
 "original_454_H7_study_basin_geographic_parent_mapping":sorted(record,key=lambda q:q["HYBAS_ID_L07"]),
 "H6_group_5fold_hash_before_new_effects":True
}
raw=(json.dumps(payload,sort_keys=True,separators=(",",":"))+"\n").encode()
out["original_H7_to_H6_parent_mapping_source_SHA256"]=hashlib.sha256(raw).hexdigest()
(ROOT/"HYBAS7_232_TO_OFFICIAL_HYBAS6_SOURCE_PARENT_EXPLICIT_REPAIR_v2.json").write_bytes(raw)
(ROOT/"HYBAS6_LEVEL6_GROUP_EXPLICIT_REPAIR_BEFORE_ECOLOGY_QA.json").write_text(json.dumps(out,indent=2,sort_keys=True))
print("HYBAS6_GEOGRAPHY_REPAIRED_SOURCE_GATE_START")
print(json.dumps(out,sort_keys=True))
print("HYBAS6_GEOGRAPHY_REPAIRED_SOURCE_GATE_END")
if not sourcepass:
    raise RuntimeError("Explicitly source-repaired Level6 parent >95% containment/coverage minimum 232/232 or >=30 unique parent /5 per fold FAILED; no new modeled geography allowed")
