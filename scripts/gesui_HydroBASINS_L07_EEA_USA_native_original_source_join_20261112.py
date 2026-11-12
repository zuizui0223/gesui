#!/usr/bin/env python3
"""SOURCE-ONLY true HydroBASINS Level7 catchment crosswalk for the original
EU EEA 2018/2022 biology sites and US EPA NRSA 2018/19 original source sites.

Independent site/visit and period keys; never fit ecological outcome
coefficients or pretend Level7 shared polygon implies upstream nutrient
delivery. Output candidate basin×period source support and native site joins.
"""
import collections,csv,hashlib,io,json,math,os,pathlib,statistics,time,urllib.request,zipfile
import shapefile
from shapely.geometry import Point,shape as shape_from_geo
from shapely.strtree import STRtree

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_hybas7_macro_source"
ROOT.mkdir(parents=True,exist_ok=True)
FILES={
 "EEA18":("results/eea2018_independent_temporal_bqe_20261108/benthic_invertebrates_2018_WWTP_source_real_EQR.csv","0f30bfbc3dac13f6a182c8a910e9899c04778322"),
 "EEA22":("results/eea2022_biology_bqe_20261108/BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv","325f199e985d9862ffee5b0edb205c61a7be3833"),
}
US_SHA="fbf53b8379d4e6564dbb51b3705b02fdd03e4a2caacad06a2314d643c8412e91"
def gitsha(raw):return hashlib.sha1(("blob "+str(len(raw))+"\0").encode()+raw).hexdigest()
points=[]
for group,(path,sha) in FILES.items():
    raw=pathlib.Path(path).read_bytes()
    if gitsha(raw)!=sha:raise RuntimeError("Original EEA source SHA mismatch "+group)
    data=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    if len(data)!={"EEA18":2634,"EEA22":2475}[group]:raise RuntimeError("Original EEA cohort changed")
    for r in data:
        points.append({"panel":group,
           "native_id":r["site_key"],
           "physical_id":r["site_key"],
           "country":r["countryCode"],
           "period":"2017-2019" if group=="EEA18" else "2020-2022",
           "lat":float(r["lat"]),"lon":float(r["lon"])})
p=pathlib.Path(os.environ["EPA_FROZEN_MANIFEST"])
v=p.read_bytes()
if hashlib.sha256(v).hexdigest()!=US_SHA:raise RuntimeError("Original NRSA 1,221 source manifest sha mismatch")
usa=json.loads(v)["primary_berw_UID_subset"]
assert len(usa)==1221
for r in usa:
    points.append({"panel":"US_NRSA18_19",
     "native_id":r["uid"],"physical_id":r["site_id"],"country":"US-"+r["state"],
     "period":"2017-2019","lat":float(r["lat"]),"lon":float(r["lon"])})
assert len(points)==6330
for p in points:
    if not (-90<=p["lat"]<=90 and -180<=p["lon"]<=180):
        raise RuntimeError("Invalid original coordinate "+p["panel"])
def download(code):
    # Official HydroSHEDS HydroBASINS repository; one Level7 polygon zip per
    # continent only. If unavailable, report source failure, no guessed joins.
    urls=[
       f"https://data.hydrosheds.org/file/hydrobasins/standard/hybas_{code}_lev07_v1c.zip",
       f"https://data.hydrosheds.org/file/HydroBASINS/standard/hybas_{code}_lev07_v1c.zip"
    ]
    attempts=[]
    for url in urls:
        for attempt in range(2):
            try:
                req=urllib.request.Request(url,headers={"Accept":"application/zip,application/octet-stream,*/*","User-Agent":"gesui-HydroBASINS-original-hydrological-source-audit/1"})
                with urllib.request.urlopen(req,timeout=75) as q:
                    b=q.read(100000000)
                    final=q.url
                if len(b)>100000000 or not zipfile.is_zipfile(io.BytesIO(b)):
                    raise RuntimeError("Not an original valid HydroBASINS shapefile ZIP")
                zipfile.ZipFile(io.BytesIO(b)).testzip()
                return {"data":b,"url":url,"final":final,"sha256":hashlib.sha256(b).hexdigest(),"nbytes":len(b)}
            except Exception as e:
                attempts.append(type(e).__name__+":"+str(e)[:120])
                time.sleep(1+attempt)
    raise RuntimeError("Original official HydroBASINS level7 "+code+" source inaccessible; do not fabricate basin IDs. "+repr(attempts))
def spatial_index(zipped,code):
    z=zipfile.ZipFile(io.BytesIO(zipped["data"]))
    names=z.namelist()
    bases={}
    for path in names:
        if path.lower().endswith((".shp",".shx",".dbf")) and "lev07" in path.lower():
            base=path.rsplit(".",1)[0]
            ext=path.rsplit(".",1)[1].lower()
            bases.setdefault(base,{})[ext]=path
    opts=[x for x in bases.values() if {"shp","shx","dbf"}.issubset(x)]
    if len(opts)!=1:raise RuntimeError("HydroBASINS official level 7 shapefile not unique "+str(names[:80]))
    e=opts[0]
    reader=shapefile.Reader(shp=io.BytesIO(z.read(e["shp"])),
                          shx=io.BytesIO(z.read(e["shx"])),
                          dbf=io.BytesIO(z.read(e["dbf"])))
    cols=[x[0] for x in reader.fields[1:]]
    if not {"HYBAS_ID","SUB_AREA","UP_AREA"}.issubset(set(cols)):
        raise RuntimeError("Unexpected HydroBASINS source headers "+str(cols))
    geoms=[]; ids=[];area=[];up=[];extra=[]
    for sr in reader.iterShapeRecords():
        rec=dict(zip(cols,sr.record))
        geom=shape_from_geo(sr.shape.__geo_interface__)
        if geom.is_empty:continue
        if not geom.is_valid:
            # Do not silently repair genuine source geometry. Log count and
            # skip invalid polygons, plus source mapping failure.
            extra.append(str(rec["HYBAS_ID"]));continue
        geoms.append(geom);ids.append(str(rec["HYBAS_ID"]))
        area.append(float(rec["SUB_AREA"]));up.append(float(rec["UP_AREA"]))
    if len(geoms)<500 or len(set(ids))!=len(ids):
        raise RuntimeError("HydroBASINS geometry extent or unique basin ID invalid for "+code)
    return STRtree(geoms),geoms,ids,area,up,{
        "geographical_region":code,
        "zip_source_url":zipped["url"],"source_sha256":zipped["sha256"],
        "source_bytes":zipped["nbytes"],"basin_shapes":len(geoms),
        "bad_geometry_count":len(extra),"source_columns":cols,
        "area_km2_min":min(area),"area_km2_median":statistics.median(area),
        "area_km2_max":max(area)}
print("Fetching authentic HydroBASINS Level7 EU and NA region polygons",flush=True)
regions={}
source=[]
for code in ("eu","na"):
    blob=download(code)
    tree,geoms,ids,area,up,summary=spatial_index(blob,code)
    regions[code]={"tree":tree,"geom":geoms,"ids":ids,"area":area,"up":up}
    source.append(summary)
    print("Original Level7 region "+code+" polygons="+str(summary["basin_shapes"])+" zipSHA="+summary["source_sha256"],flush=True)
result=collections.Counter()
assigned=[]
for r in points:
    pt=Point(r["lon"],r["lat"])
    found=[]
    for code,panel in regions.items():
        hits=panel["tree"].query(pt,predicate="intersects")
        # no ID cherry pick; all candidates are reported before any biology.
        for ind in hits:
            geom=panel["geom"][int(ind)]
            found.append((code,panel["ids"][int(ind)],bool(geom.contains(pt)),
                          panel["area"][int(ind)],panel["up"][int(ind)]))
    if len(found)!=1:
        if len(found)>1:
            inside=[z for z in found if z[2]]
            if len(inside)==1:found=inside
    if len(found)==1:
        code,hid,strict,area,up=found[0]
        assigned.append({**r,"HYBAS_ID_L07":hid,"region":code,
                         "sub_area_km2":area,"up_area_km2":up})
        result["matched_"+r["panel"]]+=1
    else:
        result["unmatched_"+r["panel"]]+=1
        if len(found)>1:result["ambiguous_"+r["panel"]]+=1
        assigned.append({**r,"HYBAS_ID_L07":None,"region":None,"sub_area_km2":None,"up_area_km2":None})
# A fixed, data-source-unrelated coverage gate. No effect signs read.
for group,n in [("EEA18",2634),("EEA22",2475),("US_NRSA18_19",1221)]:
    if result["matched_"+group]/n<.85:
        raise RuntimeError("A priori >=85% geography match gate failed for "+group+": "+str(result["matched_"+group])+"/"+str(n))
num=lambda ds:{x["HYBAS_ID_L07"] for x in assigned if x["panel"]==ds and x["HYBAS_ID_L07"] is not None}
a,b,c=(num(t) for t in ["EEA18","EEA22","US_NRSA18_19"])
site_map=collections.defaultdict(dict)
for z in assigned:
    if z["panel"] in ("EEA18","EEA22"):
        site_map[z["native_id"]][z["panel"]]=z["HYBAS_ID_L07"]
same=[v for v in site_map.values() if "EEA18" in v and "EEA22" in v]
n_stable=sum(x["EEA18"]==x["EEA22"] and x["EEA18"] is not None for x in same)
n_move=sum(x["EEA18"]!=x["EEA22"] and x["EEA18"] is not None and x["EEA22"] is not None for x in same)
out={
 "schema":"gesui_GlobalHydroBASINS_level7_EEA18_EEA22_USA18_19_SOURCE_ONLY_v1",
 "hydrobains_v1c_official_L7_zip_sources":source,
 "time_bins":["2017-2019","2020-2022"],
 "original_sources":{k:{"path":path,"GitBlobSHA1":sha} for k,(path,sha) in FILES.items()},
 "US_pre_numeric_UID_source_SHA256":US_SHA,
 "original_rows":len(points),"native_locations_by_panel":{"EEA18":2634,"EEA22":2475,"US_NRSA18_19":1221},
 "HYBAS7_geolocation_match_counters":dict(result),
 "basins_per_panel":{"EEA18":len(a),"EEA22":len(b),"USA":len(c)},
 "EEA_2018_2022_same_H7_basins":len(a&b),
 "EU_and_US_H7_overlap":len((a|b)&c),
 "EEA_repeat_monitor_sites_in_both_years":len(same),
 "EEA_repeat_sites_same_valid_H7":n_stable,
 "EEA_repeat_sites_changed_H7":n_move,
 "same_source_geo_quality_stop_threshold":.85,
 "no_biological_EQR_EPT_plant_or_butterfly_regression":True,
 "HYBAS7_is_local_subbasin_NOT_total_upstream_loaded_catchment":True,
 "river_upstream_source_must_be_routed_to_actual_monitoring_reach":True,
 "ecological_master_plant_butterfly_joint_bins_not_yet_acquired":True
}
# Only source geometry + exact monitoring ID + HYBAS ID. Original aquatic
# response values have NOT been read or exported in the new basin manifest.
manifest={"schema":"gesui_EEA_EPA_2018_2022_HYBAS7_source_geo_candidates_v1",
 "source_HydroBASINS_region_SHA256":{z["geographical_region"]:z["source_sha256"] for z in source},
 "candidate_source_site_period_rows":assigned,"biological_effects_not_examined":True}
mraw=(json.dumps(manifest,sort_keys=True,separators=(",",":"))+"\n").encode()
out["manifest_sha256"]=hashlib.sha256(mraw).hexdigest()
(ROOT/"EEA2018_EEA2022_EPA1819_NATIVE_HYBAS7_SOURCE_ONLY_GEOGRAPHY.json").write_bytes(mraw)
(ROOT/"GLOBAL_HYBAS7_EEA_EPA_SOURCE_SUPPORT_REPORT.json").write_text(json.dumps(out,indent=2,sort_keys=True))
print("HYBAS7_RIVER_SOURCE_SUPPORT_START",flush=True)
print(json.dumps(out,sort_keys=True),flush=True)
print("HYBAS7_RIVER_SOURCE_SUPPORT_END",flush=True)
