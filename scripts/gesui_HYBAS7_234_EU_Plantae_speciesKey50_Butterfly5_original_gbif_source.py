#!/usr/bin/env python3
"""Actual speciesKey-faceted GBIF plant 50-record rarefaction and FIVE butterfly
family observed-record composition in 234 source-frozen repeated EU EEA river
basins x 2017-19/2020-22; source only, NO EQR outcomes or WWTP regression.

New source cohort and GBIF 5-family keys were frozen in previous independent
GH Actions #43204064427 BEFORE any new species facet/family numerical result.
"""
import argparse, collections, concurrent.futures, hashlib, json, math, os, pathlib, time, urllib.parse, urllib.request, urllib.error

p=argparse.ArgumentParser()
p.add_argument("--mode",choices=("pilot","full"),required=True)
p.add_argument("--shard",type=int,choices=(0,1,2,3),default=0)
args=p.parse_args()
SOURCE=pathlib.Path(os.environ["H7_234_PRE_NEW_RESPONSE_SOURCE"])
SOURCE_SHA="afc524ebfa353e4d1c6ea51d208ba529a85fedb57c4d0f28d83545010bb76189"
if hashlib.sha256(SOURCE.read_bytes()).hexdigest()!=SOURCE_SHA:
    raise RuntimeError("Original 234 H7×2 source-only cohort+butterfly family source SHA changed")
j=json.loads(SOURCE.read_text())
basins=j["original_two_year_EEA_H7_subbasins"]
families=j["verified_butterfly_family_taxa"]
if len(basins)!=234 or set(families)!={"Lycaenidae","Nymphalidae","Pieridae","Papilionidae","Hesperiidae"}:
    raise RuntimeError("Pre-GBIF outcome real source basins or family definition changed")
if j["fixed_periods"]!=["2017-2019","2020-2022"]:
    raise RuntimeError("Original two periods were changed")
for x in basins:
    if x["status"]!="all_original_polygons_WKT_queryable":raise RuntimeError("Original native WKT missing")
for f,tax in families.items():
    if tax["rank"]!="FAMILY" or tax["scientificName"]!=f:raise RuntimeError("Original original taxonomic key drift")
if families["Lycaenidae"]["taxonKey"]!=5473:raise RuntimeError("Legacy butterfly key changed")
BACKBONE="d7dddbf4-2cf0-4f39-9b2a-bb099caae36c"
periods={"2017-2019":"2017,2019","2020-2022":"2020,2022"}
if args.mode=="pilot":
    # Source-only deterministic fixed prior to NEW species facet/family source value.
    b=sorted(basins,key=lambda x:hashlib.sha256(("H7_EU_SPECIES_FAMILY_PILOT_v1|"+x["HYBAS_ID_L07"]).encode()).hexdigest())[:6]
else:
    b=[x for x in basins if int(hashlib.sha256(x["HYBAS_ID_L07"].encode()).hexdigest()[:8],16)%4==args.shard]
b=sorted(b,key=lambda x:x["HYBAS_ID_L07"])
print("SOURCE stage",args.mode,"shard",args.shard,"basins",len(b),flush=True)
if not b:raise RuntimeError("Empty source geography shard")
def request(url):
    req=urllib.request.Request(url,headers={"Accept":"application/json","User-Agent":"gesui-H7-234-repeat-EQR-GBIF-trait-source-20261112/1"})
    with urllib.request.urlopen(req,timeout=75) as r:return r.read(),r.url
def params(component,period,field):
    wk=component["wkt"]
    if hashlib.sha256(wk.encode()).hexdigest()!=component["wkt_SHA256"]:
        raise RuntimeError("Source polygon geometry changed")
    q={
      "limit":"0","geometry":wk,"year":periods[period],
      "hasCoordinate":"true","hasGeospatialIssue":"false","occurrenceStatus":"PRESENT",
    }
    if field=="Plantae_speciesKey_facet":
        q["taxonKey"]="6"
        q["facet"]="speciesKey";q["facetLimit"]="1000"
    else:
        q["taxonKey"]=str(families[field]["taxonKey"])
        q["checklistKey"]=BACKBONE
    return q
def fetchcounts(q):
    url="https://api.gbif.org/v1/occurrence/search?"+urllib.parse.urlencode(q,safe="(),")
    error=None;status=None
    for i in range(4):
        try:
            raw,actual=request(url)
            d=json.loads(raw)
            n=d.get("count")
            if not isinstance(n,int) or n<0:
                raise RuntimeError("Source GBIF original indexed count invalid")
            return {"ok":True,"count":n,"response_sha256":hashlib.sha256(raw).hexdigest(),
                     "data":d,"url":url}
        except urllib.error.HTTPError as e:
            status=e.code;error=type(e).__name__+":"+str(e)[:90]
            if status not in (400,429,500,502,503,504):break
            time.sleep([2,4,8,12][i])
        except Exception as e:
            error=type(e).__name__+":"+str(e)[:120]
            time.sleep([1,2,4,7][i])
    return {"ok":False,"url":url,"HTTP_status":status,"error":error}
def getfield(basin,period,field):
    sourcecounts=[]
    receipts=[]
    counts=collections.Counter()
    complete=True;truncated=False;facet_source_coverage=None
    for component in basin["polygon_parts"]:
        base=params(component,period,field)
        if field=="Plantae_speciesKey_facet":
            offset=0;all_species={}
            while offset<=20000:
                req=dict(base);req["facetOffset"]=str(offset)
                ret=fetchcounts(req)
                if not ret["ok"]:
                    complete=False
                    receipts.append({"url":ret["url"],"status":"missing","error":ret.get("error"),
                       "HTTP_status":ret.get("HTTP_status"),"original_component_SHA256":component["wkt_SHA256"]})
                    break
                d=ret["data"]
                facets=d.get("facets") or []
                matching=[z for z in facets if str(z.get("field","")).upper().replace("_","")=="SPECIESKEY"]
                if len(matching)!=1:
                    complete=False
                    receipts.append({"url":ret["url"],"status":"unrecognized_original_speciesKey_facet",
                      "GBIF_facets_returned":[z.get("field") for z in facets],
                      "original_component_SHA256":component["wkt_SHA256"]})
                    break
                entries=matching[0].get("counts") or []
                for item in entries:
                    key=str(item.get("name",""))
                    n=item.get("count")
                    if not key or not isinstance(n,int) or n<0 or key in all_species:
                        raise RuntimeError("Original speciesKey facet has duplicate/unparsable taxon value")
                    all_species[key]=n
                receipts.append({"url":ret["url"],"status":"ok",
                   "source_record_count":ret["count"],"response_sha256":ret["response_sha256"],
                   "returned_facet_count":len(entries),"source_facetOffset":offset,
                   "original_component_SHA256":component["wkt_SHA256"]})
                if offset==0:sourcecounts.append(ret["count"])
                if len(entries)<1000:break
                offset+=1000
                if offset>=20000:truncated=True;complete=False;break
            if not complete:break
            for species,n in all_species.items():counts[species]+=n
        else:
            ret=fetchcounts(base)
            if not ret["ok"]:
                complete=False
                receipts.append({"url":ret["url"],"status":"missing","HTTP_status":ret.get("HTTP_status"),
                    "error":ret.get("error"),"original_component_SHA256":component["wkt_SHA256"]})
                break
            sourcecounts.append(ret["count"])
            receipts.append({"url":ret["url"],"status":"ok","source_record_count":ret["count"],
                    "response_sha256":ret["response_sha256"],
                    "original_component_SHA256":component["wkt_SHA256"]})
    res={"HYBAS_ID_L07":basin["HYBAS_ID_L07"],"period":period,"source_field":field,
          "original_components":len(basin["polygon_parts"]),
          "source_requests":receipts,"status":"ok" if complete else "missing",
          "source_record_index_count":sum(sourcecounts) if complete else None}
    if field=="Plantae_speciesKey_facet":
        if complete:
            N=sum(counts.values())
            reported=sum(sourcecounts)
            coverage=N/reported if reported else 0.
            if coverage>1.000001:
                raise RuntimeError("Original GBIF Plantae speciesKey facet counts exceed source original records")
            source_s50=None
            if N>=50 and coverage>=.70 and len(counts)>=10:
                def missed(n_i):
                    # Calculate probability none of species i in sample of 50
                    if N-n_i<50:return 0.
                    return math.exp(math.lgamma(N-n_i+1)-math.lgamma(N-n_i-50+1)-
                               math.lgamma(N+1)+math.lgamma(N-50+1))
                source_s50=sum(1-missed(n) for n in counts.values())
            res.update({"indexed_original_speciesKey_distinct":len(counts),
                "speciesKey_original_record_sum":N,
                "speciesKey_index_record_coverage_fraction":coverage,
                "species_50_draw_rarefied_richness":source_s50,
                "original_rarefaction_species_eligible":source_s50 is not None,
                "original_species_source_facets_complete":not truncated})
        else:res.update({"indexed_original_speciesKey_distinct":None,
            "species_50_draw_rarefied_richness":None,
            "original_species_source_facets_complete":False})
    return res
fields=["Plantae_speciesKey_facet","Lycaenidae","Nymphalidae","Pieridae","Papilionidae","Hesperiidae"]
tasks=[(x,p,f) for x in b for p in periods for f in fields]
res=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
    futures=[ex.submit(getfield,*task) for task in tasks]
    for i,fu in enumerate(concurrent.futures.as_completed(futures)):
        res.append(fu.result())
        if (i+1)%120==0:print("Original H7 source species/family tasks",i+1,"/",len(futures),flush=True)
res.sort(key=lambda x:(x["HYBAS_ID_L07"],x["period"],x["source_field"]))
if len(res)!=len(tasks):raise RuntimeError("Frozen original jobs incomplete")
by=collections.defaultdict(dict)
miss=collections.Counter()
for r in res:
    k=(r["HYBAS_ID_L07"],r["period"])
    if r["source_field"] in by[k]:raise RuntimeError("Repeated H7×period×guild source key")
    by[k][r["source_field"]]=r
    if r["status"]!="ok":miss[r["source_field"]]+=1
site_periods=len(b)*2
plant_complete=sum(d["Plantae_speciesKey_facet"]["status"]=="ok" for d in by.values())
rarefied=sum(d["Plantae_speciesKey_facet"].get("species_50_draw_rarefied_richness") is not None for d in by.values())
family_complete=sum(all(d[f]["status"]=="ok" for f in fields[1:]) for d in by.values())
joint=sum(all(d[f]["status"]=="ok" for f in fields) and
           d["Plantae_speciesKey_facet"].get("species_50_draw_rarefied_richness") is not None and
           sum(d[f]["source_record_index_count"] for f in fields[1:])>=20
           for d in by.values())
summary={
 "schema":"HYBAS7_234_EU_BOTH_EEA_SOURCE_PLANT_S50_AND_FIVE_BUTTERFLY_FAMILIES_original_v1",
 "fixed_original_source_and_taxon_manifest_sha256":SOURCE_SHA,
 "mode":args.mode,"shard":args.shard,
 "original_fixed_EU_basins_this_shard":len(b),
 "original_expected_fixed_basin_periods":site_periods,
 "original_expected_field_queries":len(tasks),
 "source_facet_plant_complete":plant_complete,
 "source_plant_rarefied_50_eligible":rarefied,
 "all_five_butterfly_families_request_complete":family_complete,
 "joint_plant_S50_and_butterfly20_genuine_source_support":joint,
 "source_request_missing_by_field":dict(miss),
 "no_station_EQR_outcome_opened":True,
 "no_wastewater_or_service_coefficient_calculated":True,
 "plant_richness_is_rareified_speciesKey_observation_not_true field flora":True,
 "butterfly5_denominator_is_five_family_index_occurrence_NOT population":True
}
folder=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/("gesui_H7_234_EU_"+args.mode+"_shard"+str(args.shard))
folder.mkdir(parents=True,exist_ok=True)
(folder/"HYBAS7_234EU_NATIVE_PLANT_S50_BUTTERFLY5_ORIGINAL_SOURCE_RESULTS.json").write_text(json.dumps(res,indent=1))
(folder/"HYBAS7_234EU_NATIVE_SOURCE_JOINT_ECOLOGY_QA.json").write_text(json.dumps(summary,indent=2,sort_keys=True))
print("HYBAS7_234_EU_SPECIES_BUTTERFLY_SOURCE_START")
print(json.dumps(summary,sort_keys=True))
print("HYBAS7_234_EU_SPECIES_BUTTERFLY_SOURCE_END")
