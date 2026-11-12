#!/usr/bin/env python3
"""GBIF SOURCE ONLY: 1219 real HydroBASINS L7 polygons × two fixed year bins ×3 taxa.

Never fit ecology, claim photos equal plant richness, alter polygons/years after
counts, or treat HTTP missing as true biological 0. Source-polygon stage0 is an
immutable separate archived artifact. Preflight 12 basins then full four shards.
"""
import argparse,collections,concurrent.futures,hashlib,json,os,pathlib,time,urllib.parse,urllib.request,urllib.error

ap=argparse.ArgumentParser()
ap.add_argument("--mode",choices=("preflight","full"),required=True)
ap.add_argument("--shard",type=int,choices=(0,1,2,3),default=0)
args=ap.parse_args()
p=pathlib.Path(os.environ["HYBAS7_WKT_PRECOUNT_SOURCE"])
expect="47f597e0a798cd6155ee9fb69ce8916ed5629951d046d114feb9f667963310e4"
original=p.read_bytes()
if hashlib.sha256(original).hexdigest()!=expect:
    raise RuntimeError("Original HydroBASINS 1219 WKT source manifest changed")
j=json.loads(original)
a=j["basins"]
if len(a)!=1219 or not j["geography_source_gate_pass"]:
    raise RuntimeError("Source full monitored Level7 basin cohort changed")
allvalid=[b for b in a if b["status"]=="all_original_polygons_WKT_queryable"]
if len(allvalid)!=1219:raise RuntimeError("Original 1219 basin polygon source cohort drifted")
if args.mode=="preflight":
    chosen=[]
    for region in ("eu","na"):
        v=[b for b in allvalid if b["region"]==region]
        # This fixes the pilot without looking at original EQR or current GBIF.
        v=sorted(v,key=lambda b:hashlib.sha256(("PREFLIGHT_NATIVE_L7_V1|"+b["HYBAS_ID_L07"]).encode()).hexdigest())
        chosen.extend(v[:6])
else:
    chosen=[b for b in allvalid if b["original_source_shard"]==args.shard]
    chosen.sort(key=lambda b:b["HYBAS_ID_L07"])
if not chosen:raise RuntimeError("Empty fixed original source subset")
TAXA={"Plantae":6,"Lycaenidae":5473,"Animalia":1}
PERIODS={"2017-2019":"2017,2019","2020-2022":"2020,2022"}
CHECKLIST="d7dddbf4-2cf0-4f39-9b2a-bb099caae36c"
def request(url,timeout=60):
    req=urllib.request.Request(url,headers={"User-Agent":"gesui-HydroBASINS-L7-original-macro-source-20261112/1","Accept":"application/json"})
    with urllib.request.urlopen(req,timeout=timeout) as rr:
        return rr.read(),rr.url
if args.mode=="preflight":
    # Verify names/ranks before any actual GBIF index counts.
    for label,key in TAXA.items():
        u="https://api.gbif.org/v1/species/"+str(key)
        raw,actual=request(u)
        v=json.loads(raw)
        scientific=str(v.get("scientificName","")).lower()
        expected_rank="FAMILY" if label=="Lycaenidae" else "KINGDOM"
        if label.lower() not in scientific or v.get("rank")!=expected_rank:
            raise RuntimeError("Original GBIF taxonomy key not correct "+label+":"+str(v)[:1000])
        print("GBIF_PRE_COUNTS_TAXON_ID_VERIFIED",label,key,expected_rank,"responseSHA256",hashlib.sha256(raw).hexdigest(),flush=True)
def api(b,period,label,key):
    counts=[]
    requests=[]
    for component in b["polygon_parts"]:
        wk=component["wkt"]
        if hashlib.sha256(wk.encode()).hexdigest()!=component["wkt_SHA256"]:
            raise RuntimeError("Original WKT source hash changed")
        query={
           "limit":"0",
           "geometry":wk,
           "year":PERIODS[period],
           "taxonKey":str(key),
           "hasCoordinate":"true",
           "hasGeospatialIssue":"false",
           "occurrenceStatus":"PRESENT"
        }
        if label=="Lycaenidae":query["checklistKey"]=CHECKLIST
        url="https://api.gbif.org/v1/occurrence/search?"+urllib.parse.urlencode(query,safe="(),")
        ok=False;err=None;httpcode=None
        for attempt in range(4):
            try:
                raw,actual=request(url,timeout=65)
                ans=json.loads(raw)
                c=ans.get("count")
                if not isinstance(c,int) or c<0:raise RuntimeError("GBIF original integer count missing or invalid")
                counts.append(c)
                requests.append({"url":url,"original_component_sha256":component["wkt_SHA256"],
                     "count":c,"status":"ok","response_sha256":hashlib.sha256(raw).hexdigest()})
                ok=True
                break
            except urllib.error.HTTPError as e:
                httpcode=e.code
                err=type(e).__name__+" "+str(e)[:180]
                # Retry on temporary server overload only; malformed geometry
                # is a source quality fail, not a missing biological species.
                if e.code not in (400,429,500,502,503,504):
                    break
                time.sleep([2,5,10,18][attempt])
            except Exception as e:
                err=type(e).__name__+" "+str(e)[:180]
                time.sleep([1,3,7,12][attempt])
        if not ok:
            requests.append({"url":url,"original_component_sha256":component["wkt_SHA256"],
               "count":None,"status":"missing","http_code":httpcode,"error":err})
            break
    if len(counts)==len(b["polygon_parts"]):
        return {"HYBAS_ID_L07":b["HYBAS_ID_L07"],"region":b["region"],
         "period":period,"taxon":label,"status":"ok","original_polygon_count":len(b["polygon_parts"]),
         "GBIF_index_record_count":sum(counts),"request_components":requests}
    return {"HYBAS_ID_L07":b["HYBAS_ID_L07"],"region":b["region"],
        "period":period,"taxon":label,"status":"missing","original_polygon_count":len(b["polygon_parts"]),
        "GBIF_index_record_count":None,"request_components":requests}
tasks=[(b,period,label,key) for b in chosen for period in PERIODS for label,key in TAXA.items()]
print("GBIF_FIXED_SOURCE_TASKS",args.mode,args.shard,len(chosen),"basins",len(tasks),"basin_period_taxon_tasks",flush=True)
results=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
    fut=[pool.submit(api,*t) for t in tasks]
    for i,f in enumerate(concurrent.futures.as_completed(fut)):
        results.append(f.result())
        if (i+1)%100==0:print("GBIF source count queries completed",i+1,"/",len(fut),flush=True)
results.sort(key=lambda x:(x["HYBAS_ID_L07"],x["period"],x["taxon"]))
if len(results)!=len(tasks):raise RuntimeError("All fixed source tasks did not complete")
ok=collections.Counter(x["taxon"] for x in results if x["status"]=="ok")
failure=collections.Counter(str(c.get("http_code",c.get("error","unknown"))) for x in results
         if x["status"]!="ok" for c in x["request_components"] if c["status"]!="ok")
keys=collections.defaultdict(dict)
for x in results:
    k=(x["HYBAS_ID_L07"],x["period"])
    if x["taxon"] in keys[k]:raise RuntimeError("Duplicate source taxonomy key")
    keys[k][x["taxon"]]=x
complete=[(hid,period) for (hid,period),tax in keys.items()
          if len(tax)==3 and all(r["status"]=="ok" for r in tax.values())]
qualified=[k for k in complete if keys[k]["Plantae"]["GBIF_index_record_count"]>=50
           and keys[k]["Lycaenidae"]["GBIF_index_record_count"]>=10
           and keys[k]["Animalia"]["GBIF_index_record_count"]>=50]
out={
 "schema":"HYBAS7_1219_GBIF_PLANTAE_LYCAENIDAE_ANIMALIA_TWO_PERIOD_INDEX_SOURCE_"+args.mode,
 "source_wkt_manifest_SHA256":expect,
 "fixed_original_H7_eligible_basins":1219,
 "mode":args.mode,"shard":args.shard,
 "fixed_source_shard_basin_count":len(chosen),
 "fixed_source_basin_periods":len(chosen)*2,
 "fixed_source_requests":len(tasks),
 "GBIF_source_successful_requests":sum(ok.values()),
 "GBIF_source_failed_requests":len(tasks)-sum(ok.values()),
 "GBIF_source_success_by_taxon":dict(ok),
 "error_code_counts":dict(failure),
 "full_all3_basin_period_source_success":len(complete),
 "three_taxa_all_present_source_count_tier":len(qualified),
 "qualified_unique_H7_count":len({z[0] for z in qualified}),
 "successful_query_source_count_fraction":sum(ok.values())/len(tasks),
 "GBIF_original_taxon_keys":TAXA,
 "fixed_periods":PERIODS,
 "WKT_geo_is_original_river_basin_not_port_radius":True,
 "GBIF_record_count_IS_NOT_species_richness_or_true_population":True,
 "no_cross_layer_association_or_ecosystem_service_estimate":True,
 "source_not_collected_404_429_or_API_error_IS_NOT_biological_zero":True,
 "source_parts_with_multiple_geometries_sum_may_duplicate_boundary_records":True,
 "global_source_gate_NOT_EVALUATED_until_all_4_shards_successful":True
}
root=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/("gesui_H7_GBIF_"+args.mode+"_shard"+str(args.shard))
root.mkdir(parents=True,exist_ok=True)
(root/"HYBAS7_GBIF_ORIGINAL_SOURCE_COUNTS_URL_SHA.json").write_text(json.dumps(results,indent=1))
(root/"HYBAS7_GBIF_ORIGINAL_SOURCE_SUPPORT_QA.json").write_text(json.dumps(out,indent=2,sort_keys=True))
print("HYBAS7_GBIF_SOURCE_COUNTS_RESULT_START",flush=True)
print(json.dumps(out,sort_keys=True),flush=True)
print("HYBAS7_GBIF_SOURCE_COUNTS_RESULT_END",flush=True)
if args.mode=="preflight" and (len(complete)<int(.8*len(chosen)*2)):
    raise RuntimeError("Original polygon/taxon source preflight capacity below 80%; do not launch full source without source-only diagnosis")
