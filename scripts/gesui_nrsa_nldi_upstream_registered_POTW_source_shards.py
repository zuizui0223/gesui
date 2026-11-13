#!/usr/bin/env python3
"""SOURCE-only authentic EPA NRSA native COMID -> USGS NLDI UT/DM -> exact
EPA ECHO permit+feature external municipal outfall. No DMR, TP/TN, bio effects.
Each shard owns deterministic original UID subset, no post-sign selection.
"""
from __future__ import annotations
import argparse,concurrent.futures,hashlib,json,os,pathlib,re,time,urllib.error,urllib.parse,urllib.request
from collections import Counter,defaultdict

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_nrsa_npdes_nldi"
ROOT.mkdir(exist_ok=True,parents=True)
BASE="https://api.water.usgs.gov/nldi/linked-data"
INPUT="NRSA_NPDES_POTW_2ERA_SOURCE_ONLY_GEO25.json"
N_SHARDS=6
def digest(b):return hashlib.sha256(b).hexdigest()
def npdes_key(permit,feature):
 p=str(permit).strip().upper();f=str(feature).strip().upper()
 # Official EPA ICIS 'perm_feature_nmbr' is 3 characters, numeric 1=>001.
 # Non-numeric or longer IDs are never altered; still compare exact ID.
 if re.fullmatch(r"[0-9]{1,3}",f):f=f.zfill(3)
 return p+f
def call(url):
 for attempt in range(4):
  try:
   req=urllib.request.Request(url,headers={"User-Agent":"gesui-USGS-NLDI-NRSA-source-study/2026-11","Accept":"application/json"})
   with urllib.request.urlopen(req,timeout=75) as r:
    raw=r.read(12*1024*1024+1)
    if len(raw)>12*1024*1024:raise RuntimeError("USGS response exceeds bound 12MB")
    obj=json.loads(raw)
    if not isinstance(obj,dict) or obj.get("type")!="FeatureCollection" or not isinstance(obj.get("features"),list):
     raise RuntimeError("Unexpected NLDI native FeatureCollection schema")
    records=[]
    for x in obj["features"]:
     prop=x.get("properties") or {}
     identity=prop.get("identifier")
     if identity is None:continue
     comp=str(prop.get("comid","")).strip()
     records.append({"id":str(identity).upper().strip(),"comid":comp,
                     "measure":prop.get("measure"),"reachcode":prop.get("reachcode")})
    return {"status":r.status,"n_source_features":len(obj["features"]),
            "source_sha256":digest(raw),"records":records}
  except Exception as e:
   if isinstance(e,urllib.error.HTTPError) and e.code in (400,404,410,422):
    return {"status":e.code,"error":str(e),"records":None}
   if attempt==3:return {"status":getattr(e,"code",None),"error":type(e).__name__+":"+str(e)[:180],"records":None}
   time.sleep((attempt+1)*2)
def network_one(site):
 cid=str(site["comid"]).strip()
 if not cid.isdigit():
  return {"uid":site["uid"],"status":"NO_NATIVE_COMID"}
 if not site.get("candidate_outfalls_25km"):
  return {"uid":site["uid"],"status":"NO_REGISTERED_SOURCE_25KM"}
 queries={}
 for nav in ("UT","DM"):
  url=f"{BASE}/comid/{cid}/navigation/{nav}/npdes?distance=50&f=json&excludeGeometry=true"
  queries[nav]=call(url)
 good=all(queries[nav].get("status")==200 and queries[nav].get("records") is not None for nav in ("UT","DM"))
 output={"uid":site["uid"],"site_id":site["site_id"],
         "comid":cid,"huc8":site["huc8"],"state":site["state"],
         "source_geometry_near25":True,
         "source_ECHO_current_registered_2026_not_actual_period_discharge":True,
         "USGS_RIVER_graph_query":{"UT":{"status":queries["UT"].get("status"),"error":queries["UT"].get("error"),"total_provider_npdes":queries["UT"].get("n_source_features"),"sha":queries["UT"].get("source_sha256")},
                                  "DM":{"status":queries["DM"].get("status"),"error":queries["DM"].get("error"),"total_provider_npdes":queries["DM"].get("n_source_features"),"sha":queries["DM"].get("source_sha256")}},
         "status":"NETWORK_SOURCE_API_OK" if good else "NETWORK_API_UNAVAILABLE"}
 if not good:return output
 upstream=defaultdict(list);downstream=defaultdict(list)
 for mode,idx in (("UT",upstream),("DM",downstream)):
  for x in queries[mode]["records"]:
   idx[x["id"]].append(x)
 records=[];source_ids=set()
 for z in site["candidate_outfalls_25km"]:
  permit=z["permit"].upper();feat=z["feature"].upper()
  key=npdes_key(permit,feat)
  U=upstream.get(key,[]);D=downstream.get(key,[])
  if not U and not D:status="NOT_IN_50KM_UP_DOWN_NETWORK"
  elif U and D:status="AMBIGUOUS_IN_BOTH_UT_AND_DM"
  elif U and any(str(u["comid"])==cid for u in U):status="SAME_REACH_UNRESOLVED_ORDER"
  elif U:status="UPSTREAM_DIFFERENT_REACH"
  elif D and any(str(d["comid"])==cid for d in D):status="SAME_REACH_UNRESOLVED_ORDER"
  else:status="DOWNSTREAM_DIFFERENT_REACH"
  # Only physically source IDs available from BOTH independent products count.
  rows={"permit":permit,"feature":feat,"nldi_exact_id":key,"nldi_route_status":status,
        "near10":float(z["geodesic_km"])<=10,
        "geodesic_km":z["geodesic_km"],
        "original_outfall_lat":z["outfall_lat"],"original_outfall_lon":z["outfall_lon"],
        "upstream_original_nldi_comids":sorted({str(q["comid"]) for q in U}),
        "downstream_original_nldi_comids":sorted({str(q["comid"]) for q in D})}
  if status!="NOT_IN_50KM_UP_DOWN_NETWORK":records.append(rows)
  if status=="UPSTREAM_DIFFERENT_REACH":source_ids.add((permit,feat))
 output["candidate_original_ECHO_registered_n"]=len(site["candidate_outfalls_25km"])
 output["matched_NLDI_exact_registered_pairs"]=records
 output["n_verified_different_reach_upstream_candidate"]=len(source_ids)
 output["n_U_geodesic10"]=sum(x["nldi_route_status"]=="UPSTREAM_DIFFERENT_REACH" and x["near10"] for x in records)
 return output
def run(inp,shard):
 src=json.loads(inp.read_text())
 if src.get("schema")!="gesui_NRSA_EPA2026_NPDES_POTW_EXTERNAL_OUTFALL_sourceONLY_2018_19_2023_24":
  raise RuntimeError("Pre-route source geographic gate schema changed")
 if not src["biodiversity_or_nutrient_values_never_read"]:raise RuntimeError("Original geo cohort not source-only")
 both=[]
 for era in ("2018_19","2023_24"):
  if not src["eras"][era]["proceed_to_network"]:raise RuntimeError("Original candidate source geography gate was not passed")
  for s in src["eras"][era]["source_records"]:
   if not s.get("candidate_outfalls_25km"):continue
   if not str(s.get("comid","")).isdigit():continue
   both.append((era,s))
 if len({s["uid"] for era,s in both})!=len(both):raise RuntimeError("Original native UID repeated across eras")
 both.sort(key=lambda x:hashlib.sha256(x[1]["uid"].encode()).hexdigest())
 chosen=[x for i,x in enumerate(both) if i%N_SHARDS==shard]
 print("NLDI_ROUTE_SOURCE",shard,"total_geocandidate_with_COMID",len(both),"shard_n",len(chosen),flush=True)
 if len(chosen)<30:raise RuntimeError("Original shard smaller than source QA minimum")
 results=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
  futures=[ex.submit(network_one,s) for era,s in chosen]
  for i,((era,s),fut) in enumerate(zip(chosen,futures),1):
   z=fut.result()
   z["era"]=era
   results.append(z)
   if i%25==0 or i==len(chosen):
    print("NLDI_PROGRESS",shard,i,len(chosen),"source_200",sum(x["status"]=="NETWORK_SOURCE_API_OK" for x in results),
          "upstream_sites",sum(x.get("n_verified_different_reach_upstream_candidate",0)>0 for x in results),flush=True)
 qa={}
 for era in ("2018_19","2023_24"):
  rr=[r for r in results if r["era"]==era]
  ok=[r for r in rr if r["status"]=="NETWORK_SOURCE_API_OK"]
  ups=[r for r in ok if r.get("n_verified_different_reach_upstream_candidate",0)>0]
  upstream_pairs={(y["permit"],y["feature"]) for r in ups for y in r.get("matched_NLDI_exact_registered_pairs",[])
                  if y["nldi_route_status"]=="UPSTREAM_DIFFERENT_REACH"}
  qa[era]={"queried_original_geo+COMID_sites":len(rr),"USGS_both_UT_DM_success_sites":len(ok),
           "USGS_UT_DM_errors_sites":len(rr)-len(ok),"site_with_confirmed_distinct_reach_upstream_registered_external_outfall":len(ups),
           "site_near10_verified_source":sum(r.get("n_U_geodesic10",0)>0 for r in ups),
           "source_states":len({r["state"] for r in ups}),
           "source_HUC8":len({r["huc8"] for r in ups}),
           "distinct_exact_outfall_features":len(upstream_pairs)}
 return {"schema":"gesui_NRSA_trueNLDI_echoregister_exact_UT_DM_shard_sourceONLY_v1",
         "shard":shard,"N_SHARDS":N_SHARDS,
         "geo_source_sha256":hashlib.sha256(inp.read_bytes()).hexdigest(),
         "original_ECHO_national_zip_sha256":src["EPA_ECHO_official_outfall_zip_sha256"],
         "usgs_nldi_UT_DM_50km_near25km_source_only":True,
         "no_DMR_or_biological_numeric_values_read":True,
         "not_actual_discharge_or_causal_exposure":True,
         "source_qa":qa,"site_route_records":results}
def main():
 a=argparse.ArgumentParser()
 a.add_argument("--geo-source",type=pathlib.Path,required=True)
 a.add_argument("--shard",type=int,required=True)
 args=a.parse_args()
 if args.shard not in range(N_SHARDS):raise RuntimeError("Shard not in frozen 0:6")
 out=run(args.geo_source,args.shard)
 p=ROOT/f"USGS_NLDI_NRSA_POTW_SOURCE_ROUTE_SHARD{args.shard}.json"
 p.write_text(json.dumps(out,indent=2,sort_keys=True))
 print("NLDI_SHARD_SOURCE_QA_RESULT",json.dumps(out["source_qa"],sort_keys=True))
if __name__=="__main__":main()
