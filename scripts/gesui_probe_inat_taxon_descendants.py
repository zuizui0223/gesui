#!/usr/bin/env python3
"""Correct taxon-level iNaturalist coverage: resolve taxon ID before counting.

Explicitly queries descendants with taxon_id instead of taxon_name (which
may match only the family-level observations on the API). Count source-only
photo opportunities for 2018/2022 European Lycaenidae vs controls.
Counts are NOT behavior, colonization, population abundance or WWTP effects.
"""
from __future__ import annotations
import argparse,datetime,hashlib,json,time,urllib.request,urllib.parse,urllib.error
from pathlib import Path

ROOT="https://api.inaturalist.org/v1"
FIELD="Lepidopterans aggregating/Roosting/mud puddling"
BBOX={"swlat":"44","swlng":"-12","nelat":"64","nelng":"26"}
FAMILIES=("Lycaenidae","Pieridae","Papilionidae")
YEARS=(2018,2022)
N=0

def request(path,p=None):
    global N
    N+=1
    if N>60:raise ValueError("Query cap exceeded")
    url=ROOT+path+("?" + urllib.parse.urlencode(p) if p else "")
    receipt={"url":url}
    try:
        with urllib.request.urlopen(urllib.request.Request(url,headers={
             "User-Agent":"gesui-lycaenid-species-descendant-source-feasibility",
             "Accept":"application/json"}),timeout=34) as r:
            raw=r.read(500000)
            j=json.loads(raw)
            receipt.update({"http":r.status,"source_sha256":hashlib.sha256(raw).hexdigest()})
            return receipt,j
    except urllib.error.HTTPError as e:
        receipt.update({"http":e.code,"error":str(e)})
    except Exception as e:
        receipt.update({"error":type(e).__name__+": "+str(e)})
    return receipt,None

def taxon_id_for_family(name):
    r,j=request("/taxa/autocomplete",{"q":name,"per_page":20})
    entries=j.get("results",[]) if isinstance(j,dict) else []
    eligible=[x for x in entries if x.get("name")==name and x.get("rank")=="family"]
    if len(eligible)!=1 or not isinstance(eligible[0].get("id"),int):
        raise ValueError(f"Could not uniquely resolve family taxon {name}: {r}, {[{k:v for k,v in x.items() if k in ('name','id','rank')} for x in entries[:8]]}")
    x=eligible[0]
    return x["id"],{"taxonomy_api":r,"name":name,"id":x["id"],
                  "rank":x["rank"],"observations_count_indexed":x.get("observations_count")}

def obs_query(taxonid,year,ex):
    params={**BBOX,"taxon_id":taxonid,
            "d1":f"{year}-01-01","d2":f"{year}-12-31",
            "photos":"true","geo":"true","geoprivacy":"open",
            "taxon_geoprivacy":"open","per_page":"2"}
    params.update(ex)
    req,j=request("/observations",params)
    out={"request":req,"total_results":j.get("total_results") if isinstance(j,dict) else None,
         "sample_taxon_rank":[],
         "sample_field_validation":[]}
    if isinstance(j,dict):
        for x in j.get("results",[])[:2]:
            out["sample_taxon_rank"].append({
               "id":x.get("id"),"name":(x.get("taxon") or {}).get("name"),
               "rank":(x.get("taxon") or {}).get("rank"),"positional_accuracy_m":x.get("positional_accuracy"),
               "coordinates_obscured":x.get("coordinates_obscured"),"date":x.get("observed_on")})
            out["sample_field_validation"].append([{
               "id":f.get("observation_field_id"),
               "value":f.get("value")} for f in (x.get("ofvs") or x.get("observation_field_values") or [])][:12])
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--out",type=Path,required=True);args=ap.parse_args()
    taxonomy={f:taxon_id_for_family(f) for f in FAMILIES}
    queries={}
    for family,(id,meta) in taxonomy.items():
        queries[family]={}
        for year in YEARS:
            q={
              "descendant_species_photo_pool":obs_query(id,year,{}),
              "field_name_yes_candidate":obs_query(id,year,{"field:"+FIELD:"Yes"}),
              "described_as_puddling_candidate":obs_query(id,year,{"q":"puddling"}),
            }
            queries[family][str(year)]=q
            print("GESUI_INAT_BUTTERFLY_DESCENDANTS_SOURCE_ONLY",json.dumps({
                "taxon":family,"taxon_id":id,"year":year,
                "counts":{k:v["total_results"] for k,v in q.items()},
                "first_ranks":{k:v["sample_taxon_rank"][:1] for k,v in q.items()},
                "http":{k:v["request"].get("http") for k,v in q.items()}},ensure_ascii=False),flush=True)
            time.sleep(0.4)
    obj={
        "schema":"gesui_iNat_source_only_2018_2022_europe_taxon_id_descendant_puddling_gate_v02",
        "recorded_at_UTC":datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "taxonomy":{k:v[1] for k,v in taxonomy.items()},
        "queries":queries,"bbox":BBOX,
        "previous_family_name_only_count_is_invalid_as_species_coverage":True,
        "no_puddling_behavior_photo_confirmed":True,
        "no_wwtp_sewage_exposure_assigned":True,
        "no_butterfly_population_changes_tested":True,
        "source_limitations":[
            "Taxon ID selects descendants but does not validate visual feeding action",
            "Field 17428 combines aggregation, roosting and mud-puddling",
            "Field search may be incomplete / differently indexed",
            "Observed date is not date uploaded",
            "Geoprivacy, spatial uncertainty, observer effort and urbanization remain confounders",
            "0 field hits cannot be equated with 0 true puddling events",
            "Public photo records are ecological opportunity, not local occurrence zeros"
        ]
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(obj,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_INAT_PUBLIC_FAMILY_DESCENDANT_COVERAGE_GATE_SUCCESS",
          json.dumps({"requests":N,"validation_not_ecology":True}),flush=True)
if __name__=="__main__":main()
