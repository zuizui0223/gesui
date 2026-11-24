#!/usr/bin/env python3
"""Source-only original H7×period actual GBIF speciesKey identity retention.

Fixes the original source gap: the 2026-11-12 plant S50 / five-butterfly-family
artifact archived only species-count summaries, NOT per-species identities.
This script uses the SAME previously signed 234 H7 geometry and five families,
same two time windows and API source filters; it saves exact source identities,
including zero/missing and response-source SHA. It does NOT join hosts or claim
field feeding, sewage exposure, pollination, or species absences.

Run in trusted networked environment with original JSON source:
  H7_234_PRE_NEW_RESPONSE_SOURCE=/.../HYBAS7_234_REPEATED_EU_EQR_ORIGINAL_WKT_TAXON_FAMILY_PRE_NEW_OUTCOME.json \
    python scripts/gesui_H7_original_fivefamily_plants_speciesKey_retention.py \
    --mode pilot --shard 0 --out /tmp/gesui_facets
Four independent original geometry shards available with --mode full.
"""
import argparse
import collections
import concurrent.futures
import hashlib
import json
import os
import pathlib
import time
import urllib.error
import urllib.parse
import urllib.request

SOURCE_SHA="afc524ebfa353e4d1c6ea51d208ba529a85fedb57c4d0f28d83545010bb76189"
GBIF_BACKBONE="d7dddbf4-2cf0-4f39-9b2a-bb099caae36c"
FAMILIES=("Lycaenidae","Nymphalidae","Pieridae","Papilionidae","Hesperiidae")
PERIODS={"2017-2019":"2017,2019","2020-2022":"2020,2022"}

def sha(b):
    return hashlib.sha256(b).hexdigest()

def read_frozen():
    raw=pathlib.Path(os.environ["H7_234_PRE_NEW_RESPONSE_SOURCE"]).read_bytes()
    if sha(raw)!=SOURCE_SHA:
        raise RuntimeError("Original signed 234 H7 geometry or taxonomy SHA drift")
    x=json.loads(raw)
    basins=x["original_two_year_EEA_H7_subbasins"]
    families=x["verified_butterfly_family_taxa"]
    if len(basins)!=234 or set(families)!=set(FAMILIES):
        raise RuntimeError("Original source family or basin population drift")
    if x["fixed_periods"]!=list(PERIODS):
        raise RuntimeError("Changed fixed source two-window chronology")
    if families["Lycaenidae"]["taxonKey"]!=5473:
        raise RuntimeError("Source Lycaenidae GBIF key changed")
    for fam in FAMILIES:
        if families[fam]["rank"]!="FAMILY" or families[fam]["scientificName"]!=fam:
            raise RuntimeError("Family source taxonomy unknown")
    return basins,families

def http_json(params):
    url="https://api.gbif.org/v1/occurrence/search?"+urllib.parse.urlencode(params,safe="(),")
    req=urllib.request.Request(url,headers={"Accept":"application/json",
        "User-Agent":"gesui-H7-frozen-original-specieskey-source-20261114/1"})
    errors=[]
    for k in range(5):
        try:
            with urllib.request.urlopen(req,timeout=80) as response:
                raw=response.read()
            return {"ok":True,"raw":raw,"uri":url,"sha256":sha(raw)}
        except urllib.error.HTTPError as e:
            errors.append(f"{e.code}:{str(e)[:60]}")
            if e.code not in (400,429,500,502,503,504):
                break
        except Exception as e:
            errors.append(type(e).__name__+":"+str(e)[:80])
        time.sleep([1,2,4,8,13][k])
    return {"ok":False,"uri":url,"errors":errors}

def task(basin,period,field,families):
    h=basin["HYBAS_ID_L07"]
    hits=collections.Counter()
    receipts=[]
    parts=basin["polygon_parts"]
    complete=True;truncated=False
    rawtotal=0
    for part in parts:
        wkt=part["wkt"]
        if sha(wkt.encode())!=part["wkt_SHA256"]:
            raise RuntimeError("Original polygon component SHA drift")
        params={"geometry":wkt,"year":PERIODS[period],"limit":"0",
            "hasCoordinate":"true","hasGeospatialIssue":"false",
            "occurrenceStatus":"PRESENT","facet":"speciesKey",
            "facetLimit":"1000"}
        if field=="Plantae":
            params["taxonKey"]="6"
        else:
            params["taxonKey"]=str(families[field]["taxonKey"])
            params["checklistKey"]=GBIF_BACKBONE
        offset=0
        parttotal=None
        part_keys=set()
        while True:
            query=dict(params,facetOffset=str(offset))
            r=http_json(query)
            if not r["ok"]:
                receipts.append({"status":"API_missing","url":r["uri"],
                    "errors":r["errors"],"component_SHA":part["wkt_SHA256"]})
                complete=False
                break
            j=json.loads(r["raw"])
            n=j.get("count")
            if not isinstance(n,int) or n<0:
                raise RuntimeError("Invalid original GBIF occurrence count")
            if parttotal is None:parttotal=n
            elif n!=parttotal:
                complete=False
                receipts.append({"status":"nonrepeatable_total_count_within_pagination",
                    "url":r["uri"],"observed":n,"expected":parttotal})
                break
            facets=[z for z in j.get("facets",[]) if
                    str(z.get("field","")).upper().replace("_","")=="SPECIESKEY"]
            if len(facets)!=1:
                complete=False
                receipts.append({"status":"unknown_GBIF_speciesKey_facet",
                    "url":r["uri"],"source_fields":[z.get("field") for z in j.get("facets",[])]})
                break
            chunk=facets[0].get("counts",[])
            for cell in chunk:
                key=str(cell.get("name",""))
                count=cell.get("count")
                if not key.isdigit() or not isinstance(count,int) or count<0:
                    raise RuntimeError("Source speciesKey/count malformed")
                if key in part_keys:
                    raise RuntimeError("Duplicate original speciesKey across paginated facets in one polygon part")
                part_keys.add(key)
                hits[key]+=count
            receipts.append({"status":"source_success","url":r["uri"],
                "source_SHA256":r["sha256"],"facetOffset":offset,
                "facet_length":len(chunk),"component_SHA":part["wkt_SHA256"]})
            if len(chunk)<1000:break
            offset+=1000
            if offset>=20000:
                complete=False;truncated=True
                receipts.append({"status":"species_facet_cap_reached",
                    "facetOffset":offset,"component_SHA":part["wkt_SHA256"]})
                break
        if parttotal is not None:rawtotal+=parttotal
        if not complete:break
    eligible=complete and len(hits)>0
    if complete and sum(hits.values())>rawtotal:
        raise RuntimeError("Species-facet counts exceed original GBIF occurrence total")
    return {"HYBAS_ID_L07":h,"period":period,"group":field,
        "source_status":"ok" if complete else "missing",
        "source_query_occurrence_count":rawtotal if complete else None,
        "source_facetted_occurrence_count":sum(hits.values()) if complete else None,
        "source_speciesKey_count":len(hits) if complete else None,
        "source_facet_coverage":sum(hits.values())/rawtotal if eligible and rawtotal>0 else None,
        "source_speciesKey_counts":dict(sorted(hits.items())) if complete else None,
        "source_request_receipts":receipts,
        "truncated":truncated,
        "index_nonobservation_not_biological_absence":True}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--mode",choices=("pilot","full"),required=True)
    p.add_argument("--shard",type=int,choices=range(4),default=0)
    p.add_argument("--out",required=True)
    args=p.parse_args()
    basins,families=read_frozen()
    if args.mode=="pilot":
        bs=sorted(basins,key=lambda x:sha(("GESUI_20261114_source_only_facet_pilot|"+
            x["HYBAS_ID_L07"]).encode()))[:8]
    else:
        bs=[x for x in basins if int(sha(x["HYBAS_ID_L07"].encode())[:8],16)%4==args.shard]
    if not bs:raise RuntimeError("Empty source shard")
    tasks=[(z,period,group,families) for z in bs
           for period in PERIODS for group in ("Plantae",)+FAMILIES]
    result=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        futs=[pool.submit(task,*q) for q in tasks]
        for fut in concurrent.futures.as_completed(futs):
            result.append(fut.result())
    result.sort(key=lambda x:(x["HYBAS_ID_L07"],x["period"],x["group"]))
    if len(result)!=len(tasks):raise RuntimeError("Source opportunity count drift")
    problems=[q for q in result if q["source_status"]!="ok"]
    totals=collections.Counter(q["group"] for q in result if q["source_status"]=="ok")
    report={"mode":args.mode,"shard":args.shard,
        "n_original_H7":len(bs),"n_target_source_facets":len(tasks),
        "n_success":len(result)-len(problems),
        "n_missing":len(problems),
        "source_Frozen_234_source_SHA":SOURCE_SHA,
        "source_status_by_group":dict(totals),
        "note":"Only GBIF indexed source species identities, not realized host interaction or wastewater effects",
        "source_complete":len(problems)==0}
    out=pathlib.Path(args.out);out.mkdir(parents=True,exist_ok=True)
    stem=f"H7_speciesKey_{args.mode}_shard{args.shard}"
    (out/(stem+".json")).write_text(json.dumps(result,sort_keys=True,ensure_ascii=False))
    (out/(stem+"_audit.json")).write_text(json.dumps(report,sort_keys=True,indent=2))
    print(json.dumps(report,sort_keys=True),flush=True)
    if problems:raise SystemExit("Source-only missingness must be resolved or kept NA before any ecological fit")

if __name__=="__main__":
    main()
