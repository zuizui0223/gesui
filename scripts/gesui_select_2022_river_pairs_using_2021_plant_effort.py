#!/usr/bin/env python3
"""2021 Plantae-only photo-effort screening of frozen 2022 EEA WWTP pairs.

The 2022 WWTP pair registry was committed before this script. The screen
does NOT query 2022 butterfly taxa, photo share, species, or puddling behavior.
Output is an independent botanical reporting-opportunity convenience cohort,
not evidence of sewage increasing organisms or of matched observer effort.
"""
from __future__ import annotations
import argparse,csv,datetime,hashlib,json,time,urllib.parse,urllib.request,urllib.error
from pathlib import Path
from collections import defaultdict,Counter
URL="https://api.inaturalist.org/v1"
MIN_PER_SITE=8
MIN_RATIO=0.20
MAX_SELECTED=32
MAX_PER_COUNTRY=8
RADIUS_KM=5
YEAR=2021
ACCURACY=500
STATUS=("upstream_registered_WWTP","no_registered_qualifying_upstream_WWTP")

def request(path,params):
    url=URL+path+"?"+urllib.parse.urlencode(params)
    receipt={"url":url}
    for attempt in range(4):
        try:
            req=urllib.request.Request(url,headers={"Accept":"application/json","User-Agent":"gesui-2021-public-plant-photo-effort-only/0.1"})
            with urllib.request.urlopen(req,timeout=26) as response:
                raw=response.read(450000)
                x=json.loads(raw)
                receipt.update({"http_status":response.status,"source_sha256":hashlib.sha256(raw).hexdigest()})
                return receipt,x
        except urllib.error.HTTPError as e:
            receipt.update({"http_status":e.code,"error":str(e)})
            if e.code not in (429,500,502,503) or attempt==3:break
        except Exception as e:
            receipt["error"]=type(e).__name__+": "+str(e)
            if attempt==3:break
        time.sleep(2.0*(attempt+1))
    return receipt,None

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-pairs",required=True,type=Path)
    ap.add_argument("--source-qc",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args()
    qc=json.loads(a.source_qc.read_text())
    assert qc["analysis_design"].startswith("2021_Plant_photo_effort") and qc["radius_km"]==RADIUS_KM
    if hashlib.sha256(a.source_pairs.read_bytes()).hexdigest()!=qc["paired_sample_sha256"]:
        raise ValueError("Source-only new 5-km candidate pairs modified after freeze")
    rows=list(csv.DictReader(a.source_pairs.open()))
    if len(rows)!=qc["country_3deg_matched_independent_pairs"]*2:
        raise ValueError("Source-only pair registry length altered")
    tax_req,tax=request("/taxa/autocomplete",{"q":"Plantae","per_page":15})
    found=[q for q in (tax or {}).get("results",[]) if q.get("name")=="Plantae" and q.get("rank")=="kingdom"]
    if len(found)!=1:raise ValueError("Could not resolve unique Plantae kingdom iNaturalist ID")
    taxid=found[0]["id"]
    all_records=[]
    for i,row in enumerate(rows):
        params={"taxon_id":taxid,"year":YEAR,
          "lat":str(float(row["lat"])),"lng":str(float(row["lon"])),"radius":RADIUS_KM,
          "photos":"true","geo":"true","geoprivacy":"open","taxon_geoprivacy":"open",
          "quality_grade":"research","captive":"false","acc_below":ACCURACY,
          "per_page":1}
        rec,x=request("/observations",params)
        total=x.get("total_results") if isinstance(x,dict) else None
        all_records.append({"pair_id":int(row["pair_id"]),"status":row["status"],
            "country":row["countryCode"],"source_site_key":row["site_key"],
            "plant_2021_photo_count":total,"request":rec})
        if i%20==0 or i==len(rows)-1:
            print("GESUI_INAT_INDEPENDENT_PRIORYEAR_PLANT_EFFORT_PROGRESS",
                json.dumps({"queried":i+1,"total":len(rows),"latest_n":total,
                "http":rec.get("http_status")},sort_keys=True),flush=True)
        time.sleep(0.9)
    group=defaultdict(dict)
    for r in all_records:group[r["pair_id"]][r["status"]]=r
    candidates=[]
    for pairid,rr in sorted(group.items()):
        if set(rr)!=set(STATUS):raise ValueError("Malformed source-only pair")
        e=rr[STATUS[0]];c=rr[STATUS[1]]
        x=e["plant_2021_photo_count"];y=c["plant_2021_photo_count"]
        good= isinstance(x,int) and isinstance(y,int) and min(x,y)>=MIN_PER_SITE and min(x,y)/max(x,y)>=MIN_RATIO
        candidates.append({"pair_id":pairid,"country":e["country"],
             "source_plant2021_photos":x,"control_plant2021_photos":y,
             "minimum_plant2021_photos":min(x,y) if isinstance(x,int) and isinstance(y,int) else None,
             "effort_ratio":min(x,y)/max(x,y) if isinstance(x,int) and isinstance(y,int) and max(x,y)>0 else None,
             "passes_2021_nonbutterfly_effort":bool(good)})
    kept=[]
    countrycap=Counter()
    for p in sorted((r for r in candidates if r["passes_2021_nonbutterfly_effort"]),
                    key=lambda q:(-q["minimum_plant2021_photos"],q["pair_id"])):
        if countrycap[p["country"]]>=MAX_PER_COUNTRY:continue
        kept.append(p);countrycap[p["country"]]+=1
        if len(kept)>=MAX_SELECTED:break
    selected_ids=set(r["pair_id"] for r in kept)
    selected=[x for x in rows if int(x["pair_id"]) in selected_ids]
    # Frozen sampling order by source-only ID, not by butterfly results.
    selected.sort(key=lambda x:(int(x["pair_id"]),x["status"]))
    a.out.parent.mkdir(parents=True,exist_ok=True)
    with a.out.with_suffix(".csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else [])
        writer.writeheader();writer.writerows(selected)
    enough=len(kept)>=12 and len(countrycap)>=3
    result={"schema":"gesui_sourcefrozen_2021_independent_plant_reporting_effort_screen_v01",
         "source_registry_SHA256":qc["paired_sample_sha256"],
         "source_registry_QC_SHA256":hashlib.sha256(a.source_qc.read_bytes()).hexdigest(),
         "original_candidate_pairs":len(group),"processed_site_queries":len(all_records),
         "year_of_botanical_effort":YEAR,"iNaturalist_plant_kingdom_taxon_id":taxid,
         "botanical_taxon_autocomplete":tax_req,
         "radius_km":RADIUS_KM,"filter":"photographed, research grade, open geo, positional_accuracy <= 500m",
         "minimum_site_plant_photos":MIN_PER_SITE,"min_pair_photo_ratio":MIN_RATIO,
         "max_selected_pairs":MAX_SELECTED,"max_per_country":MAX_PER_COUNTRY,
         "available_plant_effort_pairs":sum(r["passes_2021_nonbutterfly_effort"] for r in candidates),
         "selected_pairs":kept,
         "selected_pair_country_counts":dict(countrycap),
         "selected_n":len(kept),
         "passes_independent_effort_gate_12_pairs_3_countries":enough,
         "all_candidate_pairs":candidates,
         "all_plant_effort_requests":all_records,
         "selected_pair_CSV_sha256":hashlib.sha256(a.out.with_suffix(".csv").read_bytes()).hexdigest(),
         "selected_pair_CSV":a.out.with_suffix(".csv").name,
         "NO_2022_BUTTERFLY_RESULT_QUERIED":True,
         "no_lycaenid_proportion_or_puddling_inference":True,
         "can_bias_generalization":"Effort-rich selection conditions on Plantae observer activity, itself linked to access/urbanization/ecology",
         "source_zero_cannot_be_interpreted_as_no_sewage":True,
         "utc":datetime.datetime.now(datetime.timezone.utc).isoformat()}
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_2021_BOTANICAL_EFFORT_FROZEN_BEFORE_2022_BUTTERFLY",
          json.dumps({"candidate_pairs":len(group),
          "qualifying_2021_plant_pairs":result["available_plant_effort_pairs"],
          "selected":len(kept),"countries":dict(countrycap),
          "gate":enough,"no_2022_butterflies_queried":True},sort_keys=True),flush=True)
if __name__=="__main__":main()
