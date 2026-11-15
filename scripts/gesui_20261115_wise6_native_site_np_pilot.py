#!/usr/bin/env python3
"""Original-source WISE6 site-ID chemistry cohort availability pilot.

The pilot selection is deterministic and based on original station identifiers
and source REPORTING status only, never on water chemistry or biology values.
No WWTP-specific effect or causal mediation is estimated.
"""
import collections
import concurrent.futures
import csv
import hashlib
import io
import json
import math
import os
import pathlib
import time
import urllib.parse
import urllib.request
import urllib.error

TABLE="[WISE_SOE].[latest].[Waterbase_T_WISE6_AggregatedData]"
NP={"TN":"EEA_31615-01-7","TP":"CAS_7723-14-0"}
FIRST_454="EU234_H7_ORIGINAL_STUDY_NATIVE_BASINPERIOD_EQR_PLANT_BUTTERFLY_SOURCE_COHORT.csv"
EEA_MAP="EEA2018_EEA2022_EPA1819_NATIVE_HYBAS7_SOURCE_ONLY_GEOGRAPHY.json"
NP_SOURCE="EEA2018_EEA2022_454_H7_ORIGINAL_REPORTED_NP_SITE_FLOW_EQUIVALENT_SOURCE_ONLY.json"

def digest(raw):
    return hashlib.sha256(raw).hexdigest()
def selection():
    with (pathlib.Path(os.environ["ORIGINAL_H7_454_DIR"])/FIRST_454).open(encoding="utf-8-sig",newline="") as f:
        original=list(csv.DictReader(f))
    keymap={(r["HYBAS_ID_L07"],r["PERIOD"]):r["country"] for r in original}
    if len(original)!=454 or len(keymap)!=454 or len({k[0] for k in keymap})!=232:
        raise RuntimeError("Original 454/232 cohort drift")
    sdata=json.loads((pathlib.Path(os.environ["ORIGINAL_EEA_SITE_H7_MAP_DIR"])/EEA_MAP).read_bytes())
    sites=collections.defaultdict(list)
    for r in sdata["candidate_source_site_period_rows"]:
        k=(str(r.get("HYBAS_ID_L07")),r.get("period"))
        if k in keymap and r.get("panel") in ("EEA18","EEA22"):
            bits=r["physical_id"].split("|")
            if len(bits)!=3 or bits[0]!=r["country"] or r["physical_id"]!=r["native_id"]:
                raise RuntimeError("Original native ID mismatch")
            sites[k].append(r["physical_id"])
    if sum(map(len,sites.values()))!=3525 or len({x for arr in sites.values() for x in arr})!=2485:
        raise RuntimeError("Original 3525/2485 site source drift")
    data=json.loads((pathlib.Path(os.environ["ORIGINAL_REPORTED_NP_DIR"])/NP_SOURCE).read_bytes())
    sc=data["source_H7_period_reported_N_P_upstream_to_true_original_station"]
    reported={(r["HYBAS_ID_L07"],r["PERIOD"]):r["source_station_annual_upstream_NP_positive_n"]>0 for r in sc}
    if len(reported)!=454 or sum(reported.values())!=56 or set(reported)!=set(keymap):
        raise RuntimeError("Original positive reported source support shifted")
    rank=lambda s:digest(("GESUI_WISE6_NATIVE_PILOT_20261115|"+str(s)).encode())
    chosen=[]
    for flag,goal,per_key in ((True,20,2),(False,8,1)):
        used=set()
        ordered=sorted((k for k in reported if reported[k] is flag),key=rank)
        for k in ordered:
            if k[0] in used:continue
            known=sorted(set(sites[k]),key=rank)
            if not known:continue
            used.add(k[0])
            chosen.append({"H7":k[0],"period":k[1],"calendar_year":2018 if k[1]=="2017-2019" else 2022,
                           "country":keymap[k],"positive_upstream_annual_report_H7":flag,
                           "native_sites":known[:per_key],"n_original_EEA_sites":len(known)})
            if sum(q["positive_upstream_annual_report_H7"] is flag for q in chosen)>=goal:break
    if sum(z["positive_upstream_annual_report_H7"] for z in chosen)!=20 or len(chosen)!=28:
        raise RuntimeError("Source-only pilot selection count changed")
    return chosen

def read_site_source(site):
    bits=site.split("|")
    if len(bits)!=3 or "'" in bits[1] or "'" in bits[2]:
        raise ValueError("Unsafe native station ID")
    sql="SELECT * FROM "+TABLE+" WHERE monitoringSiteIdentifier = '"+bits[1]+"'"
    rows=[];requests=[]
    for page in range(1,5):
        url="https://discodata.eea.europa.eu/sql?"+urllib.parse.urlencode({"query":sql,"p":page,"nrOfHits":250})
        try:
            req=urllib.request.Request(url,headers={"Accept":"application/json","User-Agent":"gesui-EEA-native-site-NP-quality-pilot/20261115"})
            with urllib.request.urlopen(req,timeout=50) as resp: raw=resp.read(1500000)
            obj=json.loads(raw); r=obj.get("results") if isinstance(obj,dict) else None
            if not isinstance(r,list):return {"status":"bad_json_results","receipts":requests,"source_rows":[]}
            requests.append({"page":page,"n":len(r),"source_sha256":digest(raw),"url":url})
            rows.extend(r)
            if len(r)<250:
                return {"status":"complete","receipts":requests,"source_rows":rows}
        except urllib.error.HTTPError as e:
            return {"status":"HTTP_"+str(e.code),"receipts":requests,"source_rows":[]}
        except Exception as e:
            return {"status":"error_"+type(e).__name__,"receipts":requests,"source_rows":[]}
        time.sleep(.25)
    return {"status":"truncated_1000","receipts":requests,"source_rows":rows}

def check_one(source,site,year,code):
    if source["status"]!="complete":
        return {"source_status":source["status"],"native_site":site,"calendar_year":year,
                "nutrient":code,"source_rows_seen":len(source["source_rows"]),"verified":False}
    bits=site.split("|")
    candidates=[r for r in source["source_rows"] if
       str(r.get("monitoringSiteIdentifier"))==bits[1]
       and str(r.get("monitoringSiteIdentifierScheme"))==bits[2]
       and int(r.get("phenomenonTimeReferenceYear") or -1)==year
       and r.get("observedPropertyDeterminandCode")==NP[code]
       and r.get("parameterWaterBodyCategory")=="RW"
       and r.get("procedureAnalysedMatrix")=="W"]
    strict=[]
    for r in candidates:
        try:
            x=float(r["resultMeanValue"])
            n=int(r["resultNumberOfSamples"])
            u=r["resultUom"]
            v=r.get("resultQualityMeanBelowLOQ")
            if math.isfinite(x) and x>=0 and n>=3 and u==("mg{N}/L" if code=="TN" else "mg{P}/L") and str(v).lower() not in ("true","1"):
                strict.append((r,x,n))
        except (ValueError,TypeError,KeyError):
            continue
    return {"source_status":"complete","native_site":site,"calendar_year":year,
            "nutrient":code,"native_candidate_records":len(candidates),
            "strict_qualified_records":len(strict),
            "verified":len(strict)==1,
            "measured_mgL":strict[0][1] if len(strict)==1 else None,
            "counted_samples":strict[0][2] if len(strict)==1 else None,
            "LOQ_assumed_absence":False}

def main():
    cohort=selection()
    bysite=sorted({s for r in cohort for s in r["native_sites"]})
    print("SOURCE_SELECTED",json.dumps({"H7_period":len(cohort),"native_sites":len(bysite),
                     "positive_H7_period":20,"no_matched_positive_report_H7_period":8}),flush=True)
    results={}
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        futures={pool.submit(read_site_source,s):s for s in bysite}
        for future in concurrent.futures.as_completed(futures):
            site=futures[future]
            results[site]=future.result()
            print("STATION_SOURCE",json.dumps({"site":site,"status":results[site]["status"],
                 "rows":len(results[site]["source_rows"])},sort_keys=True),flush=True)
    summary=[];hits=[]
    for q in cohort:
        site_results=[]
        for site in q["native_sites"]:
            for nutrient in ("TN","TP"):
                v=check_one(results[site],site,q["calendar_year"],nutrient)
                site_results.append(v)
                if v["verified"]:hits.append(dict(q,**v))
        summary.append({k:v for k,v in q.items() if k!="native_sites"}|
                       {"native_sites":q["native_sites"],
                        "site_results":site_results,
                        "same_site_verified_both_TN_TP":any(
                           any(z["verified"] and z["nutrient"]=="TN" and z["native_site"]==s for z in site_results)
                           and any(z["verified"] and z["nutrient"]=="TP" and z["native_site"]==s for z in site_results)
                           for s in q["native_sites"])})
    complete=[x for x in cohort if all(results[s]["status"]=="complete" for s in x["native_sites"])]
    same=sum(x["same_site_verified_both_TN_TP"] for x in summary)
    out=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_wise6_native_pilot"
    out.mkdir(parents=True,exist_ok=True)
    (out/"SELECTED_SOURCE_OPPORTUNITIES.json").write_text(json.dumps(cohort,indent=2))
    (out/"SOURCE_SITE_QUALITY_RESULTS.json").write_text(json.dumps(summary,indent=2))
    # The stored data are only actual exact native water quality source observations.
    (out/"WQ_VERIFIED_SITE_YEAR_NP.json").write_text(json.dumps(hits,indent=2))
    report={"schema":"GESUI_WISE6_2026_FROZEN_NATIVE_SITE_YEAR_PILOT_V1",
            "cohort_H7_period":len(cohort),"native_sites":len(bysite),
            "H7_period_fully_source_retrieved":len(complete),
            "H7_period_same_native_site_TN_AND_TP_strict":same,
            "source_positive_H7_period_both":sum(x["same_site_verified_both_TN_TP"] and x["positive_upstream_annual_report_H7"] for x in summary),
            "source_nonpositive_H7_period_both":sum(x["same_site_verified_both_TN_TP"] and not x["positive_upstream_annual_report_H7"] for x in summary),
            "source_access_status":dict(collections.Counter(x["status"] for x in results.values())),
            "incomplete_source_never_zero_imputed":True,
            "no_eco_outcomes_opened_or_wwtp_causal_model_fitted":True,
            "source_selection_prechemistry":True}
    (out/"SOURCE_COVERAGE_GATE.json").write_text(json.dumps(report,indent=2,sort_keys=True))
    print("WISE6_PILOT_FINAL",json.dumps(report,sort_keys=True),flush=True)
if __name__=="__main__":main()
