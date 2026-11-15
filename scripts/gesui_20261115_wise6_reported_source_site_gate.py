#!/usr/bin/env python3
"""Exact native EEA site-year WWTP-reported-positive -> WISE6 TN/TP pilot.

Pilot source selection is entirely before any river chemistry API result.
No biological EQR, plant, butterfly or SEM response is read or fitted.
"""
import base64
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
import urllib.request
from gesui_20261115_wise6_native_site_np_pilot import read_site_source,check_one

SOURCES={
 "EEA18":("0f30bfbc3dac13f6a182c8a910e9899c04778322",
          "2018_measured_upstream_P_t_per_year","2018_measured_upstream_N_t_per_year"),
 "EEA22":("325f199e985d9862ffee5b0edb205c61a7be3833",
          "upstream_2022_EEA_P_measured_t_per_yr","upstream_2022_EEA_N_measured_t_per_yr"),
}
STUDY="EU234_H7_ORIGINAL_STUDY_NATIVE_BASINPERIOD_EQR_PLANT_BUTTERFLY_SOURCE_COHORT.csv"
SITES="EEA2018_EEA2022_EPA1819_NATIVE_HYBAS7_SOURCE_ONLY_GEOGRAPHY.json"

def sha(x):return hashlib.sha256(x).hexdigest()
def git_sha(x):
    return hashlib.sha1(("blob "+str(len(x))+"\0").encode()+x).hexdigest()
def get_original_csv(blob_sha):
    uri="https://api.github.com/repos/zuizui0223/gesui/git/blobs/"+blob_sha
    req=urllib.request.Request(uri,headers={"Accept":"application/vnd.github+json",
          "User-Agent":"gesui-20261115-original-upstream-station-identified-source/1"})
    with urllib.request.urlopen(req,timeout=45) as r: j=json.load(r)
    raw=base64.b64decode(j["content"].replace("\n",""))
    if git_sha(raw)!=blob_sha:raise RuntimeError("Source EEA site CSV Git blob SHA changed")
    records=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    if len(records) not in (2634,2475):
        raise RuntimeError("Original EEA source station counts drift")
    result={r["site_key"]:r for r in records}
    if len(result)!=len(records):raise RuntimeError("Original EEA native station ID not unique")
    return result

def original_source_selection():
    with (pathlib.Path(os.environ["ORIGINAL_H7_454_DIR"])/STUDY).open(encoding="utf-8-sig",newline="") as f:
        cohort=list(csv.DictReader(f))
    keys={(r["HYBAS_ID_L07"],r["PERIOD"]) for r in cohort}
    if len(keys)!=454 or len(cohort)!=454:raise RuntimeError("Study H7 period original 454 drift")
    obj=json.loads((pathlib.Path(os.environ["ORIGINAL_EEA_SITE_H7_MAP_DIR"])/SITES).read_bytes())
    sources={k:[] for k in keys}
    for r in obj["candidate_source_site_period_rows"]:
        k=(str(r.get("HYBAS_ID_L07")),r.get("period"))
        if k not in keys or r.get("panel") not in SOURCES:continue
        native=r["native_id"]
        if native!=r["physical_id"] or len(native.split("|"))!=3:
            raise RuntimeError("Native EEA source station scheme drift")
        sources[k].append(native)
    if sum(map(len,sources.values()))!=3525:raise RuntimeError("3525 original EEA station years drift")
    originals={}
    for panel,(sha_code,_,__) in SOURCES.items():
        originals[panel]=get_original_csv(sha_code)
    valid=[]
    for k,ids in sources.items():
        panel="EEA18" if k[1]=="2017-2019" else "EEA22"
        this=originals[panel]
        p_key,n_key=SOURCES[panel][1:]
        for site in ids:
            if site not in this:raise RuntimeError("Site absent from original EEA source")
            r=this[site]
            try:
                p=float(r[p_key]);n=float(r[n_key]);flow=float(r["DIS_AV_CMS"])
                if not all(math.isfinite(z) for z in (p,n,flow)) or p<0 or n<0 or flow<=0:
                    continue
            except (ValueError,TypeError):continue
            valid.append({"H7":k[0],"period":k[1],
             "year":2018 if panel=="EEA18" else 2022,
             "native_site":site,"country":r["countryCode"],
             "upstream_report_positive":bool(p>0 or n>0),
             "P_tpy_reported_index":p,"N_tpy_reported_index":n,
             "flow_modelled_cms":flow,
             "source_identity":"original_benthic_monitoring_site_actual_reported_upstream_source"})
    allpositive=[r for r in valid if r["upstream_report_positive"]]
    if len(allpositive)!=132 or len({(r["H7"],r["period"]) for r in allpositive})!=56:
        raise RuntimeError("Original source-positive 132 sites / 56 H7 periods not reproduced")
    distinct=len({r["native_site"] for r in valid})
    if len(valid)!=3525 or distinct!=2485:raise RuntimeError("Original valid site population drift")
    rank=lambda s:sha(("GESUI_20261115_REPORTED_EXACT_STATION_QA|"+str(s)).encode())
    selected=[]
    for year in (2018,2022):
        pool=sorted((r for r in allpositive if r["year"]==year),key=lambda r:rank(r["native_site"]))
        usedH7=set()
        for r in pool:
            if r["H7"] in usedH7:continue
            selected.append(r)
            usedH7.add(r["H7"])
            if len(usedH7)==12:break
        if len(usedH7)!=12:raise RuntimeError("Insufficient original independent source-positive H7 year groups for pilot")
    # Reference class is no positive source MATCH, not true no-WWTP controls.
    for year in (2018,2022):
        unused={r["native_site"] for r in selected}
        possible=sorted((r for r in valid if r["year"]==year and not r["upstream_report_positive"] and r["native_site"] not in unused),
                        key=lambda r:rank(r["native_site"]))
        picked=set()
        for r in possible:
            if r["H7"] in picked:continue
            selected.append(r);picked.add(r["H7"])
            if len(picked)==4:break
    if len(selected)!=32 or sum(r["upstream_report_positive"] for r in selected)!=24:
        raise RuntimeError("Fixed source pilot selection drift")
    return selected,valid,originals

def main():
    selected,valid,originals=original_source_selection()
    unique=sorted({r["native_site"] for r in selected})
    print("GATE_PRE_CHEMISTRY",json.dumps({
      "original_site_year":len(valid),"original_source_positive_site_year":sum(r["upstream_report_positive"] for r in valid),
      "original_positive_H7_period":56,"pilot_selected":len(selected),
      "selected_positive":sum(r["upstream_report_positive"] for r in selected)},sort_keys=True),flush=True)
    results={}
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
        fut={ex.submit(read_site_source,site):site for site in unique}
        for f in concurrent.futures.as_completed(fut):
            name=fut[f]
            results[name]=f.result()
            print("EXACT_SITE",json.dumps({"id":name,"status":results[name]["status"],
                  "total_source_rows":len(results[name].get("source_rows",[]))}),flush=True)
    output=[]
    for r in selected:
        both=[check_one(results[r["native_site"]],r["native_site"],r["year"],c) for c in ("TN","TP")]
        output.append({**r,"WISE6_TN":both[0],"WISE6_TP":both[1],
                       "both_valid":all(v["verified"] for v in both),
                       "water_quality_source_complete":results[r["native_site"]]["status"]=="complete"})
    source_audit={
        "schema":"GESUI_20261115_WISE6_2026_SOURCE_POSITIVE_NATIVE_SITE_PILOT",
        "original_site_year":len(valid),
        "original_positive_upstream_report_site_year":132,
        "original_positive_upstream_report_H7_period":56,
        "selected_site_years":len(selected),
        "selected_positive_site_years":24,
        "source_status":dict(collections.Counter(z["status"] for z in results.values())),
        "TN_strict_site_years":sum(v["WISE6_TN"]["verified"] for v in output),
        "TP_strict_site_years":sum(v["WISE6_TP"]["verified"] for v in output),
        "both_TNTP_strict_site_years":sum(v["both_valid"] for v in output),
        "source_positive_both_TNTP_site_years":sum(v["both_valid"] and v["upstream_report_positive"] for v in output),
        "has_independent_counterfactual_control":False,
        "original_2018_2022_EEA_blob_SHA":{k:v[0] for k,v in SOURCES.items()},
        "not_a_sewage_causal_model":True,"missing_never_zero_imputed":True}
    out=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_reported_site_WISE6"
    out.mkdir(parents=True,exist_ok=True)
    (out/"SOURCE_POSITIVE_NATIVE_SITE_SELECTION.json").write_text(json.dumps(selected,indent=2))
    (out/"STRICT_SAME_SITE_YEAR_TNTP_OUTPUT.json").write_text(json.dumps(output,indent=2))
    (out/"PILOT_COVERAGE_AUDIT.json").write_text(json.dumps(source_audit,indent=2,sort_keys=True))
    (out/"SOURCE_FETCH_RECEIPTS.json").write_text(json.dumps({s:results[s]["receipts"] for s in results},indent=2))
    print("SOURCE_CHEMISTRY_RESULT",json.dumps(source_audit,sort_keys=True),flush=True)
if __name__=="__main__":main()
