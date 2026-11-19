#!/usr/bin/env python3
"""Fetch original French ROSEAU/SILLAGE public sludge export for source QA.

Endpoint and form keys are exactly from the government's published
/public/js/exportbasededonneeBoue_min.js (SHA pinned).
No biological traits. No form guesswork. No raw government rows committed.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,re,requests
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urljoin

BASES=(
    "https://www.assainissement.developpement-durable.gouv.fr",
    "https://assainissement.developpement-durable.gouv.fr",
    "https://www2.assainissement.developpement-durable.gouv.fr",
)
OFFICIAL_PAGE_PATH="/pages/data/basededonneeboue.php"
OFFICIAL_JS_PATH="/public/js/exportbasededonneeBoue_min.js"
OFFICIAL_EXPORT_PATH="/controllers/GetDataExportBoueHandler.php"
JS_SHA="f4cf1e2fe84e0ccf608ece6c5d363426960e2b730e21900d48265eba4a78f0e5"
FORM_KEYS=("bassin","region","annee","format","tranche","departement")
PERMITTED_YEARS={str(y) for y in range(2017,2023)}
MAX_RAW=95_000_000
SOURCE_FIELD=re.compile(r"(epand|épand|commune|code.?insee|boues|mati[eè]re.?s[eè]che|tonne|anne|year|destin|sillage|parcelle)",re.I)

def digest(raw):return hashlib.sha256(raw).hexdigest()

def sniff(raw):
    text=None
    enc=None
    for encoding in ("utf-8-sig","utf-8","cp1252","iso-8859-1"):
        try:
            text=raw.decode(encoding)
            enc=encoding
            break
        except UnicodeDecodeError:continue
    if text is None:raise ValueError("Unable to decode original official CSV")
    first=text[:4000].lower()
    if "<html" in first or "aucune donn" in first or "erreur" in first and len(text)<4000:
        raise ValueError("Official export is HTML or no-record notice, not data")
    lines=text.splitlines()
    if not lines:raise ValueError("No rows from official data")
    delimiters={"|":lines[0].count("|"),";":lines[0].count(";"),",":lines[0].count(",")}
    sep=max(delimiters,key=delimiters.get)
    if delimiters[sep]<=0:
        raise ValueError("Official export contains no recognized CSV delimiter")
    parsed=csv.reader(io.StringIO(text),delimiter=sep)
    header=next(parsed)
    if len(header)<5:
        raise ValueError(f"Too few source columns: {len(header)}")
    n=0
    nonempty=collections.Counter()
    unique=collections.defaultdict(set)
    for row in parsed:
        if not any(v.strip() for v in row):continue
        n+=1
        if len(row)!=len(header):
            raise ValueError(f"Official CSV row width changed: {len(row)} vs {len(header)} at record {n}")
        for i,v in enumerate(row):
            if v.strip():
                nonempty[i]+=1
                if SOURCE_FIELD.search(header[i]) and len(unique[i])<200000:
                    unique[i].add(v.strip())
    flags=[]
    for i,h in enumerate(header):
        if SOURCE_FIELD.search(h):
            flags.append({"column":h,"position":i,"nonempty_records":nonempty[i],
                          "nonempty_ratio":(nonempty[i]/n if n else None),
                          "unique_nonempty_values":len(unique[i])})
    return {
        "encoding":enc,"delimiter":sep,"columns":header,
        "record_rows":n,"source_relevant_fields":flags,
        "csv_parse_complete":True,
    }

def retrieve(year,client,base):
    page_url=base+OFFICIAL_PAGE_PATH
    js_url=base+OFFICIAL_JS_PATH
    export_url=base+OFFICIAL_EXPORT_PATH
    if year not in PERMITTED_YEARS:
        raise ValueError("Year outside frozen 2017–22 source cohort")
    page=client.get(page_url,timeout=(12,35))
    page.raise_for_status()
    if year not in page.text or "list-bassin" not in page.text:
        raise ValueError("Official website no longer advertises year/required basin filter")
    js=client.get(js_url,timeout=(12,35))
    js.raise_for_status()
    if digest(js.content)!=JS_SHA:
        raise ValueError("Official export JavaScript changed; endpoint/controls not re-authorized")
    for key in FORM_KEYS:
        if f'formData.append("{key}"' not in js.text:
            raise ValueError("Original government form key not present "+key)
    if 'fetch("/controllers/GetDataExportBoueHandler.php"' not in js.text:
        raise ValueError("Original public JS does not authorize expected export handler")
    body={
      "bassin":"ALL",
      "region":"ALL",
      "annee":year,
      "format":"csv",
      "tranche":"ALL",
      "departement":"ALL",
    }
    res=client.post(export_url,data=body,timeout=(15,110),headers={
      "Origin":base,
      "Referer":page_url,
      "Accept":"text/csv,text/plain,*/*"
    })
    res.raise_for_status()
    raw=res.content
    if len(raw)>MAX_RAW:
        raise ValueError("Official CSV exceeds fixed source acquisition byte limit")
    if len(raw)<3500:
        raise ValueError("Original public JS documents exports shorter than 3485 bytes as no records")
    parsed=sniff(raw)
    return {
      "year":year,"official_source_url":page_url,"official_export_endpoint":export_url,
      "official_JS_url":js_url,"original_published_JS_sha256":JS_SHA,
      "request_form_public_fields":body,
      "official_HTTP_status_code":res.status_code,
      "original_content_type":res.headers.get("Content-Type",""),
      "original_file_bytes":len(raw),"original_file_SHA256":digest(raw),
      "original_website_page_SHA256":digest(page.content),
      "original_csv_metadata":parsed,
      "government_raw_personal_or_stew_rows_NOT_redistributed":True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--years",type=str,default="2018")
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    years=[x.strip() for x in a.years.split(",") if x.strip()]
    if not years or len(set(years))!=len(years) or any(x not in PERMITTED_YEARS for x in years):
        raise ValueError("Select unique years strictly within 2017–2022")
    a.out.mkdir(parents=True,exist_ok=True)
    c=requests.Session()
    c.headers.update({"User-Agent":"gesui-academic-FrROSEAU-sludge-source-first/20261119"})
    results=[]
    for year in years:
        failures=[]
        record=None
        for base in BASES:
            try:
                record=retrieve(year,c,base)
                break
            except Exception as e:
                failures.append({"official_host":base,"reason":str(e)[:350]})
                print("GESUI_FR_SLUDGE_OFFICIAL_HOST_STOP",year,base,type(e).__name__,str(e)[:180],flush=True)
        if record is not None:
            record["official_host_attempt_failures_before_success"]=failures
            results.append(record)
            print("GESUI_FR_SLUDGE_ORIGINAL_YEAR_SOURCE",year,"records",record["original_csv_metadata"]["record_rows"],
                  "columns",len(record["original_csv_metadata"]["columns"]),
                  "application_commune_columns",[f["column"] for f in record["original_csv_metadata"]["source_relevant_fields"]
                                                  if "commune" in f["column"].lower() or "pand" in f["column"].lower()],
                  flush=True)
        else:
            results.append({"year":year,"status":"SOURCE_ACCESS_OR_FORMAT_STOP","official_host_failures":failures})
            print("GESUI_FR_SLUDGE_ORIGINAL_YEAR_STOP",year,"all three public official hosts failed",flush=True)
    report={
      "schema":"gesui_20261119_FR_ROSEAU_SILLAGE_original_public_CSV_year_source_only_v1",
      "UTC":datetime.now(timezone.utc).isoformat(),
      "scope":"Public government French sludge / communes d'epandage; not individual field application",
      "requested_years":years,"source_year_results":results,
      "all_requested_years_downloaded_and_parsed":all("original_csv_metadata" in r for r in results),
      "independent_communes_of_actual_application_NOT_ESTABLISHED":True,
      "does_not_infer_2017to22_biological_exposure":True,
      "no_EQR_plant_NUTRI_butterfly_HNS_values_read":True,
      "old_NUTS2_and_226H7_POINT_DIFFUSE_SOURCE_STOP_remains":True,
    }
    (a.out/"FR_ROSEAU_BOUES_ORIGINAL_YEAR_SOURCE_CSV_QA.json").write_text(
      json.dumps(report,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("GESUI_FR_SLUDGE_CSV_SOURCE_STATUS",json.dumps({
      "year_count":len(years),
      "accessible_years":[r["year"] for r in results if "original_csv_metadata" in r],
      "stopped_years":[r["year"] for r in results if "original_csv_metadata" not in r],
      "biological_values_read":False
    }),flush=True)

if __name__=="__main__":main()
