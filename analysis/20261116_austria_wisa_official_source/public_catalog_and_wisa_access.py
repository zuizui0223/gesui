#!/usr/bin/env python3
"""Official public Austria WISA/GZÜV source availability; no biology values.

Uses documented data.europa.eu public metadata-search API and PUBLIC WISA/H2O
landing pages. Never invokes private endpoints, undocumented business logic,
authenticated sessions, or reads benthic EQR / river TP / terrestrial NDVI.
"""
from __future__ import annotations
import argparse,hashlib,json,re,time
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urlparse,urljoin,quote_plus
import requests

SCHEMA="GESUI_20261116_WISA_GZUEV_NATIONAL_2022_PUBLIC_SOURCE_ACCESS_GATE_V1"
SEARCH_URL="https://data.europa.eu/api/hub/search/search"
SEARCH_TERMS=(
   "GZÜV Messnetz Fließgewässer",
   "GZÜV-Messnetz",
   "Gewässerzustandsüberwachung",
   "GZÜV Measuring Network Flowing Waters",
   "GZUEV Fliessgewaesser Messnetz",
)
PUBLIC_HOME=(
  "https://wasser.umweltbundesamt.at/wisa-datenabfrage/",
  "https://wasser.umweltbundesamt.at/h2odb/",
  "https://www.bmluk.gv.at/themen/wasser/wisa/daten/daten.html",
)
EXPECTED_OFFICIAL_DOMAINS=("data.gv.at","bmluk.gv.at",
  "umweltbundesamt.at","data.europa.eu")
TARGET_TOKENS=("gzuev","gzuv","flie","fliess","fließ",
               "measur","messnetz","messstelle","sampling")

class Links(HTMLParser):
    def __init__(self):super().__init__();self.links=[];self.scripts=[]
    def handle_starttag(self,tag,attrs):
        o=dict(attrs)
        if tag=="a" and o.get("href"):self.links.append({"href":o["href"],"title":o.get("title","")})
        if tag=="script" and o.get("src"):self.scripts.append(o["src"])

def digest(b):return hashlib.sha256(b).hexdigest()
def compact_dataset(rec):
    if not isinstance(rec,dict):return {"type":str(type(rec))}
    fields=("id","title","name","description","publisher","source","catalog","url",
            "uri","resource","resources","distributions","distribution","accessURL",
            "downloadURL","format","license","issued","modified")
    result={}
    for key in fields:
        if key in rec:
            a=rec[key]
            if isinstance(a,str):result[key]=a[:1400]
            elif isinstance(a,list):
                result[key]=[str(q)[:900] for q in a[:8]]
            elif isinstance(a,dict):
                result[key]={k:str(v)[:360] for k,v in list(a.items())[:9]}
            else:result[key]=a
    return result

def official_source(s):
    try:host=urlparse(s).hostname or ""
    except Exception:return False
    return any(host==v or host.endswith("."+v) for v in EXPECTED_OFFICIAL_DOMAINS)

def get(session,url,**kwargs):
    try:
        resp=session.get(url,timeout=(12,55),allow_redirects=True,**kwargs)
        return {"status":resp.status_code,"final_url":resp.url,
            "content_type":resp.headers.get("content-type",""),
            "source_response_sha256":digest(resp.content),
            "source_response_bytes":len(resp.content)},resp
    except requests.RequestException as e:
        return {"error_type":type(e).__name__,"error":str(e)[:350]},None

def search_api(session):
    records=[]
    for q in SEARCH_TERMS:
        url=SEARCH_URL+"?q="+quote_plus(q)+"&filters=dataset&limit=25"
        meta,resp=get(session,url,headers={"Accept":"application/json"})
        item={"query":q,"original_documented_search_url":url,"source_http":meta}
        if resp is not None and resp.status_code==200:
            try:
                js=resp.json()
                item["original_JSON_root_keys"]=list(js)[:25] if isinstance(js,dict) else []
                if isinstance(js,dict):
                    # Public SearchAPI shape may be {"result":...} or
                    # {"results":[...]}; never infer source record counts from
                    # an unparsed HTML response or undocumented response.
                    candidates=js.get("result",js)
                    if not isinstance(candidates,dict):candidates={}
                    hits=candidates.get("results",candidates.get("datasets",[]))
                    if not isinstance(hits,list):hits=[]
                    item["official_metadata_source_hit_count_reported"]=candidates.get(
                        "count",candidates.get("total",candidates.get("totalElements")))
                    item["actual_returned_public_dataset_items"]=len(hits)
                    item["source_public_dataset_extracts"]=[compact_dataset(x) for x in hits[:12]]
                    # Named canonical station source candidate must actually
                    # occur in official dataset title/description, never
                    # invent URL from a third-party search result.
                    shortlist=[]
                    for hit in hits:
                        title=str(hit.get("title","")) if isinstance(hit,dict) else ""
                        desc=str(hit.get("description","")) if isinstance(hit,dict) else ""
                        text=(title+" "+desc).casefold()
                        if (("gz" in text and ("fluss" in text or "flie" in text or "fließ" in text or "messnetz" in text)) or
                             "measuring network flowing waters" in text):
                            shortlist.append(compact_dataset(hit))
                    item["likely_original_GZUEV_flowing_water_public_items"]=shortlist[:20]
                else:item["api_response_not_json_object"]=True
            except (ValueError,TypeError) as e:
                item["source_JSON_parse_error"]=str(e)[:200]
        records.append(item)
        time.sleep(.15)
    return records

def websites(session):
    ans=[]
    for url in PUBLIC_HOME:
        meta,resp=get(session,url,headers={"Accept":"text/html"})
        it={"public_entrypoint":url,"source_http":meta}
        if resp is not None and resp.status_code==200:
            html=resp.text
            p=Links()
            try:p.feed(html)
            except Exception:pass
            it["requires_JS_message"]=bool("You need to enable JavaScript" in html or "Bitte aktivieren Sie JavaScript" in html)
            it["public_links"]=[{"url":urljoin(resp.url,l["href"]),"title":l["title"]}
                 for l in p.links if
                 ("qualität" in (l["href"]+l["title"]).lower() or
                 "wisa" in (l["href"]+l["title"]).lower() or
                 "gzuev" in (l["href"]+l["title"]).lower() or
                 "abfrage" in (l["href"]+l["title"]).lower())][:30]
            it["public_static_script_names"]=[urljoin(resp.url,q) for q in p.scripts[:8]]
            it["actual_2022_station_year_rows_downloaded"]=0
            it["not_automatically_a_public_tabular_export"]=True
        ans.append(it)
    return ans

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    s=requests.Session()
    s.headers.update({"User-Agent":"gesui-20261116-public-wisa-source-only-academic/1.0"})
    searched=search_api(s)
    sites=websites(s)
    candidate=sum(len(q.get("likely_original_GZUEV_flowing_water_public_items",[])) for q in searched)
    result={"schema":SCHEMA,
      "official_public_metadata_search_documentation":"https://data.europa.eu/api/hub/search/",
      "official_source_metadata_searches":searched,
      "official_public_WISA_H2O_landing_results":sites,
      "source_candidate_exact_relevant_public_metadata_items":candidate,
      "source_target_year":2022,
      "target_biotic_unit":"genuine same physically mapped river station-year with valid TN+TP and macrozoobenthos",
      "original_AT145_EEA2022_outfall_source_from_PR116_not_reselected":True,
      "public_official_GZUEV_station_geodata_SHA256_not_yet_checked":True,
      "public_source_2022_TNTP_stations_verified":False,
      "public_source_2022_MZB_station_years_verified":False,
      "source_official_EmReg_station_level_aggregate_not_single_outfall_delivered_dose":True,
      "no_2022_water_quality_benthos_ecological_numbers_opened":True,
      "no_wastewater_terrestrial_causal_fitting":True}
    if not candidate:
        result["status"]="PUBLIC_GZUEV_DATASET_DISTRIBUTION_UNRESOLVED"
    else:
        result["status"]="PUBLIC_METADATA_FOUND_YEAR2022_SAMPLE_SOURCE_NOT_YET_VERIFIED"
    (a.output/"WISA_GZUEV_OFFICIAL_SOURCE_ACCESS_PROBE_20261116.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_WISA_PUBLIC_SOURCE_ACCESS_GATE="+json.dumps({
      "status":result["status"],
      "candidate_public_original_official_GZUEV_catalog_items":candidate,
      "metadata_search_response_http":[
         {"query":q["query"],"status":q["source_http"].get("status"),
          "candidate":len(q.get("likely_original_GZUEV_flowing_water_public_items",[])),
          "returned":q.get("actual_returned_public_dataset_items")}
         for q in searched],
      "WISA_public_home":[{"url":x["public_entrypoint"],
         "http":x["source_http"].get("status"),"requires_js":x.get("requires_JS_message")}
         for x in sites],
      "original_2022_MZB_and_chemistry_same_station_verified":False},ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
