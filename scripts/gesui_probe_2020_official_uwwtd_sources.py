#!/usr/bin/env python3
"""EEA UWWTD 2019–2020 historical source URL feasibility probe, no EQR outcomes.

Purpose: discover a publicly accessible *official* source archive and
document what really exists. Never treat a 2022 plant table as 2020 data.
No biological EQR numeric scores or other EEA Biology data are touched.
"""
from __future__ import annotations
import argparse, hashlib, json, re, urllib.request, urllib.error
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path

UA="gesui/2020-source-only-archive-verification (scientific source inventory)"
TIMEOUT=23
MAX_BYTES=400000
OFFICIAL={
 "EEA_2019_2020_treatment_public_spatial_GPKG_directory":"https://sdi.eea.europa.eu/webdav/datastore/public/eea_v_4258_100_m_uwwtd-tp_p_2019-2020_v08_r00",
 "EEA_2019_2020_discharge_public_spatial_GPKG_directory":"https://sdi.eea.europa.eu/webdav/datastore/public/eea_v_4258_100_m_uwwtd-disc-pts_p_2019-2020_v08_r00",
 "EEA_2019_2020_treatment_public_spatial_REST":"https://water.discomap.eea.europa.eu/arcgis/rest/services/UWWTD/SDI_TreatmentPlants_public_ETRS89_v1/MapServer?f=pjson",
 "EEA_2019_2020_discharge_public_spatial_REST":"https://water.discomap.eea.europa.eu/arcgis/rest/services/UWWTD/SDI_DischargePoints_public_ETRS89_v1/MapServer?f=pjson",
 "EEA_UWWTD_datahub":"https://www.eea.europa.eu/en/datahub/datahubitem-view/6244937d-1c2c-47f5-bdf1-33ca01ff1715",
 "EEA_historical_2019_2020_tabular_landing":"https://sdi.eea.europa.eu/data/3a69dbaf-15d8-443d-a0ff-fa2c4bb8ab0d/Tabular/",
 "EEA_historical_tabular_catalogue_candidate":"https://sdi.eea.europa.eu/catalogue/srv/api/records/3a69dbaf-15d8-443d-a0ff-fa2c4bb8ab0d",
 "EEA_2019_2020_discharge_public_metadata":"https://sdi.eea.europa.eu/catalogue/srv/api/records/fb78cfb7-dc6f-4294-95e9-e67f96b285c5",
 "EEA_2019_2020_discharge_internal_metadata_ONLY":"https://sdi.eea.europa.eu/catalogue/srv/api/records/e68d4c7c-9273-4cea-a4e6-b48abd980c83",
 "EEA_2018_historical_tabular_reference":"https://sdi.eea.europa.eu/data/dbbccc77-ec41-4b19-8ec5-093bba1ba7a8",
 "EEA_historical_2021_2020_legacy_index":"https://www.eea.europa.eu/data-and-maps/data/waterbase-uwwtd-urban-waste-water-treatment-directive-9",
 "EEA_2021_2020_reporting_schema":"https://dd.eionet.europa.eu/datasets/3485",
 "EEA_2021_2020_final_source_access_list":"https://sdi.eea.europa.eu/catalogue/srv/api/search/records/_search?from=0&size=40",
}
KEYWORDS=("2019-2020","data call 2021","waterbase","2020","uwwtd","tabular","download")

class Links(HTMLParser):
    def __init__(self): super().__init__();self.items=[]
    def handle_starttag(self,tag,attrs):
        if tag=="a":
            a=dict(attrs);self.items.append((a.get("href"),a.get("title"),a.get("download")))
def query(url, method="GET", payload=None):
    headers={"User-Agent":UA,"Accept":"application/json,text/html,application/xml,*/*"}
    if payload is not None:headers["Content-Type"]="application/json"
    req=urllib.request.Request(url,headers=headers,method=method,data=payload)
    try:
        with urllib.request.urlopen(req,timeout=TIMEOUT) as rsp:
            body=rsp.read(MAX_BYTES+1)
            clipped=len(body)>MAX_BYTES
            body=body[:MAX_BYTES]
            return {"http_status":rsp.status,"final_url":rsp.url,
                    "content_type":rsp.headers.get("Content-Type"),
                    "content_length_header":rsp.headers.get("Content-Length"),
                    "bytes_read":len(body),"truncated":clipped,
                    "body":body.decode("utf-8","replace")}
    except urllib.error.HTTPError as e:
        return {"http_status":e.code,"error":str(e),"final_url":e.url,
                "content_type":e.headers.get("Content-Type")}
    except Exception as e:return {"error":type(e).__name__+": "+str(e)}

def flatten(node,parent="",level=0):
    if level>8:return []
    found=[]
    if isinstance(node,dict):
        for k,v in node.items():
            if isinstance(v,str) and any(x in v.lower() for x in
                ("2019-2020","uwwtd-waterbase","2021 - final","sdi.eea.europa.eu/datashare/")):
                found.append({"key":parent+"."+str(k),"value":v[:450]})
            elif isinstance(v,(list,dict)):
                found+=flatten(v,parent+"."+str(k),level+1)
    elif isinstance(node,list):
        for i,v in enumerate(node[:80]):
            found+=flatten(v,parent+"["+str(i)+"]",level+1)
    return found[:80]

def excerpt(content):
    hits=[]
    for m in re.finditer(r"2019.2020|data.call.2021|uwwtd.waterbase|datashare/s/",content,re.I):
        i=m.start()
        hits.append(re.sub(r"\\s+"," ",content[max(0,i-180):i+350])[:480])
        if len(hits)>=15:break
    return hits

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=Path);a=p.parse_args()
    notes={}
    for name,url in OFFICIAL.items():
        got=query(url)
        body=got.pop("body","")
        parser=Links()
        if "<html" in body.lower():
            try:parser.feed(body)
            except Exception:pass
        x={**got,"body_SHA256_first_400KB":hashlib.sha256(body.encode()).hexdigest() if body else None,
           "search_excerpts":excerpt(body),
           "candidate_anchors":[{"href":h[:500],"title":str(t)[:120]}
                for h,t,_ in parser.items if h and any(x in h.lower()+" "+str(t).lower() for x in KEYWORDS)][:30],
           "candidate_identifiers":sorted(set(re.findall(r"eea_t_uwwtd[-\\w]*2019-2020[-\\w]*",body,re.I)))[:20]}
        try:
            j=json.loads(body)
            x["JSON_source_matches"]=flatten(j)
            if isinstance(j,dict):
                x["JSON_root_keys"]=list(j)[:30]
                x["JSON_possible_hits"]=str(j.get("hits",{}).get("total",""))[:150]
        except Exception:pass
        notes[name]=x
        print("GESUI_2020_PUBLIC_SOURCE_ACCESS_PROBE",json.dumps({"name":name,
             "http":got.get("http_status"),"error":got.get("error"),
             "final_url":got.get("final_url"),"candidate_count":len(x["candidate_anchors"]),
             "identifier_count":len(x["candidate_identifiers"]),
             "search_excerpts":x["search_excerpts"][:2]},ensure_ascii=False),flush=True)
    # Structured public metadata search is discovery only; record failure honestly.
    # A catalogue result is NOT a verified public original 2020 nutrient-mass table.
    result={
        "schema":"gesui_EEA_historical_2020_source_access_gate_v01",
        "utc":datetime.now(timezone.utc).isoformat(),
        "candidate_urls_are_from_official_public_EEA_pages":True,
        "actual_original_2020_UWWTD_ZIP_downloaded":False,
        "actual_2020_reported_P_mass_verified":False,
        "actual_2020_reported_N_mass_verified":False,
        "actual_2020_biological_EQR_opened":False,
        "source_only_2020_station_folds_frozen":False,
        "NO_2020_ECOLOGICAL_CLAIM":True,
        "official_source_probes":notes,
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_2020_SOURCE_ACCESS_AUDIT_SAVED",
          json.dumps({"sources_probed":len(notes),"EQR_opened":False}),flush=True)

if __name__=="__main__":main()
