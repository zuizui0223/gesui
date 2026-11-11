#!/usr/bin/env python3
"""Resolve ONLY public EEA 2023 final UWWTD source WebDAV GPKG original URL.

EU public metadata record 88b3b08d-3b1b-4a65-9406-365b218e09da.
Prefer exact metadata distribution URL, not speculative private Eionet path.
No numerical benthic EQR or WWTP N/P 2023 values accessed here.
"""
from __future__ import annotations
import argparse,hashlib,json,re,urllib.request,urllib.error
from datetime import datetime,timezone
from pathlib import Path

SOURCE_ID="88b3b08d-3b1b-4a65-9406-365b218e09da"
ENDPOINT=f"https://sdi.eea.europa.eu/catalogue/inspire-pds/api/records/{SOURCE_ID}?language=all"
def crawl(obj,path=()):
    if isinstance(obj,dict):
        for k,v in obj.items():yield from crawl(v,(*path,str(k)))
    elif isinstance(obj,list):
        for idx,v in enumerate(obj):yield from crawl(v,(*path,str(idx)))
    elif isinstance(obj,str) and any(k in obj.lower() for k in ("datastore","webdav","download","uwwtd-tp_p_","2022-2023","gpkg","geojson","sdiservice")):
        yield (".".join(path),obj)
def get(url):
    req=urllib.request.Request(url,headers={"Accept":"application/json","User-Agent":"gesui2023-final-UWWTD-public-only/1.0"})
    with urllib.request.urlopen(req,timeout=80) as r:
        raw=r.read(4_000_000)
        val=json.loads(raw)
        if not isinstance(val,dict):raise ValueError("EEA source metadata unavailable")
        return val,hashlib.sha256(raw).hexdigest()
def probe(url):
    out={}
    for method in ("HEAD",):
        req=urllib.request.Request(url,method=method,headers={"User-Agent":"gesui2023-FINAL-PUBLIC-source-only/1.0"})
        try:
            with urllib.request.urlopen(req,timeout=45) as r:
                out[method]={"status":r.status,"type":r.headers.get("Content-Type"),
                             "length":r.headers.get("Content-Length"),"final":r.url}
        except Exception as e:out[method]={"error":type(e).__name__+": "+str(e)[:180]}
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    document,rawsha=get(ENDPOINT)
    candidates=[]
    for key,value in crawl(document):
        # Restrict to official EEA datahub PUBLIC download paths.
        if any(x in value.lower() for x in ("/webdav/datastore/public/", "/datastore/public?path=", "datashare/public.php/webdav/")):
            for m in re.findall(r"https?://[^\s\"'<>\\,]+",value):
                v=m.rstrip(");.")
                if "sdi.eea.europa.eu" in v and v not in candidates:
                    candidates.append(v)
    result={
        "schema":"gesui_2023_original_public_FINAL_UWWTD_metadata_source_paths_and_HTTP_gate_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "original_2023_TP_public_metadata_UUID":SOURCE_ID,
        "official_public_metadata_URL":ENDPOINT,
        "original_official_metadata_JSON_SHA256":rawsha,
        "official_public_metadata_top_keys":list(document.keys())[:100],
        "source_public_EEA_2022_2023_TP_distribution_urls":candidates[:30],
        "original_official_public_download_HEAD":{q:probe(q) for q in candidates[:12]},
        "metadata_2022_2023_TP_is_NOT_proof_plant_annually_measured_NP_in_2023":True,
        "no_2016_EEA_NP_substitution_for_2023":True,
        "no_2023_EEA_river_benthic_EQR_numeric_values_read":True
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"metadata":ENDPOINT,"source_urls":result["source_public_EEA_2022_2023_TP_distribution_urls"],
                      "heads":result["original_official_public_download_HEAD"]},indent=2))
if __name__=="__main__":main()
