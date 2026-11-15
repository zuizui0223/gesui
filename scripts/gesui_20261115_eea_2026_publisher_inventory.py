#!/usr/bin/env python3
"""Inventory official 2026 EEA Water Quality direct-download files metadata only.

Never download the 4.3 GB publisher source. This checks safe WebDAV PROPFIND,
HEAD/Range response and publisher source resource availability. Source only.
"""
import hashlib
import json
import os
import pathlib
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET

URL="https://sdi.eea.europa.eu/webdav/datastore/public/eea_t_waterbase-water-quality-icm-2026_p_1900-2025_v01_r00"
ATTEMPTS=[
 ("PROPFIND_parent","PROPFIND",URL.rstrip("/")+"/"),
 ("PROPFIND_source","PROPFIND",URL),
 ("HEAD_source","HEAD",URL.rstrip("/")+"/"),
 ("GET_source_first32768","GET",URL.rstrip("/")+"/")
]
def go(name,method,url):
    h={"User-Agent":"gesui-EEA-2026-ICM-source-inventory/20261115","Accept":"*/*"}
    if method=="PROPFIND":
        h["Depth"]="1"
        h["Content-Type"]="application/xml; charset=utf-8"
        body=b'<?xml version="1.0" encoding="utf-8"?><d:propfind xmlns:d="DAV:"><d:prop><d:displayname/><d:getcontentlength/><d:getcontenttype/><d:getetag/></d:prop></d:propfind>'
    else:body=None
    if method=="GET":h["Range"]="bytes=0-32767"
    req=urllib.request.Request(url,method=method,headers=h,data=body)
    try:
        with urllib.request.urlopen(req,timeout=32) as response:
            raw=response.read(130000 if method=="PROPFIND" else (32768 if method=="GET" else 512))
            status=response.status
            headers={k.lower():v for k,v in response.headers.items() if k.lower() in ("content-type","content-range","content-length","accept-ranges","last-modified","etag","location")}
        out={"name":name,"method":method,"url":url,"status":status,
             "headers":headers,"first_bytes_sha256":hashlib.sha256(raw).hexdigest(),
             "bytes_read":len(raw)}
        if method=="PROPFIND" and raw.lstrip().startswith(b"<"):
            try:
                root=ET.fromstring(raw)
                ns="{DAV:}"
                listing=[]
                for r in root.findall(".//"+ns+"response"):
                    hnode=r.find(ns+"href")
                    props=r.find(".//"+ns+"prop")
                    listing.append({"href":hnode.text if hnode is not None else None,
                        "size":props.findtext(ns+"getcontentlength") if props is not None else None,
                        "etag":props.findtext(ns+"getetag") if props is not None else None,
                        "type":props.findtext(ns+"getcontenttype") if props is not None else None})
                out["source_directory_entries"]=listing
            except Exception as e:out["xml_error"]=type(e).__name__
        if method=="GET":out["source_prefix_type"]=raw[:100].decode("utf-8","replace")
        return out
    except urllib.error.HTTPError as e:
        raw=e.read(1024)
        return {"name":name,"url":url,"method":method,"http_status":e.code,
                "error_body_sha256":hashlib.sha256(raw).hexdigest(),"error_prefix":raw[:120].decode("utf-8","replace")}
    except Exception as e:
        return {"name":name,"url":url,"method":method,"error":type(e).__name__+":"+str(e)[:250]}

def main():
    outcomes=[go(*t) for t in ATTEMPTS]
    out={"schema":"GESUI_EEA_2026_WATER_QUALITY_ICM_ORIGINAL_DIRECT_SOURCE_FILE_INVENTORY",
         "source_publisher_url":URL,
         "no_original_data_bulk_downloaded":True,
         "no_2018_2022_missing_concentrations_imputed":True,
         "calls":outcomes}
    for r in outcomes:
        print("PUBLISHER_INVENTORY",json.dumps({k:v for k,v in r.items() if k not in ("source_directory_entries",)},sort_keys=True),flush=True)
        if r.get("source_directory_entries"):print("PUBLISHER_FILES",json.dumps(r["source_directory_entries"][:80],sort_keys=True),flush=True)
    dest=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_WQ2026_source_inventory"
    dest.mkdir(parents=True,exist_ok=True)
    (dest/"PUBLISHER_DIRECT_FILE_INVENTORY.json").write_text(json.dumps(out,sort_keys=True,indent=2))
if __name__=="__main__":main()
