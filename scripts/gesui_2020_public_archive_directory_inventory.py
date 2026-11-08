#!/usr/bin/env python3
"""Source-only 2019–2020 official EEA WWTP/outfall WebDAV archive inventory.

The share token is publicly disclosed by the EEA dataset catalog redirect;
this script obtains public, not private/internal, resources. It does NOT read
EEA Biology/EQR, fit ecology models, or mistake source absence for zero loading.
"""
from __future__ import annotations
import argparse,base64,csv,hashlib,io,json,urllib.request,urllib.parse,urllib.error
import xml.etree.ElementTree as ET
from pathlib import Path

PUB_TOKEN="sptXqwkQr5g7Bp5"
BASE="https://sdi.eea.europa.eu/datashare/public.php/webdav"
SOURCES={
  "WWTP_2020_public":{
    "path":"eea_v_4258_100_m_uwwtd-tp_p_2019-2020_v08_r00",
    "subdirs":["Tabular","GPKG","Documents"]},
  "DISCHARGE_2020_public":{
    "path":"eea_v_4258_100_m_uwwtd-disc-pts_p_2019-2020_v08_r00",
    "subdirs":["GPKG"]},
}
NS={"d":"DAV:"}
AUTH="Basic "+base64.b64encode((PUB_TOKEN+":").encode()).decode()

def request(path,method="PROPFIND",max_read=1024*1024):
    url=BASE+"/"+"/".join(urllib.parse.quote(c,safe="") for c in path.strip("/").split("/"))
    headers={"Authorization":AUTH,"User-Agent":"gesui-EEA-2020-source-only-inventory"}
    if method=="PROPFIND":headers["Depth"]="1"
    elif method=="GET":headers["Range"]="bytes=0-4095"
    req=urllib.request.Request(url,method=method,headers=headers,
                data=(b'<?xml version="1.0"?><d:propfind xmlns:d="DAV:"><d:allprop/></d:propfind>'
                      if method=="PROPFIND" else None))
    try:
        with urllib.request.urlopen(req,timeout=45) as r:
            raw=r.read(max_read)
            return {"status":r.status,"url":r.url,
                    "content_type":r.headers.get("Content-Type"),
                    "content_length_header":r.headers.get("Content-Length"),
                    "data":raw}
    except urllib.error.HTTPError as e:return {"status":e.code,"url":e.url,"error":str(e)}
    except Exception as e:return {"error":type(e).__name__+": "+str(e),"url":url}

def listing(data):
    if data.get("status")!=207:return {"http":data.get("status"),"error":data.get("error")}
    root=ET.fromstring(data["data"])
    found=[]
    for x in root.findall(".//d:response",NS):
        href=x.findtext("d:href",default="",namespaces=NS)
        size=x.findtext(".//d:getcontentlength",default=None,namespaces=NS)
        ctype=x.findtext(".//d:getcontenttype",default=None,namespaces=NS)
        etag=x.findtext(".//d:getetag",default=None,namespaces=NS)
        rtype=x.find(".//d:resourcetype/d:collection",NS)
        found.append({"href":urllib.parse.unquote(href),"size_bytes":int(size) if size and size.isdigit() else None,
                      "mime":ctype,"etag":etag,"directory":rtype is not None})
    return {"http":207,"items":found,"n":len(found)}

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--out",required=True,type=Path);a=ap.parse_args()
    result={"schema":"gesui_2020_official_public_UWWTD_source_full_directory_inventory_v01",
            "EEA_Biology_EQR_scores_accessed":False,
            "original_2020_NP_reported_mass_not_yet_verified":True,
            "official_public_share_URL":"https://sdi.eea.europa.eu/datastore/public",
            "datasets":{}}
    for label,src in SOURCES.items():
        root=src["path"]
        data={"archive_root":root,"root":listing(request(root,max_read=200000))}
        for part in src["subdirs"]:
            path=root+"/"+part
            data[part]=listing(request(path,max_read=500000))
            print("GESUI_2020_OFFICIAL_PUBLIC_DIRECTORY",json.dumps({
                "dataset":label,"subdir":part,"http":data[part]["http"],
                "files":[{"name":x["href"].rsplit("/",1)[-1],"size":x["size_bytes"],"dir":x["directory"]}
                         for x in data[part].get("items",[])][:50]},ensure_ascii=False),flush=True)
        readme=request(root+"/README.md",method="GET",max_read=8192)
        data["README_status"]=readme.get("status")
        if "data" in readme:
            data["README_public_text"]=readme["data"].decode("utf-8","replace")[:7000]
        # Print CSV schema previews, no numeric source masses or biological outcomes loaded.
        allfiles=[]
        for part in src["subdirs"]:
            allfiles.extend(data[part].get("items",[]))
        selected=[x for x in allfiles if not x["directory"] and x["href"].lower().endswith(".csv")][:20]
        previews=[]
        for item in selected:
            fullpath=urllib.parse.unquote(item["href"]).split("/webdav/")[-1].lstrip("/")
            meta=request(fullpath,method="GET",max_read=6000)
            content=meta.get("data",b"").decode("utf-8-sig","replace")
            header=content.splitlines()[0] if content else ""
            delim=max([",",";","\t"],key=header.count)
            fields=next(csv.reader([header],delimiter=delim)) if header else []
            previews.append({"relative_path":fullpath,"http":meta.get("status"),
                             "header":fields[:130],"no_data_records_retained":True})
        data["CSV_schema_previews"]=previews
        result["datasets"][label]=data
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_2020_PUBLIC_EEA_ARCHIVE_INVENTORY_DONE",json.dumps({
        "datasets":{name:{k:z[k]["n"] for k in SOURCES[name]["subdirs"] if z[k].get("n") is not None}
                    for name,z in result["datasets"].items()},"no_EQR":True}),flush=True)
if __name__=="__main__":main()
