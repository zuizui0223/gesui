#!/usr/bin/env python3
"""Inspect official 2026 EEA public share HTML/metadata for Part2 file identity.

Metadata only: HTML and JSON capped at 250 KB; do not download dataset files,
crawl private source, use credentials, or return observation values. Investigate
Nextcloud embedded manifest and publisher catalogue distribution.
"""
from __future__ import annotations
import argparse,hashlib,html,json,re
from pathlib import Path
from urllib.parse import urlparse
import requests

PUB="https://sdi.eea.europa.eu"
UUID="82e194cc-2506-411c-b341-9eb50b6ad5e0"
SHARE=PUB+"/data/"+UUID
CAT=PUB+"/catalogue/datahub/api/records/"+UUID
WEBDAV=PUB+"/webdav/datastore/public/eea_t_waterbase-water-quality-icm-2026_p_1900-2025_v01_r00"
MAX=250000
UA="gesui-ICM-2026-publisher-embedded-public-file-index-audit/20261116"

def get(url,headers=None):
    try:
        r=requests.get(url,headers={"User-Agent":UA,**(headers or {})},stream=True,timeout=(8,30),allow_redirects=False)
        b=bytearray()
        for part in r.iter_content(8192):
            b.extend(part)
            if len(b)>MAX:raise RuntimeError("metadata cap exceeded")
        return {"http_status":r.status_code,"content_type":r.headers.get("Content-Type"),"bytes":len(b),
                "response_SHA256":hashlib.sha256(b).hexdigest(),
                "location":r.headers.get("Location","")[:200]},bytes(b)
    except Exception as e:
        return {"status":"REQUEST_ERROR","error":type(e).__name__+": "+str(e)[:140]},b""

def pure_metadata():
    a,b=get(CAT)
    result={"metadata":a,"document_properties":[],"distribution_links":[],"online_file_mentions":[]}
    if "json" not in str(a.get("content_type","")).lower():return result
    try:
        obj=json.loads(b)
        def recurse(node,path):
            if isinstance(node,dict):
                for k,v in node.items():recurse(v,path+[str(k)])
            elif isinstance(node,list):
                for i,v in enumerate(node):recurse(v,path+[str(i)])
            elif isinstance(node,str):
                s=node.strip()
                if len(s)>1000:return
                low=s.lower()
                if any(t in low for t in ["distribution","resource","2026","waterbase","aggregated","webdav","zip","sqlite"]):
                    # Avoid adding third-party unrelated technical registry URLs
                    if s.startswith("https://") and any(k in s for k in ("eea.europa.eu","eionet.europa.eu")):
                        result["distribution_links"].append({"keypath":"/".join(path)[-110:],"url":s[:400]})
                    elif any(t in low for t in ["aggregateddata","part 2","zip","sqlite"]):
                        result["online_file_mentions"].append({"keypath":"/".join(path)[-110:],"value":s[:140]})
        recurse(obj,[])
        result["distribution_links"]=result["distribution_links"][:30]
        result["online_file_mentions"]=result["online_file_mentions"][:20]
        if isinstance(obj,dict):result["catalogue_keys"]=list(obj)[:30]
    except Exception as e:result["parse_error"]=str(e)[:100]
    return result

def pure_share():
    meta,b=get(SHARE,{"Accept":"text/html"})
    out={"metadata":meta}
    if "html" not in str(meta.get("content_type","")).lower():return out
    s=b.decode("utf-8",errors="replace")
    out["page_title"]=html.unescape(re.search(r"<title>(.*?)</title>",s,re.S|re.I).group(1)[:160]) if re.search(r"<title>(.*?)</title>",s,re.S|re.I) else ""
    tags=["data-filelist","data-files","data-token","data-share-token","data-sharing-token",
          "files_sharing","oc_config","_oc_webroot","publicToken","shareToken",
          "fileList","filelist","part 2","aggregateddata","zip","sqlite",
          "files-app","fileid","file-list","e2e","signature"]
    out["tag_counts"]={x:len(re.findall(re.escape(x),s,re.I)) for x in tags}
    out["script_src_names"]=[x[:180] for x in re.findall(r'<script[^>]*src=["\x27]([^"\x27]+)',s,re.I)[:35]]
    out["DOM_data_attributes"]=sorted(set(re.findall(r'data-([a-z][a-z0-9-]+)=',s)))[:100]
    out["file_ext_looking_names"]=sorted(set(re.findall(r'(?i)[\w.-]{5,140}\.(?:csv|sqlite|sqlite3|zip|7z|gz|txt|sql)',s)))[:40]
    # Inspect *only* truncated context around share indexes; redact possible
    # HTML JS server/session secrets, no raw token/cookie allowed in output.
    for marker in ["data-filelist","data-files","shareToken","publicToken","oc_appconfig","oc_webroot","files_sharing","fileList","aggregateddata"]:
        match=re.search(re.escape(marker),s,re.I)
        if match:
            snippet=s[max(match.start()-45,0):match.end()+100]
            snippet=re.sub(r'(?i)(token|nonce|requesttoken|password|csrf)(["\x27:=\s]+)([A-Za-z0-9_\-]{5,})',r"\1\2[REDACTED]",snippet)
            snippet=re.sub(r'(?i)content=["\x27][^"\x27]{40,}',"content=[REDACTED]",snippet)
            out["html_context_"+marker]=snippet[:150]
    return out

def direct():
    out={}
    for tag,url,extra in [
        ("publisher_webdav_source_redirect",WEBDAV,{}),
        ("publisher_share_to_data",SHARE+"/",{}),
        ("share_with_list_path",SHARE+"?path=%2F",{})
    ]:
        r,b=get(url,extra)
        out[tag]={**r,"not_downloaded_dataset":True}
        if "html" in str(r.get("content_type","")):
            s=b.decode("utf-8",errors="replace")
            out[tag]["embedded_filenames"]=sorted(set(re.findall(r'[\w.-]{4,100}\.(?:zip|sqlite|csv|gz|sql)',s,re.I)))[:20]
    return out

def main(args):
    a=pure_metadata()
    b=pure_share()
    c=direct()
    result={"schema":"GESUI_EEA_2026_ICM_PUBLISHER_PUBLIC_FILELIST_HTML_AND_CATALOGUE_QA_V3",
            "official_datahub":"https://www.eea.europa.eu/en/datahub/datahubitem-view/fbf3717c-cd7b-4785-933a-d0cf510542e1",
            "public_2026_dataset_UUID":UUID,
            "source_catalogue":a,
            "publisher_nextcloud_share_page":b,
            "publisher_known_public_download_redirect":c,
            "no_2026_data_table_observation_values_requested":True,
            "no_archive_downloaded":True,
            "not_a_same_station_TNTP_source_result":True}
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/"EEA_2026_PUBLIC_PUBLISHER_FILE_INDEX_DIAGNOSTIC.json").write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print("GESUI_EEA_2026_PUBLIC_FILELIST_V3="+json.dumps(result,ensure_ascii=False,sort_keys=True),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",required=True,type=Path);main(p.parse_args())
