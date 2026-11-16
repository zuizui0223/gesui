#!/usr/bin/env python3
"""Official EEA 2026 ICM Nextcloud PUBLIC SHARE file metadata diagnostic v4.

Discovery: public HTML served under /data/<datasetUUID> uses Nextcloud app
webroot /datashare, not /. Previous /public.php/dav endpoint was WRONG ROOT.
Only public WebDAV Depth 1 metadata; never data payloads, private creds,
session tokens, or raw chemistry records. Response body capped at 350kB.
"""
from __future__ import annotations
import argparse
from base64 import b64encode
import hashlib
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from urllib.parse import urlparse
import requests

PUB="https://sdi.eea.europa.eu"
UUID="82e194cc-2506-411c-b341-9eb50b6ad5e0"
SHARE=PUB+"/data/"+UUID
ROOT="/datashare"
MAX_BYTES=350000
UA="gesui-20261116-official-EEA-ICM-2026-public-WebDAV-file-index-QA/1"
PROP=b'<?xml version="1.0" encoding="utf-8"?><d:propfind xmlns:d="DAV:"><d:prop><d:displayname/><d:getcontentlength/><d:resourcetype/><d:getcontenttype/></d:prop></d:propfind>'

def req(url,method="GET",auth=None):
    hdr={"User-Agent":UA,"Accept":"application/xml,application/json,text/html,*/*"}
    if method=="PROPFIND":
        hdr["Depth"]="1"
        hdr["Content-Type"]="application/xml"
    try:
        r=requests.request(method,url,headers=hdr,data=PROP if method=="PROPFIND" else None,
            timeout=(8,25),allow_redirects=False,stream=True,auth=auth)
        parts=bytearray()
        for b in r.iter_content(8192):
            parts+=b
            if len(parts)>MAX_BYTES:raise RuntimeError("Metadata response too large, source STOP")
        out={"http_status":r.status_code,"content_type":r.headers.get("Content-Type"),
             "response_bytes":len(parts),
             "source_SHA256":hashlib.sha256(parts).hexdigest(),
             "redirect_path":urlparse(r.headers.get("Location","")).path[:180]}
        return out,bytes(parts)
    except Exception as e:
        return {"error":type(e).__name__+": "+str(e)[:140]},b""

def scan_public_html():
    status,b=req(SHARE)
    out={"meta":status,"public_html_webroot":ROOT}
    if "html" not in str(status.get("content_type","")).lower():return out,[UUID]
    s=b.decode("utf-8","replace")
    root=re.search(r'var\s+_oc_webroot\s*=\s*"([^"]+)"',s)
    out["discovered_oc_webroot"]=root.group(1) if root else None
    if out["discovered_oc_webroot"]!=ROOT:
        raise RuntimeError("Official public Nextcloud app webroot changed")
    out["filename_matches_in_html"]=sorted(set(re.findall(r'[A-Za-z0-9_.-]{5,80}\.(?:zip|sqlite|csv|gz|txt|db)',s,re.I)))[:30]
    out["data_attributes"]={}
    for label in ["url","name","placeholder-focus","owner","protected"]:
        match=re.search(r'data-'+label+r'=["\x27]([^"\x27]{0,200})',s)
        if match:
            out["data_attributes"][label]=match.group(1)
    # Record only public share-id patterns; no session/request/CSRF tokens.
    probable=set([UUID])
    pub_patterns=[
        r'public\.php/dav/files/([A-Za-z0-9_-]{5,90})',
        r'/datashare/s/([A-Za-z0-9_-]{5,90})',
        r'"shareToken":"([A-Za-z0-9_-]{5,90})"',
        r'data-sharing-token=["\x27]([A-Za-z0-9_-]{5,90})'
    ]
    for pat in pub_patterns:
        probable.update(re.findall(pat,s))
    out["public_share_candidates_n"]=len(probable)
    # Do not log share token except public UUID; needed just for DAV calls.
    return out,sorted(probable)[:3]

def parse_xml(raw):
    ns={"d":"DAV:"}
    root=ET.fromstring(raw)
    resp=root.findall(".//d:response",ns)
    entries=[]
    for rec in resp:
        name=rec.findtext(".//d:displayname",default="",namespaces=ns)
        href=rec.findtext("d:href",default="",namespaces=ns)
        sz=rec.findtext(".//d:getcontentlength",default="",namespaces=ns)
        isdir=rec.find(".//d:collection",ns) is not None
        entries.append({"name":name[:160], "size_bytes":int(sz) if sz and sz.isdigit() else None,
                        "is_directory":isdir,
                        "relative_href_path":urlparse(href).path[-170:]})
    return entries

def public_dav_method(root,token,legacy=False,with_auth=False):
    if legacy:
        url=PUB+root+"/public.php/webdav/"
    else:
        url=PUB+root+"/public.php/dav/files/"+token+"/"
    auth=(token,"") if with_auth else None
    meta,raw=req(url,method="PROPFIND",auth=auth)
    out={"endpoint_kind":"legacy" if legacy else "public_dav_v2",
         "base_application_root":root,
         "authenticated_with_public_share_token":with_auth,
         **meta,"site_data_not_downloaded":True}
    if raw and ("xml" in str(meta.get("content_type","")).lower() or raw.lstrip().startswith(b"<?xml")):
        try:
            entries=parse_xml(raw)
            out["n_source_listing_nodes"]=len(entries)
            out["public_file_entries"]=entries[:100]
        except Exception as e:out["xml_parse_error"]=type(e).__name__+": "+str(e)[:100]
    return out

def main(a):
    page,tokens=scan_public_html()
    checks=[]
    for token in tokens:
        checks.append(public_dav_method(ROOT,token,legacy=False,with_auth=False))
        checks.append(public_dav_method(ROOT,token,legacy=False,with_auth=True))
        checks.append(public_dav_method(ROOT,token,legacy=True,with_auth=True))
    out={"schema":"GESUI_20261116_OFFICIAL_2026_ICM_PUBLISHER_CORRECT_DATASHARE_WEBDAV_V4",
         "official_dataset_url":"https://www.eea.europa.eu/en/datahub/datahubitem-view/fbf3717c-cd7b-4785-933a-d0cf510542e1",
         "nextcloud_page":page,
         "correct_nested_root_webdav_probes":checks,
         "official_2026_publisher_filelist_metadata_found":any(x.get("n_source_listing_nodes",0)>1 for x in checks),
         "actual_2026_part2_filename_found":any(
             any("aggregated" in y.get("name","").lower() and not y.get("is_directory") for y in x.get("public_file_entries",[]))
             for x in checks),
         "no_exhaustive_site_chemistry_or_TNTP_value_acquired":True,
         "no_data_archive_downloaded":True,
         "no_wastewater_causal_or_HNS_EIV_model":True}
    a.output.mkdir(parents=True,exist_ok=True)
    (a.output/"OFFICIAL_2026_ICM_DATASHARE_DAV_SOURCE_FILE_INDEX.json").write_text(
        json.dumps(out,ensure_ascii=False,indent=2)+"\n")
    print("GESUI_EEA_2026_DATASHARE_V4="+json.dumps(out,ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",required=True,type=Path)
    main(ap.parse_args())
