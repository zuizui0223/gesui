#!/usr/bin/env python3
"""Metadata-only audit of the exact 2026 EEA Nextcloud publication share.
Never download full corpus. Recover public initial state and bounded DAV index.
"""
import base64, hashlib, html, json, os, pathlib, re, urllib.error, urllib.request
import xml.etree.ElementTree as ET

SHARE="82e194cc-2506-411c-b341-9eb50b6ad5e0"
HOST="https://sdi.eea.europa.eu"
URL=HOST+"/data/"+SHARE
DAV_ENDPOINTS=[
  HOST+"/public.php/dav/files/"+SHARE+"/",
  HOST+"/public.php/webdav/",
  HOST+"/webdav/datastore/public/eea_t_waterbase-water-quality-icm-2026_p_1900-2025_v01_r00/"
]
def req(url,method="GET",body=None,headers=None,limit=210000):
    h={"User-Agent":"gesui-EEA-2026-public-source-metadata-only/20261115","Accept":"text/html,application/xml;q=0.9,*/*;q=0.5"}
    h.update(headers or {})
    q=urllib.request.Request(url,method=method,headers=h,data=body)
    try:
        with urllib.request.urlopen(q,timeout=32) as response:
            data=response.read(limit)
            out={"request_url":url,"status":response.status,"final_url":response.url,
                 "response_type":response.headers.get("Content-Type"),
                 "content_length_header":response.headers.get("Content-Length"),
                 "byte_count_bounded":len(data),"response_sha256_of_bounded_bytes":hashlib.sha256(data).hexdigest(),
                 "redirects_to_login":"login" in response.url.lower()}
        return out,data
    except urllib.error.HTTPError as e:
        blob=e.read(2048)
        return {"request_url":url,"status":e.code,"error_sha256":hashlib.sha256(blob).hexdigest(),
                "error_excerpt":blob[:200].decode("utf-8","replace")},b""
    except Exception as e:
        return {"request_url":url,"error":type(e).__name__+":"+str(e)[:220]},b""
def get_html():
    result,data=req(URL)
    s=data.decode("utf-8","replace")
    result["html_title"]=re.findall(r"<title[^>]*>(.*?)</title",s,re.I|re.S)[:2]
    result["source_mentions"]={term:s.lower().count(term.lower()) for term in [
       "AggregatedData","DisaggregatedData","files","download","initial-state","shareToken","OC.","__NEXT_DATA__","42."]}
    attrs=re.findall(r'<[^>]*data-initial-state[^>]*>',s,re.I)
    result["state_tag_count"]=len(attrs)
    decoded=[]
    for tag in attrs[:25]:
        app=re.search(r'data-app=["\']([^"\']+)',tag,re.I)
        key=re.search(r'data-key=["\']([^"\']+)',tag,re.I)
        value=re.search(r'data-initial-state=["\']([^"\']+)',tag,re.I)
        payload=html.unescape(value.group(1)) if value else ""
        try:
            plain=base64.b64decode(payload).decode("utf-8","replace")
            obj=json.loads(plain)
            # Keep public file titles/path only, no credentials/cookies/internal tokens.
            names=sorted(set(re.findall(r"[^\"'\\\s/]{1,130}\.(?:zip|csv|sqlite|txt|sql|gz|7z|md)",plain,re.I)))[:35]
            decoded.append({"app":app.group(1) if app else None,
                            "key":key.group(1) if key else None,
                            "decoded_bytes":len(plain),"top_level_keys":list(obj)[:30] if isinstance(obj,dict) else None,
                            "file_like_names":names,
                            "contains_share_token":bool(re.search(r"shareToken",plain,re.I))})
        except Exception as e:
            decoded.append({"app":app.group(1) if app else None,"key":key.group(1) if key else None,"decode_error":type(e).__name__})
    result["public_initial_state"]=decoded
    result["file_like_strings_in_html"]=sorted(set(re.findall(r"[^\"'\\\s/]{1,150}\.(?:zip|csv|sqlite|txt|sql|gz|7z)",s,re.I)))[:60]
    result["app_path_strings"]=sorted(set(re.findall(r"(?:/data/|/s/|/public.php/)[A-Za-z0-9_./%-]{4,130}",s)))[:25]
    return result
def dav(url):
    body=b'<?xml version="1.0" encoding="utf-8"?><d:propfind xmlns:d="DAV:"><d:prop><d:displayname/><d:getcontentlength/><d:getcontenttype/><d:getetag/></d:prop></d:propfind>'
    result,data=req(url,method="PROPFIND",body=body,headers={
        "Depth":"1","X-Requested-With":"XMLHttpRequest","Content-Type":"application/xml; charset=utf-8"},limit=130000)
    if data.lstrip().startswith(b"<?xml") or data.lstrip().startswith(b"<d:multistatus"):
        try:
            root=ET.fromstring(data)
            entries=[]
            for r in root.findall(".//{DAV:}response"):
                p=r.find(".//{DAV:}prop")
                entries.append({"href":r.findtext("{DAV:}href"),
                                "name":p.findtext("{DAV:}displayname") if p is not None else None,
                                "size":p.findtext("{DAV:}getcontentlength") if p is not None else None})
            result["entries"]=entries[:120]
        except Exception as e:result["xml_error"]=type(e).__name__
    return result
def main():
    out={"schema":"GESUI_EEA_2026_PUBLIC_SHARE_METADATA_ONLY_V1","publication_url":URL,
         "no_bulk_download":True,"no_water_chemistry_biological_fit":True,
         "share_html":get_html(),"dav_trials":[dav(u) for u in DAV_ENDPOINTS]}
    for k in ("share_html",):
        x=out[k]
        print("SHARE_HTML",json.dumps({z:v for z,v in x.items() if z not in ("public_initial_state",)},sort_keys=True),flush=True)
        print("SHARE_STATE",json.dumps(x["public_initial_state"],ensure_ascii=False),flush=True)
    for x in out["dav_trials"]:
        print("PUBLIC_DAV",json.dumps(x,ensure_ascii=False,sort_keys=True),flush=True)
    target=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_EEA_2026_public_share_metadata"
    target.mkdir(parents=True,exist_ok=True)
    (target/"OFFICIAL_2026_PUBLIC_SHARE_METADATA_AUDIT.json").write_text(json.dumps(out,indent=2,ensure_ascii=False))
if __name__=="__main__":main()
