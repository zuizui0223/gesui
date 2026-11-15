#!/usr/bin/env python3
"""Recover file names and sizes from official PUBLIC EEA 2026 share metadata.
Public catalogue link only; never download any dataset or disclose cookies.
"""
import base64,hashlib,html,html.parser,json,os,pathlib,re,urllib.request,urllib.error,xml.etree.ElementTree as ET

ROOT="https://sdi.eea.europa.eu"
TOKEN="ZcJ37WK6oaWonFJ"  # Public token embedded in EEA's official direct-download HTML.
CANDIDATES=[
("share_ui",ROOT+"/data/82e194cc-2506-411c-b341-9eb50b6ad5e0","GET"),
("dav_get",ROOT+"/public.php/dav/files/"+TOKEN+"/","GET"),
("dav_options",ROOT+"/public.php/dav/files/"+TOKEN+"/","OPTIONS"),
("dav_propfind",ROOT+"/public.php/dav/files/"+TOKEN+"/","PROPFIND"),
("share_ocs",ROOT+"/ocs/v2.php/apps/files_sharing/api/v1/shareinfo?token="+TOKEN,"GET"),
("ocs_app",ROOT+"/index.php/apps/files_sharing/shareinfo?token="+TOKEN,"GET"),
]
class AttrParser(html.parser.HTMLParser):
    def __init__(self):
        super().__init__();self.tags=[]
    def handle_starttag(self,tag,attrs):
        data=dict(attrs)
        if (data.get("id","").startswith("initial-state") or
            "data-initial-state" in data or "data-share" in str(data).lower() or
            tag=="script" and "initial" in str(data).lower()):
            self.tags.append((tag,data))
def request(url,method):
    headers={"User-Agent":"gesui-EEA-2026-publisher-public-share-safe-file-index/20261115",
             "Accept":"application/xml,application/json,text/html,*/*;q=0.1"}
    if method=="PROPFIND":
        headers.update({"Depth":"1","X-Requested-With":"XMLHttpRequest","Content-Type":"application/xml"})
        body=b'<?xml version="1.0"?><d:propfind xmlns:d="DAV:"><d:prop><d:displayname/><d:getcontentlength/><d:getcontenttype/></d:prop></d:propfind>'
    else:body=None
    if "/ocs/" in url:headers["OCS-APIREQUEST"]="true"
    req=urllib.request.Request(url,data=body,method=method,headers=headers)
    try:
        with urllib.request.urlopen(req,timeout=28) as res:
            blob=res.read(350000)
            result={"status":res.status,"final_url":res.url,
                    "response_type":res.headers.get("Content-Type"),
                    "content_length":res.headers.get("Content-Length"),
                    "bounded_sha":hashlib.sha256(blob).hexdigest(),
                    "bytes_read":len(blob)}
        return result,blob
    except urllib.error.HTTPError as e:
        raw=e.read(4000)
        return {"status":e.code,"error_bounded_sha":hashlib.sha256(raw).hexdigest(),
                "error_first_text":raw[:300].decode("utf-8","replace")},b""
    except Exception as e:
        return {"error":type(e).__name__+":"+str(e)[:220]},b""
def decode_ui(raw):
    s=raw.decode("utf-8","replace")
    p=AttrParser();p.feed(s)
    found=[]
    for tag,attrs in p.tags[:60]:
        state=attrs.get("data-initial-state") or attrs.get("value") or ""
        data=html.unescape(state)
        try:
            val=base64.b64decode(data).decode("utf-8","replace")
            ob=json.loads(val)
            files=sorted(set(re.findall(r"[^/\"'\\\s]{1,145}\.(?:zip|sqlite|csv|gz|txt|sql|7z)",val,re.I)))[:80]
            found.append({"id":attrs.get("id"),"tag":tag,"state_decoded_size":len(val),
                 "top_level_keys":list(ob)[:30] if isinstance(ob,dict) else None,
                 "file_names":files,"file_listing_keys":sorted(set(re.findall(r'"(fileName|fileSize|fileId|name|path)"',val)))})
        except Exception as e:
            found.append({"id":attrs.get("id"),"tag":tag,
                          "data_attributes":sorted(attrs),"state_size":len(data),"failure":type(e).__name__})
    return {"initial_tags":found,
            "public_download_links":sorted(set(re.findall(r"/s/[A-Za-z0-9]+(?:/download)?",s)))[:12],
            "file_names_in_html":sorted(set(re.findall(r"[^/\"'\\\s]{1,145}\.(?:zip|sqlite|csv|gz|txt|sql|7z)",s,re.I)))[:50],
            "nextcloud_api_endpoints":sorted(set(re.findall(r"(?:(?:/index.php)?/(?:ocs/v2.php|apps/files_sharing|public.php))[^\s\"']{2,200}",s)))[:20],
            "potential_files_metadata_snippet":s[max(0,s.lower().find("initial-state")-130):s.lower().find("initial-state")+280] if "initial-state" in s.lower() else None}
def parse_xml(raw):
    try:
        root=ET.fromstring(raw)
        files=[]
        for item in root.findall(".//{DAV:}response"):
            props=item.find(".//{DAV:}prop")
            files.append({"href":item.findtext("{DAV:}href"),
                "name":props.findtext("{DAV:}displayname") if props is not None else None,
                "size":props.findtext("{DAV:}getcontentlength") if props is not None else None})
        return {"entries":files[:100]}
    except Exception:
        return {"xml_first":raw[:180].decode("utf-8","replace")}
def main():
    ret=[]
    for label,url,method in CANDIDATES:
        out,raw=request(url,method)
        out.update({"name":label,"method":method,"public_url":url})
        if label=="share_ui" and out.get("status")==200:out.update(decode_ui(raw))
        if label=="dav_propfind" and raw:out.update(parse_xml(raw))
        if label in ("share_ocs","ocs_app") and raw:out["bounded_body_excerpt"]=raw[:500].decode("utf-8","replace")
        ret.append(out)
        # Never print private cookies or huge HTML; source share route already public.
        print("PUBLISHER_SHARE_AUDIT",json.dumps(out,ensure_ascii=False,sort_keys=True)[:9000],flush=True)
    root=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_EEA_2026_public_file_index"
    root.mkdir(parents=True,exist_ok=True)
    (root/"EEA_PUBLIC_NEXTCLOUD_SHARE_INDEX.json").write_text(json.dumps({
     "schema":"GESUI_EEA_2026_PUBLIC_SHARE_FILE_LIST_V1",
     "no_full_download":True,"original_metadata_public":True,"source_trial_receipts":ret},indent=2,ensure_ascii=False))
if __name__=="__main__":main()
