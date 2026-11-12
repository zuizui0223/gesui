#!/usr/bin/env python3
"""Source-route probe only, for CC BY 4.0 HJA AS015 version 2.
Attempts official repository and the publisher's official hosted page only;
does not create or alter ecological values. Logs status and narrowly relevant
href patterns. Original 2026-11-12 run 43124861340 had EDI PASTA 403.
"""
import html.parser,json,pathlib,os,re,urllib.request,urllib.error,time
PATHS={
 "HJA_new_catalog":"https://andrewsforest.oregonstate.edu/data/datacatalog/AS015",
 "HJA_legacy_catalog":"https://andlter.forestry.oregonstate.edu/data/abstract.aspx?dbcode=AS015",
 "HJA_legacy_detail":"https://andlter.forestry.oregonstate.edu/data/abstractdetail.aspx?dbcode=AS015",
 "EDI_portal_dataviewer":"https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-and.12548.2",
 "EDI_portal_metadataviewer":"https://portal.edirepository.org/nis/metadataviewer?packageid=knb-lter-and.12548.2",
 "EDI_PASTA_alternate_names":"https://pasta.edirepository.org/package/name/eml/knb-lter-and/12548/2",
 "EDI_PASTA_primary_metadata":"https://pasta.lternet.edu/package/metadata/eml/knb-lter-and/12548/2",
 "EDI_PASTA_primary_names":"https://pasta.lternet.edu/package/name/eml/knb-lter-and/12548/2"
}
class Links(html.parser.HTMLParser):
    def __init__(self):
        super().__init__();self.links=[]
    def handle_starttag(self,t,a):
        if t.lower()!="a":return
        v=dict(a)
        href=v.get("href","")
        if any(k in href.lower() for k in ("as015",".csv","data/",".zip","pasta","lternet","dataviewer","download","entity")):
            if len(href)<1200:self.links.append(href)
result={}
for label,url in PATHS.items():
    try:
        req=urllib.request.Request(url,headers={
          "User-Agent":"Mozilla/5.0 (compatible; AcademicSourceQA/1.0)",
          "Accept":"text/html,text/plain,application/xml,application/json,text/csv,*/*"})
        with urllib.request.urlopen(req,timeout=25) as r:
            data=r.read(800000)
            code=r.status;final=r.url
        text=data.decode("utf8","replace")
        ls=Links()
        try:ls.feed(text)
        except Exception:pass
        result[label]={"code":code,"bytes_read":len(data),"final_url":final,
          "links_matching":list(dict.fromkeys(ls.links))[:120],
          "body_head":text[:160].replace("\n"," ")}
    except Exception as e:
        result[label]={"error":type(e).__name__,"detail":str(e)[:150]}
    print(label,":",result[label]["code"] if "code" in result[label] else result[label],flush=True)
    if "links_matching" in result[label]:
        print("matched-links",str(result[label]["links_matching"][:25])[:3500],flush=True)
root=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_AS015_route_probe"
root.mkdir(exist_ok=True)
(root/"AS015_OFFICIAL_SOURCE_ROUTE_PROBE.json").write_text(json.dumps(result,indent=2))
