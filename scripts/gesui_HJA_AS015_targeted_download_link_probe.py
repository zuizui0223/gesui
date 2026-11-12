#!/usr/bin/env python3
"""Official HJA AS015 download link discovery. No ecological data analysis."""
import urllib.request,urllib.parse,html.parser,re,xml.etree.ElementTree as ET,hashlib
targets={
"new":"https://andrewsforest.oregonstate.edu/data/datacatalog/AS015",
"legacy":"https://andlter.forestry.oregonstate.edu/data/abstract.aspx?dbcode=AS015",
"legacy_ASCII_from_link":"https://andlter.forestry.oregonstate.edu/data/register/default.aspx?datapage=ascii_data&entnum=1&dbcode=AS015",
"official_EML_XML_from_link":"https://andrewsforest.oregonstate.edu/sites/default/files/lter/fsdbdata/eml_archive/AS015_2.xml",
}
class A(html.parser.HTMLParser):
    def __init__(self):super().__init__();self.a=[];self.prev=None
    def handle_starttag(self,tag,attrs):
        if tag.lower()=="a":self.prev=dict(attrs).get("href","");self.a.append((self.prev,""))
    def handle_data(self,data):
        if self.prev and self.a:
            h,txt=self.a[-1];self.a[-1]=(h,txt+data.strip())
    def handle_endtag(self,tag):
        if tag.lower()=="a":self.prev=None
for lab,url in targets.items():
    try:
        req=urllib.request.Request(url,headers={"User-Agent":"Mozilla/5.0", "Accept":"text/html,application/xml,text/csv,*/*"})
        with urllib.request.urlopen(req,timeout=35) as r:
            raw=r.read(2000000);ft=r.headers.get("Content-Type","");actual=r.url;status=r.status
        txt=raw.decode("utf-8-sig","replace")
        print("\nTARGET",lab,"HTTP",status,"len",len(raw),"mime",ft,"resolved",actual,flush=True)
        print("body-snippet",txt[:220].replace("\n"," "))
        if lab=="official_EML_XML_from_link":
            print("XML source digest",hashlib.sha256(raw).hexdigest())
            root=ET.fromstring(raw)
            for el in root.iter():
                if any(x in el.tag.lower() for x in ("objectname","url","onlinename","entityname","packageid")):
                    s=(el.text or "").strip()
                    if s:print("XML",el.tag,repr(s[:400]))
            continue
        a=A();a.feed(txt)
        candidates=[(href,label) for href,label in a.a if
              any(k in (href+" "+label).lower() for k in
                    ("as015","ascii_data","pasta","lternet","download data","download","csv","entnum","entityid"))]
        print("CANDIDATE_HREFS",repr(candidates[-55:])[:14000])
        for word in ("AS01501_v1.csv","Download Data","Data Downloads","ASCII","pasta.lternet","/package/data","datapage=ascii_data"):
            loc=txt.lower().find(word.lower())
            if loc>=0:print("CONTEXT",word,repr(txt[max(0,loc-250):loc+470])[:700])
    except Exception as e:print("ERROR",lab,type(e).__name__,str(e)[:350],flush=True)
