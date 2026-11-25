#!/usr/bin/env python3
"""Retrieve publisher-declared SI.1 specimen data, without ecological fitting."""
import csv,hashlib,io,json,pathlib,re,sys,urllib.parse,zipfile
import requests,openpyxl
from bs4 import BeautifulSoup
DOI="10.1017/S0266467426100510"
ARTICLE="https://www.cambridge.org/core/journals/journal-of-tropical-ecology/article/choice-of-carrion-bait-for-monitoring-tropical-butterflies/2F44CCD3E080BDA4F06009AD4F236294"
COMP=DOI+".sm001"
OUT=pathlib.Path(sys.argv[1]);OUT.mkdir(parents=True,exist_ok=True)
S=requests.Session()
S.headers["User-Agent"]="Academic butterfly ecological source integrity audit"
allowed={"doi.org","api.crossref.org","cambridge.org","www.cambridge.org","static.cambridge.org","resolve.cambridge.org","assets.cambridge.org","content.cambridge.org"}
def get(url,limit=12*1024*1024):
    original=url
    for _ in range(10):
        assert urllib.parse.urlparse(url).scheme=="https" and urllib.parse.urlparse(url).hostname in allowed,("unreviewed origin",url)
        with S.get(url,stream=True,allow_redirects=False,timeout=(15,45)) as r:
            if r.is_redirect:
                url=urllib.parse.urljoin(url,r.headers["Location"]);continue
            r.raise_for_status();b=bytearray()
            for chunk in r.iter_content(1<<16):
                b.extend(chunk);assert len(b)<=limit
            return bytes(b),{"requested":original,"resolved":url,"status":r.status_code,"bytes":len(b),"sha256":hashlib.sha256(b).hexdigest(),"type":r.headers.get("Content-Type"),"disposition":r.headers.get("Content-Disposition")}
    raise RuntimeError("redirect loop")
def describe(b):
    if zipfile.is_zipfile(io.BytesIO(b)):
        with zipfile.ZipFile(io.BytesIO(b)) as z:
            if "xl/workbook.xml" in z.namelist():
                w=openpyxl.load_workbook(io.BytesIO(b),read_only=True,data_only=True)
                data=[{"sheet":s.title,"rows":s.max_row,"columns":s.max_column,"header_candidates":[[str(x)[:200] if x is not None else None for x in r] for r in s.iter_rows(min_row=1,max_row=min(3,s.max_row),values_only=True)]} for s in w]
                w.close();return "xlsx",data
            return "zip", [{"member":x.filename,"bytes":x.file_size} for x in z.infolist()]
    if b.startswith(b"\xd0\xcf\x11\xe0"):return "xls",None
    try:
        s=b.decode("utf-8-sig")
        if not re.search(r"<(?:html|!doctype|head|body|script)\b",s[:4000],re.I) and "\n" in s:
            first=s.splitlines()[:3]
            dialect=csv.Sniffer().sniff("\n".join(first),delimiters=",;\t")
            return "csv",{"header_candidates":list(csv.reader(first,dialect)),"lines":len(s.splitlines())}
    except Exception:pass
    return None,None
r={"status":"SOURCE_STOP","component":COMP,"attempts":[],"source_only":True}
todo=["https://doi.org/"+COMP]
try:
    b,m=get("https://api.crossref.org/works/"+COMP);j=json.loads(b)["message"];assert j["DOI"].lower()==COMP.lower()
    (OUT/"crossref_component.json").write_bytes(b);r["crossref"]=m
    for u in [j.get("URL"),j.get("resource",{}).get("primary",{}).get("URL"),*[x.get("URL") for x in j.get("link",[])]]:
        if u and u not in todo and urllib.parse.urlparse(u).hostname in allowed:todo.append(u)
except Exception as ex:r["crossref_error"]=str(ex)[:400]
todo.append(ARTICLE)
for u in todo:
    try:
        b,m=get(u);fmt,schema=describe(b)
        if fmt:
            (OUT/("SI1."+fmt)).write_bytes(b);r.update(status="SOURCE_ACQUIRED_SCHEMA_REVIEW",source=m,format=fmt,schema=schema);break
        soup=BeautifulSoup(b,"html.parser");links=[]
        for tag in soup.find_all(["a","link","meta"],href=True):
            h=urllib.parse.urljoin(m["resolved"],tag.get("href",""))
            context=" ".join(tag.stripped_strings)+" "+str(tag.attrs)
            # Only publisher-declared supplement #1 or direct SI.1 URL, never guess a file identifier.
            if re.search(r"S0266467426100510(?:sup0*1|[_\.-]sm0*1)|supplementary[-_]?material[-_]?1(?:\b|/)",h,re.I) or (re.search(r"\bSI[._ ]?1\b|supplementary material 1\b",context,re.I) and re.search(r"\.(csv|xls[x]?|zip)(?:\?|$)",h,re.I)):
                if urllib.parse.urlparse(h).hostname in allowed and h not in todo:links.append(h)
        for h in sorted(set(links)):todo.append(h)
        r["attempts"].append({"meta":m,"declared_links":links})
        (OUT/("landing_"+m["sha256"][:12]+".html")).write_bytes(b)
    except Exception as ex:r["attempts"].append({"url":u,"error":type(ex).__name__+": "+str(ex)[:700]})
(OUT/"receipt.json").write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n")
print("DIRECT BAIT SOURCE "+json.dumps(r,ensure_ascii=False),flush=True)
assert r["status"]=="SOURCE_ACQUIRED_SCHEMA_REVIEW",r["status"]
