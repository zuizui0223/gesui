#!/usr/bin/env python3
"""Source-only audit of official French ROSEAU/SILLAGE sludge-spreading exporter.

GET-only public discovery. Does not fabricate an export endpoint, login, fill
forms, or read butterfly / plant ecological values. Preserves raw response
hashes but not credentials, authentication tokens, or cookies.
"""
from __future__ import annotations
import argparse,hashlib,json,re,time
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urljoin,urlparse
import requests
from bs4 import BeautifulSoup

PAGES=[
  "https://www.assainissement.developpement-durable.gouv.fr/pages/data/basededonneeboue.php",
  "https://assainissement.developpement-durable.gouv.fr/pages/data/basededonneeboue.php",
  "https://www2.assainissement.developpement-durable.gouv.fr/pages/data/basededonneeboue.php",
]
MAX_HTML_BYTES=1200000
MAX_JS_BYTES=900000
MAX_JS_ASSETS=10
REDACT=re.compile(r"(token|csrf|auth|password|key|secret|session|signature|nonce)",re.I)
CANDIDATE=re.compile(r"(boues?|sludge|epandage|épandage|export|download|télécharg|csv|geojson)",re.I)
WANTED_YEARS=set(map(str,range(2017,2023)))

def sha(raw:bytes)->str:return hashlib.sha256(raw).hexdigest()

def form_info(soup,base):
    forms=[]
    all_selects=[]
    for f in soup.find_all("form"):
        method=(f.get("method") or "GET").upper()
        action=urljoin(base,f.get("action") or base)
        controls=[]
        for e in f.find_all(["select","input","button","textarea"]):
            name=e.get("name") or e.get("id") or ""
            typ=e.get("type") or ("select" if e.name=="select" else e.name)
            options=[{"value":str(z.get("value") or ""),"label":z.get_text(" ",strip=True)[:130]}
                     for z in e.find_all("option")] if e.name=="select" else []
            controls.append({"kind":e.name,"name":name,"input_type":typ,
                             "required":e.has_attr("required"),"option_count":len(options),
                             "options":options[:500],
                             "default_value_redacted":bool(REDACT.search(name))})
        forms.append({"method":method,"action":action,"same_origin":urlparse(action).netloc==urlparse(base).netloc,
                      "controls":controls})
    for e in soup.find_all("select"):
        n=e.get("name") or e.get("id") or ""
        opts=[{"value":str(z.get("value") or ""),"label":z.get_text(" ",strip=True)[:130]}
              for z in e.find_all("option")]
        all_selects.append({"name":n,"options":opts[:600],"n_options":len(opts)})
    year_opts=set()
    for e in all_selects:
        if "anne" in e["name"].lower() or "year" in e["name"].lower() or any(
            re.search(r"^20[0-3][0-9]$",v["label"]) for v in e["options"]):
            for z in e["options"]:
                for val in (z["value"],z["label"]):
                    if val in WANTED_YEARS:year_opts.add(val)
    return {"forms":forms,"selects":all_selects,"visible_2017_to_2022_options":sorted(year_opts)}

def safe_js_excerpt(text):
    excerpts=[]
    for match in CANDIDATE.finditer(text):
        chunk=text[max(0,match.start()-110):match.end()+180].replace("\n"," ")
        # Redact any tokens, usernames, passwords or signed URLs by not showing literal dynamic values.
        chunk=re.sub(r"(token|csrf|api[_-]?key|password|secret)\s*[:=]\s*['\"][^'\"]{1,200}['\"]",
                     r"\1=[REDACTED]",chunk,flags=re.I)
        excerpts.append(chunk[:300])
        if len(excerpts)>=15:break
    return excerpts

def fetch(url,session):
    req=session.get(url,timeout=40,allow_redirects=True,headers={
        "User-Agent":"gesui-Fr-biosolids-official-source-audit/2026-11-19",
        "Accept":"text/html,application/javascript,text/javascript,*/*"
    })
    if len(req.content)>MAX_HTML_BYTES and "javascript" not in req.headers.get("Content-Type","").lower():
        raise ValueError("Official page response unexpectedly large")
    return req

def probe(outdir:Path):
    session=requests.Session()
    records=[]
    first_valid=None
    for url in PAGES:
        try:
            r=fetch(url,session)
            body=r.content
            sou=BeautifulSoup(body,"html.parser")
            fields=form_info(sou,r.url)
            record={
              "requested_url":url,"status_code":r.status_code,
              "final_url":r.url,"redirect_count":len(r.history),
              "content_type":r.headers.get("Content-Type"),
              "response_bytes":len(body),"response_sha256":sha(body),
              "html_title":sou.title.get_text(" ",strip=True)[:150] if sou.title else "",
              "page_contains_boues":bool(CANDIDATE.search(sou.get_text(" ",strip=True))),
              **fields
            }
            scripts=[]
            js_inline=[]
            for i,el in enumerate(sou.find_all("script")):
                src=el.get("src")
                if src:
                    scripts.append(urljoin(r.url,src))
                else:
                    inline=el.get_text(" ",strip=False)
                    if CANDIDATE.search(inline):
                        js_inline.append({"index":i,"sha256":sha(inline.encode()),
                                          "hits":safe_js_excerpt(inline)})
            record["scripts_same_origin"]=[x for x in scripts if urlparse(x).netloc==urlparse(r.url).netloc]
            record["scripts_off_origin_count"]=sum(urlparse(x).netloc!=urlparse(r.url).netloc for x in scripts)
            record["inline_script_export_clues"]=js_inline[:15]
            links=[]
            for x in sou.find_all("a",href=True):
                href=urljoin(r.url,x["href"])
                if CANDIDATE.search(x.get_text(" ",strip=True)+" "+href):
                    links.append({"text":x.get_text(" ",strip=True)[:100],
                                  "href":href.split("?")[0][:300],
                                  "same_origin":urlparse(href).netloc==urlparse(r.url).netloc})
            record["export_like_links"]=links[:70]
            records.append(record)
            if r.ok and record["page_contains_boues"] and first_valid is None:
                first_valid=(record,session)
        except Exception as exc:
            records.append({"requested_url":url,"probe_error":str(exc)[:400]})
    assets=[]
    if first_valid:
        page,sess=first_valid
        # Public on-host JS only, no opaque API calls.
        for src in page.get("scripts_same_origin",[])[:MAX_JS_ASSETS]:
            try:
                rr=fetch(src,sess)
                raw=rr.content
                if len(raw)>MAX_JS_BYTES:
                    assets.append({"script_url":src.split("?")[0],"status_code":rr.status_code,
                                   "bytes":len(raw),"too_large_to_inspect":True})
                    continue
                decoded=raw.decode("utf-8","replace")
                if src.split("?")[0].endswith("/public/js/exportbasededonneeBoue_min.js"):
                    (outdir/"OFFICIAL_PUBLIC_ROSEAU_BOUES_EXPORT_SCRIPT.js").write_bytes(raw)
                assets.append({"script_url":src.split("?")[0],"status_code":rr.status_code,
                               "bytes":len(raw),"sha256":sha(raw),
                               "export_related_clues":safe_js_excerpt(decoded)})
                time.sleep(.1)
            except Exception as e:
                assets.append({"script_url":src.split("?")[0],"error":str(e)[:200]})
    good=[z for z in records if z.get("status_code")==200 and z.get("page_contains_boues")]
    all_years=set().union(*[set(z.get("visible_2017_to_2022_options",[])) for z in good]) if good else set()
    all_form_names=set().union(*[set(e["name"] for f in z.get("forms",[]) for e in f["controls"]) for z in good]) if good else set()
    return {
       "schema":"gesui_french_sillage_ROSEAU_sludge_spreading_public_form_GET_ONLY_source_probe_v1",
       "UTC":datetime.now(timezone.utc).isoformat(),
       "official_page_probes":records,
       "same_origin_public_script_audit":assets,
       "source_gate":{
         "official_page_accessible":bool(good),
         "form_with_post_or_get_action_detected":any(z["forms"] for z in good),
         "form_fields_found":sorted(all_form_names),
         "2017_to_2022_years_explicit_in_server_HTML":sorted(all_years),
         "all_six_years_visible_server_HTML":WANTED_YEARS.issubset(all_years),
         "actual_CSV_or_GeoJSON_downloaded":False,
         "application_commune_codes_verified":False,
         "actual_spreading_event_verified":False,
         "source_status":"FORM_DISCOVERY_ONLY" if good else "OFFICIAL_PAGE_ACCESS_STOP",
         "do_not_open_EQR_NUTRI_HNS":True,
       },
       "warnings":[
         "Dynamic form may load options via JavaScript; absence in server HTML does not establish source data absent",
         "Commune d'epandage (municipality listed for spreading) is not the field/parcel polygon",
         "No POST/download performed; form controls must be validated before approved source export",
         "Absence from reported application municipalities never proves no sludge usage"
       ]
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    report=probe(a.out)
    (a.out/"FR_ROSEAU_SILLAGE_OFFICIAL_PUBLIC_FORM_SOURCE_PROBE.json").write_text(
        json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    g=report["source_gate"]
    print("GESUI_FR_BIOSOLIDS_OFFICIAL_FORM_SOURCE_PROBE",json.dumps(g,ensure_ascii=False),flush=True)
    for z in report["official_page_probes"]:
        print("SOURCE_PAGE",z.get("requested_url"),z.get("status_code"),z.get("final_url"),
              "forms",len(z.get("forms",[])),"years",z.get("visible_2017_to_2022_options"),flush=True)

if __name__=="__main__":main()
