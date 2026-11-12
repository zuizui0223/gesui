#!/usr/bin/env python3
"""Diagnose GBIF HTTP 400 from a fixed 12-basin PREcount polygon source
without ecological inference, access bypass, or changing original source IDs.
"""
import collections,hashlib,json,os,pathlib,urllib.request,urllib.error,urllib.parse,time
F=pathlib.Path(os.environ["GBIF_PREFLIGHT_SOURCE_ARTIFACT"])/"HYBAS7_GBIF_ORIGINAL_SOURCE_COUNTS_URL_SHA.json"
P=pathlib.Path(os.environ["GBIF_FROZEN_WKT_SOURCE_MANIFEST"])/"HYBAS7_1219_ORIGINAL_NATIVE_POLYGONS_WKT_BEFORE_NEW_GBIF_COUNTS.json"
src=P.read_bytes()
if hashlib.sha256(src).hexdigest()!="47f597e0a798cd6155ee9fb69ce8916ed5629951d046d114feb9f667963310e4":
    raise RuntimeError("WKT original immutable source checksum changed")
rows=json.loads(F.read_text())
failed=[]
unique={}
for a in rows:
    for q in a["request_components"]:
        if q["status"]=="missing":
            lookup=(a["HYBAS_ID_L07"],q["original_component_sha256"],a["taxon"])
            unique[lookup]=(a,q)
print("GBIF stage0 original source diagnostic: missing requests",len(unique),flush=True)
for k,(item,z) in list(unique.items())[:15]:
    url=z["url"];params=urllib.parse.parse_qs(urllib.parse.urlsplit(url).query)
    g=params.get("geometry",[""])[0]
    countpoints=g.count(",")
    req=urllib.request.Request(url,headers={"User-Agent":"gesui-HYBAS7-source-HTTP400-geometry-probe/1","Accept":"application/json"})
    entry={"HYBAS_ID_L07":item["HYBAS_ID_L07"],"taxon":item["taxon"],
          "region":item["region"],"period":item["period"],
          "source_original_component_sha256":z["original_component_sha256"],
          "source_WKT_length":len(g),"source_WKT_n_comma":countpoints,
          "geometry_type":g.split(" ",1)[0],
          "http_code_previous":z.get("http_code"),
          "parameters_except_WKT":{p:v for p,v in params.items() if p!="geometry"}}
    try:
        with urllib.request.urlopen(req,timeout=45) as r:
            rr=r.read(1200)
            entry["resend_HTTP_status"]=r.status
            entry["response_head_1200"]=rr.decode("utf8","replace")[:250]
    except urllib.error.HTTPError as e:
        entry["resend_HTTP_status"]=e.code
        entry["error_detail"]=e.read(4000).decode("utf8","replace")[:1800]
    except Exception as e:
        entry["resend_error"]=type(e).__name__+":"+str(e)[:160]
    failed.append(entry)
    print("SOURCE_HTTP_FAILURE",json.dumps(entry,sort_keys=True),flush=True)
path=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"HYBAS7_GBIF_ORIGINAL_SOURCE_400_TECHNICAL_DIAGNOSTIC.json"
path.write_text(json.dumps({"schema":"H7_Original_GBIF_source_HTTP400_diagnostic_v1",
 "fixed_prepilot_basins_12":True,"source_failures":failed,
 "NO_ECOLOGICAL_COEFFICIENTS":True},indent=2,sort_keys=True))
