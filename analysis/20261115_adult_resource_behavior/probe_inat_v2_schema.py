#!/usr/bin/env python3
"""One original iNat v2 observation source shape probe, no outcomes or raw ids in logs."""
import json,urllib.parse,urllib.request
BASE="https://api.inaturalist.org/v2/observations"
FIELDS="id,user.id,taxon.id,taxon.ancestor_ids,taxon.name,taxon.rank,ofvs.id,ofvs.name,ofvs.value,photos.url,photos.license_code,photos.attribution,observed_on,quality_grade"
q={"field:Butterfly Puddling Substrate":"","per_page":3,"page":1,"fields":FIELDS,"order":"asc","order_by":"id","verifiable":"any"}
url=BASE+"?"+urllib.parse.urlencode(q)
req=urllib.request.Request(url,headers={"User-Agent":"gesui-sourceonly-v2-schema-probe/1.0","Accept":"application/json"})
with urllib.request.urlopen(req,timeout=65) as s:
    raw=s.read(1000000)
d=json.loads(raw)
if not isinstance(d,dict) or not isinstance(d.get("results"),list):raise ValueError("v2 not observation response")
sample=d["results"][0] if d["results"] else {}
structure={k:(sorted(v) if isinstance(v,dict) else [sorted(a) for a in v[:2]] if isinstance(v,list) and v and isinstance(v[0],dict) else type(v).__name__) for k,v in sample.items()}
print("GESUI_INAT_V2_PROBE="+json.dumps({"total":d.get("total_results"),"resultcount":len(d["results"]),"keys":structure},sort_keys=True))
