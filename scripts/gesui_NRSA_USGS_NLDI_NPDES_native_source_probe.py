#!/usr/bin/env python3
"""Independent NLDI original API source-only compatibility probe.
No ecological value, no DMR discharge values, no wastewater causal claims.
"""
from __future__ import annotations
import hashlib,json,os,pathlib,urllib.error,urllib.parse,urllib.request,time
import pandas as pd
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_nrsa_nldi_probe"
ROOT.mkdir(parents=True,exist_ok=True)
BASE="https://api.water.usgs.gov/nldi/linked-data"
URL="https://www.epa.gov/system/files/other-files/2023-01/NRSA_1819_SiteInfo.csv"
EXPECTED="1b5b2119c881290191edeb47f1d80a7711f10440e9851574cfa86fd57be770e3"
def api(path,timeout=48):
 u=BASE+path
 try:
  req=urllib.request.Request(u,headers={"User-Agent":"gesui-NRSA-NHDPlusV2-source-QA/1.0","Accept":"application/json"})
  with urllib.request.urlopen(req,timeout=timeout) as r:
   raw=r.read(2097152)
   val=json.loads(raw)
   if isinstance(val,list):
    return {"status":r.status,"root_type":"list","n":len(val),"list_head":val[:3],
            "response_sha256":hashlib.sha256(raw).hexdigest()}
   feats=val.get("features",[]) if isinstance(val,dict) else []
   return {"status":r.status,"root_type":type(val).__name__,"keys":list(val.keys())[:15] if isinstance(val,dict) else [],
           "features_n":len(feats),"feature_head":[{"geometry_type":(z.get("geometry") or {}).get("type"),"properties":z.get("properties",{})} for z in feats[:3]],
           "response_sha256":hashlib.sha256(raw).hexdigest()}
 except Exception as e:
  return {"status":getattr(e,"code",None),"error":type(e).__name__+":"+str(e)[:250]}
def fetch_sites():
 p=ROOT/"site.csv"
 req=urllib.request.Request(URL,headers={"User-Agent":"gesui-nrsa-nldi/1.0"})
 with urllib.request.urlopen(req,timeout=120) as r,p.open("wb") as f:
  while b:=r.read(1048576):f.write(b)
 if hashlib.sha256(p.read_bytes()).hexdigest()!=EXPECTED:raise RuntimeError("Source SHA changed")
 x=pd.read_csv(p,dtype=str,keep_default_na=False,usecols=["UID","SITE_ID","COMID","DATE_COL","VISIT_NO"])
 p.unlink()
 x=x[(x.COMID.str.fullmatch(r"[0-9]+",na=False)) & x.UID.ne("")].copy()
 x["hash"]=x.UID.map(lambda v:hashlib.sha256(v.encode()).hexdigest())
 return x.sort_values("hash").head(5).to_dict("records")
def main():
 probes={}
 probes["registered_sources"]=api("?f=json")
 probes["official_provider_nrsa"]=api("/epa_nrsa?limit=2&offset=0&f=json")
 probes["official_provider_npdes"]=api("/npdes?limit=2&offset=0&f=json")
 src=fetch_sites()
 probes["NRSA_2018_sample"]=[{"native_uid_sha":hashlib.sha256(s["UID"].encode()).hexdigest(),
                             "site_id":s["SITE_ID"],"comid":s["COMID"]} for s in src]
 for idx,s in enumerate(src):
  cid=s["COMID"]
  prefix=f"/comid/{cid}/navigation"
  for op in ("UT/npdes","DM/npdes","UT/flowlines"):
   d=5 if op=="UT/flowlines" else 50
   path=f"{prefix}/{op}?f=json&distance={d}&excludeGeometry=true"
   probes[f"n{idx}_COMID_{op.replace('/','_')}"]=api(path)
  # NRSA native site ID format is tested, not assumed.
  if idx==0:probes["epa_nrsa_site_id_lookup"]=api("/epa_nrsa/"+urllib.parse.quote(s["SITE_ID"],safe="")+"?f=json")
 out=ROOT/"USGS_NLDI_NRSA_NPDES_NAVIGATION_SOURCE_PROBE.json"
 out.write_text(json.dumps(probes,indent=2,sort_keys=True,default=str))
 print("NRSA_NLDI_PROVIDER_PROBE")
 print(json.dumps(probes,sort_keys=True,default=str)[:16000])
if __name__=="__main__":main()
