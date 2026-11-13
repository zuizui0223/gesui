#!/usr/bin/env python3
"""Pre-ecology national NRSA 2018/19 + 2023/24 site vs actual EPA municipal
registered external outfall geography. Geodesic SCREEN only, no water routing,
actual outfall discharge, DMR nutrient values, or ecological outcomes.
Original fixed EPA SITE UID gates and official ECHO weekly raw ZIP + SHA.
"""
from __future__ import annotations
import csv,hashlib,json,math,os,pathlib,re,time,urllib.request,zipfile
from collections import Counter,defaultdict
import numpy as np,pandas as pd
from scipy.spatial import cKDTree

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_nrsa_wwtp_source"
ROOT.mkdir(parents=True,exist_ok=True)
SITE_URL={
 "2018_19":("https://www.epa.gov/system/files/other-files/2023-01/NRSA_1819_SiteInfo.csv","1b5b2119c881290191edeb47f1d80a7711f10440e9851574cfa86fd57be770e3"),
 "2023_24":("https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_siteinfo.csv","f7d14cf1574d12828b3aa97bd7b08b95f22aec2115753d5b75514f8eed43d013")
}
ECHO="https://echo.epa.gov/files/echodownloads/npdes_outfalls_layer.zip"
SOURCE_2018_SHA="fbf53b8379d4e6564dbb51b3705b02fdd03e4a2caacad06a2314d643c8412e91"
DATE_WINDOWS={"2018_19":("2018-01-01","2019-12-31"),"2023_24":("2023-01-01","2024-12-31")}
COLS=["EXTERNAL_PERMIT_NMBR","PERM_FEATURE_NMBR","SIC_CODES","NAICS_CODES",
 "LATITUDE83","LONGITUDE83","LATLONG_TYPE","SUB_TYPE_DESC",
 "PERMIT_EFFECTIVE_DATE","PERMIT_EXPIRATION_DATE","PERMIT_TERMINATION_DATE","STATE_CODE"]
def h(s):return hashlib.sha256(s.encode()).hexdigest()
def sha(p):
 H=hashlib.sha256()
 with p.open("rb") as f:
  while b:=f.read(1048576):H.update(b)
 return H.hexdigest()
def fetch(url,label,known_sha=None):
 p=ROOT/label
 for k in range(4):
  try:
   req=urllib.request.Request(url,headers={"User-Agent":"gesui NRSA academic source-only study"})
   with urllib.request.urlopen(req,timeout=120) as r,p.open("wb") as f:
    while b:=r.read(2**20):f.write(b)
   if p.stat().st_size<100:raise RuntimeError("Source unexpectedly small")
   dig=sha(p)
   if known_sha and dig!=known_sha:raise RuntimeError(f"Official SHA mismatch on {label}")
   return p,dig
  except Exception as e:
   p.unlink(missing_ok=True)
   if k==3:raise
   time.sleep(2**k)
def source_site_ids():
 p=pathlib.Path(os.environ["NRSA_BENTHOS_18_MANIFEST"])
 if sha(p)!=SOURCE_2018_SHA:raise RuntimeError("Original 2018 benthic pre-outcome UID manifest changed")
 a=json.loads(p.read_text())
 x={r["uid"] for r in a["primary_berw_UID_subset"]}
 if len(x)!=1221:raise RuntimeError("2018 frozen site IDs changed")
 y=json.loads(pathlib.Path(os.environ["NRSA_PERIPHYTON_2ERA_CROSSWALK"]).read_text())
 maps=y["2023_24"]["uid_bridge_manifest"]
 z=set(maps["biomass"].values())&set(maps["chla"].values())&set(maps["water"].values())
 if len(z)!=2129:raise RuntimeError("2023 original before-biomass-outcome site IDs changed")
 return {"2018_19":x,"2023_24":z}
def native_sites(era,sites):
 p,origsha=fetch(SITE_URL[era][0],era+"_siteinfo.csv",SITE_URL[era][1])
 needed=["UID","SITE_ID","LAT_DD83","LON_DD83","HUC8","COMID","PSTL_CODE","AG_ECO9","DATE_COL","VISIT_NO"]
 df=pd.read_csv(p,usecols=needed,dtype=str,keep_default_na=False)
 p.unlink()
 if df.UID.duplicated().any():raise RuntimeError("Original native site table duplicate UID")
 df=df[df.UID.isin(sites)].copy()
 if len(df)!=len(sites):raise RuntimeError("Frozen native source UIDs not fully present")
 df["lat"]=pd.to_numeric(df.LAT_DD83,errors="coerce")
 df["lon"]=pd.to_numeric(df.LON_DD83,errors="coerce")
 df=df[df.lat.between(24,50)&df.lon.between(-125,-66)].copy()
 df["COMID"]=df.COMID.str.strip()
 return df,origsha
def outfalls():
 p,dig=fetch(ECHO,"official_npdes_outfalls.zip")
 with zipfile.ZipFile(p) as z:
  names=[x for x in z.namelist() if x.rsplit("/",1)[-1].upper()=="NPDES_OUTFALLS_LAYER.CSV"]
  if len(names)!=1:raise RuntimeError("Original ECHO ZIP does not expose 1 native source CSV")
  with z.open(names[0]) as f:orig=pd.read_csv(f,usecols=COLS,dtype=str,keep_default_na=False,low_memory=False)
 p.unlink()
 if len(orig)<100000:raise RuntimeError("ECHO dataset unexpectedly small")
 total=len(orig)
 m=orig.LATLONG_TYPE.str.strip().str.lower().eq("permitted feature")
 m &=orig.SUB_TYPE_DESC.str.strip().str.lower().eq("external outfall")
 code=orig.SIC_CODES.fillna("")+";"+orig.NAICS_CODES.fillna("")
 m &=code.str.contains(r"(?<!\\d)(?:221320|4952)(?!\\d)",regex=True)
 raw=orig[m].copy()
 raw["lat"]=pd.to_numeric(raw.LATITUDE83,errors="coerce")
 raw["lon"]=pd.to_numeric(raw.LONGITUDE83,errors="coerce")
 raw=raw[raw.lat.between(24,50)&raw.lon.between(-125,-66)].copy()
 raw["permit"]=raw.EXTERNAL_PERMIT_NMBR.str.strip().str.upper()
 raw["feature"]=raw.PERM_FEATURE_NMBR.str.strip().str.upper()
 raw=raw[(raw.permit!="")&(raw.feature!="")].copy()
 raw=raw.drop_duplicates(["permit","feature","lat","lon"]).copy()
 mult=raw.groupby(["permit","feature"],sort=False).size()
 bad=set(mult[mult>1].index)
 raw=raw[[ (p,f) not in bad for p,f in zip(raw.permit,raw.feature)]].copy()
 if raw.empty:raise RuntimeError("Source strict true sewage external outfalls unexpectedly zero")
 return raw,dig,{"original_official_national_rows":total,"external_municipal_valid_coordinates":len(raw),
                "ambiguous_permit_feature_pairs_excluded":len(bad)}
def unitxyz(la,lo):
 a=np.radians(np.asarray(la,dtype=float));b=np.radians(np.asarray(lo,dtype=float))
 return np.stack([np.cos(a)*np.cos(b),np.cos(a)*np.sin(b),np.sin(a)],axis=1)
def haversine(lat,lon,lat2,lon2):
 a,b,c,d=map(math.radians,[lat,lon,lat2,lon2])
 z=math.sin((c-a)/2)**2+math.cos(a)*math.cos(c)*math.sin((d-b)/2)**2
 return 12742.0176*math.asin(math.sqrt(max(0,min(1,z))))
def plausible(f,era):
 st,en=DATE_WINDOWS[era]
 eff=pd.to_datetime(f.PERMIT_EFFECTIVE_DATE,errors="coerce")
 exp=pd.to_datetime(f.PERMIT_EXPIRATION_DATE,errors="coerce")
 term=pd.to_datetime(f.PERMIT_TERMINATION_DATE,errors="coerce")
 return eff.le(pd.Timestamp(en)) & (exp.isna()|exp.ge(pd.Timestamp(st))) & (term.isna()|term.ge(pd.Timestamp(st)))
def compare(sites,f,era):
 flag=plausible(f,era)
 ff=f[flag].copy().reset_index(drop=True)
 if ff.empty:raise RuntimeError("No plausible historical registered sources")
 tree=cKDTree(unitxyz(ff.lat,ff.lon))
 rad=2*math.sin(25/6371.0088/2)
 record=[]
 source_stats=Counter()
 for x in sites.itertuples(index=False):
  ix=tree.query_ball_point(unitxyz([x.lat],[x.lon])[0],rad+1e-10)
  cand=[]
  for j in ix:
   fac=ff.iloc[j]
   dist=haversine(x.lat,x.lon,fac.lat,fac.lon)
   if dist<=25:
    cand.append({"permit":fac.permit,"feature":fac.feature,"outfall_lat":float(fac.lat),
                 "outfall_lon":float(fac.lon),"geodesic_km":round(dist,5)})
  cand.sort(key=lambda y:(y["geodesic_km"],y["permit"],y["feature"]))
  source_stats["n_sites"]+=1
  source_stats["n_sites_with_COMID"]+=int(str(x.COMID).isdigit())
  source_stats["near5"]+=int(any(z["geodesic_km"]<=5 for z in cand))
  source_stats["near10"]+=int(any(z["geodesic_km"]<=10 for z in cand))
  source_stats["near25"]+=int(bool(cand))
  record.append({"uid":str(x.UID),"site_id":str(x.SITE_ID),"comid":str(x.COMID),
                 "huc8":str(x.HUC8),"state":str(x.PSTL_CODE),"eco9":str(x.AG_ECO9),
                 "date_col":str(x.DATE_COL),"visit_no":str(x.VISIT_NO),
                 "site_lat":float(x.lat),"site_lon":float(x.lon),
                 "candidate_outfalls_25km":cand})
 near=[r for r in record if r["candidate_outfalls_25km"]]
 source_stats["distinct_HUC8_near25"]=len({x["huc8"] for x in near})
 source_stats["distinct_states_near25"]=len({x["state"] for x in near})
 source_stats["distinct_permits_near25"]=len({z["permit"] for x in near for z in x["candidate_outfalls_25km"]})
 source_stats["distinct_features_near25"]=len({(z["permit"],z["feature"]) for x in near for z in x["candidate_outfalls_25km"]})
 return record,dict(source_stats),len(ff)
def main():
 siteids=source_site_ids()
 src,source_sha,profile=outfalls()
 result={"schema":"gesui_NRSA_EPA2026_NPDES_POTW_EXTERNAL_OUTFALL_sourceONLY_2018_19_2023_24",
         "EPA_ECHO_official_outfall_zip_sha256":source_sha,"ECHO_2026_live_report_historical_incompleteness":True,
         "source_profile":profile,"eras":{},
         "original_EPA_2018_benthos_UID_manifest_sha256":SOURCE_2018_SHA,
         "biodiversity_or_nutrient_values_never_read":True,
         "nearby_outfall_is_NOT_verified_water_reaching_site":True,
         "not_actual_DMR_discharge":True}
 for era in SITE_URL:
  ss,shaorig=native_sites(era,siteids[era])
  rows,st,nfac=compare(ss,src,era)
  result["eras"][era]={"native_site_source_sha256":shaorig,"source_sites":st,
                        "date_plausible_external_outfalls":nfac,"source_records":rows}
  print("NRSA_POTW_SOURCE",era,json.dumps(st,sort_keys=True),flush=True)
  if st["near25"]<20 or st["distinct_permits_near25"]<15 or st["distinct_states_near25"]<5:
   result["eras"][era]["proceed_to_network"]=False
  else:result["eras"][era]["proceed_to_network"]=True
 out=ROOT/"NRSA_NPDES_POTW_2ERA_SOURCE_ONLY_GEO25.json"
 out.write_text(json.dumps(result,sort_keys=True,indent=2))
 print("NRSA_POTW_SOURCE_ONLY_RESULTS",json.dumps({k:{"stats":v["source_sites"],"outfalls":v["date_plausible_external_outfalls"],"route_gate":v["proceed_to_network"]} for k,v in result["eras"].items()},sort_keys=True))
if __name__=="__main__":main()
