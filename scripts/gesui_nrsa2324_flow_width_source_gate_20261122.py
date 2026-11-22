#!/usr/bin/env python3
"""SOURCE-ONLY EPA NRSA23/24 NHDPlus modeled flow and field width QA.
No biological outcome values, nutrient concentrations or model coefficients.
"""
import argparse,collections,csv,hashlib,json,math,pathlib,statistics

SOURCE_SHA={
 "siteinfo":"f7d14cf1574d12828b3aa97bd7b08b95f22aec2115753d5b75514f8eed43d013",
 "physicalhabitat":"8a6094bb87f8a57cf45a35a1272b137e6a9327cefd158f4bee88b025eaf9b2f4"
}
def file_sha(p):
 h=hashlib.sha256()
 with open(p,"rb") as f:
  for b in iter(lambda:f.read(1048576),b""):h.update(b)
 return h.hexdigest()
def source(path):
 with path.open("r",newline="",encoding="utf-8-sig") as f:
  rd=csv.DictReader(f)
  hdr=rd.fieldnames
  by_uid={}
  for rec in rd:
   key=rec.get("UID","").strip()
   if not key:continue
   if key in by_uid:raise RuntimeError("Original duplicate source UID in "+str(path))
   by_uid[key]=rec
 return hdr,by_uid
def fnum(s):
 try:
  z=float(s)
  return z if math.isfinite(z) else None
 except (TypeError,ValueError):return None
def quant(v):
 if not v:return None
 v=sorted(v)
 return {"n":len(v),"min":v[0],"q25":v[(len(v)-1)//4],"median":statistics.median(v),
         "q75":v[(3*(len(v)-1))//4],"max":v[-1]}
def pearson(a,b):
 if len(a)<3:return None
 ma=statistics.mean(a);mb=statistics.mean(b)
 x=[z-ma for z in a];y=[z-mb for z in b]
 den=math.sqrt(sum(z*z for z in x)*sum(z*z for z in y))
 return sum(x[i]*y[i] for i in range(len(a)))/den if den>0 else None
def main():
 p=argparse.ArgumentParser()
 p.add_argument("--siteinfo",required=True,type=pathlib.Path)
 p.add_argument("--physicalhabitat",required=True,type=pathlib.Path)
 p.add_argument("--out",required=True,type=pathlib.Path)
 q=p.parse_args()
 for k,pin in SOURCE_SHA.items():
  if file_sha(getattr(q,k))!=pin:raise RuntimeError("Official source digest changed "+k)
 header,sites=source(q.siteinfo);hdr2,hab=source(q.physicalhabitat)
 required={"MEAN_FLOW","FLOW_SIZE","WSAREASQKM","URBN_NRS23","HUC8",
   "SITE_ID","VISIT_NO","DATE_COL"}
 if not required.issubset(header):raise RuntimeError("Official site table missing hydrology field")
 if not {"UID","XWIDTH","PROTOCOL","XCDENMID","XCDENBK","SITE_ID","VISIT_NO","DATE_COL"}.issubset(hdr2):
  raise RuntimeError("PHab source missing required field")
 if len(set(sites)&set(hab))!=2139:raise RuntimeError("Prior five-table 2139 native UID source changed")
 mismatch=0
 results={}
 for proto in ("WADEABLE","BOATABLE"):
  rows=[];num=collections.Counter();categories=collections.Counter()
  for uid,h in hab.items():
   s=sites.get(uid)
   if s is None:continue
   if (s["SITE_ID"],s["VISIT_NO"],s["DATE_COL"])!=(h["SITE_ID"],h["VISIT_NO"],h["DATE_COL"]):
    mismatch+=1;continue
   if h["PROTOCOL"]!=proto:continue
   flow=fnum(s["MEAN_FLOW"]);area=fnum(s["WSAREASQKM"]);width=fnum(h["XWIDTH"])
   can=fnum(h["XCDENMID"] if proto=="WADEABLE" else h["XCDENBK"])
   flow_ok=flow is not None and flow>0
   area_ok=area is not None and area>0
   width_ok=width is not None and width>0
   can_ok=can is not None and 0<=can<=100
   num["total"]+=1
   num["flow_gt0"]+=int(flow_ok);num["area_gt0"]+=int(area_ok)
   num["width_gt0"]+=int(width_ok);num["canopy_valid"]+=int(can_ok)
   num["flow_area_width_canopy_valid"]+=int(flow_ok and area_ok and width_ok and can_ok)
   categories[("flow_size",s["FLOW_SIZE"].strip())]+=1
   categories[("urban",s["URBN_NRS23"].strip())]+=1
   rows.append({"uid":uid,"flow":flow,"area":area,"width":width,
    "canopy":can,"huc8":s["HUC8"],"site_id":s["SITE_ID"],
    "valid":flow_ok and area_ok and width_ok and can_ok})
  valid=[r for r in rows if r["valid"]]
  results[proto]={
    "source_rows":dict(num),
    "source_valid_physical_sites":len({r["site_id"] for r in valid}),
    "source_valid_HUC8":len({r["huc8"] for r in valid}),
    "mean_flow_cfs":quant([r["flow"] for r in valid]),
    "watershed_km2":quant([r["area"] for r in valid]),
    "field_width_m":quant([r["width"] for r in valid]),
    "log_flow_width_Pearson_source_covariates":pearson(
      [math.log10(r["flow"]) for r in valid],[math.log10(r["width"]) for r in valid]),
    "log_area_width_Pearson_source_covariates":pearson(
      [math.log10(r["area"]) for r in valid],[math.log10(r["width"]) for r in valid]),
    "flow_size_labels":{label:value for (typ,label),value in categories.items() if typ=="flow_size"},
    "urban_labels":{label:value for (typ,label),value in categories.items() if typ=="urban"},
    "no_biology_or_nutrient_value_accessed":True
   }
 if mismatch:raise RuntimeError("Site/date/visit mismatch: "+str(mismatch))
 out={"schema":"gesui_NRSA2324_original_hydrology_width_source_only_v1",
   "source_sha":SOURCE_SHA,"field_site_date_mismatches":mismatch,
   "source_eligible_UIDs":len(sites.keys() & hab.keys()),
   "method_results":results,"no_biology_or_nutrient_numbers_read":True}
 q.out.parent.mkdir(parents=True,exist_ok=True)
 q.out.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n")
 print("NRSA_FLOW_WIDTH_SOURCE_ONLY",json.dumps(results,sort_keys=True),flush=True)
if __name__=="__main__":main()
