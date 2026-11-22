#!/usr/bin/env python3
"""Source-only official EPA RV and common flow/width-support audit.

No macroinvertebrate responses, chemistry numbers or selected effect signs read.
"""
import argparse,collections,hashlib,json,math,pathlib
from gesui_nrsa2324_flow_width_source_gate_20261122 import file_sha,source,fnum,SOURCE_SHA

def quantile(values,q):
 s=sorted(values)
 if not s:return None
 k=(len(s)-1)*q
 a=int(math.floor(k));b=int(math.ceil(k))
 return s[a]*(b-k)+s[b]*(k-a) if b!=a else s[a]
def count(rows):
 return {"visits":len(rows),
  "physical_sites":len({r["site"] for r in rows}),
  "HUC8":len({r["huc8"] for r in rows}),
  "states":len({r["state"] for r in rows})}
def main():
 ap=argparse.ArgumentParser()
 ap.add_argument("--siteinfo",type=pathlib.Path,required=True)
 ap.add_argument("--physicalhabitat",type=pathlib.Path,required=True)
 ap.add_argument("--out",type=pathlib.Path,required=True)
 a=ap.parse_args()
 for k,h in SOURCE_SHA.items():
  if file_sha(getattr(a,k))!=h:raise RuntimeError("Source SHA changed "+k)
 _,ss=source(a.siteinfo)
 _,hh=source(a.physicalhabitat)
 if len(set(ss)&set(hh))!=2139:raise RuntimeError("Source 2139 common UID changed")
 rows={}
 matched=0
 for method,proto,canopy in (("BERW","WADEABLE","XCDENMID"),("BETB","BOATABLE","XCDENBK")):
  selected=[]
  for uid,p in hh.items():
   s=ss.get(uid)
   if s is None:continue
   if (s["SITE_ID"],s["VISIT_NO"],s["DATE_COL"])!=(p["SITE_ID"],p["VISIT_NO"],p["DATE_COL"]):
    raise ValueError("Original source UID/site/visit/date mismatch")
   if p["PROTOCOL"]!=proto or s["FLOW_SIZE"].strip()!="RV":continue
   matched+=1
   flow=fnum(s["MEAN_FLOW"])
   area=fnum(s["WSAREASQKM"])
   width=fnum(p["XWIDTH"])
   can=fnum(p[canopy])
   if not(flow is not None and area is not None and width is not None and can is not None and
          flow>0 and area>0 and width>0 and 0<=can<=100):continue
   selected.append({"uid":uid,"site":s["SITE_ID"],"huc8":s["HUC8"],"state":s["PSTL_CODE"],
                    "log_flow":math.log10(flow),"log_width":math.log10(width)})
  rows[method]=selected
 quant={}
 for key,values in rows.items():
  quant[key]={"log10_flow_q05_q95":[quantile([r["log_flow"] for r in values],.05),
                                   quantile([r["log_flow"] for r in values],.95)],
              "log10_width_q05_q95":[quantile([r["log_width"] for r in values],.05),
                                    quantile([r["log_width"] for r in values],.95)]}
 fmin=max(quant[m]["log10_flow_q05_q95"][0] for m in rows)
 fmax=min(quant[m]["log10_flow_q05_q95"][1] for m in rows)
 wmin=max(quant[m]["log10_width_q05_q95"][0] for m in rows)
 wmax=min(quant[m]["log10_width_q05_q95"][1] for m in rows)
 if not(fmin<fmax and wmin<wmax):raise RuntimeError("No common source flow or width support")
 outcomes={}
 for method,values in rows.items():
  p1=values
  p2=[r for r in values if fmin<=r["log_flow"]<=fmax]
  p3=[r for r in p2 if wmin<=r["log_width"]<=wmax]
  rr={}
  for k,sub,threshold in [
   ("RV_source_only",p1,(100,70)),
   ("RV_common_flow",p2,(80,50)),
   ("RV_common_flow_and_width",p3,(60,40))]:
    ct=count(sub)
    rr[k]={**ct,"source_GO":ct["visits"]>=threshold[0] and ct["HUC8"]>=threshold[1]}
  outcomes[method]={"survey_subsets":rr,
    "original_RV_source_quartile_ranges":quant[method]}
 output={"schema":"gesui_NRSA2324_original_RV_flow_width_overlap_source_only_v1",
   "source_SHA":SOURCE_SHA,"common_UID_all_original":2139,
   "source_only_candidate_RV_protocol_rows":matched,
   "fixed_common_log10_flow_bounds":[fmin,fmax],
   "fixed_common_flow_cfs_bounds":[10**fmin,10**fmax],
   "fixed_common_log10_width_bounds":[wmin,wmax],
   "fixed_common_width_m_bounds":[10**wmin,10**wmax],
   "methods":outcomes,
   "outcome_or_water_chemistry_values_read":False,
   "RV_source_label_not_field_measured_discharge":True}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
 print("NRSA_RV_SOURCE_ONLY",json.dumps(output,sort_keys=True),flush=True)
if __name__=="__main__":main()
