#!/usr/bin/env python3
"""Source-only EPA NRSA2023 substratum and hydro-landscape field coverage QA.

NO benthic outcome values or posthoc selected biological coefficients read.
"""
from __future__ import annotations
import argparse,hashlib,json,pathlib
import numpy as np,pandas as pd
import gesui_NRSA_fourthcycle_2023_EPT_heldout_transport as prev
from gesui_nrsa2324_fourthcycle_source_only_20261122 import SHA,sha

PHYS=["PCT_FN","PCT_SAFN","LSUB_DMM","LRBS","XDEPTH","XFC_NAT","W1_HAG","XEMBED"]
HYDRO=["MEAN_FLOW","WSAREASQKM","URBN_NRS23"]
BLOCKS={"particle_grain":["PCT_FN","PCT_SAFN","LSUB_DMM"],
        "hydrological_scale":["MEAN_FLOW","WSAREASQKM","URBN_NRS23"]}
def profile(s,name):
 raw=s.astype(str).str.strip()
 if name=="URBN_NRS23":
  good=raw.isin(["Urban","Non-Urban"])
  return {"columns_present":True,"valid":int(good.sum()),
          "source_labels":raw.value_counts().head(12).to_dict(),
          "missing_or_other":int((~good).sum())}
 n=pd.to_numeric(raw,errors="coerce")
 finite=np.isfinite(n)
 if name in ["PCT_FN","PCT_SAFN"]:
  good=finite & n.between(0,100)
 elif name in ["MEAN_FLOW","WSAREASQKM","XDEPTH"]:
  good=finite & (n>0)
 elif name in ["XFC_NAT","W1_HAG","XEMBED"]:
  good=finite & (n>=0)
 else:
  good=finite
 z=n.loc[good]
 return {"columns_present":True,"valid":int(good.sum()),
     "missing_or_invalid":int((~good).sum()),
     "mean":float(z.mean()) if len(z) else None,
     "median":float(z.median()) if len(z) else None,
     "q10":float(z.quantile(.1)) if len(z) else None,
     "q90":float(z.quantile(.9)) if len(z) else None,
     "min":float(z.min()) if len(z) else None,
     "max":float(z.max()) if len(z) else None}
def main():
 a=argparse.ArgumentParser()
 for k in ("sites","water","habitat","chla","out"):
  a.add_argument("--"+k,type=pathlib.Path,required=True)
 q=a.parse_args()
 for k in ("sites","water","habitat","chla"):
  actual=sha(getattr(q,k))
  if actual!=SHA[k]:raise RuntimeError("Original EPA pinned source changed "+k+":"+actual)
 d=prev.covariates(q,{})
 if len(d)!=1165 or d.huc8.nunique()!=693 or d.site_id.nunique()!=1046:
  raise RuntimeError("2023 original prebiological source cohort changed")
 site=pd.read_csv(q.sites,dtype=str,keep_default_na=False,low_memory=False).set_index("UID")
 ph=pd.read_csv(q.habitat,dtype=str,keep_default_na=False,low_memory=False).set_index("UID")
 profiles={}
 fullsource={"siteinfo":sorted(site.columns.tolist()),
             "physicalhabitat":sorted(ph.columns.tolist())}
 for k in PHYS:
  profiles[k]=profile(ph.loc[d.index,k],k) if k in ph else {"columns_present":False}
 for k in HYDRO:
  profiles[k]=profile(site.loc[d.index,k],k) if k in site else {"columns_present":False}
 blocks={}
 for block,fields in BLOCKS.items():
  if not all(profiles[k]["columns_present"] for k in fields):
   blocks[block]={"source_go":False,"reason":"at least one required field absent from original EPA source"}
   continue
  valid=pd.Series(True,index=d.index)
  for k in fields:
   src=ph.loc[d.index,k] if k in ph else site.loc[d.index,k]
   v=src.astype(str).str.strip()
   if k=="URBN_NRS23":ok=v.isin(["Urban","Non-Urban"])
   else:
    num=pd.to_numeric(v,errors="coerce")
    ok=num.notna() & np.isfinite(num)
    if k in ("PCT_FN","PCT_SAFN"):ok=ok&num.between(0,100)
    if k in ("MEAN_FLOW","WSAREASQKM"):ok=ok&(num>0)
   valid &=ok
  huc=d.loc[valid].groupby("huc8").size()
  blocks[block]={"complete_UID":int(valid.sum()),"complete_HUC8":int(huc.size),
    "complete_fraction_UID":float(valid.mean()),
    "full_693_HUC8_fraction":float(huc.size/693),
    "source_go_preliminary":float(huc.size/693)>=.80,
    "reason":"final source GO also requires both old and never-sampled source subcohorts >=80% support"}
 output={"schema":"NRSA2023_substrate_hydro_original_source_only_20261122_v1",
   "original_EPA_SHA":{k:SHA[k] for k in ("sites","water","habitat","chla")},
   "source_original_2023_eligible_UID":len(d),"source_original_2023_HUC8":d.huc8.nunique(),
   "profiles":profiles,"candidate_blocks":blocks,
   "all_numeric_outcomes_unread":True,
   "original_original_header_field_sets":{k:[x for x in v if any(term in x.upper() for term in (
     "PCT_","LSUB","LRBS","SLOPE","GRAD","XDEPTH","W1H_","W1_HA","HAB","RBS","FIN","SED","URBN","FLOW","AREA"))][:120]
    for k,v in fullsource.items()}}
 q.out.parent.mkdir(parents=True,exist_ok=True)
 q.out.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
 print("EPA_NRSA_2023_SUBSTRATE_SOURCE_COVERAGE",json.dumps(
   {"n_UID":len(d),"HUC8":d.huc8.nunique(),
    "fields":profiles,"blocks":blocks},sort_keys=True),flush=True)
if __name__=="__main__":main()
