#!/usr/bin/env python3
"""Observed river nutrient-controlled physical habitat decomposition,
EPA NRSA2023-24, posthoc same original source cohorts. No sewage causality.
"""
import argparse,json,pathlib
import pandas as pd
import gesui_nrsa2324_tp_canopy_benthos_interaction_20261122 as orig
import gesui_nrsa2324_ambient_np_habitat_additive_posthoc_20261122 as base
V={
 "B_plus_N":base.BASE+base.N,
 "B_plus_N_width":base.BASE+base.N+["log_width"],
 "B_plus_N_canopy":base.BASE+base.N+["canopy"],
 "B_plus_N_width_canopy":base.BASE+base.N+["log_width","canopy"]
}
D=[
 ("width_after_N","B_plus_N","B_plus_N_width"),
 ("canopy_after_N","B_plus_N","B_plus_N_canopy"),
 ("width_after_N_canopy","B_plus_N_canopy","B_plus_N_width_canopy"),
 ("canopy_after_N_width","B_plus_N_width","B_plus_N_width_canopy")]
def main():
 p=argparse.ArgumentParser()
 for k in orig.SHA:p.add_argument("--"+k,required=True,type=pathlib.Path)
 p.add_argument("--out",required=True,type=pathlib.Path)
 a=p.parse_args()
 data={}
 for k,h in orig.SHA.items():
  path=getattr(a,k)
  if orig.sha(path)!=h:raise RuntimeError("Original EPA SHA mismatch "+k)
  data[k]=pd.read_csv(path,dtype=str,keep_default_na=False,low_memory=False)
 df,qa,n=orig.panel(data)
 if n!=2139 or sum(qa.values()):raise RuntimeError("Original same-UID field source changed")
 output={}
 for method in ("BERW","BETB"):
  output[method]={}
  for response in ("MMI_BENT","EPT_fraction"):
   d=orig.cohort(df,method,response)
   if len(d)<250 or d.huc8.nunique()<100:raise RuntimeError("Original source gate changed")
   y=d[response].to_numpy(float)
   fitted={k:base.oof(d,y,cols) for k,cols in V.items()}
   ref=float(fitted["B_plus_N_width_canopy"][0].mean())
   origref=base.EXPECTED[(method,response)]
   if abs(ref-origref)>max(1e-9,1e-7*abs(origref)):
    raise RuntimeError("Whole original HUC8 heldout baseline MSE mismatch")
   comps={}
   for tag,baseline,full in D:
    comps[tag]=base.gain(fitted[baseline][0],fitted[full][0])
    b=fitted[baseline][1];f=fitted[full][1]
    comps[tag]["visit_weighted_gain_pct"]=100*(b-f)/b
   output[method][response]={"UID":len(d),"SITE_ID":d.site.nunique(),"HUC8":d.huc8.nunique(),
     "original_reference_exact_reproduced":True,
     "model_HUC8_MSE":{k:float(v[0].mean()) for k,v in fitted.items()},
     "comparisons":comps}
   print("NRSA_WIDTH_VS_CANOPY",method,response,json.dumps(comps,sort_keys=True),flush=True)
 result={"schema":"gesui_NRSA2324_width_vs_canopy_adjusted_ambientNP_posthoc_v1",
   "source_UID":n,"original_shas":orig.SHA,
   "models":V,"results":output,
   "explicitly_posthoc_after_prior_interaction_and_combined_main_results":True,
   "not_sewage_derived_TP":True,"not_restoration_causal":True}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
 print("NRSA_WIDTH_CANOPY_SUMMARY",json.dumps({m:{r:{k:round(q["HUC8_equal_MSE_gain_pct"],3)
   for k,q in x["comparisons"].items()} for r,x in rows.items()} for m,rows in output.items()},sort_keys=True))
if __name__=="__main__":main()
