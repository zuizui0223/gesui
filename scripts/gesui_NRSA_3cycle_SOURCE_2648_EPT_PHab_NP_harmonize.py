#!/usr/bin/env python3
"""Blind-before-transfer validation and source harmonization of three real
original EPA NRSA survey cycles. Inputs are reconstructed from pinned
publisher CSVs by their prior approved independent GitHub scripts, and live
only inside one GitHub Actions job. Does NOT fit any new predictive models.
"""
from __future__ import annotations
import hashlib,json,os,pathlib
import numpy as np,pandas as pd

BASE=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))
OUT=BASE/"gesui_threecycle_ept_transfer"
OUT.mkdir(parents=True,exist_ok=True)
FIELDS=["source_cycle","source_UID","site_id","huc8","state","eco9","lat","lon",
        "year","stream_order","season_sin","season_cos","log_cond",
        "log_tn","log_tp","XEMBED","XCDENMID","XPCAN","XFC_OHV",
        "W1_HALL","log_width","ept_fraction"]
EXPECT={"2008_09":(586,437,48),"2013_14":(1116,706,48),"2018_19":(946,603,47)}
NAMES={
 "2008_09":("gesui0809_reconstructed_EPT","AUDIT_ONLY_2008_RAW_EPT_CHEMW_PHAB_SOURCE_PANEL.csv"),
 "2013_14":("gesui1314_independent_numeric","AUDIT_ONLY_NRSA1314_BERW_EPT_CHLA_SOURCE_PANEL.csv"),
 "2018_19":("gesui1819_triad_numeric","AUDIT_ONLY_2018_19_EPT_Chla_PHAB_SOURCE_PANEL.csv")
}
def src_cycle(name):
 folder,filename=NAMES[name]
 p=BASE/folder/filename
 if not p.is_file():raise RuntimeError("Source previously EPA reconstructed panel missing: "+str(p))
 source=pd.read_csv(p,dtype={"UID":str,"site_id":str,"huc8":str,"state":str,"eco9":str},low_memory=False)
 if "UID" not in source or source.UID.duplicated().any():
  raise RuntimeError("Original source UID repeated or missing for "+name)
 if name=="2008_09":
  # A source-verified conservative counted-fraction estimate; order unknowns
  # are not invented as EPT or non-EPT. Preserve its logical upper bound.
  df=source.rename(columns={"ept_share_raw_count_lower_bound":"ept_fraction"}).copy()
  df["ept_fraction_upper_bound"]=source["ept_share_source_upper_bound"]
 elif name=="2013_14":
  df=source.rename(columns={
   "logcond":"log_cond","logTN":"log_tn","logTP":"log_tp",
   "logwidth":"log_width"}).copy()
 elif name=="2018_19":
  flags=(source.chla_flag.fillna("").eq("") &
         source.tn_mgL_rawflag.fillna("").eq("") &
         source.tp_mgL_rawflag.fillna("").eq(""))
  df=source.loc[flags].copy()
  df=df.rename(columns={"log_afdm":"log_afdm_original"})
 else:raise RuntimeError("Unknown frozen original cycle")
 df=df.rename(columns={"UID":"source_UID"})
 df["source_cycle"]=name
 df["huc8"]=df.huc8.astype(str).str.strip()
 df["state"]=df.state.astype(str).str.strip().str.upper()
 df["eco9"]=df.eco9.astype(str).str.strip().str.upper()
 df["site_id"]=df.site_id.astype(str).str.strip()
 required=[x for x in FIELDS if x not in ("source_cycle","source_UID","site_id","huc8","state","eco9")]
 missing=[x for x in required if x not in df]
 if missing:raise RuntimeError("Original EPA cycle "+name+" lacks common six-PHab variable "+str(missing))
 # Only original previous precommitted quality criteria. Check no new rows
 # need selection; don't silently reduce sample to favor predictive signs.
 numeric=required
 bad=np.zeros(len(df),dtype=bool)
 for k in numeric:
  z=pd.to_numeric(df[k],errors="coerce")
  bad|=~np.isfinite(z.to_numpy(float))
  df[k]=z
 if bad.any():
  raise RuntimeError("Original frozen cycle "+name+" contains "+str(int(bad.sum()))+
                     " rows missing the previously source/quality-fixed common fields")
 if df[["huc8","state","eco9","site_id"]].eq("").any().any():
  raise RuntimeError("Missing native site/HUC/state/region identifier in "+name)
 if not df.huc8.str.fullmatch(r"\d{8}").all():
  invalid=df.loc[~df.huc8.str.fullmatch(r"\d{8}"),"huc8"].unique()[:8]
  raise RuntimeError("HUC8 source labels not exactly eight digits "+name+":"+str(invalid))
 if (df.ept_fraction<0).any() or (df.ept_fraction>1).any():
  raise RuntimeError("Original EPA EPT share outside biological 0–1 range")
 low,high={"2008_09":(2008,2009),"2013_14":(2013,2014),"2018_19":(2018,2019)}[name]
 if not df.year.isin((low,high)).all():
  raise RuntimeError("Observed survey years changed "+name+":"+str(df.year.value_counts().to_dict()))
 a,b,c=EXPECT[name]
 if len(df)!=a or df.huc8.nunique()!=b or df.state.nunique()!=c:
  raise RuntimeError("Original pre-source-frozen n/HUC8/states altered: "+
    name+":"+str((len(df),df.huc8.nunique(),df.state.nunique(),(a,b,c))))
 if name=="2008_09":
  if df.ept_fraction_upper_bound.isna().any() or not (df.ept_fraction_upper_bound>=df.ept_fraction).all():
   raise RuntimeError("Source 2008 original unknown ORDER upper bound invalid")
 info={
  "cycle":name,"source_n_UIDs":len(df),"source_original_HUC8":df.huc8.nunique(),
  "source_states":df.state.nunique(),"source_unique_SITE_IDs":df.site_id.nunique(),
  "source_year_counts":{str(int(k)):int(v) for k,v in df.year.value_counts().items()},
  "ecoregions":sorted(df.eco9.unique().tolist()),
  "source_huc8_sorted_sha256":hashlib.sha256("\n".join(sorted(set(df.huc8))).encode()).hexdigest(),
  "source_aggregation_huc8_not_pseudo_independent_UID":True,
  "method":"BERWW original TOTAL300 count+TAXA_ID→ORDER" if name=="2008_09" else "BERW original TOTAL300 count",
  "same_panel_selected_pre_transfer_on_original_EPA_source_chem_chla_quality":name!="2008_09"}
 df=df.sort_values("source_UID").reset_index(drop=True)
 return df,info
def main():
 panels={};evidence={}
 for cycle in EXPECT:
  d,info=src_cycle(cycle)
  panels[cycle]=d
  evidence[cycle]=info
 all_groups={k:set(df.huc8) for k,df in panels.items()}
 keys=list(panels)
 overlap={}
 for i in range(len(keys)):
  for j in range(i+1,len(keys)):
   a,b=keys[i],keys[j]
   overlap[a+"__"+b]={"n_shared_HUC8":len(all_groups[a]&all_groups[b]),
     "only_"+a:len(all_groups[a]-all_groups[b]),
     "only_"+b:len(all_groups[b]-all_groups[a])}
 universe=set.union(*all_groups.values())
 triple=len(set.intersection(*all_groups.values()))
 dataset=pd.concat([panels[x] for x in keys],ignore_index=True)
 if dataset.duplicated(["source_cycle","source_UID"]).any() or len(dataset)!=2648:
  raise RuntimeError("Wrong real 3 national cycle source sample/duplicates")
 dataset=dataset.sort_values(["source_cycle","source_UID"]).reset_index(drop=True)
 cols=FIELDS+["ept_fraction_upper_bound"]
 # Upper bound applies to 2008 only and is NA for the later cycles.
 if "ept_fraction_upper_bound" not in dataset:raise RuntimeError("Source 2008 unknown ORDER bound lost")
 dataset[cols].to_csv(OUT/"AUDIT_ONLY_NRSA_THREE_REAL_CYCLE_EPT_NP_PHAB.csv",index=False)
 desc={
  "schema":"gesui_NR​​SA_3source_cycles_REAL_BERW_EPT_ambientNP_six_PHAB_source_only_v1",
  "source_cycle_QA":evidence,
  "original_overlap_nHUC8":overlap,
  "distinct_union_HUC8":len(universe),"n_threeway_shared_HUC8":triple,
  "original_survey_UID_total":len(dataset),
  "mean_raw_EPT_by_cycle_NOT_ECOLOGICAL_EFFECT":{k:float(v.ept_fraction.mean()) for k,v in panels.items()},
  "original_category_levels":{k:{"states":sorted(d.state.unique()),"eco9":sorted(d.eco9.unique())} for k,d in panels.items()},
  "no_external_model_or_effect_coefficients_read":True,
  "original_EP​​A_periphyton_absence_in_2008":True,
  "not_actual_sewage_effluent_N_P_mass_or_ecosystem_services":True}
 path=OUT/"NRSA_THREECYCLE_ORIGINAL_SOURCES_AND_HUC8_OVERLAP_QA.json"
 path.write_text(json.dumps(desc,indent=2,sort_keys=True,ensure_ascii=False))
 print("THREECYCLE_SOURCE_ONLY_QA",json.dumps({
  "cycle_n_HUC8_state":{k:[v["source_n_UIDs"],v["source_original_HUC8"],v["source_states"]] for k,v in evidence.items()},
  "HUC8_overlap":overlap,"union":len(universe),"triple":triple,
  "eco9_category_sets":{k:sorted(d.eco9.unique()) for k,d in panels.items()}},sort_keys=True))
if __name__=="__main__":main()
