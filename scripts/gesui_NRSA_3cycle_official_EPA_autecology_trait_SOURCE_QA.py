#!/usr/bin/env python3
"""Official EPA 2008/2013 original benthic autecology trait source; 2018
original benthic count taxonomic identity. STRICT trait mapping by native
taxon ID + exactly normalized reported TARGET_TAXON name + ORDER when
available. No other cross-cycle species/taxon guess, no fake missing traits,
no ecological trait regressions in this source-only audit.
"""
from __future__ import annotations
import hashlib,json,os,pathlib,re,time,urllib.request
import numpy as np,pandas as pd
T=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))
R=T/"gesui_threecycle_traits"
R.mkdir(parents=True,exist_ok=True)
native=pd.read_csv(T/"gesui_threecycle_ept_transfer"/"AUDIT_ONLY_NRSA_THREE_REAL_CYCLE_EPT_NP_PHAB.csv",
 dtype={"source_UID":str,"source_cycle":str,"huc8":str})
if len(native)!=2648:raise RuntimeError("Source fixed 3 era original n2648 missing")
counts=pd.read_csv(T/"gesui_threecycle_rarefied_taxa"/"AUDIT_ONLY_ORIGINAL_EPA_3CYCLE_TAXON_COUNTS.csv",
 dtype={"cycle":str,"UID":str,"key":str,"order":str},keep_default_na=False)
if counts.duplicated(["cycle","UID","key"]).any():raise RuntimeError("Source original counts repeated TAXA_ID")
FILES={
 "taxa2008":("https://www.epa.gov/sites/default/files/2016-06/nrsa_0809_benttaxa.csv","6165d5a7ec8d4ee3e85149255fd967f6268d870d49c57de6302ca4e124673b26"),
 "taxa2013":("https://www.epa.gov/sites/default/files/2019-04/nrsa1314_benttaxa_02012019.csv",None),
 "count2018":("https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_benthic_macroinvertebrate_count_-_data.csv","4e9b3a0798b25f67c7a96d4800b94f86544bbd623abaf119358e067b158324a7")
}
def fetch(key):
 url,expect=FILES[key];p=R/(key+".csv")
 for z in range(4):
  try:
   req=urllib.request.Request(url,headers={"User-Agent":"gesui-NRSA-original-benthic-autecology-traits-qa/20261114"})
   with urllib.request.urlopen(req,timeout=120) as q,p.open("wb") as w:
    while b:=q.read(1024*1024):w.write(b)
   H=hashlib.sha256(p.read_bytes()).hexdigest()
   if expect is not None and H!=expect:
    raise RuntimeError("Previously signed EPA original source SHA changed "+key)
   break
  except Exception:
   p.unlink(missing_ok=True)
   if z==3:raise
   time.sleep(2**z)
 df=pd.read_csv(p,dtype=str,keep_default_na=False,encoding="latin1",low_memory=False)
 p.unlink()
 df.columns=df.columns.str.upper().str.strip()
 return df,{"official_original_url":url,"publisher_byte_sha256":H,"row_n":len(df),"actual_source_columns":list(df)}
def normalize(s):
 return s.astype(str).str.strip().str.upper().str.replace(r"\s+"," ",regex=True)
traits={}
evidence={}
for key in ("taxa2008","taxa2013"):
 tab,receipt=fetch(key)
 for col in ("TAXA_ID","ORDER"):
  if col not in tab:raise RuntimeError("EPA genuine taxon list missing "+col)
 tab["id_key"]=normalize(tab["TAXA_ID"])
 tab["ord_key"]=normalize(tab["ORDER"])
 candidates=[x for x in ("TARGET_TAXON","TARGET_TAXON_NAME","TAXON_NAME","TAXA_NAME") if x in tab]
 namecol=candidates[0] if candidates else None
 if namecol is not None:tab["name_key"]=normalize(tab[namecol])
 else:tab["name_key"]=""
 receipt["taxon_name_field"]=namecol
 receipt["native_unique_taxonIDs"]=tab.id_key[tab.id_key.ne("")].nunique()
 receipt["duplicate_rows_same_original_taxonID"]=int(tab.id_key[tab.id_key.ne("")].duplicated().sum())
 receipt["source_unique_ID_multi_order"]=int((tab[tab.id_key.ne("")].groupby("id_key").ord_key.nunique()>1).sum())
 receipt["source_unique_ID_multi_name"]=int((tab[tab.id_key.ne("")].groupby("id_key").name_key.nunique()>1).sum())
 for col in ("FFG","HABIT","PTV","VOLTINISM","FFG_WSA","HABIT_WSA","PTV_WSA"):
  if col in tab:
   v=normalize(tab[col])
   receipt["actual_"+col+"_counts"]=v.value_counts().head(28).to_dict()
   receipt["n_unique_original_ID_with_"+col+"_nonempty"]=int(tab.loc[v.ne(""),"id_key"].nunique())
  else:receipt["actual_"+col+"_STATUS"]="NOT_IN_ORIGINAL_SOURCE"
 evidence[key]=receipt
 traits[key]=tab
# External survey's original count may include exact taxon name and order,
# but the data dictionary must be examined before claiming cross-cycle IDs.
B18,receipt=fetch("count2018")
B18=B18[B18.SAMPLE_TYPE.str.strip().str.upper().eq("BERW")&
        B18.UID.isin(set(native.loc[native.source_cycle.eq("2018_19"),"source_UID"]))].copy()
NAMES=[c for c in ("TARGET_TAXON","TARGET_TAXON_NAME","TAXON_NAME","TAXA_NAME") if c in B18]
namecol=NAMES[0] if NAMES else None
receipt["original_2018_count_name_field"]=namecol
receipt["original_2018_BERW_rows_on_source_fixed_946_UIDs"]=len(B18)
receipt["original_2018_target_source_taxonID_distinct"]=B18.TAXA_ID.str.strip().nunique()
receipt["real_count_taxa_name_field_MISSING_blocks_crosscycle_source_join"]=not bool(namecol)
if namecol:
 B18["id_key"]=normalize(B18.TAXA_ID)
 B18["name_key"]=normalize(B18[namecol])
 B18["ord_key"]=normalize(B18.ORDER)
else:
 B18["id_key"]=normalize(B18.TAXA_ID)
 B18["name_key"]=""
 B18["ord_key"]=normalize(B18.ORDER)
multi=(B18.groupby("id_key").name_key.nunique()>1)
receipt["same_source_2018_taxonID_multiple_target_names"]=int(multi.sum())
evidence["count2018"]=receipt
# Trait coverage among actual EPT/nonEPT specimen counts, not merely
# occurrence records. 2018 only accepted when same original native ID,
# exact published name AND taxonomic ORDER match old original EPA list.
FIELDS=["FFG","HABIT","PTV","VOLTINISM"]
trait_maps={}
for name,tab in traits.items():
 gr=tab[tab.id_key.ne("")].copy()
 agg=gr.groupby("id_key")
 consistent=(agg.ord_key.nunique().eq(1)&agg.name_key.nunique().eq(1))
 for field in FIELDS:
  if field in gr:
   gr[field+"_norm"]=normalize(gr[field])
   consistent=consistent&gr.groupby("id_key")[field+"_norm"].nunique().le(1)
 z=gr[gr.id_key.isin(set(consistent[consistent].index))].drop_duplicates("id_key").set_index("id_key")
 trait_maps[name]=z
 evidence[name]["n_ID_with_noncontradictory_order_name_and_FFG_HABIT_PTV_VOLTINISM"]=len(z)
# 2018 original source identity may not have name field. Don't claim
# join success without names even though taxonID numbering may look alike.
mapping={
 "2008_09":trait_maps["taxa2008"],
 "2013_14":trait_maps["taxa2013"]}
stats={}
for era in ("2008_09","2013_14","2018_19"):
 site=counts[counts.cycle.eq(era)].copy()
 site["key"]=normalize(site.key)
 site["order"]=normalize(site.order)
 site["count"]=pd.to_numeric(site["count"],errors="raise")
 if era!="2018_19":
  ref=mapping[era]
  site=site.merge(ref[["ord_key","name_key"]+[k+"_norm" for k in FIELDS if k+"_norm" in ref]],
           left_on="key",right_index=True,how="left",validate="many_to_one")
  good_order=site["order"].eq(site["ord_key"])&site["ord_key"].ne("")
  known_name=site["name_key"].notna()&site["name_key"].ne("")
  good_identity=good_order&known_name
 else:
  if namecol:
   b=B18[B18.TOTAL300.ne("")].copy()
   b["key"]=normalize(b.TAXA_ID)
   # Never choose conflicting name when same ID has multiple names.
   ids=b.groupby("key").name_key.nunique()
   stable=ids[ids.eq(1)].index
   b=b[b.key.isin(stable)].drop_duplicates("key").set_index("key")
   ident=b[["name_key","ord_key"]]
   site=site.merge(ident,left_on="key",right_index=True,how="left",validate="many_to_one")
   previous=[]
   for oldname in ("taxa2013","taxa2008"):
    ref=trait_maps[oldname]
    z=site.merge(ref[["name_key","ord_key"]+[k+"_norm" for k in FIELDS if k+"_norm" in ref]],
     left_on="key",right_index=True,how="left",validate="many_to_one",suffixes=("_18","_prior"))
    z["trait_verified"]=z.name_key_18.eq(z.name_key_prior)&z.ord_key_18.eq(z.ord_key_prior)&z.name_key_18.notna()&z.name_key_18.ne("")
    previous.append((oldname,z))
   use=previous[0][1]
   for k in FIELDS:
    col=k+"_norm"
    # Only nonempty source trait from 2013 if identity verified; when
    # missing, choose 2008 strictly matching ID+name+order. Contradictory
    # double matches are not resolved by privileging the later source.
    z0=previous[0][1];z1=previous[1][1]
    aa=z0[col] if col in z0 else pd.Series("",index=z0.index)
    bb=z1[col] if col in z1 else pd.Series("",index=z1.index)
    aok=z0.trait_verified & aa.fillna("").ne("")
    bok=z1.trait_verified & bb.fillna("").ne("")
    mismatch=aok&bok & aa.ne(bb)
    v=aa.where(aok,bb.where(bok,""))
    v=v.mask(mismatch,"")
    site[col]=v.values
   good_identity=previous[0][1].trait_verified|previous[1][1].trait_verified
   known_name=good_identity
  else:
   for field in FIELDS:site[field+"_norm"]=""
   good_identity=pd.Series(False,index=site.index)
   known_name=good_identity
 # Include only original orders EPT for EPT trait source coverage.
 ept=site.order.isin(("EPHEMEROPTERA","PLECOPTERA","TRICHOPTERA"))
 denom=float(site.loc[ept,"count"].sum())
 result={"native_positive_taxon_ID_rows":len(site),
   "known_name_and_order_matching_prior_original_source_fraction_of_EPT_individuals":
    float(site.loc[ept&good_identity,"count"].sum()/denom) if denom else 0,
   "source_counted_EPT_individuals_total":denom,
   "nonmatching_original_source_ID_or_name_or_ORDER_count_rows":int((~good_identity).sum())}
 for trait in FIELDS:
  key=trait+"_norm"
  can=good_identity & site[key].fillna("").ne("")
  result["true_original_"+trait+"_verified_EPT_counted_individual_coverage"]=float(site.loc[ept&can,"count"].sum()/denom) if denom else 0
  result["original_"+trait+"_known_category_counts"]=site.loc[ept&can,key].value_counts().head(20).to_dict()
 stats[era]=result
 print("NRSA_3CYCLE_OFFICIAL_BENTHIC_TRAITS_SOURCE_COVERAGE",era,json.dumps(result,sort_keys=True),flush=True)
# Conservative source gate for genuinely comparable 3-era FFG, PTV, HABIT.
passes={trait:all(stats[era]["true_original_"+trait+"_verified_EPT_counted_individual_coverage"]>=.70
                 for era in stats) for trait in ("FFG","HABIT","PTV")}
OUT={
 "schema":"gesui_NR​​SA_three_orig_EPA_autecology_trait_ID_name_ORDER_coverage_SOURCE_ONLY_v1",
 "original_EPA_trait_source_receipts":evidence,
 "original_2018_taxonID_only_not_sufficient_to_prove_crosscycle_trait_identity":True,
 "population_source_2648_native_UID_from_prior_3cycle":True,
 "original_3era_trait_EPT_count_weighted_identity_coverage":stats,
 "predeclared_trait_EPT_70pct_individual_coverage_gate":passes,
 "do_not_assign_traits_to_unknown_taxon_order":True,
 "taxon_trait_missingness_is_NOT_zero_resistance_or_zero_FFG":True,
 "no_trait_regressions_or_selection_on_positive_effects":True,
 "no_actual_sewage_source_or_emergence_services":True}
(R/"NRSA_THREECYCLE_REAL_EPA_TRAITS_IDENTITY_AND_COVERAGE_SOURCE_ONLY.json").write_text(
 json.dumps(OUT,sort_keys=True,indent=2))
print("NRSA_THREECYCLE_OFFICIAL_TRAITS_SOURCE_GATE",json.dumps(passes,sort_keys=True),flush=True)
