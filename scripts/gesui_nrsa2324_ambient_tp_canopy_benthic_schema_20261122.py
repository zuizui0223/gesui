#!/usr/bin/env python3
"""Original EPA NRSA 2023-24 water TP + canopy + benthic source HEADERS/UID QA.

NO biological numerical responses, no stream TP values read or fitted;
only source headers, identifiers, original metadata and source availability.
"""
import argparse,collections,csv,hashlib,json,pathlib
EXPECTED={
 "siteinfo":"f7d14cf1574d12828b3aa97bd7b08b95f22aec2115753d5b75514f8eed43d013",
 "waterchem":"b34d3f9c4b4bdc763cb5e53a856bbab329f1646d248f223f78e9404bc4560d71",
 "physicalhabitat":"8a6094bb87f8a57cf45a35a1272b137e6a9327cefd158f4bee88b025eaf9b2f4",
 "benthicmmi":"5ef607e8e4726ed8c3d24b9d15fe351e7300784859654b812fa595c09a4009dc",
 "benthiccount":"5290d679f7eb444e97c8d45096beb4531b450186cef190580f41f4b9430b0805"
}
def sha(p):
 h=hashlib.sha256()
 with open(p,"rb") as f:
  for chunk in iter(lambda:f.read(1048576),b""):h.update(chunk)
 return h.hexdigest()
def inspect(path,kind):
 uids=set();meta={};total=0
 with open(path,"r",encoding="utf-8-sig",newline="") as f:
  rd=csv.DictReader(f)
  cols=rd.fieldnames
  if "UID" not in cols:raise RuntimeError("Official EPA native UID missing "+kind)
  for r in rd:
   total+=1;uid=str(r["UID"]).strip()
   if not uid:continue
   uids.add(uid)
   if kind=="siteinfo":
    entry={k:r.get(k,"").strip() for k in ("SITE_ID","VISIT_NO","DATE_COL","HUC8","PSTL_CODE","AG_ECO9")}
   elif kind=="benthiccount":
    entry={k:r.get(k,"").strip() for k in ("SITE_ID","VISIT_NO","DATE_COL","SAMPLE_TYPE")}
   elif kind=="physicalhabitat":
    entry={k:r.get(k,"").strip() for k in ("SITE_ID","VISIT_NO","DATE_COL","PROTOCOL")}
   else:
    entry={k:r.get(k,"").strip() for k in ("SITE_ID","VISIT_NO","DATE_COL")}
   if uid in meta and meta[uid]!=entry:raise RuntimeError("Same-UID metadata disagreement within "+kind)
   meta[uid]=entry
 return {"rows":total,"unique_uid":len(uids),"columns":cols},uids,meta
def main():
 ap=argparse.ArgumentParser()
 for k in EXPECTED:ap.add_argument("--"+k,type=pathlib.Path,required=True)
 ap.add_argument("--out",type=pathlib.Path,required=True)
 a=ap.parse_args();reports={};uids={};metadata={}
 for name,hash_value in EXPECTED.items():
  path=getattr(a,name);actual=sha(path)
  if actual!=hash_value:raise RuntimeError("EPA source SHA changed "+name+":"+actual)
  reports[name],uids[name],metadata[name]=inspect(path,name)
 base=set.intersection(*uids.values())
 d=collections.Counter()
 for uid in base:
  v=metadata["siteinfo"][uid];b=metadata["benthiccount"][uid];h=metadata["physicalhabitat"][uid]
  if v["SITE_ID"]!=b["SITE_ID"] or v["SITE_ID"]!=h["SITE_ID"] or v["DATE_COL"]!=b["DATE_COL"] or v["DATE_COL"]!=h["DATE_COL"]:
   d["same_UID_site_date_mismatch"]+=1;continue
  d[b["SAMPLE_TYPE"]+"|"+h["PROTOCOL"]]+=1
 watercols=reports["waterchem"]["columns"]
 benthiccols=reports["benthiccount"]["columns"]
 mmicols=reports["benthicmmi"]["columns"]
 physicalcols=reports["physicalhabitat"]["columns"]
 output={
  "schema":"NRSA_2023_24_AMBIENT_TP_CANOPY_BENTHOS_SOURCE_SCHEMA_GATE_v1",
  "source_sha256":EXPECTED,"sources":reports,
  "all_five_shared_original_UID":len(base),"original_method_counts":dict(d),
  "official_water_tp_fields_found":{n:n in watercols for n in ["PTL_RESULT","PTL_UNITS","PTL_NARS_FLAG","NTL_RESULT","NTL_UNITS","COND_RESULT"]},
  "field_canopy_columns":{n:n in physicalcols for n in ["XCDENMID","XCDENBK","XEMBED","XWIDTH","PROTOCOL"]},
  "benthic_count_columns":{n:n in benthiccols for n in ["ORDER","TOTAL300","SAMPLE_TYPE","TAXA_ID"]},
  "published_benthic_MMI_candidate_fields":[c for c in mmicols if "MMI" in c.upper() or "EPT" in c.upper()],
  "site_numeric_outcomes_read":False,
  "not_wastewater_dose":True,
  "no_biological_coefficient_estimated":True}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
 print("SOURCE_ONLY_AMBIENT_TP_CANOPY_BENTHOS_SCHEMA")
 print(json.dumps({"shared_uid":len(base),"methods":dict(d),"benthic_source_fields":output["benthic_count_columns"],
   "MMI_fields":output["published_benthic_MMI_candidate_fields"],
   "tp_fields":output["official_water_tp_fields_found"],
   "canopy_fields":output["field_canopy_columns"]},sort_keys=True))
if __name__=="__main__":main()
