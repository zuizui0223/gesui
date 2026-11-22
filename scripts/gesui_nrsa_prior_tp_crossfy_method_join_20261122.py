#!/usr/bin/env python3
"""Merge three official EPA DMR FY source-opportunity reports by NRSA native UID.

Join ORIGINAL NRSA benthic-method and canopy metadata only (no biological outcomes).
Do not invent concentration doses, unexposed controls or cross-protocol effects.
"""
import argparse,csv,datetime,glob,hashlib,json,pathlib

BENTHIC_SHA="5290d679f7eb444e97c8d45096beb4531b450186cef190580f41f4b9430b0805"
PHAB_SHA="8a6094bb87f8a57cf45a35a1272b137e6a9327cefd158f4bee88b025eaf9b2f4"

def sha(path):
 h=hashlib.sha256()
 with open(path,"rb") as f:
  for b in iter(lambda:f.read(1048576),b""):h.update(b)
 return h.hexdigest()
def dt(x):
 s=str(x).strip()
 for f in ("%Y-%m-%d","%m/%d/%Y","%m/%d/%Y %H:%M:%S"):
  try:return datetime.datetime.strptime(s,f).date().isoformat()
  except ValueError:pass
 raise ValueError("unparsable original survey date "+repr(x))
def load_metadata(p,cols):
 result={}
 with open(p,newline="",encoding="utf-8-sig") as f:
  r=csv.DictReader(f)
  if not set(cols+["UID"]).issubset(r.fieldnames):raise ValueError("Missing EPA metadata columns")
  for x in r:
   uid=str(x["UID"]).strip()
   if not uid:continue
   value={c:str(x[c]).strip() for c in cols}
   if uid in result and result[uid]!=value:raise ValueError("Conflicting metadata within UID "+uid)
   result[uid]=value
 return result
def canopy_ok(v):
 try:
  v=float(v);return 0<=v<=100
 except (ValueError,TypeError):return False
def counts(rows):
 return {"visits":len(rows),"physical_sites":len({r["site_id"] for r in rows}),
         "HUC8":len({r["huc8"] for r in rows}),
         "states":len({r["state"] for r in rows}),
         "monthly_average_reporting_outfalls":len({h for r in rows for h in r["source_hashes"]})}
def merge_records(reports):
 if {x["fiscal"] for x in reports}!={2023,2024,2025}:raise ValueError("Three distinct fiscal reports required")
 by_uid={}
 for rep in reports:
  if rep.get("schema")!="gesui_prior_window_tp_source_only_v2":raise ValueError("Expected v2 source opportunity")
  opportunities=rep["source_only_site_opportunities"]
  if len(opportunities)!=195 or len({x["uid"] for x in opportunities})!=195:raise ValueError("Source 195 UID census changed")
  if by_uid and set(by_uid)!={x["uid"] for x in opportunities}:raise ValueError("Different UID set across FYs")
  for x in opportunities:
   uid=x["uid"];meta={k:x[k] for k in ("site_id","huc8","state","survey_date")}
   if uid not in by_uid:
    by_uid[uid]={"uid":uid,**meta,"windows":{}}
    for w in ("30","90","365"):
     by_uid[uid]["windows"][w]={"any_concentration":False,"any_mg_l":False,
      "monthly_avg_C2_MK_mg_l":False,"source_hashes":set()}
   y=by_uid[uid]
   if {k:y[k] for k in meta}!=meta:raise ValueError("UID metadata changed among fiscal scans")
   for w,src in x["windows"].items():
    dest=y["windows"][w]
    for k in ("any_concentration","any_mg_l","monthly_avg_C2_MK_mg_l"):
     dest[k] |= bool(src[k])
    dest["source_hashes"].update(src["monthly_average_exact_outfall_hashes"])
 return by_uid
def main():
 a=argparse.ArgumentParser()
 a.add_argument("--fy-dir",required=True,type=pathlib.Path)
 a.add_argument("--benthic",required=True,type=pathlib.Path)
 a.add_argument("--phab",required=True,type=pathlib.Path)
 a.add_argument("--out",required=True,type=pathlib.Path)
 z=a.parse_args()
 if sha(z.benthic)!=BENTHIC_SHA or sha(z.phab)!=PHAB_SHA:raise ValueError("Official EPA source original SHA mismatch")
 reports=[json.loads(x.read_text()) for x in sorted(z.fy_dir.glob("prior_tp_fy*.json"))]
 if len(reports)!=3:raise ValueError("Three original FY source reports required")
 merged=merge_records(reports)
 b=load_metadata(z.benthic,["SITE_ID","UNIQUE_ID","VISIT_NO","DATE_COL","SAMPLE_TYPE"])
 p=load_metadata(z.phab,["SITE_ID","UNIQUE_ID","VISIT_NO","DATE_COL","PROTOCOL","XCDENBK","XCDENMID"])
 diag={"unmatched_uid":0,"site_or_date_mismatch":0,"not_eligible_method":0,"missing_bank_canopy":0,"missing_midcanopy":0}
 allrows=[]
 for uid,x in merged.items():
  if uid not in b or uid not in p:
   diag["unmatched_uid"]+=1;continue
  a1=b[uid];a2=p[uid]
  if not (x["site_id"]==a1["SITE_ID"]==a2["SITE_ID"] and
          a1["UNIQUE_ID"]==a2["UNIQUE_ID"] and a1["VISIT_NO"]==a2["VISIT_NO"] and
          x["survey_date"]==dt(a1["DATE_COL"])==dt(a2["DATE_COL"])):
   diag["site_or_date_mismatch"]+=1;continue
  w=a1["SAMPLE_TYPE"]=="BERW" and a2["PROTOCOL"]=="WADEABLE" and canopy_ok(a2["XCDENMID"])
  bo=a1["SAMPLE_TYPE"]=="BETB" and a2["PROTOCOL"]=="BOATABLE" and canopy_ok(a2["XCDENBK"])
  if not (w or bo):
   diag["not_eligible_method"]+=1
   diag["missing_midcanopy"]+=int(a1["SAMPLE_TYPE"]=="BERW" and not canopy_ok(a2["XCDENMID"]))
   diag["missing_bank_canopy"]+=int(a1["SAMPLE_TYPE"]=="BETB" and not canopy_ok(a2["XCDENBK"]))
   continue
  for window,op in x["windows"].items():
   allrows.append({"uid":uid,"site_id":x["site_id"],"huc8":x["huc8"],"state":x["state"],
                   "survey_year":x["survey_date"][:4],"method":"BERW" if w else "BETB",
                   "window":window, "any_concentration":op["any_concentration"],
                   "any_mg_l":op["any_mg_l"],
                   "monthly_avg_C2_MK_mg_l":op["monthly_avg_C2_MK_mg_l"],
                   "source_hashes":sorted(op["source_hashes"])})
 results={}
 for win in ("30","90","365"):
  result={}
  for protocol in ("ALL","BERW","BETB"):
   part=[r for r in allrows if r["window"]==win and (protocol=="ALL" or r["method"]==protocol)]
   result[protocol]={"total_eligible":counts(part),
       "any_concentration":counts([r for r in part if r["any_concentration"]]),
       "mg_l_any_stat":counts([r for r in part if r["any_mg_l"]]),
       "mg_l_C2_MK_monthly_average":counts([r for r in part if r["monthly_avg_C2_MK_mg_l"]])}
  results[win]=result
 # Need a true delivered river exposure, not just concentration at an outlet.
 # Therefore stop before the 2x2 exposure x canopy interaction.
 result={"schema":"gesui_NRSA_prewindow_TP_method_canopy_source_union_v1",
   "official_NRSA_SHA":{"benthiccount":BENTHIC_SHA,"phab":PHAB_SHA},
   "original_fiscal_archives":{str(r["fiscal"]):r["archive_sha256"] for r in reports},
   "source_population":195,"by_window":results,"qa":diag,
   "biological_outcomes_read":False,"delivered_stream_TP_dose_identified":False,
   "unexposed_controls_identified":False,"exposure_x_canopy_4cell_possible_to_evaluate":False,
   "go_to_wastewater_moderation_regression":False,
   "interpretation":"Source opportunity only. No in-stream dose or no-discharge controls; never sum upstream outlet concentrations. Method strata are separate."}
 z.out.parent.mkdir(parents=True,exist_ok=True)
 z.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
 print("NRSA_CROSS_FY_METHOD_CANOPY_SOURCE_UNION_START")
 print(json.dumps(result,sort_keys=True))
 print("NRSA_CROSS_FY_METHOD_CANOPY_SOURCE_UNION_END")

if __name__=="__main__":main()
