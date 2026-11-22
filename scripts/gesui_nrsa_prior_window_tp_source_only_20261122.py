#!/usr/bin/env python3
"""Source-only prior-window original EPA DMR TP audit. Never fits biological outcomes."""
import argparse,collections,csv,datetime,hashlib,json,math,pathlib,zipfile
PARSER=argparse.ArgumentParser()
for key in ("fiscal","archive","geo","topology","out"): PARSER.add_argument("--"+key,required=True)
def norm(s):
 s=str(s).strip().upper()
 return s.zfill(3) if s.isdigit() and len(s)<=3 else s
def date(s):
 s=str(s).strip()
 if not s:return None
 for fmt in ("%m/%d/%Y","%Y-%m-%d","%m/%d/%Y %H:%M:%S"):
  try:return datetime.datetime.strptime(s[:19],fmt).date()
  except ValueError:pass
 return None
def main():
 a=PARSER.parse_args();fy=int(a.fiscal)
 geo=json.loads(pathlib.Path(a.geo).read_text());top=json.loads(pathlib.Path(a.topology).read_text())
 pairs=top["eras"]["2023_24"]["verified_source_only_upstream_candidate_pairs"]
 sites={r["uid"]:r for r in geo["eras"]["2023_24"]["source_records"]}
 source=collections.defaultdict(set); meta={}
 for p in pairs:
  uid=p["uid"]; d=date(sites[uid]["date_col"])
  if d is None:raise ValueError("Invalid original NRSA collection date")
  source[(p["permit"].strip().upper(),norm(p["feature"]))].add(uid)
  meta[uid]={"date":d,"huc":p["huc8"],"state":p["state"]}
 if len(meta)!=195 or len(source)!=657:raise ValueError("Frozen 195/657 source gate changed")
 # One outfall may have repeated DMR rows and multiple reporting periods.
 events=collections.defaultdict(list);unit_counts=collections.Counter();qc=collections.Counter()
 with zipfile.ZipFile(a.archive) as z:
  names=[n for n in z.namelist() if n.lower().split("/")[-1] in
         (f"npdes_dmrs_fy{fy}.csv",f"npdes_dmr_fy{fy}.csv")]
  if len(names)!=1:raise ValueError("Official archive missing expected DMR CSV")
  with z.open(names[0]) as f:
   import io
   reader=csv.DictReader(io.TextIOWrapper(f,encoding="utf-8-sig",errors="replace",newline=""))
   required={"EXTERNAL_PERMIT_NMBR","PERM_FEATURE_NMBR","PERM_FEATURE_TYPE_CODE","MONITORING_LOCATION_CODE","MONITORING_PERIOD_END_DATE","PARAMETER_CODE","DMR_VALUE_NMBR","NODI_CODE","VALUE_TYPE_CODE","DMR_UNIT_DESC","DMR_VALUE_QUALIFIER_CODE","DMR_VALUE_ID","STATISTICAL_BASE_CODE"}
   if not required.issubset(reader.fieldnames):raise ValueError("DMR original columns changed")
   seen=set()
   for row in reader:
    qc["total_rows"]+=1
    key=(row["EXTERNAL_PERMIT_NMBR"].strip().upper(),norm(row["PERM_FEATURE_NMBR"]))
    if key not in source:continue
    if row["PERM_FEATURE_TYPE_CODE"].strip().upper()!="EXO" or row["MONITORING_LOCATION_CODE"].strip().upper() not in {"1","2","EG"}:continue
    if row["PARAMETER_CODE"].strip()!="00665":continue
    d=date(row["MONITORING_PERIOD_END_DATE"])
    if d is None or d.year not in (2022,2023,2024):continue
    if row["NODI_CODE"].strip():continue
    try:v=float(row["DMR_VALUE_NMBR"])
    except ValueError:continue
    if not math.isfinite(v):continue
    ident=row["DMR_VALUE_ID"].strip()
    if ident:
     if ident in seen:continue
     seen.add(ident)
    typ=row["VALUE_TYPE_CODE"].strip().upper();unit=row["DMR_UNIT_DESC"].strip().upper()
    qualifier=row["DMR_VALUE_QUALIFIER_CODE"].strip().upper()
    stat=row["STATISTICAL_BASE_CODE"].strip().upper()
    qc["exact_outfall_tp_numeric_rows"]+=1
    unit_counts[(typ,unit,stat,qualifier)]+=1
    if not typ.startswith("C"):continue
    if qualifier not in ("","="):qc["qualified_excluded"]+=1;continue
    # Deliberately do NOT collapse mg/L, ug/L, daily max and monthly avg.
    for uid in source[key]:
     lag=(meta[uid]["date"]-d).days
     if lag>=0 and lag<=365:
      events[(uid,key)].append((lag,v,unit,stat,typ))
 counts={}
 for days in (30,90,365):
  qualifying=collections.defaultdict(list)
  for (uid,key),rows in events.items():
   for lag,v,unit,stat,typ in rows:
    if lag<=days:qualifying[uid].append((key,unit,stat,typ,v))
  uids=set(qualifying)
  comparable={uid for uid,rows in qualifying.items() if any(unit in {"MG/L","MG/L AS P","MILLIGRAMS/LITER","MILLIGRAMS PER LITER"} and v>=0 for key,unit,stat,typ,v in rows)}
  counts[str(days)]={"any_concentration_report_site_visits":len(uids),
    "distinct_huc8":len({meta[u]["huc"] for u in uids}),
    "distinct_states":len({meta[u]["state"] for u in uids}),
    "mg_per_l_candidate_site_visits":len(comparable),
    "mg_per_l_candidate_huc8":len({meta[u]["huc"] for u in comparable}),
    "note":"mg/L opportunity only; statistical bases and outfall aggregation not yet harmonized"}
 result={"schema":"gesui_prior_window_tp_source_only_v1","fiscal":fy,
  "archive_sha256":hashlib.sha256(pathlib.Path(a.archive).read_bytes()).hexdigest(),
  "original_candidate_visits":195,"original_upstream_outfalls":657,
  "counts":counts,"original_unit_stat_qualifier_profiles":[{"value_type":t,"unit":u,"statistical_base":s,"qualifier":q,"n":n} for (t,u,s,q),n in sorted(unit_counts.items())],
  "qc":dict(qc),"no_biological_outcomes":True,"no_delivered_dose_or_controls":True,
  "warning":"Fiscal archives overlap calendar-year visits; merge UID×DMR_VALUE_ID across fiscal years before final site counts"}
 pathlib.Path(a.out).write_text(json.dumps(result,indent=2,sort_keys=True))
 print(json.dumps({"fiscal":fy,"counts":counts,"qc":dict(qc)},sort_keys=True))
if __name__=="__main__":main()
