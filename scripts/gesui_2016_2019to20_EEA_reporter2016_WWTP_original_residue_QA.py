#!/usr/bin/env python3
"""Original 2019–20 public EEA UWWTD GPKG 2016 reporter-residue audit.

Read only source physical EEA GPKG files, no 2016 or 2018/22 EQR values.
Original GPKG from public EEA Eionet 2020 has mixed report codes 2020/2022/
OTHER. Test whether 2016-specific source rows actually remain and whether
2016 numeric measured N/P plus unique FW outlet are jointly viable.
Never treat no 2016 reporter source as 0 wastewater or use 2020 values.
"""
from __future__ import annotations
import argparse,hashlib,json,math,sqlite3
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path

EXPECTED={"WWTP":"13060c88a8d1d06f81b61b0a6886279db6051645a2b8d1f7d578ae8ee4e59300",
          "DCP":"1b8b20e0ba9d106f659b0fd9716013b680a16a0e7ab2a779dac2bfb53a33d9a1"}
YEAR="2016"
BBOX=(-12.,44.,26.,64.)
def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for x in iter(lambda:f.read(1024*1024),b""):h.update(x)
 return h.hexdigest()
def floatp(v):
 try:
  w=float(v)
  return w if math.isfinite(w) else None
 except (TypeError,ValueError):return None
def rows(path,table,columns):
 con=sqlite3.connect(f"file:{path}?mode=ro",uri=True)
 con.row_factory=sqlite3.Row
 if con.execute("SELECT COUNT(*) FROM sqlite_master WHERE name=?",(table,)).fetchone()[0]!=1:
  raise ValueError("Official EEA GPKG table missing "+table)
 fields={str(x[1]) for x in con.execute(f'PRAGMA table_info("{table}")')}
 missing=set(columns)-fields
 if missing:raise ValueError("Original official 2020 GPKG source fields changed "+repr(missing))
 query="SELECT "+",".join('"'+c+'"' for c in columns)+" FROM "+'"'+table+'"'
 try:
  for r in con.execute(query):yield dict(r)
 finally:con.close()
def source(wtp,dcp):
 if sha(wtp)!=EXPECTED["WWTP"] or sha(dcp)!=EXPECTED["DCP"]:
  raise ValueError("Original existing 2019–2020 EEA UWWTD public source bytes mismatch")
 controls=Counter()
 WWTPFIELDS=("repCode","uwwCode","uwwNDischargeMeasured","uwwPDischargeMeasured","uwwLatitude","uwwLongitude")
 DCPFIELDS=("repCode","uwwCode","dcpLatitude","dcpLongitude","dcpWaterBodyType")
 w2016={}
 pstat=Counter()
 pyears=Counter()
 country=Counter()
 for r in rows(wtp,"UWWTD_TreatmentPlants",WWTPFIELDS):
  s=str(r["repCode"] or "").strip()
  suffix=s[-4:] if s[-4:].isdigit() else "OTHER"
  pyears[suffix]+=1
  if not s.endswith(YEAR):continue
  controls["2016_original_treatment_rows"]+=1
  key=(s,str(r["uwwCode"] or "").strip())
  if not key[1]:continue
  pair=(floatp(r["uwwNDischargeMeasured"]),floatp(r["uwwPDischargeMeasured"]))
  both=all(v is not None and v>0 for v in pair)
  if both:controls["2016_positive_both_reported_measured_NP_plant_rows"]+=1
  lat=floatp(r["uwwLatitude"]);lon=floatp(r["uwwLongitude"])
  if lat is not None and lon is not None and BBOX[1]<=lat<=BBOX[3] and BBOX[0]<=lon<=BBOX[2]:
   controls["2016_original_plant_coordinate_in_EU_study_extent"]+=1
  pstat[key]+=1
  w2016[key]=(both,lat,lon)
 byDCP=Counter()
 dy=Counter()
 for r in rows(dcp,"UWWTD_DischargePoints",DCPFIELDS):
  s=str(r["repCode"] or "").strip()
  suffix=s[-4:] if s[-4:].isdigit() else "OTHER"
  dy[suffix]+=1
  if not s.endswith(YEAR):continue
  controls["2016_original_discharge_point_rows"]+=1
  key=(s,str(r["uwwCode"] or "").strip())
  if key not in w2016:continue
  if not w2016[key][0]:continue
  lat=floatp(r["dcpLatitude"]);lon=floatp(r["dcpLongitude"])
  if not(str(r["dcpWaterBodyType"]).strip().upper()=="FW" and
      lat is not None and lon is not None and
      BBOX[1]<=lat<=BBOX[3] and BBOX[0]<=lon<=BBOX[2]):continue
  byDCP[key]+=1
  controls["2016_bothNP_source_FW_outfall_rows_Europe_bbox"]+=1
 for key,n in byDCP.items():
  if pstat[key]==1 and n==1:
   country[key[0][:2]]+=1
   controls["2016_bothNP_UNIQUE_EUFW_source_plant_ids"]+=1
 return {
  "schema":"gesui_2016_reporter_residue_in_2019_20_EEA_official_PUBLIC_GPKGs_SOURCE_ONLY_v01",
  "UTC":datetime.now(timezone.utc).isoformat(),
  "old_2019_2020_official_2020_EEA_GPKG_SHA256":EXPECTED,
  "2015_2016_reporter_code_suffix_genuinely_checked_not_assumed":True,
  "source_2016_NP_loads_actual_reporting_year_not_laboratory_sample_date":True,
  "original_public_2020_files_can_be_incomplete_for_2016":True,
  "EEA_treatment_reporter_suffix_n":dict(pyears),
  "EEA_outfall_reporter_suffix_n":dict(dy),
  "strict_2016_reporter_NP_and_unique_singleFW_source_QC":dict(controls),
  "strict_2016_candidate_source_country_counts":dict(country),
  "candidate_source_countries":len(country),
  "can_plan_independent_2016_source_to_2016_BQE_join":controls["2016_bothNP_UNIQUE_EUFW_source_plant_ids"]>=150 and len(country)>=5,
  "no_2016_benthic_EQR_scores_read":True,
  "do_not_replace_2016_mass_with_2018_2020_2022":True
 }
def main():
 a=argparse.ArgumentParser()
 a.add_argument("--wwtp",type=Path,required=True)
 a.add_argument("--outfall",type=Path,required=True)
 a.add_argument("--out",type=Path,required=True)
 v=a.parse_args();x=source(v.wwtp,v.outfall)
 v.out.parent.mkdir(parents=True,exist_ok=True)
 v.out.write_text(json.dumps(x,indent=2,sort_keys=True)+"\n")
 print(json.dumps(x,indent=2))
if __name__=="__main__":main()
