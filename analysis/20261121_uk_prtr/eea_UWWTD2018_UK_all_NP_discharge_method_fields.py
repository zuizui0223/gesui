#!/usr/bin/env python3
"""UWWTD2018 authentic UK5 original N/P discharge METHOD field completeness.

Existing source-locked EEA 2018 UK5 has 1799 plants and 1809 discharge
points; measured uwwNDischargeMeasured and uwwPDischargeMeasured >0 have
no UK data. This does NOT mean other CALCULATED or ESTIMATED reported
P/N discharge fields are empty. This strictly SOURCE-only audits EVERY
original UWWTPS.csv column concerning source N/P discharge, preserving
the original reporter, field definitions and nonzero/empty counts.

Do NOT equate calculated/estimated reporter mass to independent measured
outlet P/N kg/y or guaranteed downstream delivery; report method/source.
No source facility rows or original mass quantities are republished.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,collections,math,urllib.request,zipfile,tempfile
from pathlib import Path
import eea_UWWTD2018_UK_outlet_source_gate as source

REPORTER="UK5"
EXPECTED=1799
MAX_COLS=400
def numeric(v):
    try:
        n=float(str(v or "").strip().replace(",",""))
        return n if math.isfinite(n) else None
    except (ValueError,TypeError):return None

def make(args):
    out={"schema":"GESUI_20261121_EEA2018_UWWTD_UK_ALL_NP_REPORT_METHOD_SOURCE_CENSUS_V1",
         "source_only":True,"status":"SOURCE_STOP",
         "no_facility_discharge_mass_or_physical_outlet_rows_published":True,
         "no_BMWP_or_ChemPop_ecological_model":True,
         "not_direct_source_resolved_PN_delivered_to_river":True}
    try:
        with tempfile.TemporaryDirectory(prefix="gesui_EEA2018_UK_PN_") as tmp:
            p=Path(tmp)/"eea.zip"
            req=urllib.request.Request(source.URL,headers={"User-Agent":"gesui-EEA2018-UK-source-complete-discharge-NP-field-audit"})
            sha=hashlib.sha256();n=0
            with urllib.request.urlopen(req,timeout=90) as resp,p.open("wb") as dest:
                if resp.status!=200:raise RuntimeError("Original EEA 2018 ZIP not 200")
                while True:
                    x=resp.read(131072)
                    if not x:break
                    n+=len(x)
                    if n>source.MAX_BYTES:raise RuntimeError("Original archive exceeds size bound")
                    dest.write(x);sha.update(x)
            with zipfile.ZipFile(p) as z:
                name=source.BASE+"UWWTPS.csv"
                raw=z.read(name)
                if hashlib.sha256(raw).hexdigest()!=source.EXPECT["UWWTPS.csv"]:
                    raise RuntimeError("Original immutable EEA UWWTPS scientific source SHA drift")
            reader=csv.DictReader(io.StringIO(raw.decode("utf-8-sig","replace")))
            cols=list(reader.fieldnames or [])
            if not cols or len(cols)>MAX_COLS:raise RuntimeError("Unexpected original historical UWWTPS columns")
            rows=[r for r in reader if str(r.get("repCode","")).strip()==REPORTER]
            if len(rows)!=EXPECTED:raise RuntimeError("Original UK5 2018 UWWTPS source rows changed")
            colnames=[c for c in cols if ("discharge" in c.lower() or
                     c.lower().startswith("uwwn") or c.lower().startswith("uwwp") or
                     "phosph" in c.lower() or "nitrogen" in c.lower())]
            report={}
            for col in colnames:
                v=[str(r.get(col,"") or "").strip() for r in rows]
                present=[q for q in v if q and q.lower() not in ("na","null","nan")]
                positive=[numeric(q) for q in present]
                pos=sum(x is not None and x>0 for x in positive)
                zero=sum(x is not None and x==0 for x in positive)
                str_labels=collections.Counter(q for q in present if numeric(q) is None)
                report[col]={
                    "present_original_UK5_field_rows":len(present),
                    "distinct_original_UK5_field_values":len(set(present)),
                    "numeric_positive_original_UK5_rows":pos,
                    "numeric_explicit_zero_original_UK5_rows":zero,
                    "numeric_negative_original_UK5_rows":sum(x is not None and x<0 for x in positive),
                    "nonnumeric_or_flag_original_rows":sum(x is None for x in positive),
                    "nonnumeric_labels_top8":str_labels.most_common(8) if len(str_labels)<100 else [],
                    "this_is_column_source_availability_NOT_mass_units_verified":True,
                }
            out.update({
                "status":"EEA_UWWTD2018_UK_PN_ALL_ORIGINAL_REPORT_SOURCE_COLUMNS_AUDITED",
                "original_zip_SHA256":sha.hexdigest(),
                "original_2018_UWWTPS_uncompressed_SHA256":source.EXPECT["UWWTPS.csv"],
                "original_2018_UK5_plant_rows":len(rows),
                "original_UWWTPS_csv_column_total":len(cols),
                "discharge_nutrient_relevant_original_column_names":colnames,
                "source_field_coverage_by_original_name":report,
                "no_source_quantity_or_facility_ID_or_original_locations_emitted":True
            })
    except Exception as e:
        out.update({"status":"EEA_UWWTD2018_UK_SOURCE_FIELD_GATE_STOP",
                    "source_error":type(e).__name__+": "+str(e)[:600]})
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,sort_keys=True,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_EEA2018_UK_NP_OTHER_METHOD_SOURCE_CENSUS="+json.dumps({
       "status":out["status"],"error":out.get("source_error"),
       "original_column_names":out.get("discharge_nutrient_relevant_original_column_names"),
       "original_field_source_coverage":out.get("source_field_coverage_by_original_name")
    },sort_keys=True,ensure_ascii=False,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    make(p.parse_args())
