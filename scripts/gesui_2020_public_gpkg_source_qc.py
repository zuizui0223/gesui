#!/usr/bin/env python3
"""Hash-locked EEA 2019–2020 public GPKG WWTP/outfall source-only QA.

Only the official public UWWTD plant/outfall 2019–2020 GPKGs are read.
Never open EEA biological EQR or assume the report period = an individual
2019/2020 measured sampling year. Original datasets are not rehosted.
"""
from __future__ import annotations
import argparse,base64,hashlib,json,math,sqlite3,urllib.request
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import quote

TOKEN="sptXqwkQr5g7Bp5"
ROOT="https://sdi.eea.europa.eu/datashare/public.php/webdav/"
SOURCES={
 "WWTP":{
  "relative":"eea_v_4258_100_m_uwwtd-tp_p_2019-2020_v08_r00/GPKG/UWWTD_TreatmentPlants.gpkg",
  "expected_bytes":7159808},
 "DCP":{
  "relative":"eea_v_4258_100_m_uwwtd-disc-pts_p_2019-2020_v08_r00/GPKG/UWWTD_DischargePoints.gpkg",
  "expected_bytes":7073792},
}

def download(name,target):
    src=SOURCES[name]
    u=ROOT+"/".join(quote(x,safe="") for x in src["relative"].split("/"))
    auth="Basic "+base64.b64encode((TOKEN+":").encode()).decode()
    req=urllib.request.Request(u,headers={
        "Authorization":auth,"User-Agent":"gesui-EEA2020-WWTP-GPKG-provenance"})
    h=hashlib.sha256();size=0
    with urllib.request.urlopen(req,timeout=100) as r,target.open("wb") as f:
        if r.status !=200:raise ValueError("Official public WebDAV download not OK")
        for b in iter(lambda:r.read(1024*1024),b""):
            f.write(b);h.update(b);size+=len(b)
    if size!=src["expected_bytes"]:
        raise ValueError(f"EEA historical GPKG byte length changed: {name} {size}")
    with target.open("rb") as f:
        if f.read(16)!=b"SQLite format 3\x00":
            raise ValueError("EEA public GPKG not SQLite GeoPackage")
    return {"public_original_path":u,"source_SHA256":h.hexdigest(),
            "download_bytes":size,"public_reporting_period":"2019-2020",
            "public_version":"v08_r00"}

def quote_field(field):
    return '"'+str(field).replace('"','""')+'"'

def inspect(path):
    db=sqlite3.connect(f"file:{path}?mode=ro",uri=True)
    db.row_factory=sqlite3.Row
    try:
        tables=[str(z["name"]) for z in db.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'")]
        feature=[str(z["table_name"]) for z in db.execute(
            "SELECT table_name FROM gpkg_contents WHERE data_type='features'")]
        allinfo={}
        for tab in feature:
            cols=[{"name":x[1],"type":x[2]} for x in db.execute(
                "PRAGMA table_info("+quote_field(tab)+")")]
            n=db.execute("SELECT COUNT(*) FROM "+quote_field(tab)).fetchone()[0]
            allinfo[tab]={"features_n":n,"column_info":cols}
        return db,tables,feature,allinfo
    except:
        db.close();raise

def firstcol(cols,*choices):
    lookup={c.lower():c for c in cols}
    for x in choices:
        if x.lower() in lookup:return lookup[x.lower()]
    return None

def numeric(x):
    try:
        if x is None or str(x).strip()=="":return None
        z=float(x)
        if math.isfinite(z):return z
    except (TypeError,ValueError):pass
    return None

def extract(db,tname,cols,kind):
    rep=firstcol(cols,"repCode")
    key=firstcol(cols,"uwwCode")
    p=firstcol(cols,"uwwPDischargeMeasured")
    n=firstcol(cols,"uwwNDischargeMeasured")
    lat=firstcol(cols,"dcpLatitude","uwwLatitude")
    lon=firstcol(cols,"dcpLongitude","uwwLongitude")
    water=firstcol(cols,"dcpWaterBodyType")
    extra=[c for c in cols if any(v in c.lower() for v in
      ("report","period","measured","discharge","waterbody","beginlife","endoflife","state","coordinate","latitude","longitude"))][:80]
    required=[rep,key]
    if not all(required):
        return {"error":"Official source lacks exact repCode+uwwCode join",
                "key_candidates":required,"candidate_data_fields":extra},None
    select=list(dict.fromkeys([c for c in [rep,key,p,n,lat,lon,water] if c]))
    query="SELECT "+",".join(map(quote_field,select))+" FROM "+quote_field(tname)
    raw=[]
    for row in db.execute(query):
        r=dict(zip(select,row))
        skey=str(r[rep]).strip()+"|"+str(r[key]).strip() if r[rep] is not None and r[key] is not None else None
        if skey and not skey.startswith("None") and not skey.endswith("|"):
            r["source_key"]=skey
        else:r["source_key"]=None
        raw.append(r)
    cnt=Counter(x["source_key"] for x in raw if x["source_key"] is not None)
    unique=[x for x in raw if x["source_key"] and cnt[x["source_key"]]==1]
    out={"features_n":len(raw),"distinct_valid_full_keys":len(cnt),
         "unique_joinable_full_keys":len(unique),
         "duplicate_full_key_groups":sum(v>1 for v in cnt.values()),
         "joining_fields":[rep,key],"model_year_reference_field_candidates":extra,
         "observed_measured_P_field":p,"observed_measured_N_field":n,
         "latitude_field":lat,"longitude_field":lon,"waterbody_field":water}
    rep_values=[str(x[rep]).strip() for x in unique if x[rep] is not None]
    out["official_reporter_code_year_suffix_counts"]={str(k):int(v) for k,v in
          Counter((v[-4:] if v[-4:].isdigit() else "OTHER") for v in rep_values).most_common(15)}
    out["official_reporter_code_year_suffix_2020_n"]=sum(v.endswith("2020") for v in rep_values)
    if p and n:
        ppos=sum(numeric(x[p]) is not None and numeric(x[p])>0 for x in unique)
        npos=sum(numeric(x[n]) is not None and numeric(x[n])>0 for x in unique)
        both=sum(numeric(x[p]) is not None and numeric(x[p])>0 and
                 numeric(x[n]) is not None and numeric(x[n])>0 for x in unique)
        out.update({"unique_plant_positive_reported_P":ppos,
                    "unique_plant_positive_reported_N":npos,
                    "unique_plant_both_positive_reported_P_N":both,
                    "reported_values_not_yet_verified_same_calendar_year":True})
    if water:
        out["waterbody_type_counts"]={str(k):int(v) for k,v in
              Counter(str(x[water]) for x in unique).most_common(30)}
    if lat and lon:
        bbox=sum(numeric(x[lat]) is not None and numeric(x[lon]) is not None and
                 44<=numeric(x[lat])<=64 and -12<=numeric(x[lon])<=26 for x in unique)
        out["unique_source_coords_in_original_Europe_study_bbox"]=bbox
    return out,unique

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--temporary-dir",type=Path,required=True);args=ap.parse_args()
    args.temporary_dir.mkdir(parents=True,exist_ok=True)
    report={"schema":"gesui_EEA_public_2019_2020_GPKG_reported_NP_source_only_v01",
            "created_utc":datetime.now(timezone.utc).isoformat(),
            "no_2020_biological_EQR_opened":True,
            "same_calendar_2020_reported_mass_year_not_yet_verified":True,
            "datasets":{}}
    rows_by_kind={}
    for k in SOURCES:
        path=args.temporary_dir/(k+".gpkg")
        provenance=download(k,path)
        db,tables,feature,info=inspect(path)
        data={"source":provenance,"tables":tables,"feature_layers":info}
        if len(feature)!=1:
            data["error"]="Not exactly one source feature layer; do not infer join"
        else:
            tab=feature[0]
            cols=[x["name"] for x in info[tab]["column_info"]]
            qc,rows=extract(db,tab,cols,k)
            data["quality"]=qc
            rows_by_kind[k]=rows
            if rows is not None and k=="WWTP" and qc.get("observed_measured_P_field") and qc.get("observed_measured_N_field"):
                possible=[x for x in rows
                          if numeric(x[qc["observed_measured_P_field"]]) is not None and
                             numeric(x[qc["observed_measured_P_field"]])>0 and
                             numeric(x[qc["observed_measured_N_field"]]) is not None and
                             numeric(x[qc["observed_measured_N_field"]])>0]
                report["WWTP_measured_positive_P_and_N_unique_source_keys_count"]=len(possible)
                # Do not persist 2020 mass without explicit original calendar-year and unit checks.
        db.close()
        report["datasets"][k]=data
        print("GESUI_2020_ORIGINAL_SOURCE_GPKG_QC",json.dumps({
            "dataset":k,"SHA256":provenance["source_SHA256"],
            "file_size":provenance["download_bytes"],"layers":feature,
            "source_quality":data.get("quality",{})},sort_keys=True),flush=True)
    # Source-only 1:1 plant and discharge point linkage, exact original IDs.
    # The public spatial product includes active/correct-coordinate source rows;
    # not all treatment sources or receiving paths are included.
    if rows_by_kind.get("WWTP") is not None and rows_by_kind.get("DCP") is not None:
        plants=rows_by_kind["WWTP"];outlets=rows_by_kind["DCP"]
        seen={x["source_key"]:x for x in outlets}
        matched=[]
        for p in plants:
            pv=numeric(p.get("uwwPDischargeMeasured"))
            nv=numeric(p.get("uwwNDischargeMeasured"))
            d=seen.get(p["source_key"])
            if not(pv is not None and pv>0 and nv is not None and nv>0 and d):
                continue
            # This spatial product mixes reporter cycles! To measure an
            # actual 2020 reporter cohort, never silently include the
            # records whose official repCode suffix is 2022/other years.
            if not str(p.get("repCode","")).strip().endswith("_2020"):
                continue
            if not str(d.get("repCode","")).strip().endswith("_2020"):
                continue
            if str(d.get("dcpWaterBodyType","")).strip().upper()!="FW":
                continue
            lat=numeric(d.get("dcpLatitude"));lon=numeric(d.get("dcpLongitude"))
            if lat is None or lon is None or not(44<=lat<=64 and -12<=lon<=26):
                continue
            matched.append({
                "source_key":p["source_key"],
                "repCode":p["repCode"],"uwwCode":p["uwwCode"],
                "uwwNDischargeMeasured_t_per_year":nv,
                "uwwPDischargeMeasured_t_per_year":pv,
                "dcpLatitude":lat,"dcpLongitude":lon,
                "dcpWaterBodyType":"FW",
            })
        import csv
        source_csv=args.out.parent/"2020_PUBLIC_CYCLE_BOTH_MEASURED_NP_UNIQUE_FW_OUTFALLS.csv"
        source_csv.parent.mkdir(parents=True,exist_ok=True)
        cols=["source_key","repCode","uwwCode","uwwNDischargeMeasured_t_per_year",
              "uwwPDischargeMeasured_t_per_year","dcpLatitude","dcpLongitude","dcpWaterBodyType"]
        with source_csv.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(matched)
        report["official_both_measured_positive_NP_unique_FW_outlets_in_study_bbox"]=len(matched)
        report["official_2020_reporter_suffix_enforced_in_candidate_CSV"]=True
        report["official_2020_reporter_code_yyyy_ref_not_individual_measurement_date"]=True
        report["source_candidate_file_sha256"]=hashlib.sha256(source_csv.read_bytes()).hexdigest()
        report["source_candidate_file_path"]=source_csv.name
        report["source_candidate_original_2020_reporting_period_not_precise_measurement_date"]=True
        print("GESUI_2020_PUBLIC_BOTH_NP_MEASURED_FRESHWATER_OUTFALL_SOURCE_ROWS",
              json.dumps({"eligible_count":len(matched),"original_year_suffix_counts":
                dict(Counter(str(x["repCode"])[-4:] for x in matched).most_common(15)),
                "source_file_SHA256":report["source_candidate_file_sha256"]}),flush=True)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GESUI_2020_SOURCE_ONLY_HISTORICAL_GPKG_ORIGINAL_BYTES_VERIFIED",flush=True)

if __name__=="__main__":main()
