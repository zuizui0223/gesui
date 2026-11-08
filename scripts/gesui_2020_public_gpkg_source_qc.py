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
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GESUI_2020_SOURCE_ONLY_HISTORICAL_GPKG_ORIGINAL_BYTES_VERIFIED",flush=True)

if __name__=="__main__":main()
