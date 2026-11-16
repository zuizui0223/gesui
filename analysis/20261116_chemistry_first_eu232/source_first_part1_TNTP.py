#!/usr/bin/env python3
"""Full EEA2026 original Part1: chemistry-first genuine EU232 H7 source gate.

Post-Part2 STOP separately pre-frozen original WQ station-year frame, not a
benthic-station selection. Stream official 96,597,294 individual raw sample
rows from 1.639GB DEFLATE64 ZIP and compute strict same-site/year N+P means.
Never fit a butterfly/plant/EQR model or claim WWTP effluent source delivery.
No original source full ZIP/CSV is rehosted.
"""
from __future__ import annotations

import argparse,collections,csv,hashlib,io,json,math,subprocess,sys,zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261116_chemistry_first_eu232"))
from source_first_part2_TNTP import load_frozen_geography,MIN_SUPPORT
sys.path.insert(0,str(ROOT/"analysis/20261116_reported_np_source"))
from acquire_eea_2026_part1_exact_site import download_source,CODES,UNITS

ORIGINAL_ZIP_SHA="ff779521749b40b952cce35a4671240c3b6c47d87f32a8fef9240a4ac503a082"
ORIGINAL_CSV_ROWCOUNT=96597294
SCHEMA="GESUI_20261116_CHEMISTRY_FIRST_PART1_WQ_EU232_STRICT_SITE_YEAR_TNTP_V1"

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def mean_src(obj):
    if len(obj)<3:return None
    return float(sum(obj.values())/len(obj))

def full_scan(zip_path,target):
    by=collections.defaultdict(lambda:{"TN":{},"TP":{}})
    conflict=set()
    count=collections.Counter()
    years=collections.Counter()
    with zipfile.ZipFile(zip_path) as z:
        valid=[f.filename for f in z.infolist() if f.filename.lower().endswith(".csv")]
        if len(valid)!=1:
            raise RuntimeError("Original official source Part1 CSV member not unique")
        member=valid[0]
        method=z.getinfo(member).compress_type
        if method!=9:
            raise RuntimeError("Expected original EEA Part1 DEFLATE64 method9 source changed")
        proc=subprocess.Popen(["7z","x","-so",str(zip_path),member],
                              stdout=subprocess.PIPE,stderr=subprocess.PIPE,bufsize=1024*1024)
        total=0
        try:
            with proc.stdout:
                f=io.TextIOWrapper(proc.stdout,encoding="utf-8-sig",newline="",errors="strict")
                reader=csv.reader(f)
                header=next(reader)
                required=["countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme",
                 "parameterWaterBodyCategory","observedPropertyDeterminandCode",
                 "procedureAnalysedMatrix","resultUom","phenomenonTimeSamplingDate",
                 "sampleIdentifier","resultObservedValue","resultQualityObservedValueBelowLOQ"]
                if not set(required).issubset(header):
                    raise RuntimeError("Original 2026 Part1 source schema changed")
                i={name:header.index(name) for name in required}
                for row in reader:
                    total+=1
                    if total%10000000==0:
                        print("GESUI_CHEMISTRY_FIRST_PART1_ROWS_SCANNED="+str(total),flush=True)
                    date=row[i["phenomenonTimeSamplingDate"]].strip()
                    if len(date)<4 or not date[:4].isdigit():continue
                    year=int(date[:4])
                    if year not in (2018,2022):continue
                    years[str(year)]+=1
                    key=(row[i["countryCode"]],
                         row[i["monitoringSiteIdentifier"]],
                         row[i["monitoringSiteIdentifierScheme"]],year)
                    if key not in target:continue
                    count["original_site_year_source_any_parameter"]+=1
                    if row[i["parameterWaterBodyCategory"]]!="RW" or row[i["procedureAnalysedMatrix"]]!="W":
                        count["not_river_water"]+=1;continue
                    code=row[i["observedPropertyDeterminandCode"]]
                    if code not in CODES.values():continue
                    nutrient="TN" if code==CODES["TN"] else "TP"
                    count[nutrient+"_original_site_year_measured_rows"]+=1
                    if row[i["resultUom"]]!=UNITS[nutrient]:
                        count[nutrient+"_wrong_unit"]+=1;continue
                    if row[i["resultQualityObservedValueBelowLOQ"]].strip().lower() in ("true","1"):
                        count[nutrient+"_below_LOQ"]+=1;continue
                    try:value=float(row[i["resultObservedValue"]])
                    except (TypeError,ValueError):
                        count[nutrient+"_missing_numeric"]+=1;continue
                    if not math.isfinite(value) or value<0:
                        count[nutrient+"_invalid_numeric"]+=1;continue
                    sid=row[i["sampleIdentifier"]].strip()
                    sample=("id",sid) if sid else ("source_date_value",date,format(value,".12g"))
                    previous=by[key][nutrient].get(sample)
                    if previous is not None:
                        if abs(previous-value)>1e-10:
                            count[nutrient+"_conflicting_duplicate_source_sample"]+=1
                            conflict.add((key,nutrient))
                        else:
                            count[nutrient+"_exact_duplicate_source_sample"]+=1
                        continue
                    by[key][nutrient][sample]=value
                    count[nutrient+"_valid_distinct_samples"]+=1
            rc=proc.wait(timeout=45)
            if rc!=0:
                msg=proc.stderr.read(5000).decode("utf-8","replace") if proc.stderr else ""
                raise RuntimeError("EEA original DEFLATE64 sample streaming failed rc="+str(rc)+" "+msg)
        finally:
            if proc.poll() is None:
                proc.kill()
                proc.wait(timeout=15)
    if total!=ORIGINAL_CSV_ROWCOUNT:
        raise RuntimeError("Official original WQ Part1 full-source row-count mismatch "+str(total))
    return by,conflict,{"official_original_part1_rows":total,
             "part1_rows_2018_and_2022":dict(years),
             "source_selected_and_quality_counts":dict(count),
             "original_CSV_required_fields_verified":required,
             "original_zip_compression_method":method,
             "conflicting_siteyear_nutrient_pairs_excluded":len(conflict)}

def summary(rows):
    h7per={(r["h"],r["period"]) for r in rows}
    return {"H7period":len(h7per),
            "distinct_H7":len({r["h"] for r in rows}),
            "true_H6":len({r["h6"] for r in rows}),
            "country_groups":len({r["country"] for r in rows}),
            "strict_valid_original_WQ_siteyears":len(rows)}

def to_basins(rows):
    groups=collections.defaultdict(list)
    for r in rows:
        groups[(r["h"],r["h6"],r["period"],r["country"])].append(r)
    ans=[]
    for (h,h6,period,country),g in sorted(groups.items()):
        def med(field):
            vals=sorted(r[field] for r in g)
            n=len(vals)
            return float(vals[n//2]) if n%2 else float((vals[n//2-1]+vals[n//2])/2)
        ans.append({"h":h,"h6":h6,"period":period,"country":country,
                    "n_source_valid_WQ_stations":len(g),
                    "original_part1_strict_TN_mgN_L_site_median":med("TN_mgN_L"),
                    "original_part1_strict_TP_mgP_L_site_median":med("TP_mgP_L"),
                    "ambient_TNTP_not_sewage_derived_allocation":True})
    return ans

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    frozen,targets=load_frozen_geography(a.geography)
    if frozen.get("source_candidate_CSV_sha256")!="c5cde7bd8134b254c5e95b06849292fa260fd6f01832e6b6d2e1994a903c5910":
        raise RuntimeError("Water-first station-years source freeze changed since Part2-only gate")
    tmp=a.output/"OFFICIAL_2026_PART1_TEMP_DO_NOT_REHOST.zip"
    receipt=download_source(tmp)
    if receipt["original_part1_SHA256"]!=ORIGINAL_ZIP_SHA:
        raise RuntimeError("Original official WISE6 Part1 ZIP SHA mismatch")
    by,conflict,qa=full_scan(tmp,targets)
    tmp.unlink()
    count=collections.Counter()
    allvalid=[]
    for key,orig in sorted(targets.items()):
        vals=by.get(key,{"TN":{},"TP":{}})
        n=vals["TN"];p=vals["TP"]
        tn_ok=len(n)>=3 and (key,"TN") not in conflict
        tp_ok=len(p)>=3 and (key,"TP") not in conflict
        if tn_ok:count["part1_strict_TN_siteyear"]+=1
        if tp_ok:count["part1_strict_TP_siteyear"]+=1
        if not(tn_ok and tp_ok):continue
        count["part1_same_siteyear_TNTP"]+=1
        allvalid.append({
            "country":orig["source_WQ_country"],
            "monitoringSiteIdentifier":orig["source_WQ_monitoringSiteIdentifier"],
            "monitoringSiteIdentifierScheme":orig["source_WQ_monitoringSiteIdentifierScheme"],
            "year":key[3],
            "h":orig["official_original_HYBAS_ID_L07"],
            "h6":orig["original_true_HYBAS_ID_L06"],
            "period":orig["original_GBIF_source_H7_period"],
            "WQ_country_matches_original_GBIF_basin_country":
               orig["WQ_station_country_equals_original_GBIF_basin_majority_country"].lower()=="true",
            "TN_mgN_L":mean_src(n),
            "TP_mgP_L":mean_src(p),
            "valid_TN_samples":len(n),"valid_TP_samples":len(p)})
    primary=[r for r in allvalid if r["WQ_country_matches_original_GBIF_basin_country"]]
    primaryrows=to_basins(primary)
    if len(primaryrows)!=summary(primary)["H7period"]:
        raise RuntimeError("Original WQ primary source H7period repeated across countries")
    if primaryrows:
        out=a.output/"CHEMISTRY_FIRST_EU232_PART1_H7_PERIOD_AMBIENT_TNTP.csv"
        with out.open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=list(primaryrows[0]))
            w.writeheader();w.writerows(primaryrows)
        psha=sha(out)
    else:psha=None
    source={
       "schema":SCHEMA,
       "official_original_EEA2026_Part1_ZIP_SHA256":ORIGINAL_ZIP_SHA,
       "original_download_source":receipt,
       "full_CSV_source_and_sample_QA":qa,
       "frozen_H7_WQ_station_year_CSV_SHA256":frozen["source_candidate_CSV_sha256"],
       "original_selected_WQ_station_year_source_candidates":len(targets),
       "source_siteyear_trait_quality_counts":dict(count),
       "full_any_country_qualified":summary(allvalid),
       "primary_original_GBIF_basin_country_matched_qualified":summary(primary),
       "original_published_2026_Part2_only_prior_GATE":"STOP_SOURCE_CHEMISTRY_PANEL_SUPPORT",
       "precommitted_identical_source_minimum":MIN_SUPPORT,
       "source_derived_plant_butterfly_correspondence_not_opened":True,
       "strict_3_distinct_source_samples_per_nutrient_and_station_year":True,
       "source_original_full_WQ_not_rehosted":True,
       "original_2018_or_2022_WQ_NOT_same_time_as_GBIF_3year_records":True,
       "no_sewage_source_apportionment_or_EEA_EQR_site_identity":True,
       "part1_primary_H7_period_CSV_SHA256":psha}
    ok=all(source["primary_original_GBIF_basin_country_matched_qualified"][key]>=value
           for key,value in MIN_SUPPORT.items())
    source["source_support_gate_pass"]=bool(ok)
    source["status"]="PASS_PART1_SOURCE_SUPPORT_ONLY_NO_ECOLOGICAL_CLAIM" if ok else "STOP_PART1_SOURCE_CHEMISTRY_SUPPORT"
    (a.output/"CHEMISTRY_FIRST_EU232_PART1_SOURCE_GATE.json").write_text(json.dumps(source,indent=2,ensure_ascii=False)+"\n")
    print("GESUI_CHEMISTRY_FIRST_EU232_PART1="+json.dumps({
       "status":source["status"],"source_gate_pass":bool(ok),
       "original_Part1_row_count":qa["official_original_part1_rows"],
       "source_siteyears":len(allvalid),"primary_siteyears":len(primary),
       "source_primary":source["primary_original_GBIF_basin_country_matched_qualified"],
       "source_all_countries":source["full_any_country_qualified"],
       "source_ZIP_SHA256":ORIGINAL_ZIP_SHA,
       "NO_ECOLOGICAL_MODEL_FITTED":True},ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--geography",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
