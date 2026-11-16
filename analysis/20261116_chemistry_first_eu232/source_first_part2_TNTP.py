#!/usr/bin/env python3
"""Source-only strict official Part2 chemistry-first EU original-232 H7 support audit.

Read selected original 2026 annual water-quality station keys, then the full
genuine public Part2 (6,450,192 row) original source, zero ecological fields.
Map year 2018 to original early GBIF 2017-19 and 2022 to GBIF 2020-22.
Use exact reported source TN & TP code, river RW, water W, canonical mg/L,
>=3 source samples, LOQ valid, single annual record per nutrient+station+year.
Do not fabricate values for missing stations or for neither-report sites.

A hydrobasin×period nutrient index is the median across qualified ORIGINAL
WQ stations in that watershed and year. It is ambient river chemistry, not
facility-derived WWTP effluent and not on a butterfly host-plant microhabitat.
Only source cohort opportunity is assessed; ecological outcomes NOT OPENED.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,math,sys,zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261116_reported_np_source"))
from acquire_eea_2026_part2_exact_site import (
    download_part2, parse_strict,CODES,EXPECTED_ZIP_BYTES)
from source_first_h7_wq_geography import PROTOCOL,YEARS,MIN_SUPPORT

ZIP_SHA="315396800d76e46db415164851fb4f99765b7ec2654828fc42782d902e679a7a"
EXPECTED_SOURCE_ROWS=6450192

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def load_frozen_geography(folder):
    gate=json.loads((folder/"CHEMISTRY_FIRST_232H7_WQ_GEOGRAPHY_GATE.json").read_text())
    file=folder/"ORIGINAL_2026_EEA_WQ_STATIONS_IN_ORIGINAL_232_H7_SOURCE_ONLY.csv"
    if (gate.get("schema")!=PROTOCOL or
        gate.get("source_candidate_CSV_sha256")!=sha(file) or
        gate.get("source_stage")!="GEOMETRY_ONLY_OUTCOMES_NOT_OPENED"):
        raise RuntimeError("Unchanged original EEA2026 water-geography source freeze missing")
    with file.open(newline="",encoding="utf-8-sig") as f:rows=list(csv.DictReader(f))
    if len(rows)!=gate["full_source_any_country"]["original_WQ_station_year_candidates"]:
        raise RuntimeError("Original WQ source station-year candidate list changed")
    keys={}
    for rec in rows:
        y=int(rec["source_WQ_calendar_year"])
        if y not in YEARS or rec["original_GBIF_source_H7_period"]!=YEARS[y]:
            raise RuntimeError("Original geographic year mapping drift")
        country=rec["source_WQ_country"]
        sid=rec["source_WQ_monitoringSiteIdentifier"]
        scheme=rec["source_WQ_monitoringSiteIdentifierScheme"]
        key=(country,sid,scheme,y)
        if key in keys:raise RuntimeError("Selected original WQ station-year duplicate")
        keys[key]=rec
    return gate,keys


def scan_original(zip_path,keys):
    with zipfile.ZipFile(zip_path) as z:
        names=[x.filename for x in z.infolist() if x.filename.lower().endswith(".csv")]
        if len(names)!=1:raise RuntimeError("Original Part2 archive member is not unique")
        with z.open(names[0]) as orig:
            stream=io.TextIOWrapper(orig,encoding="utf-8-sig",newline="",errors="strict")
            head=stream.readline()
            delimiter=";" if head.count(";")>head.count(",") else ","
            fields=next(csv.reader([head],delimiter=delimiter))
            needed=["countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme",
                    "phenomenonTimeReferenceYear","parameterWaterBodyCategory",
                    "procedureAnalysedMatrix","observedPropertyDeterminandCode",
                    "resultUom","resultNumberOfSamples","resultMeanValue","resultQualityMeanBelowLOQ"]
            if not all(x in fields for x in needed):
                raise RuntimeError("Published Part2 source schema no longer includes required true-water-chemistry fields")
            reader=csv.DictReader(stream,fieldnames=fields,delimiter=delimiter)
            raw_n=0
            status=collections.Counter()
            yy=collections.Counter()
            qualifier=collections.defaultdict(lambda:{"TN":[],"TP":[]})
            for row in reader:
                raw_n+=1
                try:year=int(str(row["phenomenonTimeReferenceYear"]).strip())
                except (ValueError,TypeError):continue
                if year not in YEARS:continue
                yy[str(year)]+=1
                key=(row["countryCode"],row["monitoringSiteIdentifier"],
                     row["monitoringSiteIdentifierScheme"],year)
                if key not in keys:continue
                status["source_selected_station_year_all_properties"]+=1
                if row["parameterWaterBodyCategory"]!="RW" or row["procedureAnalysedMatrix"]!="W":
                    status["invalid_waterbody_or_matrix"]+=1
                    continue
                code=row["observedPropertyDeterminandCode"]
                if code not in CODES.values():continue
                nutrient="TN" if code==CODES["TN"] else "TP"
                status[nutrient+"_exact_code_source"]+=1
                valid,reason=parse_strict(row,nutrient)
                if reason:
                    status[nutrient+"_"+reason]+=1
                    continue
                qualifier[key][nutrient].append(valid)
                status[nutrient+"_qualified_original_annual_records"]+=1
            if raw_n!=EXPECTED_SOURCE_ROWS:
                raise RuntimeError("Original official WQ Part2 full-source row count differs "+str(raw_n))
    return qualifier,{"official_part2_original_csv_rows":raw_n,
        "source_2018_2022_raw_rows":dict(yy),
        "source_validity_diagnostic_counts":dict(status),
        "verified_original_csv_columns":fields}

def counts(rows):
    return {"H7period":len({(r["h"],r["period"]) for r in rows}),
            "distinct_H7":len({r["h"] for r in rows}),
            "true_H6":len({r["h6"] for r in rows}),
            "country_groups":len({r["country"] for r in rows}),
            "qualified_original_WQ_station_years":sum(r["qualified_original_WQ_station_year_n"] for r in rows)}

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    geo,candidates=load_frozen_geography(a.geography)
    tmp=a.output/"TEMP_GENUINE_OFFICIAL_PART2_WQ_SOURCE.zip"
    meta=download_part2(tmp)
    if meta["official_2026_original_part2_ZIP_SHA256"]!=ZIP_SHA:
        raise RuntimeError("Official 2026 Part2 source SHA drift")
    values,rawQA=scan_original(tmp,candidates)
    tmp.unlink()
    status=collections.Counter()
    qualified=[]
    for key,r in sorted(candidates.items()):
        q=values.get(key,{"TN":[],"TP":[]})
        n=len(q["TN"]);p=len(q["TP"])
        if n==1:status["TN_single_valid_site_year"]+=1
        if p==1:status["TP_single_valid_site_year"]+=1
        if n>1 or p>1:status["ambiguous_annual_replicates"]+=1
        if n!=1 or p!=1:continue
        if not math.isfinite(q["TN"][0]["mean_mgL"]) or not math.isfinite(q["TP"][0]["mean_mgL"]):
            raise RuntimeError("Original chemistry validated nonfinite but passed")
        qualified.append({
           "country":r["source_WQ_country"],
           "WQ_monitoringSiteIdentifier":r["source_WQ_monitoringSiteIdentifier"],
           "WQ_monitoringSiteIdentifierScheme":r["source_WQ_monitoringSiteIdentifierScheme"],
           "year":int(r["source_WQ_calendar_year"]),
           "h":r["official_original_HYBAS_ID_L07"],
           "h6":r["original_true_HYBAS_ID_L06"],
           "period":r["original_GBIF_source_H7_period"],
           "WQ_country_matches_original_GBIF_basin_country":
               r["WQ_station_country_equals_original_GBIF_basin_majority_country"].lower()=="true",
           "TN_mgN_L":float(q["TN"][0]["mean_mgL"]),
           "TP_mgP_L":float(q["TP"][0]["mean_mgL"]),
           "TN_original_sample_n":int(q["TN"][0]["samples"]),
           "TP_original_sample_n":int(q["TP"][0]["samples"]),
           "same_H7_NOT_same_site_with_GBIF_or_EEA_EQR":True,
           "no_sewage_source_apportionment":True})
    primary=[r for r in qualified if r["WQ_country_matches_original_GBIF_basin_country"]]
    def to_basins(rows):
        d=collections.defaultdict(list)
        for r in rows:
            d[(r["h"],r["h6"],r["period"],r["country"])].append(r)
        med=[]
        for (h,h6,period,country),group in sorted(d.items()):
            valTN=sorted(x["TN_mgN_L"] for x in group)
            valTP=sorted(x["TP_mgP_L"] for x in group)
            def median(v):
                n=len(v)
                return float(v[n//2]) if n%2 else float((v[n//2-1]+v[n//2])/2)
            med.append({
                "h":h,"h6":h6,"period":period,
                "country":country,"n_source_valid_WQ_stations":len(group),
                "original_2026_part2_strict_TN_mgN_L_site_median":median(valTN),
                "original_2026_part2_strict_TP_mgP_L_site_median":median(valTP),
                "source_cross_basin_exposure_is_not_plant_root_zone_chemistry":True})
        return med
    allbasin=to_basins(qualified)
    primarybasin=to_basins(primary)
    # Cross-border catchments could feature more than one WQ country; the
    # PRIMARY country-matched version must be one row per original H7-period.
    pkey=[(x["h"],x["period"]) for x in primarybasin]
    if len(pkey)!=len(set(pkey)):
        raise RuntimeError("Primary original WQ country-matched H7period duplicates")
    if len({(r["h"],r["period"]) for r in qualified})>454:
        raise RuntimeError("Source valid chemistry outside fixed original 454 H7periods")
    summary={
       "schema":"GESUI_20261116_CHEMISTRY_FIRST_EU232_PART2_SOURCE_GATE_V1",
       "status":"SOURCE_ONLY_CHEMISTRY_COHORT_VERIFIED",
       "verified_source_part2_meta":meta,
       "full_part2_source_scan":rawQA,
       "geometry_source_pre_frozen_sha":geo["source_candidate_CSV_sha256"],
       "original_native_H7":232,"original_H7_period":454,"original_H6":90,
       "original_source_years":[2018,2022],
       "original_WQ_selected_station_years":len(candidates),
       "source_station_year_trait_support":dict(status),
       "full_any_WQ_country_siteyear":counts(allbasin),
       "primary_exact_original_GBIF_country_siteyear":counts(primarybasin),
       "original_strict_qualified_source_siteyear_n":len(qualified),
       "strict_primary_siteyear_n":len(primary),
       "predeclared_support_thresholds":MIN_SUPPORT,
       "geography_eligibility_country_must_match_original_basin":True,
       "all_supported_TN_TP_same_station_year_not_separate_single_nutrients":True,
       "source_no_report_does_NOT_prove_sewage_free":True,
       "chemistry_is_observed_ambient_not_untreated_WWTP_flux":True,
       "no_biotic_outcome_or_ecological_effect_fitted":True,
       "no_original_thirdparty_raw_WQ_full_data_rehosted":True}
    pcounts=summary["primary_exact_original_GBIF_country_siteyear"]
    source_pass=all(pcounts[k]>=v for k,v in MIN_SUPPORT.items())
    summary["source_gate_pass"]=bool(source_pass)
    summary["status"]="PASS_SOURCE_SUPPORT_ONLY_NO_ECOLOGICAL_CLAIM" if source_pass else "STOP_SOURCE_CHEMISTRY_PANEL_SUPPORT"
    # Save only compact source-derived rows in workflow artifact. Do not
    # commit the third-party publisher's full original station-level records.
    out=a.output/"CHEMISTRY_FIRST_EU232_H7_PERIOD_AMBIENT_TNTP_PART2.csv"
    if primarybasin:
        with out.open("w",encoding="utf-8",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(primarybasin[0]));w.writeheader();w.writerows(primarybasin)
        summary["source_derived_main_chemistry_rows_sha256"]=sha(out)
    else:summary["source_derived_main_chemistry_rows_sha256"]=None
    (a.output/"CHEMISTRY_FIRST_EU232_PART2_SOURCE_GATE.json").write_text(
        json.dumps(summary,indent=2,ensure_ascii=False)+"\n")
    print("GESUI_CHEMISTRY_FIRST_EU232_PART2="+json.dumps({
        "status":summary["status"],
        "source_gate_pass":source_pass,
        "original_source_rows":rawQA["official_part2_original_csv_rows"],
        "genuine_source_siteyears":len(qualified),
        "primary_siteyears":len(primary),
        "source_primary_support":pcounts,
        "source_any_country_support":summary["full_any_WQ_country_siteyear"],
        "original_2026_ZIP_SHA256":ZIP_SHA,
        "NO_ECOLOGICAL_MODEL_FITTED":True},ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--geography",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
