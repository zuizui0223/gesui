#!/usr/bin/env python3
"""Source-only EEA 2026 Part2 TNTP in nearby, same-native-H7 WQ stations.

Original source H7 GIS/nearest candidates were frozen BEFORE this chemistry
source retrieval. Unaltered original 2026 official 170,422,309B Part2 annual
CSV ZIP SHA256 fixed; scan every 6,450,192 rows. Exact WQ original ID,
scheme, country, calendar 2018 or 2022, RW/W and canonical TN and TP with
>=3 samples, correct elemental units and quality/LOQ. Never use observed
value's magnitude to select nearby site, only a predeclared quality gate.

Primary: frozen nearest WQ site among candidate IDs within <=1 km regardless
of chemistry availability. Secondary: ANY separately identified station with
valid both at same year, only an upper bound on source availability. Both
remain NOT EEA Biology same-site results, NOT same river reach, and never
identify fraction of measured N from upstream wastewater sources.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,sys,zipfile
from collections import Counter,defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261116_reported_np_source"))
from acquire_eea_2026_part2_exact_site import (
    download_part2,parse_int,parse_strict,CODES
)

PART2_SHA="315396800d76e46db415164851fb4f99765b7ec2654828fc42782d902e679a7a"
GEO_SCHEMA="GESUI_20261116_EEA2026_ORIGINAL_NATIVE_H7_GEOMETRIC_WQ_NEAR_SITE_OPPORTUNITY_V1"
GEO_SHA="b9b810f4db271f0fa6b39fecf01784041b948944b3c3f296a393766388e86b44"
UNION_SCHEMA="GESUI_20261116_EEA2026_PART1_PART2_ORIGINAL_3525_NATIVE_SITE_YEAR_TNTP_UNION_GATE_V1"
N_BIO_SITEYEAR=3525
N_NEAR_CANDIDATES=575
N_WQ_DISTINCT=429
MIN_SOURCE={"H7period":30,"H7":20,"H6":10,"country":3}

def one(base,name):
    paths=list(base.rglob(name))
    if len(paths)!=1:raise RuntimeError("Source only candidate artifact missing/duplicated "+name)
    return paths[0]
def yes(text):
    if str(text).strip().lower() not in ("true","false"):
        raise RuntimeError("Unexpected nonboolean original source "+str(text))
    return str(text).strip().lower()=="true"
def csvrows(path):
    with path.open(newline="",encoding="utf-8-sig") as f:return list(csv.DictReader(f))

def load_candidates(args):
    audit=json.loads(one(args.geo,"EEA2026_OFFICIAL_H7_GEOGRAPHIC_BIOLOGY_NEAR_WQ_SOURCE_GATE.json").read_text())
    if (audit.get("schema")!=GEO_SCHEMA or
        audit.get("source_official_HydroBASINS_Level7",{}).get("sha256")!=GEO_SHA or
        audit.get("source_candidate_original_bio_year_times_near_WQ_ids_1000m")!=N_NEAR_CANDIDATES or
        audit.get("source_verified_same_NATIVE_H7_opportunity_by_distance",{}).get("1000m",{}).get("distinct_official_WQ_station_IDs")!=N_WQ_DISTINCT):
        raise RuntimeError("Original pre-chemistry geographic station selection changed")
    candidates=csvrows(one(args.geo,"EEA2026_ORIGINAL_H7_GEOGRAPHIC_NEAR_WQ_CANDIDATE_SITE_YEAR_KEYS_SOURCE_ONLY.csv"))
    if len(candidates)!=N_NEAR_CANDIDATES:raise RuntimeError("Original 575 near same native H7 source candidates changed")
    events=defaultdict(list)
    target=set()
    for row in candidates:
        if not yes(row["both_original_station_coordinates_in_same_native_H7_polygon"]):
            raise RuntimeError("Source-selected WQ station outside original H7 polygon")
        if row["original_biology_country"]!=row["official_WQ_country"]:
            raise RuntimeError("Country mismatch in source original WQ geo")
        d=float(row["straight_geodesic_distance_m"])
        if not 0<=d<=1000:raise RuntimeError("Source unexpected outside 1km")
        event=(row["original_biology_country"],row["original_biology_station_identifier"],
               row["original_biology_station_identifier_scheme"],
               int(row["original_biology_calendar_year"]))
        wk=(row["official_WQ_country"],row["official_WQ_station_identifier"],
            row["official_WQ_station_identifier_scheme"],event[3])
        if wk[:3]==event[:3]:raise RuntimeError("Original WQ site same as Biology site, this must be alternative")
        events[event].append((d,wk,row))
        target.add(wk)
    if len({x[:3] for x in target})!=N_WQ_DISTINCT:raise RuntimeError("429 WQ source keys changed")
    if len(events)!=402:raise RuntimeError("Original 402 Biology siteyears with alternative same-H7 WQ station changed")
    if sum(len(x) for x in events.values())!=N_NEAR_CANDIDATES:
        raise RuntimeError("Candidate duplicate source count")
    for k,rows in events.items():
        if len({x[1] for x in rows})!=len(rows):
            raise RuntimeError("Duplicate candidate for original Biology site-year "+repr(k))
        rows.sort(key=lambda x:(x[0],x[1]))
    return events,target,audit

def old_biology_reporting(args,events):
    qa=json.loads(one(args.union,"EEA2026_PART1_PART2_ORIGINAL_3525_SITE_YEAR_UNION_SOURCE_GATE.json").read_text())
    if (qa.get("schema")!=UNION_SCHEMA or
        qa.get("genuine_EEA_biology_native_station_years")!=N_BIO_SITEYEAR or
        qa.get("original_station_year_source_positive")!=132 or
        qa.get("strict_report_positive_both_source_support",{}).get("H7period")!=9 or
        qa.get("source_gate_pass") is not False):
        raise RuntimeError("Original exact same-station/year TNTP report-positive STOP has changed")
    rows=csvrows(one(args.union,"EEA2026_ONLY_3525_NATIVE_SITE_YEAR_SOURCE_UNION.csv"))
    if len(rows)!=3525:raise RuntimeError("Original 3525 reporting station-year records changed")
    data={}
    for x in rows:
        k=(x["country"],x["original_monitoringSiteIdentifier"],
           x["original_monitoringSiteIdentifierScheme"],int(x["year"]))
        if k in data:raise RuntimeError("Duplicate original Biology site-year")
        data[k]={"h":str(x["h"]),"h6":str(x["h6"]),"period":str(x["period"]),
                 "country":x["country"],
                 "report_positive":yes(x["original_reported_upstream_positive"])}
    if len(data)!=3525 or sum(x["report_positive"] for x in data.values())!=132:
        raise RuntimeError("Original source positive reporting 132 siteyears drift")
    for k,items in events.items():
        row=items[0][2]
        if k not in data:raise RuntimeError("Original source Biology station-year not in genuine EQR panel")
        base=data[k]
        if (base["h"],base["h6"],base["period"])!=(
            row["original_biology_HYBAS_ID_L07"],
            row["original_true_HYBAS_ID_L06"],row["original_biology_period"]):
            raise RuntimeError("Original EQR site-year HydroBASINS ID or period inconsistent with WQ geometry")
    return data

def scan_part2(zipfile_path,target):
    with zipfile.ZipFile(zipfile_path) as z:
        file=[x for x in z.infolist() if x.filename.lower().endswith(".csv")]
        if len(file)!=1:raise RuntimeError("Not one original official Part2 CSV")
        hits=defaultdict(lambda:{"TN":[],"TP":[]})
        stage=Counter();issue=Counter();yearsource=Counter()
        with z.open(file[0]) as raw:
            text=io.TextIOWrapper(raw,encoding="utf-8-sig",newline="",errors="strict")
            header_line=text.readline()
            delimiter=";" if header_line.count(";")>header_line.count(",") else ","
            header=next(csv.reader([header_line],delimiter=delimiter))
            required=["countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme",
                "phenomenonTimeReferenceYear","parameterWaterBodyCategory",
                "observedPropertyDeterminandCode","procedureAnalysedMatrix",
                "resultUom","resultNumberOfSamples","resultQualityMeanBelowLOQ","resultMeanValue"]
            if any(f not in header for f in required):
                raise RuntimeError("Official Part2 required original source structure drift")
            reader=csv.DictReader(text,fieldnames=header,delimiter=delimiter)
            n=0
            for row in reader:
                n+=1
                yr=parse_int(row.get("phenomenonTimeReferenceYear"))
                if yr not in (2018,2022):continue
                k=(row.get("countryCode"),row.get("monitoringSiteIdentifier"),
                   row.get("monitoringSiteIdentifierScheme"),yr)
                if k not in target:continue
                stage["candidate_WQ_site_year_any_determinand_rows"]+=1
                if row.get("parameterWaterBodyCategory")!="RW" or row.get("procedureAnalysedMatrix")!="W":
                    issue["not_RW_river_or_W_water"]+=1;continue
                if row.get("observedPropertyDeterminandCode")==CODES["TN"]:nu="TN"
                elif row.get("observedPropertyDeterminandCode")==CODES["TP"]:nu="TP"
                else:continue
                stage["matched_source_"+nu+"_record"]+=1
                good,why=parse_strict(row,nu)
                if why:
                    issue[nu+"_"+why]+=1
                    continue
                # Not export concentration; record only eligible publisher row.
                hits[k][nu].append({"valid_original_annual_chemistry_source":True,
                     "original_samples":int(good["samples"])})
        if n!=6450192:
            raise RuntimeError("Official 2026 Part2 unmodified total 6450192 annual source rows changed: "+str(n))
        return hits,{"part2_official_2026_all_annual_rows":n,
                    "eligible_original_target_WQ_keys":len(target),
                    "source_matched_rows_by_stage":dict(stage),
                    "source_rejected_rows_by_quality":dict(issue),
                    "not_2026_original_source_rows_treated_as_missing_chemistry":True}

def support(rows):
    npositive=[x for x in rows if x["report_positive"]]
    pop=[(x["h"],x["period"]) for x in npositive]
    return {"station_years":len(rows),"reported_positive_station_years":len(npositive),
       "all_H7_periods":len({(x["h"],x["period"]) for x in rows}),
       "positive_H7period":len(set(pop)),
       "positive_H7":len({x["h"] for x in npositive}),
       "positive_H6":len({x["h6"] for x in npositive}),
       "positive_country":len({x["country"] for x in npositive}),
       "source_support_if_near_sites_treated_as_proxies_only":bool(
           len(set(pop))>=MIN_SOURCE["H7period"]
           and len({x["h"] for x in npositive})>=MIN_SOURCE["H7"]
           and len({x["h6"] for x in npositive})>=MIN_SOURCE["H6"]
           and len({x["country"] for x in npositive})>=MIN_SOURCE["country"])}

def main(args):
    args.output.mkdir(parents=True,exist_ok=True)
    events,target,geo=load_candidates(args)
    reported=old_biology_reporting(args,events)
    part2=args.output/"TEMP_EEA2026_WISE6_P2_original_AggregatedData.zip"
    receipt=download_part2(part2)
    if receipt.get("official_2026_original_part2_ZIP_SHA256")!=PART2_SHA:
        raise RuntimeError("Original official publisher 170MB Part2 ZIP SHA changed")
    chem,scan=scan_part2(part2,target)
    part2.unlink()
    nearest=[];anyvalid=[];any_excl_alias=[]
    summary_status=Counter()
    unique_candidates_annual=sum(bool(chem.get(k,{}).get("TN")) or bool(chem.get(k,{}).get("TP")) for k in target)
    for key,c in events.items():
        base=reported[key]
        nearestall=c[0]
        def strict(keyw):
            a=chem.get(keyw,{"TN":[],"TP":[]})
            return len(a["TN"])==1 and len(a["TP"])==1
        s={"h":base["h"],"h6":base["h6"],"period":base["period"],
           "country":base["country"],"report_positive":base["report_positive"]}
        if strict(nearestall[1]):
            nearest.append(s)
        if any(strict(wkey) for _,wkey,_ in c):
            anyvalid.append(s)
        if any(strict(wkey) and meters>25 for meters,wkey,_ in c):
            any_excl_alias.append(s)
        if not strict(nearestall[1]) and any(strict(wkey) for _,wkey,_ in c):
            summary_status["alternative_later_than_nearest_wq_station_only"]+=1
        if any(len(chem.get(wkey,{}).get("TN",[]))==1 for _,wkey,_ in c):
            summary_status["any_nearWQ_valid_TN_at_2018or2022"]+=1
        if any(len(chem.get(wkey,{}).get("TP",[]))==1 for _,wkey,_ in c):
            summary_status["any_nearWQ_valid_TP_at_2018or2022"]+=1
        if any(len(chem.get(wkey,{}).get("TN",[]))>1 or len(chem.get(wkey,{}).get("TP",[]))>1 for _,wkey,_ in c):
            summary_status["ambiguous_duplicate_sameWQ_siteyear_nutrient"]+=1
    result={
       "schema":"GESUI_20261116_EEA2026_NEARBY_OFFICIAL_SAME_H7_PART2_TNTP_SOURCE_GATE_V1",
       "source_pre_chemistry_candidate_provenance":"https://github.com/zuizui0223/gesui/actions/runs/43804933915",
       "original_true_H7_geography_source_SHA256":GEO_SHA,
       "official_publisher_Part2_source_receipt":receipt,
       "original_source_Part2_row_scan":scan,
       "biology_original_siteyear_with_at_least_one_nearby_official_WQ_station":len(events),
       "original_distinct_nearby_WQ_station_IDs":len({x[:3] for x in target}),
       "all_candidate_WQ_station_year_keys":len(target),
       "all_Biology_WQ_year_candidate_pairs":N_NEAR_CANDIDATES,
       "Part2_original_2026_candidate_station_years_with_either_valid_TN_or_TP":unique_candidates_annual,
       "source_nearby_site_year_eligibility":{
          "nearest_WQ_ID_frozen_before_chemistry":support(nearest),
          "ANY_nearby_WQ_ID_upperbound_source_only":support(anyvalid),
          "ANY_nearby_WQ_ID_morethan25m_alias_excluded":support(any_excl_alias)},
       "individual_nutrient_and_dup_source_quality_counts":dict(summary_status),
       "same_original_EEA_Biology_native_station_year_TNTP_original_GATE_STOP_preserved":True,
       "selected_new_WQ_different_monitoringSiteIdentifier_or_scheme":True,
       "selected_new_WQ_same_original_true_H7_polygon_proved":True,
       "selected_new_WQ_same_true_HydroRIVERS_segment_NOT_YET_VERIFIED":True,
       "source_TNTP_concentrations_are_AMBIENT_river_not_WWTP_source_fraction":True,
       "source_report_not_matched_DOES_NOT_mean_unexposed":True,
       "source_years_are_2018_2022_not_full_biology_3year_windows":True,
       "original_candidate_site_selection_did_NOT_use_chemical_values":True,
       "no_original_EEA_EQR_or_EIV_or_HNS_outcomes_opened":True,
       "no_ecological_effect_or_mediation_fitted":True}
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/"EEA2026_NEARBY_SAMEH7_PART2_ANNUAL_TNTP_SOURCE_GATE.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_EEA2026_PART2_SAME_H7_NEARBY_TNTP_GATE="+json.dumps(
        {k:v for k,v in result.items() if k not in ("original_source_Part2_row_scan",)},
        ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)
    print("STRICT_SAME_BIOLOGICAL_STATION_YEAR_SOURCE_GATE=STILL_STOP",flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--geo",type=Path,required=True)
    ap.add_argument("--union",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    main(ap.parse_args())
