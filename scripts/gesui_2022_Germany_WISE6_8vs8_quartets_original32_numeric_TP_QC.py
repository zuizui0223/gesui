#!/usr/bin/env python3
"""Original EEA2026 WISE6 river Total-P 2022 numeric quality at exact 32
source-frozen Germany monitoring stations forming 8×(source-vs-zero) quartets.

Site selection and 2° matching were committed BEFORE any new numeric WISE6
TP values at these 13 newly opened sites. No site replacement or threshold
optimization, no individual butterflies or actual treatment effect.
"""
from __future__ import annotations
import argparse,hashlib,json,math,statistics,urllib.request,zipfile
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

ZIP_URL="https://discodata.eea.europa.eu/download/WISE_SOE/latest/Waterbase_T_WISE6_DisaggregatedData"
ORIG_ZIP_SHA="2a2ef2a43a16ec068b7974a80ae5fca49e44697660a46ce3b91fb1813abe62d2"
ORIG_ZIP_SIZE=2024284753
ORIG_MEMBER="waterbase_t_wise6_disaggregateddata.csv"
ORIG_CSV_SIZE=36043700357
ORIG_ROWS=102569057
ORIG_SOURCE_QUARTETS_BLOB="8bd67be25968311ca16be4026888daa783fbe69f"
TP_CODE="CAS_7723-14-0"
TP_UNIT="mg{P}/L"
DATE_MIN=4
MONTH_MIN=3
MAX_BELOW_LOQ=.5
MIN_COMMON_PAIR_MONTHS=2
MIN_COMPLETE_QUARTETS=6
MIN_BLOCKS=4
RECORDS_LIMIT=150000000

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()

def blob(path):
    x=path.read_bytes()
    return hashlib.sha1(b"blob "+str(len(x)).encode()+b"\0"+x).hexdigest()

def safe_str(b):return b.strip().strip(b'"').decode("utf8",errors="replace")

def original_download(path):
    req=urllib.request.Request(ZIP_URL,headers={
      "User-Agent":"gesui-2022-WISE6-nationwide-Germany-8plus8-quartets-32-station-originalTP/1.0",
      "Accept":"application/zip,*/*"})
    h=hashlib.sha256();size=0
    with urllib.request.urlopen(req,timeout=180) as r,path.open("wb") as out:
        if r.status!=200:raise ValueError("Original official WISE6 ZIP HTTP not 200")
        for part in iter(lambda:r.read(1048576),b""):
            size+=len(part)
            if size>ORIG_ZIP_SIZE:raise ValueError("Original WISE6 ZIP changed size")
            out.write(part);h.update(part)
    if size!=ORIG_ZIP_SIZE or h.hexdigest()!=ORIG_ZIP_SHA:
        raise ValueError("Original WISE6 original numeric ZIP SHA changed")
    return {"original_URL":ZIP_URL,"compressed_size":size,"original_SHA256":h.hexdigest()}

def original_site_set(j):
    if j["predeclared_8_quartet_sourceonly_quality_GATE"] is not True or j["original_WISE6_2022_total_TP_numeric_sample_values_OPENED"]:
        raise ValueError("Original 8+8 source-only quartet gate was not outcome-blind")
    rows=j["selected_original_8_source_vs_clean_zero_monitor_quartets"]
    if len(rows)!=8:raise ValueError("Original 8 selected matched quartets changed")
    keys=[]
    for r in rows:
        for arm in ("source","registered_source_zero"):
            for end in ("upstream","downstream"):
                x=r[arm]
                if x["countryCode"]!="DE":raise ValueError("Non-Germany site in source-only quartet")
                keys.append(x[end+"_station_SHA256"])
    if len(set(keys))!=32:raise ValueError("A preselected original WISE6 station reused across different quartets")
    return rows,set(keys)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-frozen-quartets",required=True,type=Path)
    ap.add_argument("--tmpdir",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args()
    if blob(a.source_frozen_quartets)!=ORIG_SOURCE_QUARTETS_BLOB:
        raise ValueError("Frozen country+2degree 8:8 original quartet SHA changed")
    original=json.loads(a.source_frozen_quartets.read_text())
    quartets,sites=original_site_set(original)
    obs={k:{"original_2022_sample_rows":0,"below_LOQ":0,"blank_LOQ_flag":0,
            "invalid_flag":0,"missing_result":0,"invalid_numeric":0,
            "valid_by_date":defaultdict(list)} for k in sites}
    a.tmpdir.mkdir(parents=True,exist_ok=True)
    archive=a.tmpdir/"EEA_2026_original_WISE6_2022_32_original_stations.zip"
    receipt=original_download(archive)
    qc=Counter()
    with zipfile.ZipFile(archive) as z:
        if z.namelist()!=[ORIG_MEMBER] or z.getinfo(ORIG_MEMBER).file_size!=ORIG_CSV_SIZE:
            raise ValueError("Original public 36GB WISE6 CSV ZIP contents changed")
        with z.open(ORIG_MEMBER) as inp:
            hdr=inp.readline().decode("utf-8-sig").strip().split(";")
            expected=["countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme",
                      "parameterWaterBodyCategory","observedPropertyDeterminandCode",
                      "observedPropertyDeterminandLabel","procedureAnalysedMatrix","resultUom",
                      "phenomenonTimeSamplingDate","phenomenonTimeReferenceYear","sampleIdentifier",
                      "resultObservedValue","resultQualityObservedValueBelowLOQ"]
            if hdr[:13]!=expected:raise ValueError("Original EEA WISE6 numeric columns changed")
            for line in inp:
                qc["original_public_WISE6_2026_individual_rows_streamed"]+=1
                if qc["original_public_WISE6_2026_individual_rows_streamed"]>RECORDS_LIMIT:
                    raise ValueError("Unexpected original whole source row count")
                if not(line.startswith(b"DE;") or line.startswith(b'"DE";')):
                    continue
                x=line.split(b";",14)
                if len(x)<14:continue
                if safe_str(x[9])!="2022" or not safe_str(x[8]).startswith("2022-"):
                    continue
                if safe_str(x[3])!="RW" or safe_str(x[5])!="Total phosphorus" or safe_str(x[6])!="W":
                    continue
                if safe_str(x[4])!=TP_CODE or safe_str(x[7])!=TP_UNIT:
                    raise ValueError("Official Germany 2022 water TP code or unit source drift")
                qc["original_DE_2022_river_TP_water_samples_all_sites"]+=1
                h=hashlib.sha256(("DE|"+safe_str(x[2])+"|"+safe_str(x[1])).encode()).hexdigest()
                if h not in sites:continue
                datum=obs[h];datum["original_2022_sample_rows"]+=1
                flag=safe_str(x[12]).lower()
                if flag in ("true","1"):
                    datum["below_LOQ"]+=1;continue
                if flag not in ("false","0",""):
                    datum["invalid_flag"]+=1;continue
                if not flag:datum["blank_LOQ_flag"]+=1
                sval=safe_str(x[11])
                if not sval:datum["missing_result"]+=1;continue
                try:v=float(sval)
                except (ValueError,TypeError):
                    datum["invalid_numeric"]+=1;continue
                if not math.isfinite(v) or v<0:
                    datum["invalid_numeric"]+=1;continue
                datum["valid_by_date"][safe_str(x[8])[:10]].append(v)
    if qc["original_public_WISE6_2026_individual_rows_streamed"]!=ORIG_ROWS:
        raise ValueError("Original 102569057 measured water source row count changed")
    source_qual={}
    for k,v in obs.items():
        dates=sorted(v["valid_by_date"])
        months=sorted(set(d[:7] for d in dates))
        fr=v["below_LOQ"]/v["original_2022_sample_rows"] if v["original_2022_sample_rows"] else 1.
        good=len(dates)>=DATE_MIN and len(months)>=MONTH_MIN and fr<=MAX_BELOW_LOQ
        daily=[statistics.median(v["valid_by_date"][d]) for d in dates]
        source_qual[k]={
          "original_country_scheme_site_SHA256":k,
          "original_2022_all_sample_rows":v["original_2022_sample_rows"],
          "original_2022_qualified_numeric_days":len(dates),
          "original_2022_qualified_numeric_months":months,
          "original_2022_belowLOQ_sample_rows":v["below_LOQ"],
          "original_2022_belowLOQ_fraction":fr,
          "original_2022_unknown_blank_LOQflag_n":v["blank_LOQ_flag"],
          "original_2022_invalid_flag_n":v["invalid_flag"],
          "original_2022_numeric_missing_n":v["missing_result"],
          "original_2022_invalid_numeric_n":v["invalid_numeric"],
          "original_2022_annual_siteTP_median_of_daily_medians_mgP_L":
              statistics.median(daily) if good else None,
          "original_site_numeric_4day_3month_LOQ_gate":good,
        }
    source_groups=[]
    for i,q in enumerate(quartets,1):
        arms={};validquartet=True; reasons=[]
        for arm in ("source","registered_source_zero"):
            s=q[arm]
            up=source_qual[s["upstream_station_SHA256"]]
            down=source_qual[s["downstream_station_SHA256"]]
            months=set(up["original_2022_qualified_numeric_months"])&set(down["original_2022_qualified_numeric_months"])
            g=bool(up["original_site_numeric_4day_3month_LOQ_gate"] and
                   down["original_site_numeric_4day_3month_LOQ_gate"] and len(months)>=MIN_COMMON_PAIR_MONTHS)
            if not g:
                validquartet=False
                reasons.append(arm+("_site_quantQC_or_overlap_months_below2"))
            arms[arm]={
              "source_only_original_intervening_WWTP_P_t_year":
                  s["real_exact_source_between"]["strictly_between_original_P_reported_t_year"],
              "source_original_upstream_station_hash":s["upstream_station_SHA256"],
              "source_original_downstream_station_hash":s["downstream_station_SHA256"],
              "source_full_reach_path_km":s["full_reach_path_to_downstream_reach_start_km"],
              "common_original_2022_quantified_months":sorted(months),
              "common_original_quantified_months_n":len(months),
              "original_source_2station_pair_numeric_QC_and_2month_overlap":g,
            }
        source_groups.append({
          "original_source_quartet_index":i,
          "original_source_UPSTREAM_country2deg_block":q["source"]["upstream_country2degree_block"],
          "source_only_quartet_geometry_matching_cost":q["value"],
          "all_original_four_station_TP_numeric_source_eligible":validquartet,
          "reason_numeric_or_months_unavailable":reasons,
          "original_source_and_control_2station_arms":arms,
          "original_source_2022_TP_differences_NOT_COMPUTED_YET":True,
        })
    included=[q for q in source_groups if q["all_original_four_station_TP_numeric_source_eligible"]]
    blocks=set(q["original_source_UPSTREAM_country2deg_block"] for q in included)
    gate=len(included)>=MIN_COMPLETE_QUARTETS and len(blocks)>=MIN_BLOCKS
    result={
       "schema":"gesui_DE_2022_WISE6_original_8WWTPsource_vs8zero_matched_river_monitor_quartet_32site_numeric_QC_v01",
       "UTC":datetime.now(timezone.utc).isoformat(),
       "POSTHOC_after_Lycaenidae_2022_photo_and_2022_reach_P_source_results":True,
       "original_source_frozen_8X8_quartets_Git_blob_SHA1":ORIG_SOURCE_QUARTETS_BLOB,
       "original_source_only_quartet_file_SHA256":sha(a.source_frozen_quartets),
       "original_WISE6_full_public_2026_source_zip":receipt,
       "exact_original_32_monitoring_sites_preselected_without_TP_numeric_data":32,
       "actual_original_2022_WISE6_DE_riverTP_source_row_counts":dict(qc),
       "source_original_32_monitoring_sites_2022_numeric_TP_QC":list(source_qual.values()),
       "original_source_8matched_quartet_numeric_completion":source_groups,
       "source_quartets_valid_4station_4days_3months_2commonmonths_LOQ":len(included),
       "source_valid_quartets_original_country2degree_source_blocks":len(blocks),
       "original_2022_quartet_numeric_data_eligibility_GATE":gate,
       "original_8selected_matched_source_versus_control_2station_quartet_still_unmodified":True,
       "original_2022_river_TP_original_units":"mg{P}/L",
       "no_butterfly_photos_or_species_counts_used":True,
       "no_2022_TP_matched_gradients_computed_in_this_stage":True,
       "no_causal_source_effect_claim":True,
       "original_2022_site_values_accessed_AFTER_original_8:8_source_only_quartet_freeze":True,
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,sort_keys=True,ensure_ascii=False,indent=2)+"\n")
    print("GESUI_GERMANY_2022_WISE6_ORIGINAL_8SOURCE8CONTROL_TP_NUMERIC_QUALITY",
      json.dumps({"32sites":len(source_qual),
        "qualified_siteyears":sum(x["original_site_numeric_4day_3month_LOQ_gate"] for x in source_qual.values()),
        "complete_sourcezero_comparator_quartets":len(included),
        "spatial_blocks":len(blocks),"gate":gate,
        "numeric_LOQ_below_rows":sum(x["original_2022_belowLOQ_sample_rows"] for x in source_qual.values()),
        "no_quartet_differences_computed":True},sort_keys=True),flush=True)
if __name__=="__main__":main()
