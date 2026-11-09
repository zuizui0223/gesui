#!/usr/bin/env python3
"""Conditional 2023 WISE6 real TP matched four-arm same-month gradient.

Runs ONLY after original fixed 2022 eight WWTP/control four-arm quartets'
2023 *sampling-date-only* source gate is Git committed and PASSED.
Streams original EEA2026 2.024GB/36GB ZIP, numeric TP only at the exact
source-frozen 32 original monitoring site hashes. Does NOT assume official
2022 reported WWTP mass is a measured 2023 source load. Four sites per
quartet and equal quartet weighting; not a causal pollution estimate.
"""
from __future__ import annotations
import argparse,hashlib,itertools,json,math,statistics,urllib.request,zipfile
from collections import defaultdict,Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np

URL="https://discodata.eea.europa.eu/download/WISE_SOE/latest/Waterbase_T_WISE6_DisaggregatedData"
ZIP_SHA="2a2ef2a43a16ec068b7974a80ae5fca49e44697660a46ce3b91fb1813abe62d2"
ZIP_BYTES=2024284753
CSV_BYTES=36043700357
ORIGINAL_ROWS=102569057
YEAR="2023"
TP_CODE="CAS_7723-14-0"
TP_UNIT="mg{P}/L"
ORIG32=32
MIN_DAYS=4
MIN_MONTHS=3
MAX_BELOW_LOQ_FRACTION=0.50
MIN_SHARED_MONTHS=3
MIN_QUARTETS=5
MIN_BLOCKS=4
MIN_SHARED_EXACT_DAYS=3
MIN_SAMEDAY_QUARTETS=4
MIN_SAMEDAY_BLOCKS=3
LOG_OFFSET=0.001
N_BOOT=9999
RANDOM_SEED=20261109

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda:f.read(1048576),b""):h.update(chunk)
    return h.hexdigest()
def val(buf):return buf.strip().strip(b'"').decode("utf8",errors="replace")
def sitehash(c,scheme,site):
    return hashlib.sha256((c+"|"+scheme+"|"+site).encode()).hexdigest()
def download(path):
    req=urllib.request.Request(URL,headers={
        "User-Agent":"gesui-Germany-2023-original-quarter-WISE6-TP-after-date-source-gate/1.0",
        "Accept":"application/zip,*/*"})
    size=0;h=hashlib.sha256()
    with urllib.request.urlopen(req,timeout=200) as r,path.open("wb") as f:
        if r.status!=200:raise ValueError("Official EEA2026 original sample source unavailable")
        for b in iter(lambda:r.read(1048576),b""):
            size+=len(b)
            if size>ZIP_BYTES:raise ValueError("Original WISE6 source expanded beyond original ZIP")
            h.update(b);f.write(b)
    if size!=ZIP_BYTES or h.hexdigest()!=ZIP_SHA:
        raise ValueError("Original EEA2026 WISE6 individual ZIP SHA changed")
    return {"original_2026_public_WISE6_original_URL":URL,
            "verified_source_SHA256":h.hexdigest(),
            "original_public_compressed_zip_bytes":size}

def scan_numeric(path,ids):
    out={h:{"raw_rows":0,"below_LOQ":0,"invalid_flag":0,
            "status_rejected":0,"missing_numeric":0,"invalid_numeric":0,
            "missing_LOQflag":0,"status_original":Counter(),
            "daily":defaultdict(list)} for h in ids}
    qc=Counter()
    with zipfile.ZipFile(path) as z:
        if z.namelist()!=["waterbase_t_wise6_disaggregateddata.csv"] or z.infolist()[0].file_size!=CSV_BYTES:
            raise ValueError("Original WISE6 36GB individual data source member changed")
        with z.open(z.infolist()[0]) as f:
            header=f.readline().decode("utf-8-sig").rstrip("\r\n").split(";")
            expected=("countryCode monitoringSiteIdentifier monitoringSiteIdentifierScheme "
                     "parameterWaterBodyCategory observedPropertyDeterminandCode "
                     "observedPropertyDeterminandLabel procedureAnalysedMatrix resultUom "
                     "phenomenonTimeSamplingDate phenomenonTimeReferenceYear sampleIdentifier "
                     "resultObservedValue resultQualityObservedValueBelowLOQ").split()
            if header[:13]!=expected or header[21]!="resultObservationStatus":
                raise ValueError("Original EEA2026 WISE6 2023 water-TP numeric columns changed")
            for line in f:
                qc["original_source_all_WISE6_individual_rows"]+=1
                if not line.startswith((b"DE;",b'"DE";')):continue
                preview=line.split(b";",11)
                if len(preview)<12:continue
                if (val(preview[9])!=YEAR or val(preview[3])!="RW" or
                    val(preview[4])!=TP_CODE or val(preview[5])!="Total phosphorus" or
                    val(preview[6])!="W" or val(preview[7])!=TP_UNIT):
                    continue
                qc["original_2023_DE_RW_TP_sample_rows_all_sites"]+=1
                h=sitehash("DE",val(preview[2]),val(preview[1]))
                if h not in out:continue
                q=out[h];q["raw_rows"]+=1
                qc["original_2023_TP_rows_at_fixed_32_2022_sites"]+=1
                d=val(preview[8])[:10]
                fields=line.split(b";",22)
                if len(fields)<23:
                    q["invalid_numeric"]+=1;continue
                status=val(fields[21])
                q["status_original"][status]+=1
                if status not in ("","A"):
                    q["status_rejected"]+=1
                    continue
                loq=val(fields[12]).casefold()
                if loq in ("true","1"):
                    q["below_LOQ"]+=1;continue
                if loq not in ("false","0",""):
                    q["invalid_flag"]+=1;continue
                if loq=="":q["missing_LOQflag"]+=1
                v=val(fields[11])
                try:x=float(v)
                except (TypeError,ValueError):
                    q["missing_numeric"]+=1;continue
                if not math.isfinite(x) or x<0 or len(d)!=10 or d[:4]!=YEAR:
                    q["invalid_numeric"]+=1;continue
                q["daily"][d].append(x)
    if qc["original_source_all_WISE6_individual_rows"]!=ORIGINAL_ROWS:
        raise ValueError("Original WISE6 all-row denominator changed")
    return out,dict(qc)

def contrast(arms,site_daily,months_or_days,year):
    # May be called only where all four original frozen arms have sample data.
    if len(months_or_days)==0:return []
    def by_month(day_data):
        groups=defaultdict(list)
        for day,number in day_data.items():
            groups[day[:7]].append(number)
        return {k:statistics.median(v) for k,v in groups.items()}
    byarm={role:({date:site_daily[h][date] for date in site_daily[h]}
                 if len(months_or_days[0])==10 else by_month(site_daily[h]))
           for role,h in arms.items()}
    out=[]
    for m in months_or_days:
        nums={role:byarm[role].get(m) for role in arms}
        if any(x is None for x in nums.values()):
            raise ValueError("Numeric four-arm contrast lacks a shared date/month "+m)
        def lg(d,u):return math.log(nums[d]+LOG_OFFSET)-math.log(nums[u]+LOG_OFFSET)
        src=lg("source_downstream","source_upstream")
        ctrl=lg("registered_zero_downstream","registered_zero_upstream")
        out.append({"original_2023_sampling_period":m,
            "source_upstream_mgP_L":nums["source_upstream"],
            "source_downstream_mgP_L":nums["source_downstream"],
            "registered_zero_upstream_mgP_L":nums["registered_zero_upstream"],
            "registered_zero_downstream_mgP_L":nums["registered_zero_downstream"],
            "source_pair_log_gradient":src,
            "registered_zero_pair_log_gradient":ctrl,
            "source_minus_reference_logTP_gradient":src-ctrl,
            "not_same_day_when_YYYY_MM":len(m)==7,
            "reported_2022_WWTP_annualP_is_not_2023_dose":True})
    return out

def stats(quartets):
    q=np.asarray([x["source_minus_registeredzero_four_arm_logTP_gradient_mean"] for x in quartets],float)
    n=len(q)
    if n==0:raise ValueError("No eligible four-arm quartets")
    rng=np.random.default_rng(RANDOM_SEED)
    boots=np.mean(rng.choice(q,size=(N_BOOT,n),replace=True),axis=1)
    actual=float(np.mean(q))
    sign_flip=np.asarray([np.mean(q*np.asarray(bits,float)) for bits
                in itertools.product((-1,1),repeat=n)],float)
    loo=[float(np.mean(np.delete(q,i))) for i in range(n)] if n>1 else []
    return {
       "source_original_independent_2022_selected_2023_retested_quartets":n,
       "same_country_original_2degree_blocks":len({x["source_2022_original_country2degree_block"] for x in quartets}),
       "quartet_equal_mean_source_minus_sourcezero_logTP_gradient":actual,
       "quartet_equal_median_source_minus_sourcezero_logTP_gradient":float(np.median(q)),
       "quartet_positive_n":int(np.sum(q>0)),"quartet_negative_n":int(np.sum(q<0)),
       "bootstrap_descriptive_9999_quartet_equal95pct_interval":[float(v) for v in np.quantile(boots,[.025,.975])],
       "exact_quartet_signflip_reference_two_sided_fraction":float(np.mean(np.abs(sign_flip)>=abs(actual)-1e-13)),
       "quartet_leave_one_out_mean_range":[min(loo),max(loo)] if loo else None,
       "not_randomized_wastewater_treatment_pvalue":True,
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--original-2023-date-source",required=True,type=Path)
    ap.add_argument("--original-2022-samemonth-numeric",required=True,type=Path)
    ap.add_argument("--tempdir",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args()
    source=json.loads(a.original_2023_date_source.read_text())
    if (source["schema"]!="gesui_2023_original_WISE6_DE_TP_at_2022_frozen_WWTP_source_control_quartets_date_only_gate_v01" or
        not source["never_parsed_2023_resultObservedValue_LOQ_or_WWTP_chemical_load"] or
        source["original_2022_source_sites_n"]!=32):
        raise ValueError("Not precommitted original source-only 2023 dates and 32 site source gate")
    # STOP without downloading or decoding chemistry if original fixed coverage insufficient.
    if not source["strict_2023_frozen_4arm_same_MONTH_numeric_gate_PASS"]:
        report={"status":"STOP_PRENUMERIC_2023_MONTH_SOURCE_GATE_FAILED",
            "precommitted_source_SHA256":sha(a.original_2023_date_source),
            "2023_numeric_sample_rows_read":0,
            "source_calendar_quartets_ge3_months":source["strict_2023_same_month_quartets_n"]}
        a.out.parent.mkdir(parents=True,exist_ok=True)
        a.out.write_text(json.dumps(report,indent=2)+"\n")
        print("GESUI_2023_32STATION_PRENUMERIC_SAMPLE_DATE_GATE_STOP",flush=True)
        return
    prev=json.loads(a.original_2022_samemonth_numeric.read_text())
    if (not prev["source_4site_numeric_same_month_geographic_gate_pass"] or
        prev["source_4site_numeric_valid_ge3_month_quartets_n"]!=7):
        raise ValueError("Original 2022 true WISE TP same-month 7-quartet response not immutable")
    expected={h:val for h,val in source["source_only_original32_site_2023_sampling_date_metadata"].items()}
    if len(expected)!=32:raise ValueError("Unexpected 2023 32station source count")
    a.tempdir.mkdir(parents=True,exist_ok=True)
    zipp=a.tempdir/"2023_EEA2026_Original_WISE6_numeric.zip"
    original=download(zipp)
    measured,sourceQC=scan_numeric(zipp,set(expected))
    bysite={};siteResults=[]
    for h,meta in expected.items():
        d=measured[h]
        if d["raw_rows"]!=meta["original_2023_total_TP_sample_date_metadata_rows"]:
            raise ValueError("2023 numeric raw source row count differs from pre-value site-year metadata "+h)
        if not set(d["daily"]).issubset(set(meta["original_2023_distinct_qualified_source_days"])):
            raise ValueError("2023 numeric daily dates not in original source-only date list")
        quantified_day={date:statistics.median(vals) for date,vals in d["daily"].items()}
        months={s[:7] for s in quantified_day}
        loq_fraction=(d["below_LOQ"]/d["raw_rows"] if d["raw_rows"] else 1.)
        good=(len(quantified_day)>=MIN_DAYS and len(months)>=MIN_MONTHS and
              loq_fraction<=MAX_BELOW_LOQ_FRACTION)
        bysite[h]=quantified_day if good else {}
        siteResults.append({
         "source_original_2022_station_hash":h,
         "original_source_TP_rows_2023":d["raw_rows"],
         "quantified_independent_days_2023":len(quantified_day),
         "quantified_independent_months_2023":len(months),
         "source_2023_below_LOQ_rows":d["below_LOQ"],
         "source_2023_below_LOQ_fraction":loq_fraction,
         "source_2023_other_status_rejected":d["status_rejected"],
         "original_2023_resultObservationStatus_codes":dict(d["status_original"]),
         "original_unknown_LOQflag_rows":d["missing_LOQflag"],
         "original_2023_numeric_site_QC":good,
         "original_2023_source_site_year_median_mgP_L":
               statistics.median(quantified_day.values()) if good else None,
        })
    quarts=[];goodmonthly=[];gooddaily=[]
    for q in source["original8_unreplaced_quartets_2023_fourarm_sample_time_coverage"]:
        arms=q["original_WISE_source_only_fixed_4site_hashes"]
        if len(arms)!=4:raise ValueError("Original four-arm design corrupted")
        dates=[set(bysite[h]) for h in arms.values()]
        common_days=sorted(set.intersection(*dates))
        common_months=sorted(set.intersection(*[{d[:7] for d in days} for days in dates]))
        okmonth=len(common_months)>=MIN_SHARED_MONTHS
        okday=len(common_days)>=MIN_SHARED_EXACT_DAYS
        if okmonth and not q["source_2023_four_arm_ge3_same_months_gate"]:
            raise ValueError("2023 numeric same-month support unexpectedly exceeds original date-only support")
        periods=contrast(arms,bysite,common_months,YEAR) if okmonth else []
        perday=contrast(arms,bysite,common_days,YEAR) if okday else []
        row={
          "source_original_2022_quartet_index":q["original_2022_quartet_index"],
          "source_2022_original_country2degree_block":q["original_2022_same_river_country2degree_source_block"],
          "original_2022_EEA_reported_P_interposed_tpy":q["original_2022_source_interposed_reported_P_tonnes_per_year"],
          "original_station_hashes":arms,
          "2023_fourarm_shared_numeric_months":common_months,
          "2023_fourarm_shared_numeric_days":common_days,
          "2023_numeric_fourarm_month_gate":okmonth,
          "2023_numeric_fourarm_sameday_gate":okday,
          "2023_samemonth_numeric_source_vs_zero_logTP_gradients":periods,
          "2023_sameday_numeric_source_vs_zero_logTP_gradients":perday,
          "source_minus_registeredzero_four_arm_logTP_gradient_mean":
              float(np.mean([v["source_minus_reference_logTP_gradient"] for v in periods])) if periods else None,
        }
        quarts.append(row)
        if okmonth:goodmonthly.append(row)
        if okday:gooddaily.append(row)
    qblocks={r["source_2022_original_country2degree_block"] for r in goodmonthly}
    gate=len(goodmonthly)>=MIN_QUARTETS and len(qblocks)>=MIN_BLOCKS
    dblocks={r["source_2022_original_country2degree_block"] for r in gooddaily}
    dgate=(source["strict_2023_frozen_4arm_same_DAY_numeric_gate_PASS"] and
           len(gooddaily)>=MIN_SAMEDAY_QUARTETS and len(dblocks)>=MIN_SAMEDAY_BLOCKS)
    result={
      "schema":"gesui_WISE6_2023_numeric_TP_original2022_WWTP_Germany_8fixed_fourarms_same_month_replication_v01",
      "original_WISE6_source_download":original,
      "original_2023_sampling_date_source_QC_SHA256":sha(a.original_2023_date_source),
      "original_2022_source_year_matched_TP_file_SHA256":sha(a.original_2022_samemonth_numeric),
      "original_2023_WISE6_full_source_numeric_QC":sourceQC,
      "original_2023_32site_numeric_QC":siteResults,
      "all_original8_2023_numeric_fourarm_quartet_source_comparisons":quarts,
      "2023_numeric_same_month_min5_quartets_min4_blocks_GATE":gate,
      "2023_numeric_same_day_min4_quartets_min3_blocks_GATE":dgate,
      "2023_numeric_same_month_quartets_n":len(goodmonthly),
      "2023_numeric_same_month_original_blocks_n":len(qblocks),
      "2023_numeric_same_day_quartets_n":len(gooddaily),
      "source_year_2022_reported_WWTP_P_NOT_proven_2023_treatment":True,
      "no_real_2023_WWTP_effluent_mass_or_original_Lycaenidae_abundance_measured":True,
      "POSTHOC_to_known_2022_same_month_positive_quartet_result":True,
      "timestamp_utc":datetime.now(timezone.utc).isoformat()
    }
    if gate:
        result["2023_equal_quartet_numeric_same_month_summary"]=stats(goodmonthly)
        original_byq={q["original_quartet_index"]:q for q in prev["all_original_8_quartet_source_month_numeric_QC"]}
        common=[]
        for q in goodmonthly:
            orig=original_byq[q["source_original_2022_quartet_index"]]
            historical=orig["mean_matched_source_minus_control_monthly_logTP_gradient"]
            if historical is not None and orig["numeric_month_gate_this_quartet"]:
                common.append({
                  "quartet_index":q["source_original_2022_quartet_index"],
                  "original_2022_samemonth_logTP_gradient":historical,
                  "2023_samemonth_logTP_gradient":q["source_minus_registeredzero_four_arm_logTP_gradient_mean"],
                  "year23_minus22_difference":q["source_minus_registeredzero_four_arm_logTP_gradient_mean"]-historical})
        result["original22vs23_identical_quartets_only_temporal_persistence"]={
          "matched_original_quartet_n":len(common),
          "2022_original_year_mean_matchedsame_quartet":float(np.mean([p["original_2022_samemonth_logTP_gradient"] for p in common])) if common else None,
          "2023_mean_matchedsame_quartet":float(np.mean([p["2023_samemonth_logTP_gradient"] for p in common])) if common else None,
          "paired_change_2023_minus_2022_mean":float(np.mean([p["year23_minus22_difference"] for p in common])) if common else None,
          "per_original_quartet":common,
          "not_actual_2022_vs_2023_wastewater_treatment_change":True}
    if dgate:
        same_day_r=[{**q,"source_minus_registeredzero_four_arm_logTP_gradient_mean":float(np.mean([t["source_minus_reference_logTP_gradient"] for t in q["2023_sameday_numeric_source_vs_zero_logTP_gradients"]))} for q in gooddaily]
        result["2023_sameday_realTP_source_vs_zero_summary"]=stats(same_day_r)
    else:
        result["2023_sameday_realTP_source_vs_zero_summary"]=None
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GESUI_2023_WISE6_REAL_TP_ORIGINAL2022_QUARTETS_SAMEMONTH_NUMERIC",
       json.dumps({"original_8_quartets":len(quarts),
        "2023_quantified_monthly_quartets":len(goodmonthly),
        "2023_month_gate":gate,
        "2023_same_day_gate":dgate,
        "result_2023":result.get("2023_equal_quartet_numeric_same_month_summary"),
        "same_original_2022_2023_quartets_yearpair":result.get("original22vs23_identical_quartets_only_temporal_persistence",{}).get("matched_original_quartet_n"),
        "no_causal_wastewater_P":True},sort_keys=True),flush=True)
if __name__=="__main__":main()
