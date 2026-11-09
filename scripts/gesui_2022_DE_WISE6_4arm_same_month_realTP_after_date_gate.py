#!/usr/bin/env python3
"""Source-gated 2022 EEA WISE6 real TP SAME MONTH matched quartet difference.

This numeric original 36GB streaming source is opened ONLY after original
8×4 monitoring station SOURCE-ONLY exact-day/month sample dates were
Git-committed. Strict same-day gate FAILS; do not fit a same-day effect.
Same-month four-arm 7/8 quartets/5 spatial blocks PASSED. 2022 original
source/interposed reported P was previously frozen BEFORE original annual
TP values. All new monthly quantitative results remain POSTHOC overall.

This is neither a randomized intervention nor adult butterfly puddling or
host-quality/abundance research.
"""
from __future__ import annotations
import argparse,hashlib,json,math,statistics,urllib.request,zipfile,itertools
from collections import defaultdict,Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np

ORIG_URL="https://discodata.eea.europa.eu/download/WISE_SOE/latest/Waterbase_T_WISE6_DisaggregatedData"
ORIG_ZIP_SHA="2a2ef2a43a16ec068b7974a80ae5fca49e44697660a46ce3b91fb1813abe62d2"
ORIG_ZIP_N=2024284753
ORIG_MEMBER="waterbase_t_wise6_disaggregateddata.csv"
ORIG_MEMBER_N=36043700357
ORIG_ROWS=102569057
ORIG_QUARTET_BLOB="8bd67be25968311ca16be4026888daa783fbe69f"
ORIG_NUMERIC32_SHA="28b679184cbebf6376f1c65d7c160f7511f057ff22787e98ff3fa265b617b261"
TP_CODE="CAS_7723-14-0"
TP_UNIT="mg{P}/L"
SAME_DAY_MIN_QUARTETS=4
SAME_DAY_MIN_BLOCKS=3
SAME_DAY_MIN_DAYS=3
SAME_MONTH_MIN_MONTHS=3
SAME_MONTH_MIN_QUARTETS=5
SAME_MONTH_MIN_BLOCKS=4
OFFSET=0.001
BOOT=9999
SEED=20261109

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()
def sitehash(c,scheme,site):
    return hashlib.sha256((c+"|"+scheme+"|"+site).encode()).hexdigest()
def strval(b):return b.strip().strip(b'"').decode("utf8",errors="replace")
def download(path):
    request=urllib.request.Request(ORIG_URL,headers={
      "User-Agent":"gesui-2022-DE-sourcelocked-WISE6-4arm-sameMONTH-TP-numeric-QC/1.0",
      "Accept":"application/zip,*/*"})
    h=hashlib.sha256();n=0
    with urllib.request.urlopen(request,timeout=180) as r,path.open("wb") as o:
        if r.status!=200:raise ValueError("Original EEA2026 source unavailable")
        for buf in iter(lambda:r.read(1024*1024),b""):
            n+=len(buf)
            if n>ORIG_ZIP_N:raise ValueError("Source original public archive exceeds version-locked size")
            h.update(buf);o.write(buf)
    if n!=ORIG_ZIP_N or h.hexdigest()!=ORIG_ZIP_SHA:
        raise ValueError("Original 2026 WISE6 published ZIP checksum changed")
    return {"original_source_url":ORIG_URL,"original_zip_bytes":n,"source_sha256":h.hexdigest()}

def original_site_numeric_quality(path,allhashes):
    sites={h:{"source_rawTP_rows":0,"belowLOQ":0,"invalidFlag":0,"missingResult":0,
              "invalidNumeric":0,"blankLOQflag":0,"status_codes":Counter(),
              "daily":defaultdict(list)} for h in allhashes}
    q=Counter()
    with zipfile.ZipFile(path) as z:
        if z.namelist()!=[ORIG_MEMBER]:raise ValueError("Unrecognized original WISE6 ZIP member")
        if z.getinfo(ORIG_MEMBER).file_size!=ORIG_MEMBER_N:
            raise ValueError("Unrecognized original WISE6 36GB uncompressed source version")
        with z.open(ORIG_MEMBER) as f:
            expected=("countryCode monitoringSiteIdentifier monitoringSiteIdentifierScheme parameterWaterBodyCategory observedPropertyDeterminandCode observedPropertyDeterminandLabel procedureAnalysedMatrix resultUom phenomenonTimeSamplingDate phenomenonTimeReferenceYear sampleIdentifier resultObservedValue resultQualityObservedValueBelowLOQ").split()
            cols=f.readline().decode("utf-8-sig").rstrip("\r\n").split(";")
            if cols[:13]!=expected or cols[21]!="resultObservationStatus":raise ValueError("WISE6 original numerical/quality column schema drift")
            for line in f:
                q["original_source_WISE6_rows"]+=1
                if not line.startswith((b"DE;",b'"DE";')):continue
                head=line.split(b";",11)
                if len(head)<12:continue
                if strval(head[9])!="2022" or strval(head[3])!="RW" or strval(head[5])!="Total phosphorus" or strval(head[6])!="W":continue
                if strval(head[4])!=TP_CODE or strval(head[7])!=TP_UNIT:
                    raise ValueError("Original 2022 DE river Total P code or measured unit drift")
                q["2022_original_DE_RW_TP_rows"]+=1
                h=sitehash("DE",strval(head[2]),strval(head[1]))
                if h not in sites:continue
                r=sites[h];r["source_rawTP_rows"]+=1
                fields=line.split(b";",22)  # only original columns up to resultObservationStatus; Remarks never parsed
                if len(fields)<23:raise ValueError("Cannot read original WISE6 sample status field")
                status=strval(fields[21])
                r["status_codes"][status]+=1
                flag=strval(fields[12]).lower()
                if flag in ("true","1"):r["belowLOQ"]+=1;continue
                if flag not in ("false","0",""):r["invalidFlag"]+=1;continue
                if flag=="":r["blankLOQflag"]+=1
                raw=strval(fields[11])
                if not raw:r["missingResult"]+=1;continue
                try:v=float(raw)
                except (ValueError,TypeError):r["invalidNumeric"]+=1;continue
                if not math.isfinite(v) or v<0:r["invalidNumeric"]+=1;continue
                date=strval(fields[8])[:10]
                if len(date)!=10 or not date.startswith("2022-"):raise ValueError("Unexpected source TP date")
                r["daily"][date].append(v)
    if q["original_source_WISE6_rows"]!=ORIG_ROWS:
        raise ValueError("Original EEA2026 WISE6 source row count changed")
    return sites,dict(q)

def rollup(site_numeric,original32QC,source_day_metadata):
    qc={x["original_country_scheme_site_SHA256"]:x for x in original32QC["source_original_32_monitoring_sites_2022_numeric_TP_QC"]}
    exact=source_day_metadata["exact_site_source_32_original_identifiers_hashed_and_calendar"]
    nitems={}
    for key,raw in site_numeric.items():
        original=qc[key]; dates=sorted(raw["daily"])
        months=sorted(set(date[:7] for date in dates))
        if raw["source_rawTP_rows"]!=original["original_2022_all_sample_rows"]:
            raise ValueError("Source station sample count diverged from earlier independently frozen numeric QC")
        if dates and not set(dates).issubset(set(exact[key]["2022_original_distinct_sample_dates"])):
            raise ValueError("Numeric valid original WISE6 sample date not present in pre-value source-only dates")
        dayvals={date:statistics.median(raw["daily"][date]) for date in dates}
        frac=raw["belowLOQ"]/raw["source_rawTP_rows"] if raw["source_rawTP_rows"] else 1
        good=len(dates)>=4 and len(months)>=3 and frac<=0.5
        annual=statistics.median(dayvals.values()) if good else None
        old=original["original_2022_annual_siteTP_median_of_daily_medians_mgP_L"]
        if (old is None)!=(annual is None) or (old is not None and abs(float(old)-annual)>1e-10):
            raise ValueError("Existing original 2022 TP annual site median source does not reproduce exactly")
        monthsvalues=defaultdict(list)
        for date,val in dayvals.items():monthsvalues[date[:7]].append(val)
        monthly={month:statistics.median(vals) for month,vals in monthsvalues.items()}
        if good and not original["original_site_numeric_4day_3month_LOQ_gate"]:
            raise ValueError("Original independently frozen numeric TP QC eligibility differs")
        nitems[key]={
            "full_date_site_medians_mgP_L":dayvals,
            "site_month_median_of_daily_medians_mgP_L":monthly,
            "source_rawTP_rows":raw["source_rawTP_rows"],
            "below_LOQ_rows":raw["belowLOQ"],
            "status_codes_original_WISE6_resultObservationStatus":dict(raw["status_codes"]),
            "site_source_LOQ_date_year_gate_pass":good,
            "original_numeric_annual_site_median_recomputed_exactly":True,
        }
    return nitems

def quartet_evaluation(q,sites,level):
    a=q["four_arm_original_WISE_country_scheme_site_hashes"]
    d={arm:sites[h]["full_date_site_medians_mgP_L" if level=="day" else "site_month_median_of_daily_medians_mgP_L"]
       for arm,h in a.items()}
    intersection=sorted(set.intersection(*(set(v) for v in d.values())))
    values=[]
    for t in intersection:
        if any(not sites[a[arm]]["site_source_LOQ_date_year_gate_pass"] for arm in a):
            continue
        s_up=float(d["source_upstream"][t]);s_down=float(d["source_downstream"][t])
        c_up=float(d["registered_zero_upstream"][t]);c_down=float(d["registered_zero_downstream"][t])
        src_gradient=math.log(s_down+OFFSET)-math.log(s_up+OFFSET)
        ctr_gradient=math.log(c_down+OFFSET)-math.log(c_up+OFFSET)
        values.append({"YYYY_MM_DD" if level=="day" else "YYYY_MM":t,
           "source_gradient_logTP":src_gradient,"control_gradient_logTP":ctr_gradient,
           "matched_gradient_source_minus_control":src_gradient-ctr_gradient,
           "source_upstream_TP_mgP_L":s_up,"source_downstream_TP_mgP_L":s_down,
           "control_upstream_TP_mgP_L":c_up,"control_downstream_TP_mgP_L":c_down,
           "same_day_true":level=="day",
        })
    return intersection,values

def contrast_result(quartets,level):
    if not quartets:return None
    v=np.asarray([q["mean_matched_source_minus_control_monthly_logTP_gradient"] if level=="month" else q["mean_matched_source_minus_control_daily_logTP_gradient"]
        for q in quartets],float)
    rng=np.random.default_rng(SEED)
    samples=np.mean(rng.choice(v,size=(BOOT,len(v)),replace=True),axis=1)
    observed=float(np.mean(v))
    refs=[float(np.mean(v*np.asarray(sign))) for sign in itertools.product((-1,1),repeat=len(v))]
    tail=sum(abs(x)>=abs(observed)-1e-12 for x in refs)/len(refs)
    loo=[float(np.mean(np.delete(v,i))) for i in range(len(v))] if len(v)>1 else []
    return {
      "original_independent_matched_quartets_n":len(quartets),
      "distinct_original_country2degree_source_blocks_n":len({q["source_country2degree_upstream_block"] for q in quartets}),
      "quartet_equal_weight_mean_matched_source_minus_sourcezero_logTP_gradient":observed,
      "quartet_equal_weight_median_matched_gradient":float(np.median(v)),
      "positive_contrast_quartets":int(np.sum(v>0)),
      "negative_contrast_quartets":int(np.sum(v<0)),
      "descriptive_9999_quartet_resample95pct_interval":[float(x) for x in np.quantile(samples,[.025,.975])],
      "source_only_equal_quartet_signflip_2sided_reference_fraction":tail,
      "original_quartet_leave_one_out_mean_range":[min(loo),max(loo)] if loo else None,
      "signflip_is_NOT_randomized_WWTP_treatment_p_value":True,
      "annual_2022_reported_WWTP_P_load_not_monthly_wastewater_dose":True,
      "TP_is_river_phosphorus_not_Lycaenidae_adult_puddling_or_density":True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--pre-numeric-original-dates",type=Path,required=True)
    p.add_argument("--original-32site-numeric-QC",type=Path,required=True)
    p.add_argument("--tempdir",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    src=json.loads(a.pre_numeric_original_dates.read_text())
    prev=json.loads(a.original_32site_numeric_QC.read_text())
    if src["schema"]!="gesui_2022_original_EEA_WISE6_Germany_8_source_zero_TP_quartet_exact_day_and_month_preoutcome_source_gate_v01":
        raise ValueError("New date-specific input not original pre-value source-only 32 station panel")
    if src["NO_original_TP_numeric_resultObservedValue_or_LOQ_values_read"] is not True:
        raise ValueError("Date-specific source-only gate was contaminated with numeric TP")
    if src["source_original_32site_numeric_QC_SHA256"]!=sha(a.original_32site_numeric_QC):
        raise ValueError("Original 32-site numeric QC input differs after pre-value date eligibility")
    runmonth=src["source_same_calendar_month_fourstation_gate_PASS"]
    runday=src["source_exact_fourstation_calendar_date_gate_PASS"]
    root={
      "schema":"gesui_2022_real_WISE6_German_4station_WWTPinterposed_vs_registeredzero_same_day_and_month_after_nonvalue_date_gate_v01",
      "timestamp":datetime.now(timezone.utc).isoformat(),
      "original_source_first_sample_date_gate_SHA256":sha(a.pre_numeric_original_dates),
      "original_32site_numeric_QC_SHA256":sha(a.original_32site_numeric_QC),
      "POSTHOC_2022_month_temporal_stress_after_known_positive_7quarts_annual_gradient":True,
      "source_same_day_preoutcome_gate":runday,
      "source_same_month_preoutcome_gate":runmonth,
      "no_Lycaenidae_species_photo_or_behavior_data":True,
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    if not(runday or runmonth):
        root.update({"stage":"STOP_NO_NUMERIC_TP_OPENED","reason":"Both four-arm pre-value date and month gates failed"})
        a.out.write_text(json.dumps(root,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
        print("GESUI_SOURCE_GATED_2022_GERMAN_SAMEDAY_AND_MONTH_BOTH_STOP_NO_NUMERIC_TP",flush=True)
        return
    a.tempdir.mkdir(parents=True,exist_ok=True)
    zipdst=a.tempdir/"original_EEA2026_2022_WISE6_German_monthly_TP.zip"
    receipt=download(zipdst)
    sites,qc=original_site_numeric_quality(zipdst,set(src["exact_site_source_32_original_identifiers_hashed_and_calendar"]))
    annual=rollup(sites,prev,src)
    monthq=[];dayq=[];allq=[]
    for q in src["all_original_8_no_replacement_quartets_and_date_coverage"]:
        daydates,dayvals=quartet_evaluation(q,annual,"day")
        monthdates,monthvals=quartet_evaluation(q,annual,"month")
        daygood=len(dayvals)>=SAME_DAY_MIN_DAYS and runday
        monthgood=len(monthvals)>=SAME_MONTH_MIN_MONTHS and runmonth
        r={
          "original_quartet_index":q["quartet_index"],
          "source_country2degree_upstream_block":q["source_country2degree_upstream_block"],
          "original_reported_interposed_P_tonnes_year":q["original_source_interposed_P_tonnes_per_year"],
          "original_source_only_four_site_hashes":q["four_arm_original_WISE_country_scheme_site_hashes"],
          "original_four_arm_numeric_quantified_sameday_dates":daydates,
          "original_four_arm_numeric_quantified_months":monthdates,
          "sameday_fourarm_valid_n":len(dayvals),
          "samemonth_fourarm_valid_n":len(monthvals),
          "numeric_exact_day_gate_this_quartet":daygood,
          "numeric_month_gate_this_quartet":monthgood,
          "month_specific_all4_TP_measurements_and_source_control_deltas":monthvals if runmonth else [],
          "day_specific_all4_TP_measurements_and_source_control_deltas":dayvals if runday else [],
          "mean_matched_source_minus_control_monthly_logTP_gradient":statistics.mean(x["matched_gradient_source_minus_control"] for x in monthvals) if monthgood else None,
          "mean_matched_source_minus_control_daily_logTP_gradient":statistics.mean(x["matched_gradient_source_minus_control"] for x in dayvals) if daygood else None,
        }
        allq.append(r)
        if daygood:dayq.append(r)
        if monthgood:monthq.append(r)
        print("GESUI_DE2022_4ARM_2022_ORIGINAL_SAMEMONTH_NUMERIC_TP_QUARTET",
          json.dumps({"quartet":q["quartet_index"],
              "numeric_same_day":len(dayvals),"numeric_same_month":len(monthvals),
              "month_qualified":monthgood,
              "source_minus_control_matched_monthly_logTP_mean":r["mean_matched_source_minus_control_monthly_logTP_gradient"]},sort_keys=True),flush=True)
    dayok=runday and len(dayq)>=SAME_DAY_MIN_QUARTETS and len({q["source_country2degree_upstream_block"] for q in dayq})>=SAME_DAY_MIN_BLOCKS
    monthok=runmonth and len(monthq)>=SAME_MONTH_MIN_QUARTETS and len({q["source_country2degree_upstream_block"] for q in monthq})>=SAME_MONTH_MIN_BLOCKS
    root.update({
      "source_full_EEA2026_WISE6_zip_SHA_original_numeric_2022_TP":receipt,
      "source_2022_original_whole_WISE6_TP_rows_checked":qc,
      "all_original_source_32monitoring_site_raw_LOQ_and_status_qc":{k:{
            "total_source_rows":v["source_rawTP_rows"],
            "source_below_LOQ_rows":v["below_LOQ_rows"],
            "original_resultObservationStatus_source_codes":v["status_codes_original_WISE6_resultObservationStatus"],
            "qualified_station":v["site_source_LOQ_date_year_gate_pass"],
            "original_annual_site_median_replayed":v["original_numeric_annual_site_median_recomputed_exactly"]}
          for k,v in annual.items()},
      "all_original_8_quartet_source_month_numeric_QC":allq,
      "source_4site_numeric_valid_ge3_day_quartets_n":len(dayq),
      "source_4site_numeric_valid_ge3_month_quartets_n":len(monthq),
      "source_4site_numeric_same_day_geographic_gate_pass":dayok,
      "source_4site_numeric_same_month_geographic_gate_pass":monthok,
      "NO_SAMEDAY_EFFECT_ESTIMATE_IF_BEFORE_OUTCOME_DAY_GATE_FAILED":not runday,
      "source_same_day_valid_match_results":contrast_result(dayq,"day") if dayok else None,
      "source_same_month_valid_match_results":contrast_result(monthq,"month") if monthok else None,
      "source_zero_intervening_source_is_NOT_no_other_pollution":True,
      "source_original_HydroRIVERS_conversion_not_independently_official_verified":True,
      "no_2022_WWTP_actual_daily_treatment_discharge_measured":True,
      "no_butterfly_ecological_nutrient_effect_claim":True,
    })
    root["stage"]="POSTHOC_SAMEMONTH_MEASURED_TP_RESULT" if monthok else "STOP_NUMERIC_TP_TIMING_GATE"
    a.out.write_text(json.dumps(root,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_GERMANY_2022_ORIGINAL_4ARM_SAMEMONTH_SOURCE_MINUS_SOURCEZERO_TP_DONE",
      json.dumps({"original_2022_quartets":len(allq),"same_day_source_gate":runday,
        "same_month_source_gate":runmonth,
        "same_day_numeric_supported":len(dayq),"same_month_numeric_supported":len(monthq),
        "same_month_model_permitted":monthok,
        "month_effect":root["source_same_month_valid_match_results"],
        "sample_status_codes_at_source":dict(Counter(str(k) for v in annual.values() for k in v["status_codes_original_WISE6_resultObservationStatus"]))},sort_keys=True),flush=True)
if __name__=="__main__":main()
