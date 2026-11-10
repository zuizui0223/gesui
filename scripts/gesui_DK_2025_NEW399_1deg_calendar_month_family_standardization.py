#!/usr/bin/env python3
"""POSTHOC Denmark eBMS NEW399 2025 original calendar-month × region standardization.

Replays licensed original Event/Occurrence rows and SOURCE-FROZEN 399 site
and 286 linked-date eligibility; never selects site, radius, comparison family
or taxonomic mapping to optimize outcome direction.

This is an observational recording-rate/denominator audit. One-degree areas
are NOT independent randomized WWTP treatments, nor proof of downstream dose.
"""
from __future__ import annotations
import argparse,hashlib,io,json,math,zipfile
import xml.etree.ElementTree as ET
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

import gesui_DK_EBMS_V14_2023_2024_licensed_family_count_panel_after_taxonomy as original

DIR=Path("results/ebms_denmark_2023_2024_WWTP_source_only")
SOURCE=DIR/"DENMARK_2025_NEW399_SITE_EVENT_ANY_OCC_WWTP20_SOURCE_GATE_BEFORE_FAMILY.json"
COUNTS=DIR/"DENMARK_2025_NEW399_SITE_LYCAENIDAE_FAMILY_POST_SOURCE_GATE_COUNTS.json"
SITE=DIR/"DENMARK_EBMS2023_2024_EEA2022_WWTP_SOURCEONLY_EVENT_QC.json"
FAMS=("Lycaenidae","Pieridae","Papilionidae","Nymphalidae","Hesperiidae")
GROUPS=("near5","remote10")
PRIMARY=("Lycaenidae","Pieridae","Papilionidae")
SEASONS={"all_common_calendar_months":set(range(1,13)),"May_to_Aug_common_calendar_months":set((5,6,7,8))}
SEED=20261110

def sha(data):
    return hashlib.sha256(data).hexdigest()

def get_strata(rows):
    sub=[r for r in rows if r["near"]["linked_dates"]>=2 and r["remote"]["linked_dates"]>=2]
    return sub

def standardize(strata):
    if not strata:return None
    totalweight=0
    weighted={g:dict.fromkeys(FAMS,0.) for g in GROUPS}
    for s in strata:
        dn=s["near"]["linked_dates"];dr=s["remote"]["linked_dates"]
        if not(dn>=2 and dr>=2):raise ValueError("Stratum lacks pre-set common reporting support")
        w=min(dn,dr)
        totalweight+=w
        for g in GROUPS:
            # Original region/month stratum uses "near"/"remote" display keys.
            slot="near" if g=="near5" else "remote"
            for f in FAMS:
                weighted[g][f]+=w*s[slot]["counts"][f]/s[slot]["linked_dates"]
    if totalweight<=0:return None
    for g in GROUPS:
        for f in FAMS:weighted[g][f]/=totalweight
    def share(fs,g):
        denom=sum(weighted[g][f] for f in fs)
        return weighted[g]["Lycaenidae"]/denom if denom>0 else None
    def odds(value):
        return value/(1-value) if value is not None and 0<value<1 else None
    lshare3={g:share(PRIMARY,g) for g in GROUPS}
    lshare5={g:share(FAMS,g) for g in GROUPS}
    lrate=weighted["near5"]["Lycaenidae"]/weighted["remote10"]["Lycaenidae"] if weighted["remote10"]["Lycaenidae"]>0 else None
    a3,b3=odds(lshare3["near5"]),odds(lshare3["remote10"])
    a5,b5=odds(lshare5["near5"]),odds(lshare5["remote10"])
    return {
        "matched_region_month_strata":len(strata),
        "matching_weight_sum_min_reporting_days":totalweight,
        "documented_linked_dates_in_selected_strata":{g:sum(x[g]["linked_dates"] for x in strata) for g in GROUPS},
        "standardized_each_family_recorded_per_linked_date":weighted,
        "Lycaenidae_per_linked_date_near_vs_remote_ratio":lrate,
        "Lycaenidae_three_family_share":lshare3,
        "Lycaenidae_five_family_share":lshare5,
        "Lycaenidae_three_family_share_oddsratio_near_over_remote":a3/b3 if a3 is not None and b3 is not None else None,
        "Lycaenidae_five_family_share_oddsratio_near_over_remote":a5/b5 if a5 is not None and b5 is not None else None,
        "not_independent_effluent_inference":True
    }

def percentile(v,q):
    a=sorted(v)
    return a[int((len(a)-1)*q)]

def compute_bootstrap(regions,byregion,n=4000):
    # Cluster resample entire 1deg area and all its matched original months,
    # NOT pseudo-replicated site visits or WWTP IDs. Descriptive posthoc only.
    keys=("Lycaenidae_per_linked_date_near_vs_remote_ratio",
          "Lycaenidae_three_family_share_oddsratio_near_over_remote",
          "Lycaenidae_five_family_share_oddsratio_near_over_remote")
    series={k:[] for k in keys}
    state=SEED
    for __ in range(n):
        strata=[]
        for _ in regions:
            state=(1664525*state+1013904223)&0xffffffff
            k=regions[int((state/(2**32))*len(regions))]
            strata.extend(byregion[k])
        x=standardize(strata)
        if x is None:continue
        for k in keys:
            if x[k] is not None and x[k]>0:series[k].append(x[k])
    return {k:{
        "n_valid":len(a),
        "descriptive_bootstrap95_ratio":[percentile(a,.025),percentile(a,.975)] if len(a)>100 else None,
        "positive_over_unit_fraction":sum(y>1 for y in a)/len(a) if a else None
    } for k,a in series.items()}

def run(root):
    frozen=json.loads((root/SOURCE).read_text())
    family=json.loads((root/COUNTS).read_text())
    site_all=json.loads((root/SITE).read_text())
    if frozen["frozen_before_2025_NEW399_butterfly_counts_gate"]["2025_NEW399_butterfly_family_count_outcome_permitted"] is not True:
        raise ValueError("Original independently locked 399/286 source gate failed")
    if not family["source_family_model_sample_gate"]["original_2025_family_contrast_reporting_gate_PASS"]:
        raise ValueError("Original independently frozen family response support gate failed")
    if family["source_blinded_to_2025_family_frozen_PREOUTCOME_JSON_SHA256"]!=sha((root/SOURCE).read_bytes()):
        raise ValueError("Source freeze changed after the biological response")
    old={r["original_public_site_SHA256"]:r for r in family["2025_original_399_site_accepted_family_counts"]}
    if len(old)!=399:raise ValueError("Unexpected 399 licensed result site count")
    siteindex={s["site_id_SHA256"]:s for s in site_all["original_2023_2024_repeated_event_site_source_metadata"]}
    bycoord={f"{siteindex[sid]['lat']:.5f}|{siteindex[sid]['lon']:.5f}":sid for sid in old}
    if len(bycoord)!=399 or not all(sha(s.encode())==id for s,id in bycoord.items()):
        raise ValueError("Frozen public coordinate SHA mismatch")
    _,taxmanifest,_,_,taxa,familymap,_=original.load_sources(root)
    if not taxmanifest["original_individualCount_values_NOT_READ"] or len(taxa)!=112:
        raise ValueError("Accepted original independent species-family mapping changed")
    raw=original.download()
    if sha(raw)!=original.ORIGINAL_SHA:raise ValueError("Original licensed DwCA v1.4 checksum changed")
    event_lookup={}
    orig_dates=defaultdict(set)
    recorded_dates=defaultdict(set)
    count_month=defaultdict(Counter)
    qa=Counter()
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        meta=ET.fromstring(z.read("meta.xml"))
        ec=meta.find(original.NS+"core")
        extensions=[v for v in meta.findall(original.NS+"extension") if v.attrib.get("rowType")=="http://rs.tdwg.org/dwc/terms/Occurrence"]
        if ec is None or len(extensions)!=1:raise ValueError("Unexpected original DwCA structure")
        ef,efile,es,efield=original.core_terms(ec)
        of,ofile,os,ofield=original.core_terms(extensions[0])
        ei=efield("id");oi=ofield("coreid")
        unique=set()
        for line in original.stream_csv(z,efile,es):
            qa["original_Event_rows"]+=1
            if len(line)<=max(ei,ef["eventDate"],ef["decimalLatitude"],ef["decimalLongitude"]):
                raise ValueError("Short original Event row")
            eid=line[ei].strip()
            if not eid or eid in unique:raise ValueError("Original Event ID not unique")
            unique.add(eid)
            date=line[ef["eventDate"]].strip()[:10]
            try:
                yr=int(date[:4]);m=int(date[5:7]);la=float(line[ef["decimalLatitude"]]);lo=float(line[ef["decimalLongitude"]])
            except (ValueError,TypeError):continue
            if yr!=2025 or not 1<=m<=12 or not(math.isfinite(la) and math.isfinite(lo)):continue
            sid=bycoord.get(f"{la:.5f}|{lo:.5f}")
            if sid is None:continue
            event_lookup[eid]=(sid,m,date)
            orig_dates[(sid,m)].add(date)
            qa["original_NEW399_2025_Event_rows"]+=1
        for line in original.stream_csv(z,ofile,os):
            qa["original_Occurrence_rows"]+=1
            if len(line)<=max(oi,of["individualCount"],of["occurrenceStatus"],of["scientificName"]):
                raise ValueError("Short original Occurrence row")
            evt=event_lookup.get(line[oi].strip())
            if evt is None:continue
            sid,m,date=evt
            recorded_dates[(sid,m)].add(date)
            qa["original_NEW399_2025_joined_Occurrence_rows"]+=1
            status=line[of["occurrenceStatus"]].strip()
            if status and status.lower() not in ("present","presence"):continue
            n=original.floatcount(line[of["individualCount"]])
            if n is None:continue
            scientific=line[of["scientificName"]].strip()
            if scientific not in taxa:raise ValueError("New taxon absent from originally blind GBIF source map")
            fam=familymap.get(scientific)
            if fam in FAMS:
                count_month[(sid,m)][fam]+=n
                qa["original_accepted_family_rows"]+=1
    if qa["original_Event_rows"]!=15606 or qa["original_Occurrence_rows"]!=21617:
        raise ValueError("Original licensed DwCA source is missing records")
    for sid,s in old.items():
        ovisit=sum(len(orig_dates[(sid,m)]) for m in range(1,13))
        linked=sum(len(recorded_dates[(sid,m)]) for m in range(1,13))
        if ovisit!=s["original_event_visit_days"] or linked!=s["linked_visit_days"]:
            raise ValueError("Original source-locked daily recording support differs")
        for fam in FAMS:
            n=sum(count_month[(sid,m)][fam] for m in range(1,13))
            if n!=s["family_individualCount"][fam]:
                raise ValueError(f"Original new399 accepted {fam} count differs from previous independent source audit")
    eligible=[(sid,s) for sid,s in old.items() if s["source_gate_2plus_any_Occurrence_dates"]]
    if len(eligible)!=286:raise ValueError("Original preoutcome published-Occurrence support changed")
    strata=defaultdict(lambda:defaultdict(lambda:{"linked_dates":0,"counts":Counter(),"sites":set()}))
    for sid,s in eligible:
        if s["source_band"] not in GROUPS:continue
        for month in range(1,13):
            n=len(recorded_dates[(sid,month)])
            if n==0:continue
            k=(s["original_geo_1deg_cell"],month)
            cell=strata[k][s["source_band"]]
            cell["linked_dates"]+=n
            cell["sites"].add(sid)
            for fam in FAMS:cell["counts"][fam]+=count_month[(sid,month)][fam]
    rows=[]
    for (region,month),g in sorted(strata.items()):
        rows.append({
            "original_1degree_region":region,
            "calendar_month":month,
            "near":{
                "linked_dates":g["near5"]["linked_dates"],
                "source_site_rows_reporting_in_month":len(g["near5"]["sites"]),
                "counts":{f:int(g["near5"]["counts"][f]) for f in FAMS}},
            "remote":{
                "linked_dates":g["remote10"]["linked_dates"],
                "source_site_rows_reporting_in_month":len(g["remote10"]["sites"]),
                "counts":{f:int(g["remote10"]["counts"][f]) for f in FAMS}}
        })
    report={}
    for title,season in SEASONS.items():
        sub=get_strata([x for x in rows if x["calendar_month"] in season])
        rg=sorted(set(x["original_1degree_region"] for x in sub))
        days_near=sum(x["near"]["linked_dates"] for x in sub)
        days_remote=sum(x["remote"]["linked_dates"] for x in sub)
        pass_support=len(rg)>=3 and len(sub)>=8 and min(days_near,days_remote)>=50
        byregion={r:[x for x in sub if x["original_1degree_region"]==r] for r in rg}
        est=standardize(sub) if pass_support else None
        loo={}
        if pass_support:
            for left in rg:
                fits=standardize([x for x in sub if x["original_1degree_region"]!=left])
                if fits:loo[left]={k:fits[k] for k in (
                    "Lycaenidae_per_linked_date_near_vs_remote_ratio",
                    "Lycaenidae_three_family_share_oddsratio_near_over_remote",
                    "Lycaenidae_five_family_share_oddsratio_near_over_remote")}
        report[title]={
            "prespecified_after_annual_results_before_month_values_observed_source_support":{
                "min_regions":3,"min_matched_region_months":8,"min_each_group_linked_days":50,
                "common_region_month_source_effort_gate_PASS":pass_support},
            "regions_with_matched_months_n":len(rg),
            "regions_with_matched_months":rg,
            "source_near_and_remote_matched_region_months_n":len(sub),
            "source_near_matched_linked_days":days_near,
            "source_remote_matched_linked_days":days_remote,
            "posthoc_min_date_weight_standardized_family_result":est,
            "four_thousand_region_cluster_resamples_descriptive95":compute_bootstrap(rg,byregion) if pass_support else None,
            "leave_one_original_1degree_region_out":loo,
            "matched_original_source_region_month_aggregate_counts":sub,
        }
    return {
        "schema":"gesui_DK_eBMS_2025_NEW399_calendar_month_region_observation_standardization_v01",
        "utc":datetime.now(timezone.utc).isoformat(),
        "original_licensed_Denmark_DWCA_SHA256":original.ORIGINAL_SHA,
        "prior_399_2025_site_family_count_SHA256":sha((root/COUNTS).read_bytes()),
        "preoutcome_original_399_site_source_gate_SHA256":sha((root/SOURCE).read_bytes()),
        "independent_399_site_2025_five_family_counts_and_all_2plus_linked_days_exactly_reproduced":True,
        "only_286_pre_2025_family_eligibility_source_sites_in_2025_month_outcome_analysis":True,
        "all_earlier_107_2023_2025_site_ids_excluded_from_current_cohort":True,
        "licensed_original_aggregate_Event_Occurrence_2025_QC":dict(qa),
        "date_criterion_strata_based_only_on_original_published_ANY_occurrence_not_family_count":True,
        "regional_calendar_month_common_support_inference_is_strictly_posthoc_descriptive":True,
        "2022_EEA_WWTP_distance_NOT_verified_2025_effluent_or_water_receptor_exposure":True,
        "2025_model_has_no_randomized_or_causal_independent_WWTP_treatment":True,
        "seasonal_and_geographic_source_band_family_results":report,
    }

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args()
    root=Path(__file__).resolve().parents[1]
    o=run(root)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({k:{"source_support":v["prespecified_after_annual_results_before_month_values_observed_source_support"],
                        "matched_regions":v["regions_with_matched_months_n"],
                        "region_months":v["source_near_and_remote_matched_region_months_n"],
                        "matched_linked_dates_near":v["source_near_matched_linked_days"],
                        "matched_linked_dates_remote":v["source_remote_matched_linked_days"],
                        "result":v["posthoc_min_date_weight_standardized_family_result"]}
                      for k,v in o["seasonal_and_geographic_source_band_family_results"].items()},indent=2))

if __name__=="__main__":
    main()
