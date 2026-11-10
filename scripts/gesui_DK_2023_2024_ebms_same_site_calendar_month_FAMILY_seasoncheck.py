#!/usr/bin/env python3
"""POSTHOC Danish eBMS 2023–24 same-site same-calendar-month butterfly-family audit.

Original licensed EBMS V1.4 Event and Occurrence are rejoined; raw event and
individual occurrence records stay inside the workflow. Only group / cluster
aggregate results are archived. Not a wastewater treatment experiment.
"""
from __future__ import annotations
import argparse, hashlib, io, json, math, zipfile
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from datetime import datetime,timezone
from pathlib import Path

import gesui_DK_EBMS_V14_2023_2024_licensed_family_count_panel_after_taxonomy as first

ROOT=Path("results/ebms_denmark_2023_2024_WWTP_source_only")
PANEL=ROOT/"DENMARK_EBMS_2023_2024_FAMILY_SITEYEAR_COUNTS_QC_BEFORE_WWTP_CONTRAST.json"
FAMILIES=("Lycaenidae","Pieridae","Nymphalidae","Hesperiidae","Papilionidae")
YEARS=(2023,2024)
SEASONS={"same_calendar_months_all":tuple(range(1,13)),
         "same_calendar_months_May_Aug":(5,6,7,8)}
NEAR="near5";FAR="remote10"
NS=first.NS

def sha(b):return hashlib.sha256(b).hexdigest()

def hav(a,b):
    r=math.pi/180
    dl=(a["lat"]-b["lat"])*r;dg=(a["lon"]-b["lon"])*r
    v=math.sin(dl/2)**2+math.cos(a["lat"]*r)*math.cos(b["lat"]*r)*math.sin(dg/2)**2
    return 12742*math.asin(min(1, math.sqrt(v)))

def yearsources(root):
    site,tax,cluster,sites,taxa,families,srcsha=first.load_sources(root)
    panel_path=root/PANEL
    panelraw=json.loads(panel_path.read_text())
    if panelraw["original_archive_SHA256"]!=first.ORIGINAL_SHA:
        raise ValueError("Existing licensed source archive changed")
    if panelraw["exact_blinded_taxonomy_and_source_site_file_sha256"]!=srcsha:
        raise ValueError("Site/year taxonomy/source hashes changed from original certified count panel")
    panel=panelraw["site_year_family_counts_original_public_CC_BY_NC_aggregate"]
    if len(panel)!=214:raise ValueError("Original 214-site-year panel does not match")
    prior={(r["source_rounded5_site_hash"],r["site_year"]):r for r in panel}
    if len(prior)!=214:raise ValueError("Duplicate original site-year")
    closests={}
    for loc,s in sites.items():
        if s["nearest_qualified_EEA2022_DK_FW_outfall_km"]<=5:
            nn=sorted((hav(s,out),out["outfall_id_SHA256"]) for out in site["original_2022_EEA_DK_qualifying_Fw_outfalls"])
            if abs(nn[0][0]-s["nearest_qualified_EEA2022_DK_FW_outfall_km"])>1e-8:
                raise ValueError("Near source identity differs from previous source-only spatial exposure")
            closests[loc]=nn[0][1]
    if len(set(closests.values()))!=4:raise ValueError("Expected 4 independent registered source outfall clusters")
    return sites,families,prior,closests,srcsha

def counts_by_month(raw,sites,families):
    months=defaultdict(set)
    abundance=defaultdict(Counter)
    event_map={}
    all_event_ids=set()
    n=Counter()
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        root=ET.fromstring(z.read("meta.xml"))
        e=root.find(NS+"core")
        ex=[a for a in root.findall(NS+"extension") if a.attrib.get("rowType")=="http://rs.tdwg.org/dwc/terms/Occurrence"]
        if e is None or len(ex)!=1:raise ValueError("Expected Event core and exactly one licensed Occurrence extension")
        ef,efile,es,efield=first.core_terms(e)
        of,ofile,os,ofield=first.core_terms(ex[0])
        ek=efield("id");ok=ofield("coreid")
        for name in ("scientificName","individualCount","occurrenceStatus"):
            if name not in of:raise ValueError("Required original occurrence source field unavailable")
        for row in first.stream_csv(z,efile,es):
            n["all_original_Event_rows"]+=1
            if len(row)<=max(ek,ef["eventDate"],ef["decimalLatitude"],ef["decimalLongitude"]):
                raise ValueError("Unexpected short Event core row; refuse partial joined analysis")
            evt=row[ek].strip()
            if not evt or evt in all_event_ids:raise ValueError("Duplicate or empty DwCA Event core id")
            all_event_ids.add(evt)
            date=row[ef["eventDate"]].strip()[:10]
            try:
                year=int(date[:4]);month=int(date[5:7]);day=int(date[8:10])
                lat=float(row[ef["decimalLatitude"]]);lon=float(row[ef["decimalLongitude"]])
            except (ValueError,TypeError):continue
            if not (math.isfinite(lat) and math.isfinite(lon) and 1<=month<=12 and 1<=day<=31):continue
            loc=f"{lat:.5f}|{lon:.5f}"
            if year not in YEARS or loc not in sites:continue
            if sites[loc]["site_id_SHA256"]!=sha(loc.encode()):
                raise ValueError("Source site centroid SHA mismatch")
            event_map[evt]=(loc,year,month)
            months[(loc,year,month)].add(date)
            n["qualified_Event_source_records"]+=1
        for row in first.stream_csv(z,ofile,os):
            n["all_original_Occurrence_rows"]+=1
            if len(row)<=max(ok,of["scientificName"],of["individualCount"],of["occurrenceStatus"]):
                raise ValueError("Unexpected truncated licensed occurrence source row")
            coreid=row[ok].strip()
            if coreid not in event_map:continue
            loc,year,month=event_map[coreid]
            status=row[of["occurrenceStatus"]].strip()
            if status and status.lower() not in ("present","presence"):continue
            cnt=first.floatcount(row[of["individualCount"]])
            if cnt is None:continue
            fam=families.get(row[of["scientificName"]].strip())
            if fam not in FAMILIES:continue
            abundance[(loc,year,month)][fam]+=cnt
            n["original_accepted_family_counts_rows"]+=1
    if n["all_original_Event_rows"]!=15606 or n["all_original_Occurrence_rows"]!=21617:
        raise ValueError("Original official event-occurrence source bytes changed")
    return months,abundance,n

def validate_months(months,abundance,sites,prior):
    for loc,s in sites.items():
        for yr in YEARS:
            individual=Counter()
            found=set()
            for m in range(1,13):
                individual.update(abundance[(loc,yr,m)])
                found.update(months[(loc,yr,m)])
            old=prior[(s["site_id_SHA256"],yr)]
            if len(found)!=old["distinct_original_visit_days"]:
                raise ValueError(f"Original visit-day mismatch at {s['site_id_SHA256']} in {yr}")
            for fam in FAMILIES:
                if individual[fam]!=old["family_individual_counts"][fam]:
                    raise ValueError(f"Original annual family total mismatch for {fam} at {s['site_id_SHA256']} {yr}")
    return True

def category(site):
    d=site["nearest_qualified_EEA2022_DK_FW_outfall_km"]
    return NEAR if d<=5 else FAR if d>=10 else "intermediate"

def pct(a,b):
    return 100*a/b if b else None

def group_stats(rows,years):
    x={}
    for y in years:
        current=[r for r in rows if r["year"]==y]
        visit=sum(r["days"] for r in current)
        counts={f:sum(r["counts"][f] for r in current) for f in FAMILIES}
        L=counts["Lycaenidae"];P=counts["Pieridae"];Pa=counts["Papilionidae"]
        x[str(y)]={
            "original_matched_source_site_n":len(current),
            "original_site_month_survey_days":visit,
            "family_counts":counts,
            "family_individuals_per_visit_day":{f:counts[f]/visit if visit else None for f in FAMILIES},
            "Lycaenidae_share_given_three_target_families":(L/(L+P+Pa) if L+P+Pa else None),
            "target_family_count_support":L+P+Pa,
        }
    return x

def group_ratio(x,fam):
    a=x["2023"]["family_individuals_per_visit_day"][fam]
    b=x["2024"]["family_individuals_per_visit_day"][fam]
    return b/a if a and b is not None else None

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    sites,families,prior,closests,prior_sha=yearsources(root)
    raw=first.download()
    if sha(raw)!=first.ORIGINAL_SHA:raise ValueError("Original official licensed archive changed")
    months,abundance,source_qc=counts_by_month(raw,sites,families)
    assert validate_months(months,abundance,sites,prior)
    output={
        "schema":"gesui_DK_2023_2024_same_site_same_calendar_month_Lycaenidae_vs_Families_posthoc_v01",
        "utc":datetime.now(timezone.utc).isoformat(),
        "licensed_original_Denmark_eBMS_v14_SHA256":first.ORIGINAL_SHA,
        "original_source_and_taxonomy_sha256":prior_sha,
        "original_licensed_siteyear_family_source_SHA256":sha((root/PANEL).read_bytes()),
        "all_214_siteyear_original_family_sums_and_visit_days_exactly_replayed":True,
        "source_rows_qc":dict(source_qc),
        "month_matching_was_selected_posthoc_AFTER_known_positive_Germany_and_negative_DK_annual_result":True,
        "static_original_EEA2022_outfall_is_not_year_varying_butterfly_exposure":True,
        "no_causal_butterfly_population_benefit_or_harm_or_puddling_evidence":True,
        "families":FAMILIES,
        "no_original_individual_records_or_personal_observer_data_republished":True,
        "season_results":{}
    }
    for name,subset in SEASONS.items():
        rows=[]
        overlap_diag=Counter()
        raw_month_visits_by_group=defaultdict(Counter)
        for loc,s in sites.items():
            group=category(s)
            for m in subset:
                for yr in YEARS:
                    raw_month_visits_by_group[(group,yr)][m]+=len(months[(loc,yr,m)])
            common=[m for m in subset if months[(loc,2023,m)] and months[(loc,2024,m)]]
            overlap_diag[(group,"sites_with_ge1_shared_month")]+=bool(common)
            if not common:continue
            for y in YEARS:
                nv=sum(len(months[(loc,y,m)]) for m in common)
                if nv<=0:raise ValueError("Unmatched calendar month has no observed visits")
                cc={f:sum(abundance[(loc,y,m)][f] for m in common) for f in FAMILIES}
                rows.append({"site":s["site_id_SHA256"],"group":group,"year":y,
                             "days":nv,"counts":cc,
                             "original_1degree_region":s["country_1deg_lat_lon_cell"],
                             "source_cluster":closests.get(loc),
                             "shared_calendar_months":[int(x) for x in common]})
        result={"same_calendar_month_definition":list(subset),
                "within_site_calendar_month_matched_sites_per_source_class":{
                    g:int(overlap_diag[(g,"sites_with_ge1_shared_month")]) for g in (NEAR,FAR,"intermediate")},
                "original_all_site_month_sampling_days_unmatched_by_year":{str(y):{g:{str(m):int(raw_month_visits_by_group[(g,y)][m]) for m in subset} for g in (NEAR,FAR,"intermediate")} for y in YEARS},
                "matched_group_year":{},
                "four_source_outfall_near_group_2023_2024":[],
                "source_cluster_leave_one_out_5km_Lycaenidae_2023_2024_ratios":[],
                "remote_geographic_stratum_year":{}}
        for g in (NEAR,FAR,"intermediate"):
            matched=[x for x in rows if x["group"]==g]
            result["matched_group_year"][g]=group_stats(matched,YEARS)
        near=rows_for= [r for r in rows if r["group"]==NEAR]
        remote=[r for r in rows if r["group"]==FAR]
        cluster_ids=sorted(set(x["source_cluster"] for x in near))
        result["independent_registered_source_clusters_in_matched_sites"]=len(cluster_ids)
        for groupid in cluster_ids:
            rr=[r for r in near if r["source_cluster"]==groupid]
            result["four_source_outfall_near_group_2023_2024"].append({
                "outfall_id_sha256":groupid,
                "match":group_stats(rr,YEARS),
                "not_independent_site_treatment_replicates":True})
        remote_summary=group_stats(remote,YEARS)
        rr=group_ratio(remote_summary,"Lycaenidae")
        for groupid in cluster_ids:
            smaller=[r for r in near if r["source_cluster"]!=groupid]
            st=group_stats(smaller,YEARS)
            nr=group_ratio(st,"Lycaenidae")
            result["source_cluster_leave_one_out_5km_Lycaenidae_2023_2024_ratios"].append({
                "left_out_original_source_id_sha256":groupid,
                "remaining_2024_over_2023_near_Lycaenidae_per_survey_day_ratio":nr,
                "remote_2024_over_2023_rate_ratio":rr,
                "near_vs_remote_ratio_of_2024_2023_rate_ratios":nr/rr if nr is not None and rr is not None and rr>0 else None
            })
        result["descriptive_near_vs_remote_ratio_of_2024_to_2023_visit_normalized_rates_by_family"]={}
        for f in FAMILIES:
            nr=group_ratio(result["matched_group_year"][NEAR],f)
            fr=group_ratio(result["matched_group_year"][FAR],f)
            result["descriptive_near_vs_remote_ratio_of_2024_to_2023_visit_normalized_rates_by_family"][f]={
                "near_2024_over_2023":nr,"remote_2024_over_2023":fr,
                "ratio_of_rate_ratios":nr/fr if nr is not None and fr is not None and fr else None}
        for region in sorted(set(r["original_1degree_region"] for r in rows)):
            val=[r for r in rows if r["original_1degree_region"]==region]
            if {r["group"] for r in val}>={NEAR,FAR}:
                result["remote_geographic_stratum_year"][region]={
                    NEAR:group_stats([r for r in val if r["group"]==NEAR],YEARS),
                    FAR:group_stats([r for r in val if r["group"]==FAR],YEARS),
                    "not_independent_treatment_causal_test":True}
        output["season_results"][name]=result
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(output,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"source_replay":output["all_214_siteyear_original_family_sums_and_visit_days_exactly_replayed"],
                      "month_sampling":{k:v["within_site_calendar_month_matched_sites_per_source_class"] for k,v in output["season_results"].items()},
                      "family_relative_rate_ratio":{k:v["descriptive_near_vs_remote_ratio_of_2024_to_2023_visit_normalized_rates_by_family"] for k,v in output["season_results"].items()},
                      "source_cluster_units":{k:v["independent_registered_source_clusters_in_matched_sites"] for k,v in output["season_results"].items()}},indent=2))

if __name__=="__main__":
    main()
