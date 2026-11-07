#!/usr/bin/env python3
"""Network-free independent audit of new six-German-river TP + 18-site Plantae.

Uses versioned derived TP monthly source and original GBIF speciesKey facets.
Recomputes all TP/plant metrics, rank correlations and exhaustive n=6 p.
Does NOT claim that effluent caused nutrient or botanical changes.
"""
from __future__ import annotations
import csv
import hashlib
from itertools import permutations
import json
import math
from pathlib import Path
import statistics

ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/"results/independent_deu_tp_plantae_20261107"
PARENT=ROOT/"results/independent_deu_tp_source_cohort_20261107"

def close(a,b,tol=2e-6):
    if a is None or b is None:
        assert a is None and b is None,(a,b)
    else:
        assert math.isclose(float(a),float(b),abs_tol=tol,rel_tol=tol),(a,b)

def ranks(vals):
    assert len(set(vals))==len(vals),"Unexpected tied ranks"
    out=[0]*len(vals)
    for rank,index in enumerate(sorted(range(len(vals)),key=lambda i:vals[i]),1):
        out[index]=rank
    return out

def rho(a,b):
    assert len(a)==len(b)>2
    x,y=ranks(a),ranks(b)
    n=len(x)
    return 1-6*sum((u-v)**2 for u,v in zip(x,y))/(n*(n*n-1))

def rarefaction(counts,depth):
    v=[int(c) for c in counts if int(c)>0]
    N=sum(v)
    if N<depth:return None
    def comb(n,k):
        return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)
    den=comb(N,depth)
    return sum((1-math.exp(comb(N-k,depth)-den)) if N-k>=depth else 1 for k in v)

def verify(folder=FOLDER,parent=PARENT):
    receipt=json.loads((folder/"SOURCE_RECEIPT.json").read_text())
    for fname,expected in receipt["original_artifact_source_file_sha256"].items():
        if not (folder/fname).is_file():
            # Site-level raw observations were deliberately not redistributed publicly.
            assert fname=="observed_TP_station_dates_2011_2022.csv"
            continue
        sha=hashlib.sha256((folder/fname).read_bytes()).hexdigest()
        assert sha==expected,(fname,sha,expected)
    cohort=json.loads((parent/"STRICT_STATION_INDEPENDENCE_COHORTS.json").read_text())
    original=cohort["independent_source_only_cohorts"]["TP"]
    assert len(original)==6 and len({w["hydrobasin_l12"] for w in original})==6
    assert len({s["station_id"] for w in original for s in w["measured_source_stations"]["TP"].values()})==12
    design=json.loads((folder/"FROZEN_DESIGN.json").read_text())
    assert design["source_selected_without_biology_or_TP_values"] is True
    assert design["source_cohort_sha256"]==hashlib.sha256((parent/"STRICT_STATION_INDEPENDENCE_COHORTS.json").read_bytes()).hexdigest()
    assert len(design["original_source_plant_selections"])==18
    tp=json.loads((folder/"TP_SOURCE_RESULTS.json").read_text())
    plant=json.loads((folder/"PLANT_RICHNESS_AND_COMPOSITION.json").read_text())
    joint=json.loads((folder/"JOINT_TP_PLANT_DESCRIPTIVE.json").read_text())
    raw=json.loads((folder/"PLANT_original_GBIF_species_facets.json").read_text())
    exact=json.loads((folder/"EXACT_SMALL_N_AUDIT.json").read_text())
    assert len(raw["sites"])==18 and len(tp["all_results"])==6
    assert tp["source_original_zip_MD5"]=="b596dbe9833c49426147a99f10f3fe31"
    assert tp["observed_TP_positive_values_2011_2022"]==2500
    assert tp["paired_calendar_month_rows"]==836
    monthly={}
    with (folder/"paired_monthly_TP_mg_L.csv").open(newline="") as f:
        for row in csv.DictReader(f):
            wid=int(row["outfall_id"])
            monthly.setdefault(wid,[]).append(row)
            u=float(row["upstream_TP_mg_L"])
            d=float(row["downstream_TP_mg_L"])
            assert u>0 and d>0 and 2011<=int(row["yearmonth"][:4])<=2022
            close(float(row["log_down_over_up"]),math.log(d/u))
    assert sum(map(len,monthly.values()))==836 and len(monthly)==6
    by_source={int(x["outfall_id"]):x for x in tp["all_results"]}
    for wid,rows in monthly.items():
        derived=statistics.median(math.log(float(r["downstream_TP_mg_L"])
                                            /float(r["upstream_TP_mg_L"])) for r in rows)
        close(derived,by_source[wid]["median_monthly_log_down_up"])
        close(math.exp(derived),by_source[wid]["median_monthly_down_up_ratio"])
        assert len(rows)==by_source[wid]["calendar_months_with_measured_up_down"]
        assert len(rows)>=12 and by_source[wid]["primary_source_months_gate_met"]
    fac={x["site_id"]:x for x in raw["sites"]}
    assert len(fac)==18
    species_table=set()
    with (folder/"PLANT_site_species_records.csv").open(newline="") as f:
        for row in csv.DictReader(f):
            ident=(row["site_id"],str(row["gbif_speciesKey"]))
            assert ident not in species_table
            species_table.add(ident)
            assert int(row["speciesKey_observation_count"])==int(fac[row["site_id"]]["species_key_record_counts"][ident[1]])
            assert row["gbif_window"]=="2011,2022"
    assert len(species_table)==sum(len(a["species_key_record_counts"]) for a in fac.values())
    pv={int(x["wwtp_id"]):x for x in plant["wwtp_outcomes"]}
    assert len(pv)==6
    count30=0
    count10=0
    for w in original:
        wid=w["wwtp_id"]
        sites=[fac[f"{wid}_{role}"] for role in ("u2","u1","d1")]
        rows=[x["indexed_species_record_count"] for x in sites]
        depth_ratio=max(rows)/min(rows) if min(rows) else math.inf
        p30=all(x["indexed_species_record_count"]>=30 and
                x["distinct_species_keys"]>=6 for x in sites) and depth_ratio<=10
        p10=all(x["indexed_species_record_count"]>=10 and
                x["distinct_species_keys"]>=4 for x in sites)
        assert p30==pv[wid]["strict_primary_eligible"]
        assert p10==pv[wid]["depth10_sensitivity"]["eligible_by_records_and_species"]
        if p30:count30+=1
        if p10:count10+=1
        for depth,key,supported in ((30,"depth30",p30),(10,"depth10_sensitivity",p10)):
            if not supported:continue
            rare=[rarefaction(s["species_key_record_counts"].values(),depth) for s in sites]
            vals=pv[wid][key]["rarefied_expected_plant_species_u2_u1_d1"]
            for a,b in zip(rare,vals):close(a,b)
            y=math.log1p(rare[2])-2*math.log1p(rare[1])+math.log1p(rare[0])
            close(y,pv[wid][key]["extra_change_downstream_minus_prior_upstream"])
    assert count30==3 and count10==6
    entry=sorted(joint["joined"],key=lambda x:x["wwtp_id"])
    x=[q["source_fraction"] for q in entry]
    y=[q["observed_median_log_TP_down_up"] for q in entry]
    source_rho=rho(x,y)
    close(source_rho,joint["observed_TP_gradient_vs_model_wwflow"]["rho"])
    scores=[rho(x,list(p)) for p in permutations(y)]
    exact_p=sum(abs(v)>=abs(source_rho)-1e-12 for v in scores)/len(scores)
    close(exact_p,exact["two_sided_exact_rank_permutation_p"],tol=1e-12)
    assert len(scores)==720 and sum(abs(v)>=abs(source_rho)-1e-12 for v in scores)==42
    flora=[q["plant_Y_log_rarefied10"] for q in entry]
    flora_rho=rho(flora,y)
    close(flora_rho,joint["depth10_sensitivity_TP_gradient"]["rho_observed_TP_gradient_vs_GBIF_richness_curvature"])
    assert abs(flora_rho-(-0.028571428571428574))<1e-10
    return {
      "status":"PASS_INDEPENDENT_STDLIB_6_TP_18_GBIF_AUDIT",
      "real_source_run":receipt["real_GitHub_Actions_source_run"],
      "source_file_digests_audited":len(receipt["original_artifact_source_file_sha256"]),
      "original_TP_observation_rows":tp["observed_TP_positive_values_2011_2022"],
      "reconstructed_matched_TP_month_rows":sum(map(len,monthly.values())),
      "GBIF_plant_sites":len(fac),
      "new_plant_site_species_records":len(species_table),
      "source_outfalls":len(original),
      "strict_plant_supported":count30,
      "sensitive_plant_supported":count10,
      "TP_vs_water_fraction_spearman":source_rho,
      "TP_vs_water_fraction_two_sided_exact_p":exact_p,
      "TP_vs_plants_richness_depth10_spearman":flora_rho,
      "new_external_requests":0,
    }
if __name__=="__main__":
    print(json.dumps(verify(),indent=2,sort_keys=True))
