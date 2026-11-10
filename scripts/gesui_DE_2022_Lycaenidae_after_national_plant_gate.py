#!/usr/bin/env python3
"""2022 DE river Lycaenidae photographic reporting, only AFTER source/Plantae gates pass.

This is POSTHOC relative to 2022 1degree positives and NONCAUSAL.
No query is performed if original 2019-2021 plant-photo effort or geography
gates fail for both 2km and 5km. No individual photographs/observer IDs used.
"""
from __future__ import annotations
import argparse,hashlib,itertools,json,math,time
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.parse import urlencode
from urllib.request import Request,urlopen

FRAME="results/lycaenidae_2022_DE_national_73pair_sourceonly/ORIGINAL_2022_73DE_WISE6_SOURCEZERO_PAIR_100SITE_PREPHOTO_FRAME.json"
PLANTS="results/lycaenidae_2022_DE_national_73pair_sourceonly/ORIGINAL_100STATION_2KM_5KM_2019_2021_PLANTAE_EFFORT_BEFORE_2022_BUTTERFLIES.json"
TAXA={"Lycaenidae":47923,"Pieridae":48508,"Papilionidae":47223}
START,END="2022-04-01","2022-10-31"
SRC="registered_source_interposed"
ZERO="registered_qualifying_sourcezero"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def query(site,rad,family,tid):
    args={"taxon_id":tid,"place_id":7207,"lat":round(site["lat"],5),
          "lng":round(site["lon"],5),"radius":rad,"d1":START,"d2":END,
          "photos":"true","geo":"true","geoprivacy":"open",
          "taxon_geoprivacy":"open","quality_grade":"research",
          "captive":"false","acc_below":500,"per_page":1}
    url="https://api.inaturalist.org/v1/observations?"+urlencode(args)
    error=None
    for trial in range(4):
        try:
            req=Request(url,headers={"User-Agent":"gesui-national-DE-river-butterfly-photo-after-effort-gate/0.1",
                                     "Accept":"application/json"})
            with urlopen(req,timeout=45) as resp:
                raw=resp.read()
                status=resp.status
            doc=json.loads(raw)
            n=doc["total_results"]
            if status!=200 or type(n) is not int or n<0:
                raise ValueError("Response lacked valid nonnegative total_results")
            return {"n_public_photo_reports":n,"status":200,
                    "original_API_url":url,"original_JSON_sha256":hashlib.sha256(raw).hexdigest()}
        except (HTTPError,URLError,ValueError,TimeoutError) as e:
            error=str(e)
            if trial<3:time.sleep(5*(trial+1))
    return {"n_public_photo_reports":None,"status":"FAILED",
            "error":error,"original_API_url":url}

def pairIDs(p):
    return p["upstream_original_site_id_SHA256"],p["downstream_original_site_id_SHA256"]

def score(id,rec):
    vals=[rec[f]["n_public_photo_reports"] for f in TAXA]
    return sum(vals) if all(type(x) is int for x in vals) else None

def logit(id,rec):
    lyc=rec["Lycaenidae"]["n_public_photo_reports"]
    total=score(id,rec)
    return math.log((lyc+.5)/(total-lyc+.5))

def exact_signflip(x):
    n=len(x);v=sum(x)/n
    enum=2**n
    if n>18: return {"mean_difference":v,"n":n,"tail":"not_run_gt_18_quartets"}
    count=0
    for signs in itertools.product((-1,1),repeat=n):
        s=sum(a*b for a,b in zip(signs,x))/n
        if abs(s)>=abs(v)-1e-12: count+=1
    return {"mean_log_odds_difference":v,"n":n,"positive_quartets":sum(z>0 for z in x),
            "exact_two_sided_signflip_reference_fraction":count/enum,
            "not_causal_randomized_treatment_p":True}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    root=Path(__file__).resolve().parents[1]
    frame=json.loads((root/FRAME).read_text())
    plants=json.loads((root/PLANTS).read_text())
    if plants["original_73pair_100site_source_freeze_SHA256"]!=sha(root/FRAME):
        raise ValueError("Cannot open butterfly outcomes: 73-pair source frame changed")
    if plants["no_2022_butterfly_photo_counts_queried"] is not True:
        raise ValueError("Prebutterfly plant stage was contaminated by butterfly responses")
    if plants["new_station_2019_2021_plant_photo_queries_planned"]!=136:
        raise ValueError("Unexpected original 68-site source-frame opportunity size")
    pairs={p["pair_id_SHA256"]:p for p in frame["complete_original_fixed_candidate_pairs"]}
    stations={p["original_site_id_SHA256"]:p for p in frame["original_unique_public_WISE6_monitor_site_metadata"]}
    out={"schema":"gesui_DE_national_sourcelocked_2022_Lycaenidae_photoreporting_after_priorplant_gates_v01",
         "utc":datetime.now(timezone.utc).isoformat(),
         "source_frame_SHA256":sha(root/FRAME),"2019_2021_plant_gate_SHA256":sha(root/PLANTS),
         "explicitly_posthoc_after_earlier_2022_landscape_positive":True,
         "no_individual_butterfly_abundance_or_puddling_or_plant_fitness_measured":True,
         "reported_EEA_2022_P_and_2022_photo_timing_are_not_randomized_treatment":True,
         "actual_public_2022_butterfly_photo_API_queries":0,"by_radius":{}}
    for rad in (2,5):
        gate=plants["by_radius"][f"{rad}km"]
        match_go=gate["three_degree_four_arm_design_support_GATE"]["matched_2022_butterfly_photo_outcomes_permitted"]
        national_go=gate["national_pair_design_support_GATE"]["national_2022_butterfly_photograph_outcomes_permitted"]
        record={"preoutcome_fourarm_gate_pass":match_go,
                "preoutcome_national_directional_gate_pass":national_go,
                "no_2022_butterfly_query_if_both_gates_fail":True}
        if not (match_go or national_go):
            record["status"]="STOP_PREBUTTERFLY_GATES_FAILED"
            record["query_count"]=0
            out["by_radius"][f"{rad}km"]=record
            continue
        requested=set()
        if match_go:
            for m in gate["original_3deg_source_vs_registeredzero_matched_fourarm_quartets"]:
                requested.update(m["four_station_hashes"])
        if national_go:
            for p in gate["fully_site_disjoint_country_wide_directional_pairs"]:
                requested.update(pairIDs(pairs[p["pair_id_SHA256"]]))
        receipts={}
        for i,stid in enumerate(sorted(requested)):
            site=stations[stid]
            receipts[stid]={}
            for family,tid in TAXA.items():
                receipts[stid][family]=query(site,rad,family,tid)
                out["actual_public_2022_butterfly_photo_API_queries"]+=1
                time.sleep(1.15)
            print(f"[{rad}km] 2022 3family source+outcome photo receipts {i+1}/{len(requested)}",flush=True)
        record["original_2022_public_threefamily_photo_receipts"]=receipts
        record["unique_original_source_station_photo_queries"]=len(requested)
        record["query_count"]=len(requested)*len(TAXA)
        if match_go:
            eligible=[]
            for quartet in gate["original_3deg_source_vs_registeredzero_matched_fourarm_quartets"]:
                sites=quartet["four_station_hashes"]
                if not all(s in receipts and score(s,receipts[s]) is not None and score(s,receipts[s])>=5 for s in sites):
                    continue
                source_pair=pairs[quartet["source_pair_SHA256"]]
                zero_pair=pairs[quartet["registeredzero_pair_SHA256"]]
                su,sd=pairIDs(source_pair)
                cu,cd=pairIDs(zero_pair)
                delta_src=logit(sd,receipts[sd])-logit(su,receipts[su])
                delta_zero=logit(cd,receipts[cd])-logit(cu,receipts[cu])
                eligible.append({"source_pair_SHA256":source_pair["pair_id_SHA256"],
                                 "registeredzero_pair_SHA256":zero_pair["pair_id_SHA256"],
                                 "original_3deg_block":quartet["original_3deg_block"],
                                 "fourarm_Lycaenidae_relative_photo_logodds_difference":delta_src-delta_zero,
                                 "source_down_minus_up":delta_src,"registeredzero_down_minus_up":delta_zero})
            eligible_blocks={q["original_3deg_block"] for q in eligible}
            record["fourarm_photo_supported_quartets_n"]=len(eligible)
            record["fourarm_photo_supported_3deg_regions"]=len(eligible_blocks)
            record["fourarm_2022_butterfly_photo_response_numeric_gate_pass"]=len(eligible)>=6 and len(eligible_blocks)>=3
            if record["fourarm_2022_butterfly_photo_response_numeric_gate_pass"]:
                arr=[q["fourarm_Lycaenidae_relative_photo_logodds_difference"] for q in eligible]
                record["fourarm_directional_photo_response"]=exact_signflip(arr)
                record["fourarm_directional_photo_response"]["leave_one_quartet_out_mean_range"]=[
                    min(sum(x for j,x in enumerate(arr) if j!=i)/(len(arr)-1) for i in range(len(arr))),
                    max(sum(x for j,x in enumerate(arr) if j!=i)/(len(arr)-1) for i in range(len(arr)))
                ]
                record["original_photo_supported_qualified_quartet_records"]=eligible
            else:
                record["fourarm_directional_photo_response"]="NOT_EVALUATED_PHOTO_SUPPORT_GATE_FAILED"
        if national_go:
            eligible=[]
            for p in gate["fully_site_disjoint_country_wide_directional_pairs"]:
                meta=pairs[p["pair_id_SHA256"]]
                u,d=pairIDs(meta)
                if not all(s in receipts and score(s,receipts[s]) is not None and score(s,receipts[s])>=5 for s in (u,d)):
                    continue
                eligible.append({"pair_id_SHA256":meta["pair_id_SHA256"],
                                 "class":meta["exposure_class"],
                                 "original_2deg_block":meta["same_2deg_upstream_region"],
                                 "original_3deg_block":meta["upstream_3deg_region"],
                                 "down_minus_up_Lycaenidae_photo_logit":logit(d,receipts[d])-logit(u,receipts[u])})
            bytype=Counter(x["class"] for x in eligible)
            regions=defaultdict(set)
            for x in eligible:regions[x["original_3deg_block"]].add(x["class"])
            record["national_photo_qualified_pair_counts"]=dict(bytype)
            record["national_photo_qualified_3deg_blocks_both_classes"]=sum(len(v)==2 for v in regions.values())
            accepted=bytype[SRC]>=8 and bytype[ZERO]>=8 and sum(len(v)==2 for v in regions.values())>=3
            record["national_2022_directional_photo_response_gate_pass"]=bool(accepted)
            if accepted:
                values={c:[v["down_minus_up_Lycaenidae_photo_logit"] for v in eligible if v["class"]==c] for c in (SRC,ZERO)}
                record["national_directional_photo_reporting_mean_source"]=sum(values[SRC])/len(values[SRC])
                record["national_directional_photo_reporting_mean_registeredzero"]=sum(values[ZERO])/len(values[ZERO])
                record["nationwide_source_minus_registeredzero_mean_photo_logit_gradient"]=(record["national_directional_photo_reporting_mean_source"]-record["national_directional_photo_reporting_mean_registeredzero"])
                record["national_noncausal_exploratory_independent_pairs"]=eligible
            else:
                record["national_directional_photo_reporting_mean_source_minus_zero"]="NOT_EVALUATED_PHOTO_SUPPORT_GATE_FAILED"
        record["status"]="2022_PHOTO_SOURCE_QUERIED_ONLY_AFTER_PLANT_SUPPORT_GATE"
        out["by_radius"][f"{rad}km"]=record
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"queries":out["actual_public_2022_butterfly_photo_API_queries"],
                      "radii":{r:{
                          "status":v["status"],
                          "photo_outcome_query_count":v["query_count"],
                          "fourarm_inference_gate":v.get("fourarm_2022_butterfly_photo_response_numeric_gate_pass"),
                          "national_photo_inference_gate":v.get("national_2022_directional_photo_response_gate_pass"),
                      } for r,v in out["by_radius"].items()}},indent=2))

if __name__=="__main__":
    main()
