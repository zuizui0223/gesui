#!/usr/bin/env python3
"""Build exact genuine 183 unique H7 × 40 native H6 terrestrial source panel.

Source-locked OLD 234×2×6 GBIF facets, NEW recovered 30×2×6 GBIF facets,
EIV 4011 accepted plant scores, HNS 186 accepted butterfly scores,
genuine 1219-H7 Animalia index, and exact prior HydroATLAS/W replay.

Any original old HNS or plant NUTRI trait value mismatch STOPS fit. New
H7 EEA physical survey and country×method information never fabricated.
"""
from __future__ import annotations
import argparse,hashlib,json,math,csv
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
PANEL_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
OLD_FACETS_SHA="251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4"
NEW_FACETS_SHA="e076f6c8e2d93282be19aea3c75ae8431302318f068b9d75e278e2a21c20aafb"
EIV_SHA="ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba"
HNS_SHA="00c1a36a9237a3d1a6d77334aa8baa5dc8b329a49925c2e31f033fa2d302ecdb"
SOURCE_COHORT="GESUI_TERRESTRIAL_ONLY_NEW_WITHIN_H6_ORIGINAL_GEO_COHORT_BEFORE_OUTCOMES_V1"
SOURCE_ENV="GESUI_20261115_H7_1219_NATIVE_SOURCE_W_ENVIRONMENT_REPLAY_V1"
FAMS=["Hesperiidae","Lycaenidae","Nymphalidae","Papilionidae","Pieridae"]
PERIODS=["2017-2019","2020-2022"]
METRICS=["plant_n_incidence","plant_n_record","hns_incidence","hns_record","nPlant","nButterfly",
         "nAnimal","plant_coverage","hns_coverage","plant_scored_species",
         "hns_scored_species","hns_host_completeness"]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def one(folder,name):
    if not isinstance(folder,Path):folder=Path(folder)
    f=list(folder.rglob(name))
    if len(f)!=1:raise RuntimeError("Missing or duplicate genuine source file "+name+": "+str(f))
    return f[0]
def ascsv(path,dtype=None):
    return pd.read_csv(path,encoding="utf-8-sig",dtype=dtype)
def load_h7_facets(path,expected_sha,n):
    if sha(path)!=expected_sha:raise RuntimeError("Original source JSON sha wrong "+path.name)
    raw=json.loads(path.read_text())
    if len(raw)!=n:raise RuntimeError("Original facet count wrong "+path.name)
    table={}
    for r in raw:
        k=(str(r["HYBAS_ID_L07"]),str(r["period"]),str(r["group"]))
        if k in table:raise RuntimeError("Original GBIF source key duplicated "+str(k))
        if r.get("source_status")!="ok" or not isinstance(r.get("source_speciesKey_counts"),dict):
            raise RuntimeError("Original GBIF source facet still missing "+str(k))
        if not isinstance(r.get("source_query_occurrence_count"),int):
            raise RuntimeError("Original GBIF record count invalid "+str(k))
        table[k]=r
    return table

def main(a):
    co=json.loads(one(a.cohort,"FIXED_NEW_40H6_SOURCE_COHORT_GATE.json").read_text())
    if (co.get("schema")!=SOURCE_COHORT or not co.get("source_support_gate_pass")
         or co.get("source_final_true_H6")!=40 or co.get("source_final_distinct_H7")!=183):
        raise RuntimeError("Original 183-H7/40-H6 source-only cohort not frozen")
    ids=ascsv(one(a.cohort,"FIXED_183_H7_IN_40_NATIVE_H6_SOURCE_COHORT.csv"),
              dtype={"h":str,"h6":str})
    if len(ids)!=183 or ids.h.nunique()!=183 or ids.h6.nunique()!=40:
        raise RuntimeError("Original source cohort identity matrix changed")
    orig=Path("results/20261115_sync/hns_macro/H7_454_HNS_panel.csv")
    if sha(orig)!=PANEL_SHA:raise RuntimeError("Original published 454 H7 panel sha mismatch")
    old=ascsv(orig,dtype={"h":str,"h6":str})
    if len(old)!=454 or old.h.nunique()!=232 or old.h6.nunique()!=90:
        raise RuntimeError("Original published 454 source cohort changed")
    em=json.loads(one(a.environment,"H7_1219_SOURCE_REPLAY_W_ENVIRONMENT_QA.json").read_text())
    if em.get("schema")!=SOURCE_ENV or not em.get("strict_new_model_environment_source_gate_pass"):
        raise RuntimeError("Original source 232-H7 W/crop/urban/environment replay did not pass")
    new_env=ascsv(one(a.environment,"NEW_TEN_H6_DETERMINISTIC_SOURCE_ONE_H7_EACH.csv"),
                  dtype={"native_H7":str,"new_official_H6":str})
    if len(new_env)!=10:raise RuntimeError("Original new 10 selected environmental ID missing")
    envlookup={r.native_H7:r for r in new_env.itertuples(index=False)}
    oldfac=load_h7_facets(one(a.oldfacets,"GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json"),OLD_FACETS_SHA,2808)
    newfac=load_h7_facets(one(a.newfacets,"MERGED_ORIGINAL_30H7_TWO_PERIOD_SIX_FAMILY_360_SOURCE_FACETS.json"),NEW_FACETS_SHA,360)
    if set(oldfac)&set(newfac):raise RuntimeError("Old and new source H7 facets overlap, sample leakage")
    full={**oldfac,**newfac}
    expected={(h,p,f) for h in ids.h for p in PERIODS for f in ["Plantae"]+FAMS}
    if not expected.issubset(full):raise RuntimeError("Original source H7×period×6 groups missing")
    eiv=json.loads(one(a.eiv,"EU454_FULL_EIV_SOURCE_GATE.json").read_text())
    if (eiv.get("schema")!="GESUI_EU454_FULLSOURCE_PLANT_NUTRI_TRAIT_GATE_V2_HYBRID_QA" or
        not eiv.get("hard_source_support_pass") or
        eiv.get("source_EIV",{}).get("original_source_SHA256")!=EIV_SHA):
        raise RuntimeError("Original fully audited 4011 EIV scored flora source missing")
    plants=ascsv(one(a.eiv,"EU454_PLANT_11511_SOURCE_EIV_EXACT_CROSSWALK.csv"),dtype={"speciesKey":str})
    if len(plants)!=11511 or plants.speciesKey.duplicated().any():
        raise RuntimeError("Original taxonomy species roster changed")
    plantscore={r.speciesKey:float(r.EIV_NUTRI_exact) for r in plants.itertuples(index=False)
                if pd.notna(r.EIV_NUTRI_exact)}
    if len(plantscore)!=4011:raise RuntimeError("Original 4011 exact numeric plant N scores changed")
    scores_path=Path("results/20261115_sync/GESUI_303_HOST_NUTRIENT_SCORES.csv")
    if sha(scores_path)!=HNS_SHA:raise RuntimeError("Original immutable 303 butterfly HNS roster checksum drift")
    butterflies=ascsv(scores_path,dtype={"butterflyKey":str})
    if len(butterflies)!=303 or butterflies.butterflyKey.duplicated().any():
        raise RuntimeError("Original 303 butterfly roster invalid")
    aliases=butterflies.actual_used_alias_conflict.astype(str).str.lower().isin(["true","1","yes"])
    accepted=butterflies[(butterflies.taxon_status.astype(str)=="ACCEPTED")&
             ~aliases&butterflies.host_nutrient_mean.notna()]
    bp={str(z.butterflyKey):{"N":float(z.host_nutrient_mean),
       "family":str(z.family),"completeness":float(z.trait_coverage)}
       for z in accepted.itertuples(index=False)}
    if len(bp)!=186:raise RuntimeError("Original score-eligible HNS butterfly count not 186")
    master=ascsv(one(a.master,"HYBAS7_1219x2_GBIF_PLANTAE_LYCAENIDAE_ANIMALIA_HYDROWASTE_SOURCE_ONLY_MASTER.csv"),dtype={"HYBAS_ID_L07":str})
    if len(master)!=2438 or master.duplicated(["HYBAS_ID_L07","PERIOD"]).any():
        raise RuntimeError("Original genuine 1219×2 source master not unique")
    animal={(z.HYBAS_ID_L07,z.PERIOD):z.GBIF_INDEX_ANIMALIA_2017_2022_FIXED_PERIOD_RECORD_COUNT
            for z in master.itertuples(index=False)}
    oldplant=ascsv(one(a.eiv,"EU454_H7PERIOD_NUTRI_TRAIT_MEANS_SOURCEONLY.csv"),dtype={"h":str})
    if len(oldplant)!=451 or oldplant.duplicated(["h","period"]).any():
        raise RuntimeError("Published original EIV H7period traits changed")
    oldplant_map={(z.h,z.period):z for z in oldplant.itertuples(index=False)}
    oldkey={(z.h,z.period):z for z in old.itertuples(index=False)}
    site=[]
    qa={"old_HNS_diff_incidence":0.,"old_HNS_diff_record":0.,"old_plant_N_diff_incidence":0.,
        "old_plant_N_diff_record":0.,"old_hns_coverage_diff":0.,
        "old_hns_completeness_diff":0.}
    scorefail=[]
    for fixed in ids.itertuples(index=False):
        h=str(fixed.h)
        dperiod=[]
        for p in PERIODS:
            facets={f:full[(h,p,f)] for f in ["Plantae"]+FAMS}
            pl_counts={str(k):int(v) for k,v in facets["Plantae"]["source_speciesKey_counts"].items() if int(v)>0}
            plnum=[(plantscore[k],v) for k,v in pl_counts.items() if k in plantscore]
            if len(plnum)<10 or len(plnum)/max(1,len(pl_counts))<.4:
                raise RuntimeError("Frozen original source EIV plant coverage failed for H7="+h+" period="+p)
            raw_five=sum(facets[f]["source_query_occurrence_count"] for f in FAMS)
            bf=[]
            for f in FAMS:
                for k,v in facets[f]["source_speciesKey_counts"].items():
                    if int(v)<=0:continue
                    if k in bp:
                        if bp[k]["family"]!=f:raise RuntimeError("Genuine HNS family taxonomy mismatch")
                        bf.append((k,int(v),bp[k]["N"],bp[k]["completeness"]))
            if not bf or raw_five<1:raise RuntimeError("HNS source species counts empty")
            scoredcount=sum(q[1] for q in bf)
            rawanimal=animal.get((h,p))
            if pd.isna(rawanimal) or rawanimal is None or float(rawanimal)<0:
                raise RuntimeError("Authentic 1219 original Animalia source missing; no zero-imputation H7="+h)
            rec={"h":h,"h6":str(fixed.h6),"period":p,"source":fixed.source,
                "plant_n_incidence":float(np.mean([q[0] for q in plnum])),
                "plant_n_record":float(np.average([q[0] for q in plnum],weights=[q[1] for q in plnum])),
                "hns_incidence":float(np.mean([q[2] for q in bf])),
                "hns_record":float(np.average([q[2] for q in bf],weights=[q[1] for q in bf])),
                "nPlant":float(np.log1p(facets["Plantae"]["source_query_occurrence_count"])),
                "nButterfly":float(np.log1p(raw_five)),
                "nAnimal":float(np.log1p(float(rawanimal))),
                "plant_coverage":len(plnum)/len(pl_counts),
                "hns_coverage":scoredcount/raw_five,
                "plant_scored_species":len(plnum),
                "hns_scored_species":len(bf),
                "hns_host_completeness":sum(q[1]*q[3] for q in bf)/scoredcount}
            if fixed.source=="original_old_H7":
                prior=oldkey[(h,p)]
                origp=oldplant_map[(h,p)]
                for name,newv,priorv in [
                    ("old_HNS_diff_incidence",rec["hns_incidence"],prior.primary_incidence_mean),
                    ("old_HNS_diff_record",rec["hns_record"],prior.primary_record_mean),
                    ("old_plant_N_diff_incidence",rec["plant_n_incidence"],origp.Plant_N_incidence),
                    ("old_plant_N_diff_record",rec["plant_n_record"],origp.Plant_N_record),
                    ("old_hns_coverage_diff",rec["hns_coverage"],prior.primary_record_coverage_raw),
                    ("old_hns_completeness_diff",rec["hns_host_completeness"],prior.primary_host_completeness)]:
                    delta=abs(float(newv)-float(priorv))
                    qa[name]=max(qa[name],delta)
                    if delta>1e-10:scorefail.append({"H7":h,"period":p,"variable":name,"difference":delta})
            dperiod.append(rec)
        if len(dperiod)!=2 or {v["period"] for v in dperiod}!=set(PERIODS):
            raise RuntimeError("Not both source periods")
        if fixed.source=="original_old_H7":
            old_row=oldkey[(h,PERIODS[0])]
            env={k:float(getattr(old_row,k)) for k in ("W","area","crop","urban","naturalout")}
        else:
            if h not in envlookup:raise RuntimeError("New 10 H7 not in strict original source replay")
            new=envlookup[h]
            env={k:float(getattr(new,k)) for k in ("W","area","crop","urban","naturalout")}
        q={"h":h,"h6":str(fixed.h6),"original_source_cohort":str(fixed.source),
           **env}
        for k in METRICS:
            q[k]=float(sum(r[k] for r in dperiod)/2.)
        site.append(q)
    if scorefail:raise RuntimeError("Existing original old source score index did not reproduce published results: "+json.dumps(scorefail[:12]))
    if len(site)!=183 or len({z["h"] for z in site})!=183 or len({z["h6"] for z in site})!=40:
        raise RuntimeError("Original pre-outcome 183 H7 /40 H6 cohort changed")
    if not all(np.isfinite([z[k] for z in site]).all() for k in ("W","area","crop","urban","naturalout",*METRICS)):
        raise RuntimeError("Nonfinite measured source recording or trait summary")
    # Count only source values already selected by original station ranks,
    # never response direction. Fixed frame and original H6 parent remains.
    out={
      "schema":"GESUI_ORIGINAL_NEW_183H7_40H6_TERRESTRIAL_SOURCE_COHORT_AND_TRAIT_REPLAY_V1",
      "old_native_panel_sha256":PANEL_SHA,
      "old_2808_facets_sha256":OLD_FACETS_SHA,
      "extra_360_facets_sha256":NEW_FACETS_SHA,
      "original_EIV_source_sha256":EIV_SHA,
      "original_HNS_303_roster_sha256":HNS_SHA,
      "new_10H7_environment_source_exact_232_replay":True,
      "n_distinct_original_H7":183,"n_true_H6":40,
      "n_old_H7":sum(z["original_source_cohort"]=="original_old_H7" for z in site),
      "n_new_H7":sum(z["original_source_cohort"]!="original_old_H7" for z in site),
      "source_taxon_genuine_2period_required":True,
      "original_old_326_H7period_source_trait_replay_max_abs_error":qa,
      "original_old_score_mismatch_n":len(scorefail),
      "source_1219_animalia_count_available_all":True,
      "original_third_party_EIV_full_species_values_not_republished":True,
      "realized_plant_herbivory_and_local_sewage_N_NOT_measured":True,
      "no_same_EEA_bothperiod_station_data_for_new_H7":True
    }
    a.output.mkdir(parents=True,exist_ok=True)
    p=a.output/"GENUINE_TERRESTRIAL_ONLY_183_H7_40_H6_PAIRED_H7_MEANS.csv"
    pd.DataFrame(site).to_csv(p,index=False)
    out["data_source_panel_SHA256"]=sha(p)
    (a.output/"GENUINE_TERRESTRIAL_183_H7_40_H6_SOURCE_RECONSTRUCTION_QA.json").write_text(
        json.dumps(out,indent=2,ensure_ascii=False)+"\n")
    print("GESUI_183_H7_TERRESTRIAL_SOURCE_REPLAY="+json.dumps(out,sort_keys=True,ensure_ascii=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--cohort",type=Path,required=True)
    p.add_argument("--environment",type=Path,required=True)
    p.add_argument("--oldfacets",type=Path,required=True)
    p.add_argument("--newfacets",type=Path,required=True)
    p.add_argument("--eiv",type=Path,required=True)
    p.add_argument("--master",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
