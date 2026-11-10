#!/usr/bin/env python3
"""Posthoc HYDCLS inland vs tidal and conductivity sensitivity after QA.

Same original 568 source-matched USA NWCA2021 2021 probability first-visit
wetlands and 499 quality valid actual NTL/PTL/H_ALL source response. Splits
were chosen AFTER seeing overall negative nutrient-plant Shannon relation.
Retain every direction, no sewage source/cause inferred.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from gesui_us_NWCA2021_actual_site_plant_Shannon_FQAI_TP_TN_after_source_gate import (
    FILES,select_and_join,valid_rows,fitted)

MAIN=Path("results/us_NWCA2021_plant_TNTP_20261110/EPA_NWCA2021_REAL_PLANT_SHANNON_TNTP_WETLANDS_STATE_ECOREGION_QA.json")
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(root,data_dir):
    out=json.loads((root/MAIN).read_text())
    src,_,audit=select_and_join(root,{k:data_dir/v for k,v in FILES.items()})
    qual=valid_rows(src)
    if len(qual)!=out["source_joint_TNTP_and_Shannon_nonmissing_quality_sites"] or len(qual)!=499:
        raise ValueError("Original 499/568 source cohort changed in hydrologic sensitivity")
    if Counter(qual.HYDCLS)!=Counter({"INLAND":369,"TIDAL":130}):
        raise ValueError("Original site hydrologic classes changed after source audit")
    results={}
    for name in ("INLAND","TIDAL"):
        d=qual.loc[qual.HYDCLS==name].copy()
        b=fitted(d,"HALL","log10_TP_mgL+log10_TN_mgL")
        p=fitted(d,"HALL","log10_TP_mgL")
        n=fitted(d,"HALL","log10_TN_mgL")
        if d.PSTLCODE.nunique()<8:raise ValueError("Hydrologic stratum too few independent US states")
        results[name]={
          "N_actual_source_joint_quality_sites":len(d),
          "number_US_states":int(d.PSTLCODE.nunique()),
          "number_ecoregions":int(d.AGECO9.nunique()),
          "TN_and_TP_concurrent_source":b,
          "TN_only_observational_source":n,
          "TP_only_observational_source":p,
          "unadjusted_TN_v_Shannon_rho":float(spearmanr(d.NTLRESULT,d.HALL).statistic),
          "unadjusted_TP_v_Shannon_rho":float(spearmanr(d.PTLRESULT,d.HALL).statistic),
          "observational_subpopulation_not_national_sewage_treatment":True
        }
    result={
      "schema":"gesui_US_NWCA2021_499_real_water_NP_plant_Shannon_inland_tidal_COND_after_original_main_v01",
      "UTC":datetime.now(timezone.utc).isoformat(),
      "original_499wetland_TNTP_and_plant_Shannon_model_sha256":sha(root/MAIN),
      "original_source_gate_568_unmodified":True,
      "posthoc_hydrological_class_sensitivity_AFTER_original_negative_nutrient_plant_result":True,
      "original_499_nutrient_quality_sites_reproduced":True,
      "inland_v_tidal_nutrient_Shannon":results,
      "not_measured_WWTP_origin_at_each_wetland":True
    }
    cond=pd.to_numeric(qual.CONDRESULT,errors="coerce") if "CONDRESULT" in qual else pd.Series(np.nan,index=qual.index)
    m=qual.loc[cond.gt(0)].copy()
    if len(m)>=350 and m.PSTLCODE.nunique()>=20:
        m["log10_COND_uScm"]=np.log10(pd.to_numeric(m.CONDRESULT))
        result["all_wetlands_TNTP_plus_actual_specific_conductivity"]={
            "original_chemistry_source_supported_sites_n":len(m),
            "original_states":int(m.PSTLCODE.nunique()),
            "primary_TNTP_with_actual_conductivity":fitted(m,"HALL","log10_TP_mgL+log10_TN_mgL+log10_COND_uScm"),
            "conductivity_is_water_ion_proxy_not_sewage_dose":True
        }
    else:
        result["all_wetlands_TNTP_plus_actual_specific_conductivity"]={
            "original_chemistry_source_supported_sites_n":len(m),
            "STOP_insufficient_cond_data":True
        }
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-dir",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    x=run(root,a.source_dir)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(x,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({
        "hydro":{k:{"n":v["N_actual_source_joint_quality_sites"],"state":v["number_US_states"],
                    "TNandTP":v["TN_and_TP_concurrent_source"]["coefficients"],
                    "rhoTN":v["unadjusted_TN_v_Shannon_rho"]}
                 for k,v in x["inland_v_tidal_nutrient_Shannon"].items()},
        "conductivity":x["all_wetlands_TNTP_plus_actual_specific_conductivity"]
    },indent=2))

if __name__=="__main__":main()
