#!/usr/bin/env python3
"""Feurs STEU-year two-outlet SOURCE-only gate; no biological response or causal mass balance."""
import argparse, json, math
from pathlib import Path

P_ATOMIC_MASS=30.973761998
O_ATOMIC_MASS=15.999
P_OVER_P2O5=(2*P_ATOMIC_MASS)/(2*P_ATOMIC_MASS+5*O_ATOMIC_MASS)

def audit(source):
    if source.get("schema") != "gesui_feurs_primary_source_2018_2019_2021_NO_biological_outcomes_v1":
        raise ValueError("Unapproved SOURCE schema")
    if source.get("STEU_SANDRE") != "0442094S0002":
        raise ValueError("Foreign STEU cannot be mixed into Feurs validation")
    if sorted(v["year"] for v in source["years"]) != [2018, 2019, 2021]:
        raise ValueError("SOURCE records expected for 2018, 2019 and 2021")
    out={}
    for r in source["years"]:
        y=r["year"]
        if r.get("same_year_outlet_mass_pair_observable") is not False:
            raise ValueError("No complete same-year mass pair in available SOURCE")
        if r.get("true_measured_annual_effluent_total_N_t") is not None or r.get("true_measured_annual_effluent_total_P_t") is not None:
            raise ValueError("Unsupported claim of measured annual effluent N/P loads")
        year={"station":source["STEU_SANDRE"],"source_year":y,
              "true_same_year_water_and_field_nutrient_mass_pair":False,
              "source_only_NO_biological_values_read":True}
        if r["actual_sludge_spreading_confirmed"]:
            if y!=2018:
                raise ValueError("Realized parcel application only confirmed for 2018")
            for fld in ["actual_spreading_ha","actual_spreading_dry_t_including_lime","mean_sludge_N_total_pct_of_limed_dry_matter","mean_sludge_P2O5_pct_of_limed_dry_matter"]:
                if r[fld]<=0 or not math.isfinite(r[fld]):raise ValueError("Invalid sludge SOURCE")
            ms=r["actual_spreading_dry_t_including_lime"]
            n=ms*r["mean_sludge_N_total_pct_of_limed_dry_matter"]/100
            p2o5=ms*r["mean_sludge_P2O5_pct_of_limed_dry_matter"]/100
            n_by_area=r["mean_applied_N_units_per_ha"]*r["actual_spreading_ha"]/1000
            year.update({"sludge_actual_spreading_confirmed":True,
                 "applied_biosolid_dry_mass_with_lime_t":ms,
                 "sludge_N_total_assay_mean_times_applied_dry_mass_t_SCREENING_ONLY":round(n,6),
                 "sludge_P2O5_assay_mean_times_applied_dry_mass_t_SCREENING_ONLY":round(p2o5,6),
                 "sludge_elemental_P_equivalent_t_SCREENING_ONLY":round(p2o5*P_OVER_P2O5,6),
                 "independent_sludge_N_area_times_mean_reported_units_t_SCREENING_ONLY":round(n_by_area,6),
                 "N_estimate_area_vs_assay_percent_difference":round(100*(n_by_area/n-1),3),
                 "leachate_N_P_additional_unquantified":bool(r["leachate_applied_m3_not_included_in_sludge_mass"]>0),
                 "P2O5_not_elemental_P":True,
                 "do_not_interpret_as_measured_total_mass_or_plant_uptake":True
                 })
        else:
            if y not in (2019,2021):raise ValueError("Unexpected source year")
            for fld in ["treated_water_m3","mean_effluent_NTK_mg_L","mean_effluent_total_P_mg_L"]:
                if r[fld]<=0 or not math.isfinite(r[fld]):raise ValueError("Invalid water SOURCE")
            if r.get("effluent_mean_is_flow_weighted") is not False:
                raise ValueError("Do not misclassify unweighted outlet concentration")
            factor=r["treated_water_m3"]/1_000_000  # mg/L * m³ / 10^6 = tonnes
            year.update({"sludge_actual_spreading_confirmed":False,
              "sludge_dry_produced_t_including_lime_not_field_application":r.get("sludge_dry_produced_t_including_lime"),
              "treated_water_m3":r["treated_water_m3"],
              "inlet_bypass_m3_not_included":r["inlet_bypass_m3"],
              "unweighted_mean_NTK_x_annual_flow_t_SCREENING_INDEX_NOT_N_LOAD":round(r["mean_effluent_NTK_mg_L"]*factor,6),
              "unweighted_mean_PT_x_annual_flow_t_SCREENING_INDEX_NOT_P_LOAD":round(r["mean_effluent_total_P_mg_L"]*factor,6),
              "NTK_is_not_total_N":True,
              "bypass_not_included_in_water_screening_index":True,
              "unweighted_mean_not_flow_weighted":True})
        out[str(y)]=year
    return {"schema":"gesui_Feures_year_locked_source_only_mass_feasibility_v1",
            "STEU_SANDRE":source["STEU_SANDRE"],
            "P_over_P2O5_conversion_factor":round(P_OVER_P2O5,9),
            "source_year_records":out,
            "source_years_with_measured_both_aqueous_and_field_mass":[],
            "strict_paired_export_mass_GATE":"STOP",
            "strict_biological_path_GATE":"STOP",
            "nonmatching_2018_2019_2021_years_never_combined":True,
            "report_masses_not_plant_uptake":True,
            "reported_land_2018_application_NOT_equal_to_2019_or_2021_sludge_production":True}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    result=audit(json.loads(a.source.read_text(encoding="utf-8")))
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("FEURS_DUAL_OUTLET_SOURCE_ONLY",json.dumps(result,sort_keys=True))
if __name__=="__main__":main()
