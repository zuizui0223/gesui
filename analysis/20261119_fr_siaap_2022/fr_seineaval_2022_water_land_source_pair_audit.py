#!/usr/bin/env python3
"""Compare two independently source-frozen 2022 Seine Aval nutrient outflows.

Inputs are official AESN 2022 A4 validated MONTHLY effluent kg indicators,
and SIAAP 2022 actually SPREAD campaign lot-weighted NTK/P2O5 estimates.
Both belong to Seine Aval WWTP and calendar 2022. This is an accounting
*comparison* NOT a closed 2022 source partition: sludge produced/delivered
in 2021 can be applied during 2022; treated A4 and bypass A5 separate; nitrogen
forms NTK vs NGL differ; fields far from water-receiving site; EQR/NUTRI/HNS
not observed here. Never turn absence of a source into zero.
"""
import argparse,json
from decimal import Decimal
from pathlib import Path

PIN_AESN_SHA="5fbed27273b0c2c1d4590e9286fa52e724578b5cb19e68369fd1fb8e77347af7"
AESN_RUN=44538504180
SOIL_RUN=44538061031
P_ATOMIC=Decimal("30.973761998")
O_ATOMIC=Decimal("15.999")
P_FACTOR=(2*P_ATOMIC)/(2*P_ATOMIC+5*O_ATOMIC)

def D(x):
    d=Decimal(str(x))
    if not d.is_finite() or d<0:raise ValueError("Nonfinite/negative SOURCE nutrient flux")
    return d

def combine(water,soil):
    if water.get("status")!="OFFICIAL_2022_VALIDATED_MONTHLY_A4_NGL_PT_MASS_RECONSTRUCTED":
        raise ValueError("Cannot combine unconfirmed actual monthly aqueous source")
    if water.get("expected_SHA256")!=PIN_AESN_SHA or water.get("STEU")!="037800501000" or water.get("year")!=2022:
        raise ValueError("Wrong AESN original source, ST EU or year")
    w=water.get("actual_source") or {}
    if (w.get("original_source_discovery_run")!=44536921058 or
        w.get("original_reporting_year")!=2022 or
        w.get("monthly_2022_source_indicator_count")!=84 or
        w.get("original_STEU_rows")!=528 or
        w.get("A3_influent_not_effluent") is not True or
        w.get("A5_bypass_separate_not_summed_into_A4") is not True or
        w.get("no_EQR_NUTRI_HNS_used") is not True or
        w.get("biological_effect_causal_path_GATE")!="STOP"):
        raise ValueError("Aquatic source audit provenance/gates not satisfied")
    if (soil.get("schema")!="gesui_SIAAP_SeineAval_2022_actual_applied_lot_weighted_NTK_P2O5_SOURCE_ESTIMATE_v1"
        or soil.get("actually_spread_campaign_year")!=2022 or
        soil.get("actual_applied_wet_t")!=81195.15 or
        soil.get("independent_WWTP_sources")!=1 or
        soil.get("NTK_is_not_total_N_for_NGL_comparison") is not True or
        soil.get("coupled_water_land_mass_balance_and_mediation_GO") is not False or
        soil.get("source_only_no_EQR_NUTRI_HNS") is not True):
        raise ValueError("Agronomic source application/provenance/gates not satisfied")
    monthly=w.get("monthly")
    if not isinstance(monthly,list) or len(monthly)!=12:
        raise ValueError("All twelve original A4 month reports required")
    annual=w.get("annual_source_flux_indicators") or {}
    tests={}
    for point in ("A3","A4","A5"):
        for param in ("NGL","PT"):
            key=f"{point}_{param}"
            vals=[]
            for i,m in enumerate(monthly):
                if m.get("month")!=f"{i+1:02d}":
                    raise ValueError("Original monthly data shifted")
                vals.append(D(m[key+"_kg"]))
            src=sum(vals,Decimal(0))
            if abs(src-D(annual[key+"_kg_2022"]))>Decimal("0.01"):
                raise ValueError("Independent monthly sum does not match annual source report")
            tests[key+"_kg_2022"]=src
    water_A4_NGL=tests["A4_NGL_kg_2022"]
    water_A4_P=tests["A4_PT_kg_2022"]
    water_A5_NGL=tests["A5_NGL_kg_2022"]
    water_A5_P=tests["A5_PT_kg_2022"]
    inflow_NGL=tests["A3_NGL_kg_2022"]
    inflow_P=tests["A3_PT_kg_2022"]
    land_NTK=D(soil.get("estimated_applied_Kjeldahl_N_NTK_t"))*1000
    land_P=D(soil.get("estimated_applied_elemental_P_equivalent_t"))*1000
    land_P2O5=D(soil.get("estimated_applied_P2O5_t"))*1000
    if abs(land_P/land_P2O5-P_FACTOR)>Decimal("0.00001"):
        raise ValueError("P2O5-to-elemental-P mass conversion changed")
    volume=D(annual.get("A4_treated_water_m3_2022"))
    if volume<=0 or abs(sum((D(m["A4_treated_water_m3"]) for m in monthly),Decimal(0))-volume)>Decimal("0.01"):
        raise ValueError("Volume did not reconcile with 12 monthly source records")
    # This is just an independent check against original May 2023 *self-evaluation*
    # for treated-outlet removal; no hypothesis test/causal effects or inclusion
    # of A5 or abnormal events. An imperfect crosscheck is expected.
    apparent_N_removal=(1-water_A4_NGL/inflow_NGL)*100
    apparent_P_removal=(1-water_A4_P/inflow_P)*100
    if not (Decimal("0")<=apparent_N_removal<=100 and Decimal("0")<=apparent_P_removal<=100):
        raise ValueError("Impossible apparent flux reduction")
    output={
      "schema":"gesui_SeineAval_2022_two_outlet_aqueous_A4_land_applied_NP_SOURCE_comparison_v1",
      "origin_plant_name":"PARIS SEINE AVAL",
      "origin_ministry_SANDRE":"037800501000",
      "EEA_uwwCode":"FR037800501000",
      "matched_calendar_year":2022,
      "independently_frozen_source_runs":{"A4_water":AESN_RUN,"actually_applied_sludge":SOIL_RUN},
      "original_AESN_source_SHA256":PIN_AESN_SHA,
      "aquatic_official_validated_monthly_calculated_indicator_A4":{
        "NGL_t_2022":float(water_A4_NGL/1000),
        "elemental_PT_t_2022":float(water_A4_P/1000),
        "treated_water_m3_2022":float(volume),
        "A4_is_treated_effluent_outlet_only":True,
      },
      "aquatic_additional_bypass_A5_not_merged_into_A4":{
        "NGL_t_2022":float(water_A5_NGL/1000),
        "elemental_PT_t_2022":float(water_A5_P/1000),
      },
      "land_actual_spread_2022_lot_weighted_product_NP_estimates":{
        "applied_wet_sludge_t_2022":81195.15,
        "applied_Kjeldahl_NTK_t_2022":float(land_NTK/1000),
        "applied_P2O5_equivalent_t_2022":float(land_P2O5/1000),
        "applied_elemental_P_equivalent_t_2022":float(land_P/1000),
        "2021_origin_and_delivery_stock_can_be_used_in_2022":True,
        "applied_lot_mean_not_direct_field_specific_total":True,
      },
      "same_year_source_outlet_elemental_P_magnitude_comparison":{
        "estimated_land_applied_elemental_P_over_A4_treated_PT_mass":float(land_P/water_A4_P),
        "source_mass_comparison_not_2022_material_allocation_fraction":True,
        "source_receiving_geography_differs_river_vs_13_departments":True,
      },
      "nitrogen_chemical_species_not_identical_NTK_vs_NGL":True,
      "independent_regulatory_self_evaluation_crosscheck_posthoc":{
        "A3_to_A4_apparent_NGL_reduction_percent":float(apparent_N_removal),
        "A3_to_A4_apparent_PT_reduction_percent":float(apparent_P_removal),
        "SIAAP_May2023_self_evaluation_2022_NGL_percent":61.6,
        "SIAAP_May2023_self_evaluation_2022_PT_percent":80.0,
        "not_directly_equivalent_due_reporting_and_A5_bypass_and_abnormal_conditions":True,
        "original_official_source":"https://www.drieat.ile-de-france.developpement-durable.gouv.fr/IMG/pdf/support_siaap_css_sav_vf_15mai2023.pdf",
      },
      "true_closed_2022_plant_N_P_source_mass_balance":False,
      "field_specific_exposure_and_no_spreading_controls_verified":False,
      "independent_2018_Feurs_SaintJean_2022_SeineAval_not_exchangeable_without_year":True,
      "all_outcome_biology_EQR_NUTRI_HNS_unopened":True,
      "causal_N_P_transfer_and_biotic_effect_gate":"STOP"
    }
    return output

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--aqueous",required=True,type=Path)
    p.add_argument("--land",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    d=combine(json.loads(a.aqueous.read_text()),json.loads(a.land.read_text()))
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(d,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GESUI_SEINEAVAL_2022_SAME_STATION_WATER_LAND_SOURCE_PAIR",json.dumps({
       "aquatic_A4":d["aquatic_official_validated_monthly_calculated_indicator_A4"],
       "land":d["land_actual_spread_2022_lot_weighted_product_NP_estimates"],
       "same_calendar_year_P_magnitude_ratio":d["same_year_source_outlet_elemental_P_magnitude_comparison"],
       "regulatory_crosscheck":d["independent_regulatory_self_evaluation_crosscheck_posthoc"],
       "causal_biotic_gate":d["causal_N_P_transfer_and_biotic_effect_gate"]},ensure_ascii=False))
if __name__=="__main__":main()
