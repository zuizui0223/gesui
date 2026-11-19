#!/usr/bin/env python3
"""Official 2022 Seine Aval sludge campaign: two independent table decompositions.

Compare Table 2 (BTF/BTC) against Table 3 (13 departments) and the printed
TOTAL of actually spread wet/dry mass, area, and parcel count.
All department rows derive from ONE wastewater treatment plant. Department and
parcel counts are NOT independent WWTPs, catchments or documented field exposure
locations; the public report does not give mapped field geometries in this QA.
"""
import argparse
import csv
import hashlib
import json
from decimal import Decimal,InvalidOperation
from pathlib import Path

ORIGIN_CURRENT_SANDRE="037800501000" # current Ministry station fiche, unverified specifically in 2022 data
OFFICIAL_BILAN="https://www.siaap.fr/fileadmin/user_upload/Siaap/3_Equipements/Les_usines/Seine_aval/Bilan_boues_Seine_Aval/Bilan_de_campagne_2022.pdf"
OFFICIAL_SANDRE="https://assainissement.developpement-durable.gouv.fr/pages/data/fiche-037800501000"
DEPT_HEADER=("department_code","department_name","source_campaign_year","applied_wet_sludge_t",
             "applied_dry_sludge_t","actual_applied_area_ha","actual_parcel_count")
PRODUCT_HEADER=("sludge_product","source_campaign_year","applied_wet_sludge_t",
                "applied_dry_sludge_t","actual_applied_area_ha","source_process")
EXPECTED_DEPTS=("02","18","27","28","41","45","51","60","76","77","78","80","95")
EXPECTED_PRODUCTS={
  "BTF":{"wet":Decimal("24099.10"),"dry":Decimal("11905.54"),
         "area":Decimal("3650.48"),"process":"thermally_conditioned_filtered"},
  "BTC":{"wet":Decimal("57096.05"),"dry":Decimal("23197.12"),
         "area":Decimal("7210.01"),"process":"thermally_conditioned_centrifuged"}
}
TOTALS={"wet":Decimal("81195.15"),"dry":Decimal("35102.66"),
        "area":Decimal("10860.49"),"parcels":1092}


def csvread(path,expected_fields):
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        rd=csv.DictReader(f)
        if tuple(rd.fieldnames or ())!=expected_fields:
            raise ValueError("Unexpected table schema (possible private parcel/farm data)")
        return list(rd)


def pos(x):
    try:z=Decimal(x)
    except (ValueError,TypeError,InvalidOperation) as exc:
        raise ValueError("Unknown or missing original mass/area") from exc
    if not z.is_finite() or z<=0:raise ValueError("Nonpositive or invalid source mass")
    return z


def audit(depts,products):
    if len(depts)!=13 or len(products)!=2:
        raise ValueError("Expected exactly 13 original departments and 2 product families")
    dept_receipts=[]
    observed_codes=set()
    for row in depts:
        if tuple(row)!=DEPT_HEADER or row["source_campaign_year"]!="2022":
            raise ValueError("Department schema/year invalid")
        code=row["department_code"]
        if code not in EXPECTED_DEPTS or code in observed_codes or not row["department_name"]:
            raise ValueError("Missing/unexpected/duplicated SOURCE department")
        observed_codes.add(code)
        amount={k:pos(row[col]) for k,col in (
             ("wet","applied_wet_sludge_t"),("dry","applied_dry_sludge_t"),
             ("area","actual_applied_area_ha"))}
        try:n=int(row["actual_parcel_count"])
        except (TypeError,ValueError) as exc:raise ValueError("Invalid parcel counts") from exc
        if n<=0 or str(n)!=row["actual_parcel_count"]:
            raise ValueError("Negative, missing or noninteger parcel count")
        dept_receipts.append({**amount,"parcels":n,"code":code,"name":row["department_name"]})
    if observed_codes!=set(EXPECTED_DEPTS):
        raise ValueError("Unverified expected departments")
    product_receipts=[]
    products_seen=set()
    for row in products:
        if tuple(row)!=PRODUCT_HEADER or row["source_campaign_year"]!="2022":
            raise ValueError("Product schema/year invalid")
        product=row["sludge_product"]
        if product not in EXPECTED_PRODUCTS or product in products_seen:
            raise ValueError("Unknown or duplicated product")
        products_seen.add(product)
        expected=EXPECTED_PRODUCTS[product]
        quantities={k:pos(row[col]) for k,col in (
             ("wet","applied_wet_sludge_t"),("dry","applied_dry_sludge_t"),
             ("area","actual_applied_area_ha"))}
        if quantities!={k:expected[k] for k in quantities} or row["source_process"]!=expected["process"]:
            raise ValueError("Official product subtotals changed")
        product_receipts.append(quantities)
    if products_seen!=set(EXPECTED_PRODUCTS):raise ValueError("Missing product family")
    by_dept={k:sum((r[k] for r in dept_receipts),0) for k in TOTALS}
    by_product={k:sum((r[k] for r in product_receipts),Decimal(0)) for k in ("wet","dry","area")}
    if by_dept!=TOTALS:
        raise ValueError("13 geographic SOURCE subtotals fail the official 2022 TOTAL")
    if by_product!={k:TOTALS[k] for k in by_product}:
        raise ValueError("BTF+BTC fail independently stated TOTAL")
    return {
       "schema":"gesui_SIAAP_SeineAval_realized_spreading_2022_two_table_SOURCE_audit_v1",
       "origin":"SIAAP - Seine aval",
       "official_agricultural_campaign_year":2022,
       "origin_station_SANDRE_in_current_ministry_fiche":ORIGIN_CURRENT_SANDRE,
       "SANDRE_year2022_direct_original_code_not_independently_fixed":True,
       "river_discharge_recipient":"Seine",
       "official_sources":{"agronomic_campaign_bilan":OFFICIAL_BILAN,
                            "current_station_SANDRE_fiche":OFFICIAL_SANDRE},
       "actual_2022_sludge_spreading_recap":{
          "applied_wet_t":float(TOTALS["wet"]),
          "applied_dry_t":float(TOTALS["dry"]),
          "actually_spread_area_ha":float(TOTALS["area"]),
          "original_field_parcel_count":TOTALS["parcels"],
          "department_count":len(dept_receipts),
          "original_product_types":["BTF","BTC"],
          "sludge_originating_independent_stations":1
       },
       "independent_geographic_subtotals_check_pass":True,
       "independent_product_subtotals_check_pass":True,
       "source_delivery_window_includes_prior_2021":"2021-11-08 to 2022-09-30",
       "realized_field_spreading_windows_2022":[
           "2022-03-23 to 2022-05-10", "2022-07-11 to 2022-11-14"
       ],
       "sludge_delivered_in_2021_not_necessarily_applied_in_2021":True,
       "one_plant_13_departments_1092_parcels_not_independent_stations":True,
       "actual_parcel_geometries_not_downloaded":True,
       "realized_field_N_and_P_measurements_not_linked_to_parcels":True,
       "2022_aqueous_annual_measured_TN_TP_export_not_identified":True,
       "independent_H6_and_H7_actual_field_exposure_counts_confirmed":0,
       "biology_EQR_NUTRI_HNS_accessed":False,
       "source_two_outlet_mass_pair_and_biological_causal_models":"STOP",
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--departments",type=Path,required=True)
    p.add_argument("--products",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    d=csvread(a.departments,DEPT_HEADER)
    s=csvread(a.products,PRODUCT_HEADER)
    z=audit(d,s)
    z["original_transcribed_source_SHA256"]={
        "department_table":hashlib.sha256(a.departments.read_bytes()).hexdigest(),
        "product_table":hashlib.sha256(a.products.read_bytes()).hexdigest()
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(z,indent=2,ensure_ascii=False,sort_keys=True)+"\n",encoding="utf-8")
    print("GESUI_SIAAP_SEINE_AVAL_2022_REAL_SPREADING_SOURCE",json.dumps(
        z["actual_2022_sludge_spreading_recap"],sort_keys=True))


if __name__=="__main__":
    main()
