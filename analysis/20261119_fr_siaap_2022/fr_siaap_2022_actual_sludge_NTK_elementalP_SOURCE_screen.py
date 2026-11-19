#!/usr/bin/env python3
"""Estimate 2022 actual land-applied nutrients from ORIGINAL spread product lots.

Seine Aval official 2022 bilan Table 2: 24,099.10 tMB BTF; 57,096.05 tMB
BTC, actually applied. Original agronomic ANNEX Table pp11/12 gives lot
tonnage-weighted mean applied sludge composition NTK, P2O5 in kg/t MB.
Do NOT use product produced composition or t dry matter as multiplier.
NTK is Kjeldahl N (not all nitrogen species) and P2O5 is oxide equivalent.
Rounded published composition gives estimates, not exact individual parcel N/P.
"""
import argparse,csv,hashlib,json
from decimal import Decimal
from pathlib import Path

P_ATOMIC=Decimal("30.973761998")
O_ATOMIC=Decimal("15.999")
ELEMENTAL_P_PER_P2O5=(Decimal(2)*P_ATOMIC)/(Decimal(2)*P_ATOMIC+Decimal(5)*O_ATOMIC)
ACTUAL_PRODUCTS={
 "BTF":{"wet":Decimal("24099.10"),"NTK":Decimal("9.8"),"P2O5":Decimal("47.6"),"lots":46,"page":"11"},
 "BTC":{"wet":Decimal("57096.05"),"NTK":Decimal("8.2"),"P2O5":Decimal("39.5"),"lots":46,"page":"12"}
}
COMP_COLS=("sludge_product","source_campaign_year","applied_lot_count",
    "applied_sludge_NTK_kg_per_wet_t_lot_weighted_mean",
    "applied_sludge_P2O5_kg_per_wet_t_lot_weighted_mean",
    "original_source_page","original_source_kind")
TOTAL_WET=Decimal("81195.15")
ORIGINAL_BILAN="https://www.siaap.fr/fileadmin/user_upload/Siaap/3_Equipements/Les_usines/Seine_aval/Bilan_boues_Seine_Aval/Bilan_de_campagne_2022.pdf"
def D(x):
    value=Decimal(str(x))
    if not value.is_finite() or value<=0:raise ValueError("Missing or invalid original wet mass or composition")
    return value

def estimate(products,comp):
    if len(products)!=2 or len(comp)!=2:raise ValueError("Requires two public original sludge product types")
    pd={}
    for r in products:
        if r.get("source_campaign_year")!="2022" or r.get("sludge_product") in pd or r.get("sludge_product") not in ACTUAL_PRODUCTS:
            raise ValueError("Product SOURCE year/key invalid")
        product=r["sludge_product"]
        wet=D(r["applied_wet_sludge_t"])
        if wet!=ACTUAL_PRODUCTS[product]["wet"]:
            raise ValueError("Official 2022 actually spread product wet mass mismatch")
        pd[product]=wet
    cd={}
    for r in comp:
        if tuple(r)!=COMP_COLS or r.get("source_campaign_year")!="2022":
            raise ValueError("Unexpected applied lot composition source schema")
        product=r["sludge_product"]
        if product not in ACTUAL_PRODUCTS or product in cd:
            raise ValueError("Missing or duplicated applied product lot composition")
        expected=ACTUAL_PRODUCTS[product]
        if (int(r["applied_lot_count"])!=expected["lots"] or
            r["original_source_page"]!=expected["page"] or
            r["original_source_kind"]!="actually_spread_2022_lot_weighted_analyses"):
            raise ValueError("Origin/method of lot weighted mean changed")
        n=D(r["applied_sludge_NTK_kg_per_wet_t_lot_weighted_mean"])
        p=D(r["applied_sludge_P2O5_kg_per_wet_t_lot_weighted_mean"])
        if n!=expected["NTK"] or p!=expected["P2O5"]:
            raise ValueError("Original applied not produced fertilizer composition changed")
        cd[product]={"NTK":n,"P2O5":p}
    if set(pd)!=set(ACTUAL_PRODUCTS) or set(cd)!=set(ACTUAL_PRODUCTS):
        raise ValueError("Both original applied product types needed")
    if sum(pd.values())!=TOTAL_WET:
        raise ValueError("Original applied wet total does not close")
    source=[]
    sums={"NTK_kg":Decimal(0),"P2O5_kg":Decimal(0)}
    for product in sorted(pd):
        wet=pd[product]
        ng=wet*cd[product]["NTK"]
        pg=wet*cd[product]["P2O5"]
        source.append({
          "product":product,"actual_spread_wet_t":float(wet),
          "applied_lot_weighted_NTK_kg_per_wet_t":float(cd[product]["NTK"]),
          "applied_lot_weighted_P2O5_kg_per_wet_t":float(cd[product]["P2O5"]),
          "estimated_land_NTK_kg":float(ng),
          "estimated_land_P2O5_kg":float(pg),
          "estimated_land_elemental_P_kg":float(pg*ELEMENTAL_P_PER_P2O5),
          "not_field_specific_actual_uptake":True})
        sums["NTK_kg"]+=ng
        sums["P2O5_kg"]+=pg
    n=sums["NTK_kg"]; p=sums["P2O5_kg"]
    return {
      "schema":"gesui_SIAAP_SeineAval_2022_actual_applied_lot_weighted_NTK_P2O5_SOURCE_ESTIMATE_v1",
      "original_official_agronomic_bilan":ORIGINAL_BILAN,
      "actually_spread_campaign_year":2022,
      "original_bilan_annex_printed_pages":["11","12"],
      "actual_applied_wet_t":float(TOTAL_WET),
      "applied_product_lot_weighted_analysis_source_confirmed":True,
      "products":source,
      "estimated_applied_Kjeldahl_N_NTK_t":float(n/1000),
      "estimated_applied_P2O5_t":float(p/1000),
      "estimated_applied_elemental_P_equivalent_t":float(p*ELEMENTAL_P_PER_P2O5/1000),
      "elemental_P_over_P2O5_factor":float(ELEMENTAL_P_PER_P2O5),
      "published_lot_weighted_means_rounded_to_one_decimal":True,
      "not_exact_parcel_specific_or_unrounded_mass":True,
      "NTK_is_not_total_N_for_NGL_comparison":True,
      "2021_production_or_delivery_lots_can_be_spread_in_2022":True,
      "agronomic_chemical_N_P_is_not_field_host_plant_uptake":True,
      "independent_WWTP_sources":1,
      "source_only_no_EQR_NUTRI_HNS":True,
      "coupled_water_land_mass_balance_and_mediation_GO":False,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--product",type=Path,required=True)
    p.add_argument("--composition",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    with a.product.open(newline="",encoding="utf-8") as f:
        pro=list(csv.DictReader(f))
    with a.composition.open(newline="",encoding="utf-8") as f:
        rdr=csv.DictReader(f)
        if tuple(rdr.fieldnames or ())!=COMP_COLS:raise ValueError("Unexpected source composition columns")
        comp=list(rdr)
    report=estimate(pro,comp)
    report["provenance"]={"actual_product_csv_sha256":hashlib.sha256(a.product.read_bytes()).hexdigest(),
                          "actually_applied_composition_csv_sha256":hashlib.sha256(a.composition.read_bytes()).hexdigest()}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("SIAAP_2022_LAND_ACTUAL_NTK_P2O5_ESTIMATE",json.dumps({
       "NTK_t":report["estimated_applied_Kjeldahl_N_NTK_t"],
       "P2O5_t":report["estimated_applied_P2O5_t"],
       "elemental_P_t":report["estimated_applied_elemental_P_equivalent_t"],
       "chemical_lots_not_plant_uptake":True},sort_keys=True))

if __name__=="__main__":main()
