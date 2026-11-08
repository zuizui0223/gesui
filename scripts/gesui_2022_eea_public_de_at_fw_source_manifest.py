#!/usr/bin/env python3
"""2022 EEA original measured both-NP unique freshwater WWTP outlets in Germany/Austria.

Official source-only 2022 UWWTD WWTP/DischargePoints REST; no butterfly
event or species/abundance, no chemistry at butterfly habitat or EQR.
Create shareable PUBLIC original outfall coordinates for privacy-safe
data-coordinator pre-screen of independently monitored butterfly sites.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,re
from collections import Counter
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
import gesui_2022_sebms_events_to_eea_wwtp_source_only as eea

EXPECTED_WWTP_BOTH_NP=5492
EXPECTED_DCP_ALL=24416
COUNTRY_BBOX={
 "DE":(5.0,47.0,16.0,56.0),
 "AT":(9.0,46.0,18.0,50.0),
}
EU_BBOX=(-12.0,44.0,26.0,64.0)

def hash_file(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",required=True,type=Path)
    args=p.parse_args()
    plant,prec=eea.download_layer(eea.PLANT,eea.BOTH_WHERE,eea.PLANT_FIELDS,EXPECTED_WWTP_BOTH_NP)
    point,drec=eea.download_layer(eea.OUTFALL,"1=1",eea.OUT_FIELDS,EXPECTED_DCP_ALL)
    plant=eea.normalize_key(plant)
    point=eea.normalize_key(point)
    plant=plant[plant.compound_join_key.notna()].copy()
    point=point[point.compound_join_key.notna()].copy()
    pn=plant.compound_join_key.value_counts()
    dn=point.compound_join_key.value_counts()
    qc={
        "original_EEA2022_WWTP_with_both_positive_measured_NP":EXPECTED_WWTP_BOTH_NP,
        "original_EEA2022_discharge_points_all":EXPECTED_DCP_ALL,
        "source_original_WWTP_distinct_nonempty_full_key":int(pn.size),
        "source_original_DCP_distinct_nonempty_full_key":int(dn.size),
        "source_WWTP_duplicate_key_count":int(pn.gt(1).sum()),
        "source_outfall_multi_point_key_count":int(dn.gt(1).sum()),
    }
    plant=plant[plant.compound_join_key.map(pn).eq(1)].copy()
    plant["n_official_outfalls"]=plant.compound_join_key.map(dn).fillna(0).astype(int)
    qc["with_exactly_one_official_outfall"]=int(plant.n_official_outfalls.eq(1).sum())
    eligible=plant[plant.n_official_outfalls.eq(1)].merge(
       point.drop_duplicates("compound_join_key"),on="compound_join_key",
       how="inner",suffixes=("_plant","_outfall"),validate="one_to_one")
    for c in ("dcpLatitude","dcpLongitude","uwwNDischargeMeasured","uwwPDischargeMeasured"):
        eligible[c]=pd.to_numeric(eligible[c],errors="coerce")
    if not (eligible[["uwwNDischargeMeasured","uwwPDischargeMeasured"]].to_numpy(dtype=float)>0).all():
        raise ValueError("Official EEA both-positive-source precondition not maintained")
    orig_geo_n=len(eligible)
    eligible=eligible[
      eligible.dcpWaterBodyType.astype("string").str.strip().str.upper().eq("FW") &
      eligible.dcpLatitude.between(EU_BBOX[1],EU_BBOX[3]) &
      eligible.dcpLongitude.between(EU_BBOX[0],EU_BBOX[2])].copy()
    qc["original_one_outfall_before_FW_filter"]=orig_geo_n
    qc["original_unique_direct_FW_in_full_EU_study_bbox"]=len(eligible)

    source_rep=("repCode_plant" if "repCode_plant" in eligible.columns else "repCode")
    rep=eligible[source_rep].astype("string").str.strip().str.upper()
    eligible["official_reporter_code"]=rep
    reporter_frequency=Counter(rep.fillna("MISSING").astype(str))
    qc["all_original_reporter_top_codes"]=dict(reporter_frequency.most_common(70))
    out=[];country={}
    for iso,(west,south,east,north) in COUNTRY_BBOX.items():
        in_country_geo=eligible[
            eligible.dcpLatitude.between(south,north) &
            eligible.dcpLongitude.between(west,east)].copy()
        # National reporter codes are source-only verified, not guessed from photo locations.
        national=eligible[
            eligible.official_reporter_code.str.match(r"^"+iso+r"(?:_|-|[0-9])",na=False)].copy()
        national=national[national.dcpLatitude.between(south,north) &
                         national.dcpLongitude.between(west,east)].copy()
        country[iso]={
            "original_country_reporter_code_frequency":dict(Counter(national.official_reporter_code.astype(str))),
            "original_NP_both_positive_unique_direct_FW_n":len(national),
            "within_geographic_country_envelope_n_any_reporter":len(in_country_geo),
            "unique_reporter_plant_ids":int(national.compound_join_key.nunique()),
            "country_open_bbox":[west,south,east,north],
            "no_out_of_country_reporter_used_as_source":True,
        }
        print("GESUI_DE_AT_2022_EEA_MEASURED_NP_OUTFALL_COUNTRY_SOURCE",json.dumps({
            "country":iso,**country[iso]},sort_keys=True),flush=True)
        for row in national.itertuples(index=False):
            # Index columns by name rather than original tuple field names with punctuation.
            rec=row._asdict()
            key=str(rec["compound_join_key"])
            out.append({
                "countryCode":iso,
                "original_EEA_reporter_code":str(rec["official_reporter_code"]),
                "original_EEA_plant_source_key":key,
                "original_EEA_discharge_ID":str(rec.get("dcpCode") or ""),
                "discharge_lat_WGS84":round(float(rec["dcpLatitude"]),9),
                "discharge_lon_WGS84":round(float(rec["dcpLongitude"]),9),
                "original_reported_measured_N_t_year":float(rec["uwwNDischargeMeasured"]),
                "original_reported_measured_P_t_year":float(rec["uwwPDischargeMeasured"]),
                "source_one_original_FW_outfall":True,
                "original_source_cycle":"EEA_2022_UWWTD_public",
                "individual_lab_sampling_year_unconfirmed":True
            })
    out.sort(key=lambda x:(x["countryCode"],x["original_EEA_plant_source_key"]))
    if len(set(x["original_EEA_plant_source_key"] for x in out))!=len(out):
        raise ValueError("Cross-country outfall source duplicate keys")
    args.out_dir.mkdir(parents=True,exist_ok=True)
    csvp=args.out_dir/"EEA2022_DE_AT_PUBLIC_BOTH_NP_MEASURED_ONE_FW_OUTFALLS.csv"
    with csvp.open("w",newline="") as f:
        if not out:raise ValueError("Source-public DE/AT eligible outfalls none; record source gap")
        writer=csv.DictWriter(f,fieldnames=list(out[0]))
        writer.writeheader();writer.writerows(out)
    report={
        "schema":"gesui_2022_official_EEA_reported_measured_both_NP_source_DE_AT_public_FW_outfalls_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "original_API_WWTP":"https://water.discomap.eea.europa.eu/arcgis/rest/services/WISE_UWWTD/TreatmentPlants_2022/MapServer/0",
        "original_API_DCP":"https://water.discomap.eea.europa.eu/arcgis/rest/services/WISE_UWWTD/DischargePoints_2022/MapServer/0",
        "query_filters":"original both measured N/P positive, full reporter+plant ID unique, exactly one direct freshwater FW discharge point, EU and country source envelope",
        "original_source_rest_pages_WWTP":prec,
        "original_source_rest_pages_DCP":drec,
        "public_source_inventory":qc,
        "countries":country,
        "public_discharge_coordinate_manifest_CSV":csvp.name,
        "public_discharge_coordinate_manifest_SHA256":hash_file(csvp),
        "public_discharge_coordinate_manifest_rows":len(out),
        "source_included_both_NP_mass_not_actual_local_concentration":True,
        "an_outfall_source_near_butterflies_does_not_identify_an_ecological_effect":True,
        "no_butterfly_monitoring_sites_or_species_data_opened":True,
        "no_swapping_2018_or_2020_WWTP_discharge_into_2022":True,
        "missing_reported_source_does_not_prove_sewage_free":True,
        "individual_WWTP_mass_reference_lab_date_not_independently_verified":True,
        "next_step":"Provide public source-outfall coordinates to authorized butterfly monitoring coordinators for aggregate counts of independent 2022 monitored transect sites with >=2 survey visit dates near <=2km, <=5km and remote >=10km; coordinators retain protected precise site coordinates."
    }
    (args.out_dir/"EEA2022_DE_AT_PUBLIC_WWTP_NP_SOURCE_QC.json").write_text(
        json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GESUI_DE_AT_EEA2022_ORIGINAL_BOTH_NP_FW_COORDINATE_MANIFEST_QC_COMPLETE",
          json.dumps({"count":len(out),"countries":{k:v["original_NP_both_positive_unique_direct_FW_n"]
              for k,v in country.items()},"no_butterfly_data":True}),flush=True)

if __name__=="__main__":main()
