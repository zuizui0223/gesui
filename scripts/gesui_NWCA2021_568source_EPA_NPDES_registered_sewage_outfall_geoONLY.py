#!/usr/bin/env python3
"""SOURCE-ONLY EPA 2021 NWCA 568 joint wetland sites × US EPA NPDES outfalls.

Source original wetland SITE_ID+VISIT_NO frozen prior to all plant/nutrient
responses; outcome variables H_ALL/H_NAT/FQAI/TN/TP are NEVER read here.
Current ECHO NPDES permitted-feature sewage facilities are potential spatial
sources, not 2021 observed discharge or verified hydrological connectivity.
"""
from __future__ import annotations

import argparse,hashlib,io,json,math,re,zipfile
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

D=Path("results/us_NWCA2021_plant_TNTP_20261110")
ORIGINAL=D/"EPA_NWCA2021_ORIGINAL_VEG_WATERCHEM_SITEONLY_GATE.json"
SOURCEURL="https://echo.epa.gov/files/echodownloads/npdes_outfalls_layer.zip"
SOURCE_TABLE="NPDES_OUTFALLS_LAYER.csv"
SOURCE_RADIUS=(5.,10.,25.)
PRIMARY_RADIUS=10.
FAR_RADIUS=20.
TARGET_YEAR=2021
MIN_NEAR_SITES=25
MIN_PERMITS=15
MIN_NEAR_STATES=8
MIN_MATCHED_FAR_SITES=60
MIN_STATES_BOTH=6
NPDES_FIELDS=["EXTERNAL_PERMIT_NMBR","FACILITY_NAME","SIC_CODES","NAICS_CODES",
    "PERM_FEATURE_NMBR","LATITUDE83","LONGITUDE83","LATLONG_TYPE",
    "PERMIT_EFFECTIVE_DATE","PERMIT_EXPIRATION_DATE","PERMIT_TERMINATION_DATE",
    "PERMIT_STATUS_DESC","STATE_CODE"]
NWCA_FIELDS=["SITE_ID","VISIT_NO","SITETYPE","YEAR","STUDY","WGT_TP_CORE",
    "LAT_ANALYS","LON_ANALYS","PSTL_CODE","AG_ECO9","WETCLS_GRP","HYD_CLS"]

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(2*1024*1024),b""):
            h.update(b)
    return h.hexdigest()

def is_sewage_codes(raw):
    if not isinstance(raw,str):return False
    # NAICS 221320 = sewage treatment facility; SIC 4952 = sewerage systems.
    return bool(re.search(r"(?<!\d)(?:221320|4952)(?!\d)",raw))

def xyz(lat,lon):
    a=np.deg2rad(np.asarray(lat,dtype=float))
    b=np.deg2rad(np.asarray(lon,dtype=float))
    return np.column_stack((np.cos(a)*np.cos(b),np.cos(a)*np.sin(b),np.sin(a)))

def distances_km(points,center):
    lat0,lon0=center
    la=np.deg2rad(np.asarray(points[:,0],dtype=float))
    lo=np.deg2rad(np.asarray(points[:,1],dtype=float))
    dlat=la-math.radians(lat0)
    dlon=lo-math.radians(lon0)
    w=np.sin(dlat/2)**2+np.cos(math.radians(lat0))*np.cos(la)*np.sin(dlon/2)**2
    return 12742.0176*np.arcsin(np.minimum(1,np.sqrt(np.maximum(w,0.))))

def original_nwca_source_sites(root,dir):
    original=json.loads((root/ORIGINAL).read_text())
    if original["schema"]!="gesui_US_EPA_NWCA2021_original_prob_wetland_plant_shannon_waterchem_source_only_join_GATE_v01":
        raise ValueError("Original NWCA SOURCE-ONLY gate schema changed")
    if not original["all_original_2021_plant_H_ALL_FQAI_and_water_NTL_PTL_values_NOT_READ"]:
        raise ValueError("Pre-outcome NWCA site frame was not source blind")
    for key,p in {"SITEINFO":dir/"nwca21_siteinfo-data.csv",
                  "WATERCHEM":dir/"nwca21_waterchem_wide-data.csv",
                  "VEGMETRICS":dir/"nwca21_vegmetrics-data.csv",
                  "HYDROSOURCE":dir/"nwca21_hydrology_sources_wide-data.csv"}.items():
        if sha(p)!=original["official_EPA_original_source_csv_receipts"][key]["official_source_sha256"]:
            raise ValueError("Original EPA NWCA source CSV source SHA differs")
    df=pd.read_csv(dir/"nwca21_siteinfo-data.csv",usecols=NWCA_FIELDS,low_memory=False)
    df["SITE_ID"]=df["SITE_ID"].astype("string").str.strip()
    df["VISIT_NO"]=pd.to_numeric(df["VISIT_NO"],errors="coerce")
    for col in ("YEAR","WGT_TP_CORE","LAT_ANALYS","LON_ANALYS"):
        df[col]=pd.to_numeric(df[col],errors="coerce")
    good=(df.SITETYPE.astype("string").str.upper().eq("PROB")&
          df.WGT_TP_CORE.gt(0)&df.YEAR.eq(TARGET_YEAR)&df.VISIT_NO.eq(1)&
          df.LAT_ANALYS.between(24,50)&df.LON_ANALYS.between(-126,-66)&
          df.STUDY.astype("string").str.contains("NWCA",case=False,na=False))
    df=df.loc[good].copy()
    if len(df)!=original["source2021_site_qualifying_probability_weights_2021_visit1_latlon_n"] or len(df)!=920:
        raise ValueError("Original 920 source probability wetlands no longer same")
    df["source_site_visit_SHA256"]=[
        hashlib.sha256((str(a)+"|"+str(int(b))).encode()).hexdigest()
        for a,b in zip(df.SITE_ID,df.VISIT_NO)]
    keys=set(original["source_joined_original_SITE_ID_SHA256_not_plant_or_chemistry_values"])
    if len(keys)!=568:raise ValueError("Original 568 source study site signatures changed")
    df=df.loc[df.source_site_visit_SHA256.isin(keys)].copy()
    if len(df)!=568 or set(df.source_site_visit_SHA256)!=keys or df.source_site_visit_SHA256.duplicated().any():
        raise ValueError("Original 568 source-only joint plant and water chemistry population changed")
    return df,original

def ndpes_outfall_source(path):
    if path.stat().st_size<100_000:raise ValueError("Original official EPA NPDES zip appears empty")
    with zipfile.ZipFile(path) as z:
        names=[v for v in z.namelist() if v.split("/")[-1].upper()==SOURCE_TABLE]
        if len(names)!=1:raise ValueError("Expected one official NPDES_OUTFALLS_LAYER.csv original source")
        with z.open(names[0]) as f:
            header=pd.read_csv(f,nrows=0,low_memory=False).columns.tolist()
        missing=set(NPDES_FIELDS)-set(header)
        if missing:raise ValueError("Official EPA NPDES source schema missing "+repr(missing))
        with z.open(names[0]) as f:
            raw=pd.read_csv(f,usecols=NPDES_FIELDS,low_memory=False,dtype=str,
                            keep_default_na=False,on_bad_lines="error")
    if len(raw)<10000:raise ValueError("Official USA national NPDES outfalls unexpectedly small")
    raw["LAT"]=pd.to_numeric(raw.LATITUDE83,errors="coerce")
    raw["LON"]=pd.to_numeric(raw.LONGITUDE83,errors="coerce")
    coords=raw.LAT.between(24,50)&raw.LON.between(-126,-66)
    clean=raw.loc[coords].copy()
    sewage=clean.apply(lambda v:is_sewage_codes(v.SIC_CODES) or is_sewage_codes(v.NAICS_CODES),axis=1)
    clean=clean.loc[sewage].copy()
    permitted=clean.LATLONG_TYPE.str.contains(r"PERMIT|FEATURE|OUTFALL",case=False,regex=True,na=False)
    clean["PERMITTED_FEATURE_LATLON"]=permitted
    clean["FACILITY_COORDINATE_PROXY"]=~permitted
    clean["PERMIT_ID"]=clean.EXTERNAL_PERMIT_NMBR.str.strip().str.upper()
    clean["FEATURE_ID"]=clean.PERM_FEATURE_NMBR.str.strip().str.upper()
    clean=clean.loc[(clean.PERMIT_ID!="")&(clean.FEATURE_ID!="")].copy()
    effective=pd.to_datetime(clean.PERMIT_EFFECTIVE_DATE,errors="coerce")
    expiration=pd.to_datetime(clean.PERMIT_EXPIRATION_DATE,errors="coerce")
    termination=pd.to_datetime(clean.PERMIT_TERMINATION_DATE,errors="coerce")
    clean["PERMIT_DATES_PLAUSIBLE_2021"]=(
        effective.le(pd.Timestamp("2021-12-31"))&
        (expiration.isna()|expiration.ge(pd.Timestamp("2021-01-01")))&
        (termination.isna()|termination.ge(pd.Timestamp("2021-01-01")))
    )
    # Duplicate exact source records are a reporting phenomenon, not extra
    # independently permitted discharge systems.
    duplicates=int(clean.duplicated(["PERMIT_ID","FEATURE_ID","LAT","LON"]).sum())
    clean=clean.drop_duplicates(["PERMIT_ID","FEATURE_ID","LAT","LON"]).copy()
    counts={
        "original_EPA_national_NPDES_csv_rows":int(len(raw)),
        "US_contiguous_valid_latlon_rows":int(coords.sum()),
        "sewage_NAICS221320_or_SIC4952_after_coordinates_n":int(sewage.sum()),
        "rows_after_unique_permit_feature_and_usable_ID":int(len(clean)),
        "duplicate_exact_permit_feature_geography_removed":duplicates,
        "all_NPDES_2026_current_source_latlong_type":dict(Counter(raw.LATLONG_TYPE)),
        "sewage_source_latlong_type":dict(Counter(clean.LATLONG_TYPE)),
        "sewage_valid_permitted_feature_geometry_n":int(clean.PERMITTED_FEATURE_LATLON.sum()),
        "sewage_source_permit_dates_plausible2021_n":int(clean.PERMIT_DATES_PLAUSIBLE_2021.sum()),
        "sewage_strict_both_original_feature_geometry_and_dates_n":
            int((clean.PERMITTED_FEATURE_LATLON&clean.PERMIT_DATES_PLAUSIBLE_2021).sum()),
        "sewage_registered_distinct_permits_all_usable":int(clean.PERMIT_ID.nunique()),
        "sewage_registered_distinct_permits_strict":int(clean.loc[
            clean.PERMITTED_FEATURE_LATLON&clean.PERMIT_DATES_PLAUSIBLE_2021,"PERMIT_ID"].nunique()),
        "2026_XX_current_NPDES_register_not_actual_historical_2021_NP_release":True
    }
    return clean,counts

def assess_sites(nwca,facilities,mode):
    if mode=="strict":f=facilities.loc[facilities.PERMITTED_FEATURE_LATLON&facilities.PERMIT_DATES_PLAUSIBLE_2021].copy()
    elif mode=="feature_all_current":f=facilities.loc[facilities.PERMITTED_FEATURE_LATLON].copy()
    elif mode=="all_facility_or_feature_proxy":f=facilities.copy()
    else:raise ValueError("Unexpected fixed source sensitivity")
    if len(f)==0:return [],{"site_n":len(nwca),"outfall_unique_points":0,"usable_source_permits":0}
    srcpos=f[["LAT","LON"]].to_numpy(dtype=float)
    tree=cKDTree(xyz(srcpos[:,0],srcpos[:,1]))
    chord=2*math.sin((25./6371.0088)/2)
    results=[]
    for row in nwca.itertuples(index=False):
        site=(float(row.LAT_ANALYS),float(row.LON_ANALYS))
        cand=tree.query_ball_point(xyz(np.array([site[0]]),np.array([site[1]]))[0],r=chord+1e-11)
        if cand:
            dd=distances_km(srcpos[cand],site)
            dmin=float(np.min(dd))
            lookup=[(int(idx),float(dist)) for idx,dist in zip(cand,dd)]
            near10=[idx for idx,dist in lookup if dist<=10]
            near5=[idx for idx,dist in lookup if dist<=5]
            near25=[idx for idx,dist in lookup if dist<=25]
            permit_ids=set(f.iloc[near10].PERMIT_ID) if near10 else set()
            nearest=min(lookup,key=lambda v:(v[1],str(f.iloc[v[0]].PERMIT_ID)))
            nearpermit=f.iloc[nearest[0]].PERMIT_ID
        else:dmin=None;near10=[];near5=[];near25=[];permit_ids=set();nearpermit=None
        results.append({
            "NWCA2021_source_site_visit_SHA256":row.source_site_visit_SHA256,
            "original_state":str(row.PSTL_CODE),
            "original_ecoregion":str(row.AG_ECO9),
            "original_wetland_group":str(row.WETCLS_GRP),
            "original_hydro_class":str(row.HYD_CLS),
            "nearest_outfall_geodesic_km_within25km_else_null":dmin,
            "nearest_independent_permit_hash_SHA256":hashlib.sha256(str(nearpermit).encode()).hexdigest() if nearpermit else None,
            "outfall_rows_within_5km":len(near5),
            "outfall_rows_within_10km":len(near10),
            "outfall_rows_within_25km":len(near25),
            "unique_NPDES_permit_IDs_10km":len(permit_ids),
            "independent_10km_NPDES_permits_SHA256":sorted(hashlib.sha256(str(pid).encode()).hexdigest() for pid in permit_ids),
            "source_point_not_confirmed_river_to_wetland_water_transport":True
        })
    near=[r for r in results if r["outfall_rows_within_10km"]]
    far=[r for r in results if r["nearest_outfall_geodesic_km_within25km_else_null"] is None or r["nearest_outfall_geodesic_km_within25km_else_null"]>=20]
    matched={(r["original_state"],r["original_ecoregion"]) for r in near}
    far_matched=[r for r in far if (r["original_state"],r["original_ecoregion"]) in matched]
    shared_states=sorted({r["original_state"] for r in far_matched}&{r["original_state"] for r in near})
    countperm={pid for r in near for pid in r["independent_10km_NPDES_permits_SHA256"]}
    stats={
        "site_n":len(results),
        "outfall_unique_source_coordinates":len(f),
        "strictly_sewage_permit_unique":int(f.PERMIT_ID.nunique()),
        "site_near5km_n":sum(r["outfall_rows_within_5km"]>0 for r in results),
        "site_near10km_n":len(near),
        "site_near25km_n":sum(r["outfall_rows_within_25km"]>0 for r in results),
        "site_far_geodesic20km_or_more_n":len(far),
        "states_near_10km":len({r["original_state"] for r in near}),
        "nearest_10km_administrative_unique_NPDES_permits_n":len(countperm),
        "states_with_near_and_far_from_same_state_and_ecoregion":len(shared_states),
        "matched_far_sites_from_same_state_and_ecoregion_as_a_near_site_n":len(far_matched),
        "source_near10_state_counts":dict(Counter(r["original_state"] for r in near)),
        "source_near10_ecoregion_counts":dict(Counter(r["original_ecoregion"] for r in near)),
        "source_far20_state_counts":dict(Counter(r["original_state"] for r in far)),
        "source_near10_vs_far20_shared_state_names":shared_states,
        "site_point_proximity_to_outfall_is_NOT_hydrological_exposure":True
    }
    return results,stats

def run(root,nwca_dir,echo_zip):
    sites,orig=original_nwca_source_sites(root,nwca_dir)
    facilities,facility_qc=ndpes_outfall_source(echo_zip)
    result={}
    details={}
    for mode in ("strict","feature_all_current","all_facility_or_feature_proxy"):
        a,s=assess_sites(sites,facilities,mode)
        result[mode]=s
        details[mode]=a
    strict=result["strict"]
    passed=(strict["site_near10km_n"]>=MIN_NEAR_SITES and
            strict["nearest_10km_administrative_unique_NPDES_permits_n"]>=MIN_PERMITS and
            strict["states_near_10km"]>=MIN_NEAR_STATES and
            strict["matched_far_sites_from_same_state_and_ecoregion_as_a_near_site_n"]>=MIN_MATCHED_FAR_SITES and
            strict["states_with_near_and_far_from_same_state_and_ecoregion"]>=MIN_STATES_BOTH)
    return {
        "schema":"gesui_US_EPA_NWCA2021_568source_NPDES2026_reported_sewage_outfall_GEO_2021permitDates_sourceONLY_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "original_preoutcome_EPA_NWCA2021_568joint_SITE_ID_SHA256":sha(root/ORIGINAL),
        "source_EPA_ECHO_NPDES_2026_live_weekly_original_ZIP_SHA256":sha(echo_zip),
        "source_ECHO_official_ZIP_URL":SOURCEURL,
        "original_2021_NWCA_plant_chemistry_568_site_rows_from_preoutcome_gate":len(sites),
        "current_official_ECHO_registered_2026_NPDES_source_profile":facility_qc,
        "predeclared_nearest_sewage_outfall_geodesic_source_controls":result,
        "all_original_568_NWCA_source_site_to_sewage_source_distance_receipts_without_plant_or_TNTP_values":details,
        "source_precommitted_GEO_eligibility_gate":{
            "min_2021_wetland_original_study_sites_within_strict_sewage_outfall10km":MIN_NEAR_SITES,
            "min_distinct_NPDES_administrative_sewage_permits":MIN_PERMITS,
            "min_US_states_near_source":MIN_NEAR_STATES,
            "min_far_sites_same_original_state_ecoregion":MIN_MATCHED_FAR_SITES,
            "min_states_with_near_and_far_in_original_state_ecoregion":MIN_STATES_BOTH,
            "permitted_2021_outfall_distance_association_source_support_PASS":bool(passed)
        },
        "source_NPDES_current_weekly_2026_export_CANNOT_confirm_2021_WWTP_discharge":True,
        "source_NPDES_latlong_type_facility_not_actual_outfall_recorded_separately":True,
        "no_original_NWCA_H_ALL_FQAI_or_TNTP_outcomes_accessed_in_source_stage":True,
        "NO_hydrological_source_to_wetland_path_confirmed":True,
        "NO_nutrient_mass_or_WWTP_induced_biodiversity_claims":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--nwca-source-dir",type=Path,required=True)
    p.add_argument("--npdes-zip",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    res=run(root,a.nwca_source_dir,a.npdes_zip)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(res,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "source_total_ECHO":res["current_official_ECHO_registered_2026_NPDES_source_profile"],
        "source_distance_by_sewage_outfall_coordinate_age":res["predeclared_nearest_sewage_outfall_geodesic_source_controls"],
        "precommitted_biological_sample_gate":res["source_precommitted_GEO_eligibility_gate"]
    },indent=2))
if __name__=="__main__":
    main()
