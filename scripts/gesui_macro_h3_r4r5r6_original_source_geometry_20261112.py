#!/usr/bin/env python3
"""SOURCE ONLY: audit original Europe/US river points vs 2015 global WPI port points
on exactly the same H3 r5 spatial support, with r4/r6 ONLY coverage sensitivity.

Not a HydroBASINS basin assignment or an ecological coefficient. Port 2011-24
GBIF plant aggregate and 2015 model WWp CANNOT be matched to 2018/2022
river observations as if they were same-period taxon counts/WWTP effluent.
"""
import csv,hashlib,json,os,pathlib,time,urllib.request,collections,statistics
import h3

REPO="zuizui0223/gesui"
GIT_SHA={
 "EEA2018":"0f30bfbc3dac13f6a182c8a910e9899c04778322",
 "EEA2022":"325f199e985d9862ffee5b0edb205c61a7be3833",
 "Port2015":"c0868d0839e7e68e9503b87df52477160cd6783e"
}
ORIGINAL_PORT_GIT_COMMIT="89a874a96c4751792ba9baa47eb101ab536bee8d"
US_FROZEN_MANIFEST_SHA256="fbf53b8379d4e6564dbb51b3705b02fdd03e4a2caacad06a2314d643c8412e91"
FILES={
 "EEA2018":"results/eea2018_independent_temporal_bqe_20261108/benthic_invertebrates_2018_WWTP_source_real_EQR.csv",
 "EEA2022":"results/eea2022_biology_bqe_20261108/BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv",
 "Port2015":"results/global_all_Plantae_WWp2015_20261110/ORIGINAL_ALL1055_42COUNTRY_WWp2015_BEFORE_NEW808_GBIF_SOURCE.json"
}
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_global_common_scale_geometry"
ROOT.mkdir(parents=True,exist_ok=True)
def gitblob(raw):
    return hashlib.sha1(("blob "+str(len(raw))+"\0").encode()+raw).hexdigest()
def fetch_original(name):
    expected=GIT_SHA[name]
    if name.startswith("EEA"):
        path=pathlib.Path(FILES[name]);raw=path.read_bytes()
    else:
        url="https://raw.githubusercontent.com/"+REPO+"/"+ORIGINAL_PORT_GIT_COMMIT+"/"+FILES[name]
        req=urllib.request.Request(url,headers={"User-Agent":"gesui-global-common-support-source-QA/1","Accept":"application/json"})
        err=None
        for k in range(4):
            try:
                with urllib.request.urlopen(req,timeout=80) as x:raw=x.read(4000000)
                break
            except Exception as e:
                err=e;time.sleep(1+2*k)
        else:raise RuntimeError("Cannot access exact previously frozen WPI port source "+repr(err))
    if gitblob(raw)!=expected:raise RuntimeError("Original source modified for "+name)
    return raw
def csv_rows(raw):
    return list(csv.DictReader(raw.decode("utf-8-sig").splitlines()))
eu18=csv_rows(fetch_original("EEA2018"))
eu22=csv_rows(fetch_original("EEA2022"))
ports=json.loads(fetch_original("Port2015"))["source_all_1055_public_country_port_site_metadata_no_new_GBIF_response"]
if len(eu18)!=2634 or len(eu22)!=2475 or len(ports)!=1055:
    raise RuntimeError("Exact original EEA site or WPI port cohort size changed")
path=pathlib.Path(os.environ["EPA_FROZEN_MANIFEST"])
raw=path.read_bytes()
if hashlib.sha256(raw).hexdigest()!=US_FROZEN_MANIFEST_SHA256:
    raise RuntimeError("Original EPA US native UID pre-numeric source is not frozen original")
us=json.loads(raw)["primary_berw_UID_subset"]
if len(us)!=1221 or len(set(u["uid"] for u in us))!=1221:
    raise RuntimeError("US EPA native UID pre-numeric source cohort unexpectedly different")
# Data dictionaries carry ONLY source geography and coarse timing. Ecological
# values/GBIF flora are not read/used in this QA.
def record(name,key,lat,lon,place,time,sitelabel):
    try:
        lat=float(lat);lon=float(lon)
        if not (-90<=lat<=90 and -180<=lon<=180):
            raise ValueError("Outside global lat/lon")
        if not (-85<=lat<=85):
            raise ValueError("H3 pentagon/near poles out of safe reporting coverage")
        cells={str(res):h3.latlng_to_cell(lat,lon,res) for res in [4,5,6]}
    except (TypeError,ValueError) as e:
        return {"dataset":name,"id":str(key),"valid":False,"reason":str(e)}
    return {"dataset":name,"id":str(key),"valid":True,
            "lat":lat,"lon":lon,"place":str(place),"period":time,
            "physical_site":str(sitelabel),"cells":cells}
records=[]
records.extend(record("EEA18","2018:"+str(x["site_key"]),x["lat"],x["lon"],x["countryCode"],"2017-2019",x["site_key"]) for x in eu18)
records.extend(record("EEA22","2022:"+str(x["site_key"]),x["lat"],x["lon"],x["countryCode"],"2020-2022",x["site_key"]) for x in eu22)
records.extend(record("US_NRSA18_19",x["uid"],x["lat"],x["lon"],x["state"],"2017-2019",x["site_id"]) for x in us)
records.extend(record("WPI_PORT_2015_GEOGRAPHY_ONLY",x["original_port_id"],x["original_port_lat"],x["original_port_lon"],x["original_source_country"],"2015_STATIC_WW_MODEL_2011_2024_GBIF_AGGREGATED_NOT_ALIGNED",x["original_port_id"]) for x in ports)
accepted=[r for r in records if r["valid"]]
missing=[r for r in records if not r["valid"]]
source_counts=collections.Counter(x["dataset"] for x in accepted)
expected={"EEA18":2634,"EEA22":2475,"US_NRSA18_19":1221,"WPI_PORT_2015_GEOGRAPHY_ONLY":1055}
if len(accepted)!=sum(expected.values()) or any(source_counts[t]!=expected[t] for t in expected):
    raise RuntimeError("At least one original source coordinate invalid or original source count drifted; source STOP")
unique18={r["physical_site"] for r in accepted if r["dataset"]=="EEA18"}
unique22={r["physical_site"] for r in accepted if r["dataset"]=="EEA22"}
overlap18_22=unique18&unique22
# SAME station identity if source key matches, regardless of cell.
if not len(overlap18_22)>1000:raise RuntimeError("Original 2018/22 EEA overlap audit unexpectedly changed")
results={
 "schema":"GESUI_SAME_GEOGRAPHY_H3_SOURCE_ONLY_20261112_v1",
 "source_types_are_NOT_EQUIVALENT":{
    "EEA18":"2018 original actual EU WWTP reported upstream source + EQR source; early time bin",
    "EEA22":"2022 original actual EU WWTP reported upstream source + EQR source; late time bin",
    "US_NRSA18_19":"2018-19 original same-site EPA in-stream TN/TP and EPT source; early time bin only",
    "WPI_PORT_2015_GEOGRAPHY_ONLY":"2015 wastewater generated (not effluent or measured N/P); source coastal port buffer has 2011-2024 GBIF aggregate and NO comparable year-resolved plant data"
 },
 "exact_source_git_blobs":GIT_SHA,
 "original_port_commit":ORIGINAL_PORT_GIT_COMMIT,
 "US_original_pre_result_source_manifest_sha256":US_FROZEN_MANIFEST_SHA256,
 "original_input_row_counts":expected,"source_valid_rows":dict(source_counts),
 "source_latlon_invalid_rows":len(missing),
 "EEA_2018_2022_same_original_monitoring_sites":len(overlap18_22),
 "US_source_UID_visits":len(us),
 "US_distinct_physical_SITE_ID":len({x["site_id"] for x in us}),
 "H3_cells_are_NOT_upstream_connected_watersheds":True,
 "hydrobains_level7_remains_required_for_ecological_inference":True,
 "no_single_basin_year_multi_taxon_result":True,"no_regression_or_ecological_effect_calculated":True,
 "levels":{}
}
types=list(expected)
for res in [4,5,6]:
    cells={ds:{x["cells"][str(res)] for x in accepted if x["dataset"]==ds} for ds in types}
    keys={ds:collections.defaultdict(set) for ds in types}
    for x in accepted:
        keys[x["dataset"]][x["cells"][str(res)]].add(x["physical_site"])
    union=set().union(*cells.values())
    areas=[h3.cell_area(c,unit="km^2") for c in union]
    intersections={}
    for i,a in enumerate(types):
        for b in types[i+1:]:
            both=cells[a]&cells[b]
            intersections[a+"_AND_"+b]={
               "same_H3_cell_count":len(both),
               "of_"+a+"_cell_count_fraction":len(both)/len(cells[a]),
               "of_"+b+"_cell_count_fraction":len(both)/len(cells[b])
            }
    eu_same=keys["EEA18"]
    eu_late=keys["EEA22"]
    results["levels"][str(res)]={
       "source_unique_H3_cells":{k:len(v) for k,v in cells.items()},
       "all_original_candidate_H3_cells":len(union),
       "cell_area_km2_min":min(areas),"cell_area_km2_median":statistics.median(areas),"cell_area_km2_max":max(areas),
       "cross_dataset_coincident_cells_only_not_stations_or_exposure":intersections,
       "EEA_2018_2022_coincident_H3_cells":len(cells["EEA18"]&cells["EEA22"]),
       "EEA_2018_2022_exact_station_keys_in_source":len(overlap18_22),
       "matched_numerical_taxonomic_butterfly_and_plant_macro_unit_existing":0,
       "H3_r_is_spatial_coverage_not_observed_river_to_port_connection":True
    }
# Persist aggregated data and a tiny example row with nomenclature only.
outfile=ROOT/"SOURCE_ONLY_ORIGINAL_EEA_EPA_WPI_H3_R4_R5_R6_SPATIAL_CROSSWALK_QA.json"
outfile.write_text(json.dumps(results,indent=2,sort_keys=True))
print("MACRO_H3_SOURCE_GEOMETRY_RESULT_START",flush=True)
print(json.dumps(results,sort_keys=True),flush=True)
print("MACRO_H3_SOURCE_GEOMETRY_RESULT_END",flush=True)
