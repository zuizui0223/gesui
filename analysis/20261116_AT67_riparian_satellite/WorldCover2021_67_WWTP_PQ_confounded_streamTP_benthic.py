#!/usr/bin/env python3
"""2021 ESA WorldCover original land-use locally at 67 original 2022 Austrian
river TP+BQE sites. Adversarial control of landcover geography on original
2022 upstream WWTP annual P/Q → measured stream TP and benthic EQR.

Prior source-S2 10site natural-class NDVI smoke STOP remains. WorldCover
is a *2021 covariate*, not a downstream biomass or vegetation response.
Never assert causal wastewater transport or year2022 plant production.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,math,sys
from pathlib import Path
import numpy as np
import rasterio
from rasterio.errors import RasterioIOError
from rasterio.windows import from_bounds,Window

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261116_austria_wisa_official_source"))
from AT68_2022_real_WWTP_P_to_ambient_TP_to_benthic_EQR import cohort,P_CONVERSION
from AT68_equal_Q_updown_mass_adversarial_qa import model_cv,comparison

SCHEMA="GESUI_20261116_AT67_2021_WORLD_COVER_LOCAL_LAND_USE_AS_2022_WWTP_WATER_CONFOUNDER_V1"
WC_ROOT="https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/"
EXCLUDE_SITE="FW54110117"
ORIG_REGION={"15_3":9,"15_4":14,"15_5":21,"16_4":13,"16_5":10}
R_M=1000.
R_NEAR_M=250.
SOURCE_MIN_FULL_PIXELS=20000
SOURCE_MIN_INNER_PIXELS=1000
MIN_VALID_STATIONS=60
CLASSES=(10,20,30,40,50,60,80,90)

def worldcover_url(la,lo):
    tile=("N" if la>=0 else "S")+f"{abs(la):02d}"+("E" if lo>=0 else "W")+f"{abs(lo):03d}"
    return WC_ROOT+f"ESA_WorldCover_10m_2021_v200_{tile}_Map.tif",tile

def source_site_landcover(lat,lon):
    # Original WorldCover2021 tiles are 3°; check every tile intersecting
    # the 1km circle, so no silently truncated source sample at a tile edge.
    dlat=R_M/111132.;dlon=R_M/(111320*math.cos(math.radians(lat)))
    minla,maxla=lat-dlat,lat+dlat
    minlo,maxlo=lon-dlon,lon+dlon
    las=list(range(int(math.floor(minla/3)*3),int(math.floor(maxla/3)*3)+1,3))
    los=list(range(int(math.floor(minlo/3)*3),int(math.floor(maxlo/3)*3)+1,3))
    one=collections.Counter();outer=collections.Counter()
    validfull=0;validnear=0;failed=[]
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
           CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
           GDAL_HTTP_MAX_RETRY="2",GDAL_HTTP_RETRY_DELAY="1",
           GDAL_CACHEMAX=64,VSI_CACHE="TRUE",VSI_CACHE_SIZE=1000000):
      for la in las:
        for lo in los:
            url,tile=worldcover_url(la,lo)
            try:
                with rasterio.open("/vsicurl/"+url) as source:
                    if source.crs is None or source.count!=1 or not (
                       "4326" in str(source.crs)):
                        raise RuntimeError("Original ESA WorldCover2021 tile not WGS84 one-band")
                    x0=max(minlo,source.bounds.left);x1=min(maxlo,source.bounds.right)
                    y0=max(minla,source.bounds.bottom);y1=min(maxla,source.bounds.top)
                    if x1<=x0 or y1<=y0:continue
                    window=from_bounds(x0,y0,x1,y1,source.transform)
                    window=Window(int(math.floor(window.col_off)),int(math.floor(window.row_off)),
                        int(math.ceil(window.width))+2,int(math.ceil(window.height))+2)
                    data=source.read(1,window=window,boundless=True,fill_value=0)
                    trans=source.window_transform(window)
                    rows,cols=np.indices(data.shape)
                    X=trans.c+(cols+.5)*trans.a+(rows+.5)*trans.b
                    Y=trans.f+(cols+.5)*trans.d+(rows+.5)*trans.e
                    d2=(((X-lon)*111320*math.cos(math.radians(lat)))**2
                        + ((Y-lat)*111132)**2)
                    valid=data>0
                    inside=(d2<=R_M**2)&valid
                    near=(d2<=R_NEAR_M**2)&valid
                    n=int(inside.sum());k=int(near.sum())
                    validfull+=n;validnear+=k
                    for cls in CLASSES:
                        one[cls]+=int((inside&(data==cls)).sum())
                        outer[cls]+=int((near&(data==cls)).sum())
            except (RasterioIOError,RuntimeError,OSError) as ex:
                failed.append({"source_tile":tile,"error_type":type(ex).__name__,
                    "source_error":str(ex)[:160]})
    if failed or validfull<SOURCE_MIN_FULL_PIXELS or validnear<SOURCE_MIN_INNER_PIXELS:
        return None,{"error_type":"original_source_tiles_or_real_coverage_missing",
            "original_full_circle_valid_cells":validfull,
            "original_250m_near_valid_cells":validnear,"errors":failed}
    fraction={f"wc1km_{key}_fraction":one[key]/validfull for key in CLASSES}
    fraction.update({f"wc250m_{key}_fraction":outer[key]/validnear for key in CLASSES})
    fraction["wc1km_natural_fraction"]=sum(one[k] for k in (10,20,30,90))/validfull
    fraction["wc250m_natural_fraction"]=sum(outer[k] for k in (10,20,30,90))/validnear
    fraction["wc2021_original1km_10m_valid_pixels"]=validfull
    fraction["wc2021_original250m_10m_valid_pixels"]=validnear
    return fraction,None

def main(args):
    args.output.mkdir(parents=True,exist_ok=True)
    source=cohort()
    keep=[q for q in source if q["station"]!=EXCLUDE_SITE]
    if len(keep)!=67 or collections.Counter(x["block"] for x in keep)!=ORIG_REGION:
        raise RuntimeError("Exact 67 2022 Austrian source field site/region geometric QA frame drift")
    good=[];fail=[]
    for i,s in enumerate(sorted(keep,key=lambda z:z["station"]),1):
        result,why=source_site_landcover(s["lat"],s["lon"])
        if result:
            s.update(result)
            s["logPdown_flow_proxy"]=math.log1p(
                 P_CONVERSION*s["Pdown_t_yr"]/math.expm1(s["logQ"]))
            good.append(s)
        else:
            fail.append({"original_source_PKINT":"AT"+s["station"],
                         "original_source_region":s["block"],**why})
        print("GESUI_AT67_SOURCE_WC_LANDCOVER_LOCATION="+str(i)+"/67 "+s["station"]+
             " ok="+str(bool(result)),flush=True)
    blocks=collections.Counter(q["block"] for q in good)
    passed=len(good)>=MIN_VALID_STATIONS and len(blocks)==5 and min(blocks.values())>=5
    output={"schema":SCHEMA,
      "original_2022_Austria_genuine_benthic_2022_and_annual_total_P_sites":68,
      "one_prior_field_coordinate_anomaly_removed_before_landcover":EXCLUDE_SITE,
      "original_study_TP_EQR_WWTP_source_sites":67,
      "source_original_WorldCover_2021_product_v200":"ESA_WorldCover_10m_2021_v200",
      "source_original_WorldCover_2021_license":"CC-BY-4.0",
      "source_original_WorldCover_2021_source_root":WC_ROOT,
      "source_circular_distance_from_true_water_station_m":R_M,
      "source_near_river_circle_m":R_NEAR_M,
      "source_valid_1km_and_250m_landcover_source_stations":len(good),
      "source_valid_sites_by_original_3deg_geo_region":dict(blocks),
      "source_missing_geometry_cases":fail,
      "source_landcover_coverage_GO":bool(passed),
      "source_missing_no_2022_water_EQR_or_NDVI_outcome_imputation":True,
      "source_2021_landcover_is_not_2022_vegetation_response":True,
      "source_S2_2022_natural_land_near_outer_smoke_5_of_10_STILL_STOP":True,
      "source_landcover_adjusted_2022_same_site_model_not_causal":True}
    if passed:
        B=["logQ","logArea","lat","lon"]
        LC=["wc1km_40_fraction","wc1km_50_fraction","wc1km_10_fraction","wc1km_30_fraction"]
        water={
            "B_geographic_hydrology":B,
            "B_plus_local_2021_landcover":B+LC,
            "landcover_plus_UP_source_P_div_flow":B+LC+["logPup_flow_proxy"],
            "landcover_plus_DOWN_source_P_div_flow_placebo":B+LC+["logPdown_flow_proxy"],
            "landcover_plus_both_directions":B+LC+["logPup_flow_proxy","logPdown_flow_proxy"]}
        eco={
            "B_geographic_hydrology":B,
            "B_plus_local_2021_landcover":B+LC,
            "landcover_plus_real_2022_river_TP":B+LC+["logTP"],
            "landcover_plus_UP_source_P_div_flow":B+LC+["logPup_flow_proxy"],
            "landcover_plus_DOWN_source_P_div_flow_placebo":B+LC+["logPdown_flow_proxy"],
            "landcover_plus_real_TP_and_UP_source_P_div_flow":
                B+LC+["logTP","logPup_flow_proxy"]}
        first=model_cv(good,"logTP",water)
        second=model_cv(good,"eqr",eco)
        comp={
            "real_2022_stream_lnTP":{
              "landcover_over_original_background":comparison(first,
                "B_geographic_hydrology","B_plus_local_2021_landcover"),
              "upstream_P_Q_after_landcover":comparison(first,
                "B_plus_local_2021_landcover","landcover_plus_UP_source_P_div_flow"),
              "downstream_P_Q_after_landcover":comparison(first,
                "B_plus_local_2021_landcover","landcover_plus_DOWN_source_P_div_flow_placebo")},
            "real_2022_benthic_EQR":{
              "landcover_over_original_background":comparison(second,
                "B_geographic_hydrology","B_plus_local_2021_landcover"),
              "real_river_TP_after_landcover":comparison(second,
                "B_plus_local_2021_landcover","landcover_plus_real_2022_river_TP"),
              "upstream_source_P_Q_after_landcover":comparison(second,
                "B_plus_local_2021_landcover","landcover_plus_UP_source_P_div_flow"),
              "downstream_source_P_Q_after_landcover":comparison(second,
                "B_plus_local_2021_landcover","landcover_plus_DOWN_source_P_div_flow_placebo"),
              "upstream_source_P_Q_after_landcover_and_real_river_TP":comparison(second,
                "landcover_plus_real_2022_river_TP","landcover_plus_real_TP_and_UP_source_P_div_flow")}}
        output["observational_2022_source_site_sample_n"]=len(good)
        output["source_5_region_2022_WorldCover_adjusted_MSE_predictions"]=comp
        output["outcome_source_2022_stream_lnTP_all_models"]=first
        output["outcome_source_2022_real_benthic_EQR_all_models"]=second
        output["status"]="PASS_POSTHOC_SOURCE_2021_LANDCOVER_CONFOUNDING_AUDIT_NONCAUSAL"
    else:
        output["status"]="STOP_SOURCE_2021_WORLDCOVER_2022_WATER_CONFOUNDING_GEOGRAPHY"
    (args.output/"AUSTRIA67_2022_SOURCE_WWTP_P_REAL_TP_EQR_2021_WORLDCOVER_LANDUSE_QA.json").write_text(
        json.dumps(output,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    summary={k:v["source_5block_equal_MSE_added_gain_pct"]
             for group,vals in output.get("source_5_region_2022_WorldCover_adjusted_MSE_predictions",{}).items()
             for k,v in vals.items()}
    print("GESUI_AUSTRIA67_WORLDCOVER2021_CONFUNDING_ADVERSARIAL_QA="+json.dumps({
      "status":output["status"],"original67_source":len(keep),
      "original_landcover_valid":len(good),
      "original5_geo_region_counts":dict(blocks),
      "background_vs_source_regional_P_Q_and_TP_gain":summary,
      "not_a_2022_vegetation_response_or_WWTP_causal_estimate":True
    },sort_keys=True,ensure_ascii=False,allow_nan=False),flush=True)

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    main(parser.parse_args())
