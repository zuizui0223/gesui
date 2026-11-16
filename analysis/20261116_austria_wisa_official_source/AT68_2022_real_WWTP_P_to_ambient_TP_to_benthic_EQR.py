#!/usr/bin/env python3
"""Same physical site Austria 2022: ORIGINAL annual WWTP P source ↔ actually
measured river TP ↔ macroinvertebrate normalized EQR, with downstream source
negative control and proper 5 geographic-block out-of-block predictions.

Causal sewage→water→biotic mediation is NOT identified. Source annual P t/yr
is hydrographically modeled aggregation of eligible EEA2022 positive-report
FW outlets within 50km, with missing reporting NOT absence. Actual stream TP
annual mean is total unfiltered phosphorus from original BMLUK national 2022
FW1, not amount of WWTP P delivered. Only same-site 2022 benthic EQR.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,json,math,re,sys,time,xml.etree.ElementTree as ET,zipfile
from pathlib import Path
import numpy as np
import requests
from scipy.stats import spearmanr,pearsonr

BMLUK_URL="https://www.bmluk.gv.at/dam/jcr%3A24629723-d7da-465f-9c1f-048783b885a5/Wasserg%C3%BCte-Jahresbericht_2020-2022_Ergebnistabellen.zip"
BMLUK_ZIP_SHA="6a0db0d0b8360c11f580333e5eabc3479c68fe8a5d623b87ef61f0185e5bd31a"
BMLUK_ZIP_BYTES=601117
BMLUK_FILE="FW1_Fließgewässer_Chemie_2020_2022.xlsx"
BMLUK_WORKBOOK_SHA="17de22544e9d043f83a928154847fb3fcd938894ce7007cfd06c000444b6eeb2"
BENTHIC_PATH=Path("results/eea2022_biology_bqe_20261108/BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv")
BENTHIC_GITBLOB="325f199e985d9862ffee5b0edb205c61a7be3833"
SCHEMA="GESUI_AT2022_SAME68_WWTP_REPORTED_P_TO_ACTUAL_TP_AND_EQR_GEOBLOCK_V1"
REGIONS={"15_3":9,"15_4":15,"15_5":21,"16_4":13,"16_5":10}
ALPHA=10.
P_CONVERSION=1e9/(365*24*3600*1000) # annual tonnes P / Q (m³/s) -> mg P/L if all delivered
PRIMARY_N=68
ORIG_UP_POS=53
ORIG_DN_POS=54
NS={"m":"http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
RELNS={"r":"http://schemas.openxmlformats.org/officeDocument/2006/relationships"}

def sha(data):return hashlib.sha256(data).hexdigest()
def blob(data):return hashlib.sha1(("blob "+str(len(data))+"\0").encode()+data).hexdigest()

def original_workbook_bytes():
    res=requests.get(BMLUK_URL,timeout=(15,85),headers={"User-Agent":"gesui-original-Austrian-2022-P-source-to-biota/20261116"})
    if res.status_code!=200 or not res.url.startswith("https://www.bmluk.gv.at/"):
        raise RuntimeError("Actual BMLUK source 2022 report not published from government host")
    if len(res.content)!=BMLUK_ZIP_BYTES or sha(res.content)!=BMLUK_ZIP_SHA:
        raise RuntimeError("Original official BMLUK 2020-22 ZIP source byte integrity failure")
    with zipfile.ZipFile(__import__("io").BytesIO(res.content)) as z:
        if BMLUK_FILE not in z.namelist() or z.testzip() is not None:
            raise RuntimeError("Official annual river FW1 source table missing/damaged")
        workbook=z.read(BMLUK_FILE)
    if sha(workbook)!=BMLUK_WORKBOOK_SHA:
        raise RuntimeError("Original FW1 source annual 2022 workbook SHA mismatch")
    return workbook

def source_FW1_100_same_station_TP(workbook):
    # Parse OOXML source as a ZIP with Python's XML stdlib, without
    # editing/recalculating the published national workbook.
    with zipfile.ZipFile(__import__("io").BytesIO(workbook)) as z:
        roots=ET.fromstring(z.read("xl/workbook.xml"))
        shs=roots.find("m:sheets",NS)
        if shs is None:raise RuntimeError("Official workbook has no original sheet list")
        names=[(s.attrib.get("name"),s.attrib.get("{"+RELNS["r"]+"}id")) for s in shs]
        if len([a for a,b in names if a=="FW_Chemie_Mittelwerte_2022"])!=1:
            raise RuntimeError("Original FW1 annual 2022 sheet disappeared")
        rid=dict(names)["FW_Chemie_Mittelwerte_2022"]
        rel=ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        candidates=[e.attrib["Target"] for e in rel if e.attrib.get("Id")==rid]
        if len(candidates)!=1:raise RuntimeError("Original source spreadsheet 2022 relationship non-unique")
        name=candidates[0].lstrip("/")
        if name.startswith("xl/"):sheetname=name
        else:sheetname="xl/"+name
        if sheetname not in z.namelist():
            raise RuntimeError("Original 2022 FW1 worksheet target not actually in OOXML source")
        shared=[]
        if "xl/sharedStrings.xml" in z.namelist():
            root=ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in root.findall("m:si",NS):
                shared.append("".join(e.text or "" for e in si.iter() if e.tag.endswith("}t")))
        sheet=ET.fromstring(z.read(sheetname))
        content={}
        for row in sheet.findall(".//m:sheetData/m:row",NS):
            rr=int(row.attrib["r"])
            if rr not in set(range(1,104)):continue
            data={}
            for c in row.findall("m:c",NS):
                coordinate=c.attrib.get("r","")
                m=re.match(r"([A-Z]+)(\d+)",coordinate)
                if not m:continue
                col=m.group(1)
                typ=c.attrib.get("t")
                v=c.find("m:v",NS)
                if v is not None and v.text is not None:
                    if typ=="s":data[col]=shared[int(v.text)]
                    elif typ in ("str","inlineStr","e"):data[col]=str(v.text)
                    else:
                        try:data[col]=float(v.text)
                        except ValueError:data[col]=v.text
                elif typ=="inlineStr":
                    data[col]="".join(t.text or "" for t in c.iter() if t.tag.endswith("}t"))
            content[rr]=data
    header=content.get(3,{})
    if header.get("A")!="Messstelle" or "Gesamtphosphor-P" not in str(header.get("R")) or "unfiltriert" not in str(header.get("R")):
        raise RuntimeError("Original national FW1 source 2022 R-column is NOT unfiltered total-P annual mean")
    out={}
    for row in range(4,104):
        q=content.get(row,{})
        site=q.get("A")
        if not isinstance(site,str) or not site.startswith("FW"):
            raise RuntimeError("Original FW1 annual station source row ID missing "+str(row))
        if site in out:raise RuntimeError("Original 2022 FW1 river station ID duplicated")
        p=q.get("R")
        if not isinstance(p,(float,int)) or not math.isfinite(p) or p<=0:
            raise RuntimeError("Original FW1 measured total-P annual mean invalid "+str((site,p)))
        out[site]=float(p)
    if len(out)!=100:
        raise RuntimeError("Original published FW1 total-P river stations are not 100")
    return out

def source_benthic():
    raw=BENTHIC_PATH.read_bytes()
    if blob(raw)!=BENTHIC_GITBLOB:
        raise RuntimeError("Original EEA 2022 benthic EQR source Git blob mismatch")
    with BENTHIC_PATH.open(encoding="utf-8-sig",newline="") as f:
        source=list(csv.DictReader(f))
    if len(source)!=2475:
        raise RuntimeError("Original 2022 Europe measured benthic 2475 station cohort missing")
    a=[z for z in source if z["countryCode"]=="AT"]
    ids=[z["site_key"] for z in a]
    if len(a)!=68 or len(set(ids))!=68:
        raise RuntimeError("Original same-station observed Austrian EQR68 source drift")
    if collections.Counter(z["region_3deg_block"] for z in a)!=REGIONS:
        raise RuntimeError("Original genuine five region-block cohort changed")
    if {z["procedureClassificationSystem"] for z in a}!={"MI-RW-AT"}:
        raise RuntimeError("Original country 2022 benthic source method changed")
    return a

def cohort():
    p=source_FW1_100_same_station_TP(original_workbook_bytes())
    b=source_benthic()
    o=[]
    for r in b:
        parts=r["site_key"].split("|")
        if (len(parts)!=3 or parts[0]!="AT" or parts[2]!="euMonitoringSiteCode"
             or not parts[1].startswith("ATFW")):
            raise RuntimeError("Original monitored EEA site not exact national same station code")
        station=parts[1][2:]
        if station not in p:raise RuntimeError("Original actually observed national 2022 TP not at one BQE site "+station)
        try:
            y=float(r["normalised_EQR"]);q=float(r["DIS_AV_CMS"]);area=float(r["UPLAND_SKM"])
            up=float(r["upstream_2022_EEA_P_measured_t_per_yr"])
            down=float(r["downstream_2022_EEA_P_t_per_yr_placebo"])
            lat=float(r["lat"]);lon=float(r["lon"])
        except (ValueError,TypeError) as e:
            raise RuntimeError("Original measured EQR and source geography missing "+station) from e
        if (not all(map(math.isfinite,[y,q,area,up,down,lat,lon])) or
            y<0 or y>1.2 or q<=0 or area<=0 or up<0 or down<0):
            raise RuntimeError("Original EQR, reported P or river geometry invalid "+station)
        o.append({"station":station,"block":r["region_3deg_block"],
             "eqr":y,"TP_mgP_L":p[station],"logTP":math.log(p[station]),
             "logQ":math.log1p(q),"logArea":math.log1p(area),
             "lat":lat,"lon":lon,"Pup_t_yr":up,"Pdown_t_yr":down,
             "logPup":math.log1p(up),"logPdown":math.log1p(down),
             "logPup_flow_proxy":math.log1p(P_CONVERSION*up/q)})
    if len(o)!=68 or len({x["station"] for x in o})!=68:
        raise RuntimeError("Source 68 true source ecology join cardinality changed")
    if sum(x["Pup_t_yr"]>0 for x in o)!=ORIG_UP_POS or sum(x["Pdown_t_yr"]>0 for x in o)!=ORIG_DN_POS:
        raise RuntimeError("Original up/down-positive EEA source support n53/n54 changed")
    if abs(sum(x["eqr"] for x in o)-42.74)>1e-8:
        raise RuntimeError("Original published 2022 Austrian 68 benthic EQR source total changed")
    return o

def full_ols(d,y_key,base,extra):
    order=base+[extra]
    Y=np.array([r[y_key] for r in d],float)
    X=np.array([[r[k] for k in order] for r in d],float)
    mu=X.mean(axis=0);scale=X.std(axis=0,ddof=0)
    if (scale<1e-12).any():
        raise RuntimeError("Full-source regression column has zero variance")
    X=(X-mu)/scale
    xx=np.column_stack([np.ones(len(d)),X])
    coef=np.linalg.lstsq(xx,Y,rcond=None)[0]
    resid=Y-xx@coef
    return {"beta_focal_predictor_per_source_SD":float(coef[-1]),
        "in_sample_R2_adjusted_for_source_covariates":float(1-resid@resid/((Y-Y.mean())@(Y-Y.mean()))),
        "n":len(Y),"n_numeric_covariates_including_focal":len(order)}

def ridge_predict(train,test,features,y_key):
    X=np.array([[v[k] for k in features] for v in train],float)
    Xtest=np.array([[v[k] for k in features] for v in test],float)
    mu=X.mean(axis=0);sd=X.std(axis=0,ddof=0)
    sd[sd<1e-12]=1.
    x=(X-mu)/sd;xt=(Xtest-mu)/sd
    Y=np.array([v[y_key] for v in train],float)
    XX=np.column_stack((np.ones(len(train)),x))
    XT=np.column_stack((np.ones(len(test)),xt))
    pen=np.diag([0.]+[ALPHA]*len(features))
    b=np.linalg.solve(XX.T@XX+pen,XX.T@Y)
    return XT@b

def cv(d,y,models):
    regions=sorted(REGIONS)
    fold={}
    for group in regions:
        train=[z for z in d if z["block"]!=group]
        test=[z for z in d if z["block"]==group]
        if len(test)!=REGIONS[group]:raise RuntimeError("Unchanged source geography changed")
        yy=np.array([z[y] for z in test],float)
        fold[group]={"source_test_sites":len(test)}
        for model,features in models.items():
            pred=ridge_predict(train,test,features,y)
            fold[group][model]=float(np.mean((yy-pred)**2))
    results={}
    for key in models:
        block_mse=[fold[g][key] for g in regions]
        results[key]={"source_H6_proxy_3degree_region_equal_heldout_MSE":float(np.mean(block_mse)),
                      "source_each_region_test_MSE":{g:fold[g][key] for g in regions}}
    for model,ref in [
        ("M1_upstream_P","M0"),
        ("M2_downstream_placebo_P","M0"),
        ("M3_upstream_P_plus_downstream","M0"),
        ("M1_reported_upstream_P","M0"),
        ("M2_downstream_reported_P","M0"),
        ("M3_measured_river_TP","M0"),
        ("M4_stream_TP_plus_upstream_P","M3_measured_river_TP"),
        ("M5_stream_TP_plus_downstream_placebo","M3_measured_river_TP"),
        ("M6_upstream_source_flow_proxy","M0"),
        ("M4_source_P_flow_proxy","M0"),
        ("M5_source_P_flow_proxy_after_TP","M3_measured_river_TP")
    ]:
        if model not in models or ref not in models:continue
        old=results[ref]["source_H6_proxy_3degree_region_equal_heldout_MSE"]
        new=results[model]["source_H6_proxy_3degree_region_equal_heldout_MSE"]
        diffs={g:100*(fold[g][ref]-fold[g][model])/fold[g][ref] for g in regions}
        results[model]["source_added_over_"+ref+"_MSE_gain_pct"]=float(100*(old-new)/old)
        results[model]["source_5_geographic_regions_added_gain_pct"]=diffs
        results[model]["source_regions_with_positive_gain"]=sum(x>0 for x in diffs.values())
    return {"outcome":y,"original_n":len(d),"original_n_distinct_geographic_regions":len(regions),
       "fixed_Ridge_alpha":ALPHA,"region_group_sizes":REGIONS,"models":results}

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    d=cohort()
    baseline=["logQ","logArea","lat","lon"]
    water={
       "M0":baseline,
       "M1_upstream_P":baseline+["logPup"],
       "M2_downstream_placebo_P":baseline+["logPdown"],
       "M3_upstream_P_plus_downstream":baseline+["logPup","logPdown"],
       "M4_source_P_flow_proxy":baseline+["logPup_flow_proxy"]}
    eco={
       "M0":baseline,
       "M1_reported_upstream_P":baseline+["logPup"],
       "M2_downstream_reported_P":baseline+["logPdown"],
       "M3_measured_river_TP":baseline+["logTP"],
       "M4_stream_TP_plus_upstream_P":baseline+["logTP","logPup"],
       "M5_stream_TP_plus_downstream_placebo":baseline+["logTP","logPdown"],
       "M6_upstream_source_flow_proxy":baseline+["logPup_flow_proxy"],
       "M5_source_P_flow_proxy_after_TP":baseline+["logTP","logPup_flow_proxy"]}
    waterres=cv(d,"logTP",water)
    eqrres=cv(d,"eqr",eco)
    tp=np.array([z["TP_mgP_L"] for z in d],float)
    eq=np.array([z["eqr"] for z in d],float)
    pup=np.array([z["Pup_t_yr"] for z in d],float)
    pdown=np.array([z["Pdown_t_yr"] for z in d],float)
    from scipy.stats import pearsonr,spearmanr
    stats={
      "observed_2022_stream_TP_and_EQR_Spearman":float(spearmanr(tp,eq).statistic),
      "observed_2022_stream_lnTP_and_EQR_Pearson":float(pearsonr(np.log(tp),eq).statistic),
      "upstream_reported_annual_P_vs_measured_stream_TP_Spearman":float(spearmanr(pup,tp).statistic),
      "downstream_2022_placebo_P_vs_measured_stream_TP_Spearman":float(spearmanr(pdown,tp).statistic),
      "upstream_reported_annual_P_vs_EQR_Spearman":float(spearmanr(pup,eq).statistic),
      "upstream_vs_downstream_reported_P_Spearman":float(spearmanr(pup,pdown).statistic),
      "source_reported_upstream_positive_count":int((pup>0).sum()),
      "source_reported_downstream_positive_count":int((pdown>0).sum())}
    if not (-.591 <stats["observed_2022_stream_TP_and_EQR_Spearman"] < -.581):
        raise RuntimeError("Original Austrian 68 actual TP EQR independent source association could not be reproduced")
    result={"schema":SCHEMA,
      "status":"AUTHENTIC_OBSERVATIONAL_SOURCE_WATER_NUTRIENT_BENTHIC_DIRECT_CHAIN_FITTED_NO_CAUSAL_CLAIM",
      "original_sources":{
        "genuine_2020_2022_BMLUK_Fw1_zip_SHA256":BMLUK_ZIP_SHA,
        "FW1_2022_annual_mean_TP_xlsx_SHA256":BMLUK_WORKBOOK_SHA,
        "original_EEA2022_measured_benthic_file_gitblob":BENTHIC_GITBLOB,
        "original_source_2022_TRUE_annual_reported_P_50km_upstream_graph":True,
        "original_source_zero_means_no_matched_positive_both_measure_report_not_unpolluted":True},
      "source_2022_same_station_n":len(d),
      "source_2022_P_source_positive":int((pup>0).sum()),
      "original_five_spatial_3degree_region_groups":REGIONS,
      "observational_raw_diagnostics":stats,
      "source_total_P_2022_annual_mean_range_mgP_L":[float(np.min(tp)),float(np.max(tp))],
      "source_OLS_adjusted_descriptive_betas":{
        "upstream_reported_P_to_actual_2022_water_logTP":full_ols(d,"logTP",baseline,"logPup"),
        "downstream_placebo_to_actual_2022_water_logTP":full_ols(d,"logTP",baseline,"logPdown"),
        "upstream_flow_normalized_P_to_actual_2022_water_logTP":full_ols(d,"logTP",baseline,"logPup_flow_proxy"),
        "upstream_reported_P_to_2022_benthic_EQR":full_ols(d,"eqr",baseline,"logPup"),
        "actual_stream_TP_to_2022_benthic_EQR":full_ols(d,"eqr",baseline,"logTP"),
        "upstream_reported_P_after_actual_stream_TP_to_2022_benthic_EQR":full_ols(d,"eqr",baseline+["logTP"],"logPup")},
      "regional_spatial_holdout_for_real_stream_2022_logTP":waterres,
      "regional_spatial_holdout_for_real_benthic_2022_EQR":eqrres,
      "original_5_region_blocks_too_few_for_confirmatory_cluster_test":True,
      "all_original_2022_at_stations_only_and_actual_stream_2022_TP":True,
      "source_effect_noncausal_reporting_and_landuse_confounded":True,
      "downstream_2022_reported_mass_placebo_must_be_compared":True,
      "source_site_2022_annual_mean_P_not_measured_WWTP_contribution":True,
      "source_2022_benthic_EQR_not_macroinvertebrate_species_richness_biomass":True,
      "no_corresponding_local_2022_riparian_plant_or_butterfly_data_or_EcosystemService":True,
      "no_2022_wastewater_to_land_or_foodweb_effect_claimed":True}
    path=a.output/"AUSTRIA68_2022_REPORTED_WWTP_P_ACTUAL_STREAM_TP_EQR_REACHCHAIN.json"
    path.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    key={
       "n":len(d),"positive_upstream_P":sum(z["Pup_t_yr"]>0 for z in d),
       "spearman":stats,
       "TP_models":{k:v.get("source_added_over_M0_MSE_gain_pct")
                   for k,v in waterres["models"].items()},
       "EQ_models":{k:{"baseline_gain":v.get("source_added_over_M0_MSE_gain_pct"),
                      "incremental_after_TP":v.get("source_added_over_M3_measured_river_TP_MSE_gain_pct"),
                      "positive_geographic_regions":v.get("source_regions_with_positive_gain")}
                     for k,v in eqrres["models"].items()},
       "source_WWTP_P_to_waterTP_adjusted_beta":result["source_OLS_adjusted_descriptive_betas"]["upstream_reported_P_to_actual_2022_water_logTP"],
       "source_P_EQR_after_ambient_TP_beta":result["source_OLS_adjusted_descriptive_betas"]["upstream_reported_P_after_actual_stream_TP_to_2022_benthic_EQR"],
       "not_source_identified_causal":True}
    print("GESUI_AT68_2022_WWTP_P_MEASURED_TP_BENTHIC_CHAIN="+json.dumps(key,ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
