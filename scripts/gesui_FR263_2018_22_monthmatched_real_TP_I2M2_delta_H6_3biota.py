#!/usr/bin/env python3
"""Same ORIGINAL 263 observed 2018/2022 Peters TP/EEA I2M2 field station
pairs, 73 L7 basins, 24 L6 hydrologic groups.

New seasonal measurement sensitivity: equal weight to log median TP in
calendar months sampled in BOTH 2018 and 2022 (>=3 common months, all 263
passed independently SOURCE-ONLY before new outcome values were opened).

Compare ORIGINAL raw annual ΔTP versus month-aligned ΔTP on exact SAME
original physical site rows, outcomes, covariates, H6 spatial folds. All
posthoc, no effluent attribution, observed ecosystem service or causal lag.
"""
import collections,csv,datetime,hashlib,io,json,math,os,pathlib,statistics,urllib.request,zipfile
import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_FR263_2018_22_TP_MONTH_BALANCED"
ROOT.mkdir(parents=True,exist_ok=True)
MONTH_SHA="a02a9f724fb1b1671a9f248ed169620bc25db8b9a5914be566f6728e25b2c4f1"
ZIP_SHA="26f2af11f4972987ee6f7c3f497ad902820f145851f71919d4f45b589e1d809e"
ZIP_MD5="b596dbe9833c49426147a99f10f3fe31"
ZIP_BYTES=25003722
def readjson(name,sha=None):
    p=pathlib.Path(os.environ[name]);b=p.read_bytes()
    if sha is not None and hashlib.sha256(b).hexdigest()!=sha:raise RuntimeError("Original source SHA differs "+name)
    return json.loads(b)
m=readjson("PRE_OUTCOME_COMMON_MONTH_WQ_SOURCE",MONTH_SHA)
if m["schema"]!="GESUI_FR263_ORIGINAL_2018_2022_TP_COMMON_CALENDAR_MONTH_ELIGIBILITY_SOURCE_ONLY_v1":
    raise RuntimeError("Original month opportunities source data type differs")
source=m["original_263_physically_identical_2018_2022_WQ_and_EEA_source_matched_months"]
if len(source)!=263 or any(not z["source_eligible_at_least3_common_calendar_months"] for z in source):
    raise RuntimeError("Original month-matched eligibility not all 263 physical original source sites")
if len({z["original_HYBAS_ID_L07"] for z in source})!=73:
    raise RuntimeError("Source original 73 Level7 hydrologic basins drift")
if len({z["original_H6_parent"] for z in source})!=24:
    raise RuntimeError("Source original 24 parent hydrological regions drift")
oldreport=readjson("ORIGINAL_263_SAME_SITE_RAW_ANNUAL_TP_EQR_STUDY")
if oldreport["schema"]!="GESUI_FR_TRUE_263_SAME_SITE_I2M2_MEASURED_TP_2018_22_DELTA_73H7_24H6_3BIOTA_v1":
    raise RuntimeError("Original full annual measured TP/I2M2 exact source result schema drift")
def data_csv(env):
    with open(os.environ[env],newline="",encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))
prior_sites=data_csv("ORIGINAL_263_REAL_TP_I2M2_SITE_DATA")
prior_H7=data_csv("ORIGINAL_73_H7_TP_I2M2_GBIF_STUDY_DATA")
if len(prior_sites)!=263 or len(prior_H7)!=73:raise RuntimeError("Previous original TP and EQR biology frame not 263/73")
stationmap={z["Peters_identical_WQ_ID_2018_2022"]:z for z in prior_sites}
if len(stationmap)!=263 or set(stationmap)!={r["WQ_ID"] for r in source}:
    raise RuntimeError("Exact original source physical WQ ID mismatch with old TP/I2M2 site dataset")
h7map={z["HYBAS_ID_L07"]:z for z in prior_H7}
if len(h7map)!=73:raise RuntimeError("Original 73 H7 repeated-year output keys not unique")

# Publisher authenticated source, access only preselected 263 official station
# monthly date/value CSVs inside the archive.
url="https://ndownloader.figshare.com/files/70896288"
request=urllib.request.Request(url,headers={"User-Agent":"gesui-FR2018-22-month-matched-river-TP-I2M2-source/1","Accept":"application/zip,application/octet-stream,*/*"})
b=bytearray()
with urllib.request.urlopen(request,timeout=100) as rr:
    while True:
        chunk=rr.read(131072)
        if not chunk:break
        b.extend(chunk)
        if len(b)>27_000_000:raise RuntimeError("Original official Peters TP.2 archive expanded unexpectedly")
raw=bytes(b)
if len(raw)!=ZIP_BYTES or hashlib.md5(raw).hexdigest()!=ZIP_MD5 or hashlib.sha256(raw).hexdigest()!=ZIP_SHA:
    raise RuntimeError("Original month-aligned TP sample data source archive not publisher original")
archive=zipfile.ZipFile(io.BytesIO(raw))

month_sites=[]
for row in source:
    k=row["WQ_ID"]
    fname="TP.2_-_Harmonised/"+k+".csv"
    rr=csv.DictReader(io.StringIO(archive.read(fname).decode("utf-8-sig"),newline=""))
    if not {"date","value"}.issubset(rr.fieldnames or []):
        raise RuntimeError("Original source per-station measured TP date/value fields changed")
    vals={2018:collections.defaultdict(list),2022:collections.defaultdict(list)}
    all_annual={2018:[],2022:[]}
    for q in rr:
        try:
            day=datetime.date.fromisoformat(str(q["date"]).strip()[:10])
            val=float(q["value"])
        except (TypeError,ValueError):
            continue
        if day.year in (2018,2022) and math.isfinite(val) and val>0:
            vals[day.year][day.month].append(val)
            all_annual[day.year].append(val)
    old=stationmap[k]
    for y in (2018,2022):
        prev=float(old[f"real_original_{y}_measured_TP_median_mg_L"])
        got=float(np.median(all_annual[y]))
        if not math.isclose(prev,got,rel_tol=1e-9,abs_tol=1e-10):
            raise RuntimeError(f"Original annual TP per-station measured median independently failed {k}/{y}")
        original_num=int(row[f"source_{y}_original_positive_readings_n"])
        if len(all_annual[y])!=original_num:
            raise RuntimeError("Before-new-outcome month sampling original source count changed")
    shared=sorted(set(vals[2018])&set(vals[2022]))
    if shared!=row["source_common_original_calendar_months"] or len(shared)<3:
        raise RuntimeError("Original source-frozen >=3 common month membership not reproduced")
    monthly_med={y:{m:float(np.median(vals[y][m])) for m in shared} for y in (2018,2022)}
    delta=float(np.mean([math.log(monthly_med[2022][m])-math.log(monthly_med[2018][m]) for m in shared]))
    mean18=float(np.mean([math.log(monthly_med[2018][m]) for m in shared]))
    mean22=float(np.mean([math.log(monthly_med[2022][m]) for m in shared]))
    annual_delta=float(old["real_original_delta_log_TP_2022_minus_2018"])
    if not math.isclose(annual_delta,math.log(np.median(all_annual[2022]))-math.log(np.median(all_annual[2018])),rel_tol=1e-10,abs_tol=1e-12):
        raise RuntimeError("Original annual delta TP from original raw station source mismatch")
    if old["HYBAS_ID_L07"]!=row["original_HYBAS_ID_L07"] or old["EEA_identical_native_site_key_2018_2022"]!=row["original_EEA_identical_site_key"]:
        raise RuntimeError("Original source TP/EQR physical station or H7 group identity differs")
    matched18=sum(len(vals[2018][m]) for m in shared)
    matched22=sum(len(vals[2022][m]) for m in shared)
    if matched18!=row["source_2018_n_samples_in_common_months"] or matched22!=row["source_2022_n_samples_in_common_months"]:
        raise RuntimeError("Original source-only month-opportunity count changed")
    month_sites.append({
       "Peters_original_WQ_ID":k,
       "EEA_identical_station_key":row["original_EEA_identical_site_key"],
       "HYBAS_ID_L07":old["HYBAS_ID_L07"],
       "HYBAS_ID_L06":old["HYBAS_ID_L06"],
       "source_original_calendar_months_shared":",".join(map(str,shared)),
       "n_original_shared_calendar_months":len(shared),
       "n_original_2018_TP_dates_common_months":matched18,
       "n_original_2022_TP_dates_common_months":matched22,
       "actual_2018_common_month_equal_logTP":mean18,
       "actual_2022_common_month_equal_logTP":mean22,
       "actual_2018_to_2022_MONTHMATCHED_delta_logTP":delta,
       "actual_2018_to_2022_original_annual_delta_logTP":annual_delta,
       "delta_original_month_common_log1p_obs_effort":math.log1p(matched22)-math.log1p(matched18),
       "actual_2018_to_2022_same_site_delta_I2M2_EQR":float(old["real_original_delta_I2M2_EQR_2022_minus_2018"]),
       "original_I2M2_country_method_station_repeated_unmodified":True
    })
if len(month_sites)!=263:
    raise RuntimeError("Cannot construct complete original 263 physical 2018/2022 WQ month subset")
byh=collections.defaultdict(list)
for r in month_sites:byh[r["HYBAS_ID_L07"]].append(r)
if set(byh)!=set(h7map):raise RuntimeError("Original matched same 73 study H7 not preserved")
rows=[]
for h in sorted(byh):
    ss=byh[h]
    old=h7map[h]
    av=lambda k:float(np.mean([float(z[k]) for z in ss]))
    rec=dict(old)
    if len(ss)!=int(old["original_same_physical_EEA_I2M2_and_WQ_TP_station_n"]):
        raise RuntimeError("Month source primary model silently changed original physical station counts")
    rec["observed_original_common_month_logTP_delta_H7"]=av("actual_2018_to_2022_MONTHMATCHED_delta_logTP")
    rec["observed_original_annual_logTP_delta_H7_same_sites"]=av("actual_2018_to_2022_original_annual_delta_logTP")
    rec["observed_original_common_month_equal_logTP_2018"]=av("actual_2018_common_month_equal_logTP")
    rec["observed_original_common_month_equal_logTP_2022"]=av("actual_2022_common_month_equal_logTP")
    rec["original_COMMON_month_delta_log1p_TP_sample_effort"]=av("delta_original_month_common_log1p_obs_effort")
    rec["n_original_mean_shared_months_per_station"]=av("n_original_shared_calendar_months")
    if not math.isclose(float(old["mean_site_delta_log_original_measured_TP_2022_minus_2018"]),rec["observed_original_annual_logTP_delta_H7_same_sites"],rel_tol=1e-10,abs_tol=1e-12):
        raise RuntimeError("Previously recorded annual TP source per-H7 median differs on same source eligible sites")
    rec["observed_delta_EQR_same_site_H7"]=float(old["mean_site_delta_original_EQR_2022_minus_2018"])
    rows.append(rec)
assert len(rows)==73
fold=np.array([int(r["original_H6_source_fold"]) for r in rows])
h6=np.array([r["HYBAS_ID_L06"] for r in rows])
G,inv=np.unique(h6,return_inverse=True)
if len(G)!=24 or min(np.bincount(fold,minlength=5))<10:raise RuntimeError("Original 73 H7/24 H6 original source fold changed")
Y=np.array([r["observed_delta_EQR_same_site_H7"] for r in rows])
Tm=np.array([r["observed_original_common_month_logTP_delta_H7"] for r in rows])
Ta=np.array([r["observed_original_annual_logTP_delta_H7_same_sites"] for r in rows])
Xkeys=[
 "log_original_H7_local_area_km2","source_original_upstream_cropland_pct",
 "source_original_upstream_urban_pct","log_original_H7_outlet_natural_DIS_m3s",
 "mean_original_site_modelled_RIVER_Q_log1p",
 "log1p_original_same_site_2018_22_pair_n",
 "original_COMMON_month_delta_log1p_TP_sample_effort"]
X=np.array([[float(r[k]) for k in Xkeys] for r in rows])
Xyear=np.column_stack([X,[r["observed_original_common_month_equal_logTP_2018"] for r in rows]])
def pred(X,y,groups):
    p=np.empty(len(y))
    for name in sorted(set(groups)):
        tr=groups!=name;te=~tr
        if not np.any(tr) or not np.any(te):raise RuntimeError("Missing original hydrologic fold")
        A=X[tr];B=X[te];mu=A.mean(axis=0);sd=A.std(axis=0)
        sd=np.where(sd>1e-12,sd,1.)
        model=Ridge(alpha=10.,fit_intercept=True,solver="svd")
        model.fit((A-mu)/sd,y[tr]);p[te]=model.predict((B-mu)/sd)
    if not np.all(np.isfinite(p)):raise RuntimeError("Measured TP month-matched original parent holdout predictions invalid")
    return (y-p)**2
rng=np.random.default_rng(20261113)
boots=rng.integers(0,len(G),(2000,len(G)))
def gain(a,b):
    mse0=float(a.mean());mse1=float(b.mean())
    z0=np.bincount(inv,weights=a,minlength=len(G))[boots].sum(axis=1)
    z1=np.bincount(inv,weights=b,minlength=len(G))[boots].sum(axis=1)
    d=100*(z0-z1)/z0
    return {"original_baseline_MSE":mse0,
            "original_additional_MSE":mse1,
            "signed_gain_pct":100*(mse0-mse1)/mse0,
            "source_H6_bootstrap_2000_descriptive95_pct":np.percentile(d,[2.5,97.5]).tolist(),
            "source_H6_bootstrap_improving_fraction":float(np.mean(d>0))}
es0=pred(X,Y,fold)
esm=pred(np.column_stack([X,Tm]),Y,fold)
esa=pred(np.column_stack([X,Ta]),Y,fold)
e0LO=pred(X,Y,h6)
emLO=pred(np.column_stack([X,Tm]),Y,h6)
eaLO=pred(np.column_stack([X,Ta]),Y,h6)
esYear=pred(Xyear,Y,fold)
esmYear=pred(np.column_stack([Xyear,Tm]),Y,fold)
# Year-specific TP measurement abundance source has previously been used in
# 2018/2022 annual analyses; compute matched-month effort only now. The raw
# annual delta model below shares B0 and same 73 source H7 with month model.
def cluster_ols(x,y,feature_name):
    n=len(y)
    V=np.column_stack([np.ones(n),x]);mu=V[:,1:].mean(axis=0)
    sd=V[:,1:].std(axis=0);sd=np.where(sd>1e-12,sd,1)
    V[:,1:]=(V[:,1:]-mu)/sd
    coeff=np.linalg.lstsq(V,y,rcond=None)[0]
    resid=y-V@coeff
    invA=np.linalg.pinv(V.T@V)
    meat=np.zeros((V.shape[1],V.shape[1]))
    for i in range(len(G)):
        sel=inv==i;v=V[sel].T@resid[sel];meat+=np.outer(v,v)
    rk=np.linalg.matrix_rank(V)
    M=invA@meat@invA*(len(G)/(len(G)-1))*(n-1)/(n-rk)
    se=float(np.sqrt(max(M[-1,-1],0)))
    return {"term":feature_name,"standardized_delta_EQR_per_1SD_monthTP_delta":float(coeff[-1]),
        "source_H6_parent_cluster_descriptive95":[float(coeff[-1]-1.96*se),float(coeff[-1]+1.96*se)],
        "source_1SD_ΔlnTP":float(sd[-1]),"rank":int(rk)}
stationTm=np.array([r["actual_2018_to_2022_MONTHMATCHED_delta_logTP"] for r in month_sites])
stationE=np.array([r["actual_2018_to_2022_same_site_delta_I2M2_EQR"] for r in month_sites])
within=collections.defaultdict(list)
for s in month_sites:within[s["HYBAS_ID_L07"]].append(s)
a=[];b=[]
for H,ss in within.items():
    if len(ss)<2:continue
    xx=np.array([x["actual_2018_to_2022_MONTHMATCHED_delta_logTP"] for x in ss])
    yy=np.array([x["actual_2018_to_2022_same_site_delta_I2M2_EQR"] for x in ss])
    xx-=xx.mean(); yy-=yy.mean()
    a.append(float(xx@yy));b.append(float(xx@xx))
within_slope=float(sum(a)/sum(b))
# Whole-H7 source-index plant/butterfly delta, NOT physically co-located 2022
# field plants/pollinators.
terrestrial={}
for label,outcome,extra in [
 ("GBIF_Plantae_source_speciesKey_S50_change","delta_GBIF_original_plant_S50",
  ["delta_GBIF_plant_record_opportunity_log1p","delta_GBIF_Animalia_record_opportunity_log1p"]),
 ("GBIF_Lycaenidae_fivefamily_index_logodds_change","delta_GBIF_original_Lycaenidae_logodds",
  ["delta_GBIF_fivefamily_record_opportunity_log1p","delta_GBIF_Animalia_record_opportunity_log1p"])
]:
    if any(z[outcome] in (None,"","None") for z in rows):
        terrestrial[label]={"status":"source_missing_record_field_stop"}
        continue
    yy=np.array([float(z[outcome]) for z in rows])
    xx=np.array([[float(z[k]) for k in Xkeys+extra] for z in rows])
    e0=pred(xx,yy,fold);e1=pred(np.column_stack([xx,Tm]),yy,fold)
    l0=pred(xx,yy,h6);l1=pred(np.column_stack([xx,Tm]),yy,h6)
    terrestrial[label]={
      "n_original_H7_2GBIF_periods":len(yy),
      "source_H7_TPchange_GBIF_indexchange_Spearman":float(spearmanr(Tm,yy).statistic),
      "source_H6_group5fold_after_GBIF_opportunity_controls":gain(e0,e1),
      "source_24_H6_leaveparentout_after_GBIF_opportunity_controls":gain(l0,l1),
      "original_plant_butterfly_values_are_GBIF_index_not_services":True}
out={
 "schema":"GESUI_FR263_REAL_TP_2018_22_MONTHMATCHED_I2M2_EQR_73H7_24H6_v1",
 "source_original_common_calendar_month_WQ263_source_SHA256":MONTH_SHA,
 "source_original_author_TP2_zip_SHA256":ZIP_SHA,
 "source_original_2018_2022_same_physical_EEA_WQ_263":len(month_sites),
 "source_original_same_H7_basins":len(rows),
 "source_original_true_H6_parent_region_n":len(G),
 "source_original_five_H6_geoblock_H7_counts":np.bincount(fold,minlength=5).tolist(),
 "source_monthmatched_minimum_calendar_months":3,
 "source_monthmatched_original_shared_month_range":[min(r["n_original_shared_calendar_months"] for r in month_sites),max(r["n_original_shared_calendar_months"] for r in month_sites)],
 "source_original_original_annual_TP2018_2022_medians_ALL_263_original_sites_reproduced":True,
 "new_monthmatchedTP_vs_original_rawannualTP_H7_rank_Spearman":float(spearmanr(Tm,Ta).statistic),
 "new_monthmatchedTP_vs_original_rawannualTP_H7_Pearson":float(np.corrcoef(Tm,Ta)[0,1]),
 "original_H7_monthmatched_delta_logTP_vs_delta_I2M2_EQR_Spearman":float(spearmanr(Tm,Y).statistic),
 "original_H7_rawannual_delta_logTP_vs_delta_I2M2_EQR_Spearman_same_rows":float(spearmanr(Ta,Y).statistic),
 "source_sitelevel_monthmatched_delta_TP_vs_I2M2_EQR_Spearman_descriptive":float(spearmanr(stationTm,stationE).statistic),
 "source_within_H7_site_demeaned_monthmatched_deltaTP_to_deltaEQ_residual_slope":within_slope,
 "source_24_H6_clustered_monthdeltaTP_EQR_adjusted_OLS":cluster_ols(np.column_stack([X,Tm]),Y,"equal_common_calendar_months_delta_lnTP"),
 "new_monthmatch_primary_h6_geo5fold":gain(es0,esm),
 "same_source_population_rawannual_TPdelta_and_same_month_effort_baseline_h6_geo5fold":gain(es0,esa),
 "posthoc_new_monthmatch_h6_leaveoneparentout":gain(e0LO,emLO),
 "posthoc_same_source_rawannual_TPdelta_h6_leaveoneparentout":gain(e0LO,eaLO),
 "sensitivity_B0_plus_original_2018_monthmatched_baselineTP_then_month_delta":gain(esYear,esmYear),
 "same_H7_2017_19_2020_22_GBIF_plant_butterfly_change":terrestrial,
 "warnings":[
   "New month-aligned 2018/2022 TP is measured ambient river chemistry, NOT WWTP-specific P and no I2M2 causality is identified.",
   "The same 263 stations and 73 H7 units had a previous null raw annual TP change, the month-aligned test is explicitly POSTHOC and not external biological replication.",
   "2018 and2022 common calendar months reduce seasonal composition differences, but not climatic conditions, sample day mismatches or aquatic organism lag/threshold responses.",
   "Source 2018 and2022 I2M2 EEA EQR columns may differ in annual reporting or normalization even with identical method code and physical station ID.",
   "GBIF 2017–19/2020–22 plant rarefied speciesKey and Lycaenidae 5family record ratio are survey indices for whole H7, not direct field vegetation or insect ecosystem services.",
   "Four reported 2022 matched upstream WWTP P sources among 273 original TP/EQR stations, original source reporting sparse; no source-positive/untreated treatment contrast."
 ]}
with (ROOT/"FR263_SAME_SITE_I2M2_2018_22_MONTHMATCHED_TP_H6_BENTHIC_PLANT_BUTTERFLY_RESULTS.json").open("w") as f:
    json.dump(out,f,sort_keys=True,indent=2,allow_nan=False)
with (ROOT/"FR263_SOURCE_PHYSICAL_WQ_EEA_MATCHED_CALENDAR_MONTH_TP_EQR.csv").open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(month_sites[0]));w.writeheader();w.writerows(month_sites)
with (ROOT/"FR73_NATIVE_H7_MONTHMATCHED_TP_I2M2_GBIF_DELTAS.csv").open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print("FR263_2018_22_TP_MATCHED_MONTH_MACRO_START")
print(json.dumps({k:v for k,v in out.items() if k!="warnings"},sort_keys=True,allow_nan=False),flush=True)
print("FR263_2018_22_TP_MATCHED_MONTH_MACRO_END")
