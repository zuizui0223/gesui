#!/usr/bin/env python3
"""True historical observed 2018/2022 river TP delta versus exact same EEA
French I2M2 benthic EQR station delta, n263 same WQ+EEA native station IDs
across 73 H7/24 pre-frozen H6 hydrological region blocks.

Original TP.2 v5 source-archive MD5 verified. Dates/biological IDs and
geographic folds were SOURCE-FROZEN before any numerical 2018 TP/ΔEQR
model was opened. This is exploratory POSTHOC relative to known 2022
negative crosssection, NOT sewage-source causality or measured pollination.
"""
import collections,csv,datetime,hashlib,io,json,math,os,pathlib,statistics,urllib.request,zipfile
import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_FR_2018_2022_TRUE_TP_I2M2_WITHINSITE"
ROOT.mkdir(parents=True,exist_ok=True)
SOURCE_TP_ELIGIBLE_SHA="2190ac905326f4f2f4f278eb0debb4a7bfda7c4c71347d4941b54cc19d3de264"
H6_SOURCE_FOLD_SHA="294eb69b76cd4cf40f1022a5e83278529b1d2bbceac9a38931004610a23eb052"
ATLAS_SHA="ce6780e01d2b98ed8e852574adb827cb3125be792a164e8be9e2ebb917b3a686"
ZIP_ID=70896288
TP_ZIP_MD5="b596dbe9833c49426147a99f10f3fe31"
TP_ZIP_SHA="26f2af11f4972987ee6f7c3f497ad902820f145851f71919d4f45b589e1d809e"
TP_ZIP_SIZE=25003722
EEA={
 2018:("results/eea2018_independent_temporal_bqe_20261108/benthic_invertebrates_2018_WWTP_source_real_EQR.csv","0f30bfbc3dac13f6a182c8a910e9899c04778322"),
 2022:("results/eea2022_biology_bqe_20261108/BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv","325f199e985d9862ffee5b0edb205c61a7be3833")
}
WQ_DERIVED="results/eea2022_measured_nutrient_load_river_20261107/TP_2022_reported_wwtp_load_and_measured_river_station.csv"
WQ_DERIVED_BLOB="a75a4cf140de64fd6038079a3fac6756a461f3f8"
def fromjson(env,sha):
    b=pathlib.Path(os.environ[env]).read_bytes()
    if hashlib.sha256(b).hexdigest()!=sha:
        raise RuntimeError("Previously source-frozen file SHA invalid "+env)
    return json.loads(b)
source=fromjson("ORIGINAL_DATEWISE_2018_2022_TP_SOURCE_ELIGIBLE",SOURCE_TP_ELIGIBLE_SHA)
fold_source=fromjson("ORIGINAL_74H7_H6_PRE_OUTCOME_SOURCE_FOLDS",H6_SOURCE_FOLD_SHA)
atlas=fromjson("ORIGINAL_HYDROATLAS_L07_1219_SOURCE",ATLAS_SHA)["original_1219_H7_covariate_by_real_ID"]
source_sites=source["original_source_only_exact_same_EEA2018_2022_biological_site_and_annual_measured_TP_support"]
sites=[x for x in source_sites if x["source_2018_eligible"] and x["source_2022_eligible"]]
if len(source_sites)!=263 or len(sites)!=263 or len({x["original_WQ_ID"] for x in sites})!=263:
    raise RuntimeError("Original pre-numeric source 263 sites with both 2018 and 2022 real TP eligibility changed")
foldinfo=fold_source["matched_original_74_H7_fold_source"]
fold_by={x["H7"]:int(x["fold"]) for x in foldinfo}
parents={x["H7"]:x["H6"] for x in foldinfo}
if len(foldinfo)!=74 or len({parents[x["HYBAS_ID_L07"]] for x in sites})!=24:
    raise RuntimeError("True original H6 source only fold 74→24 group geography changed")
def blobsha(b):return hashlib.sha1(("blob "+str(len(b))+"\0").encode()+b).hexdigest()
def loadcsv(path,sha,expected):
    raw=pathlib.Path(path).read_bytes()
    if blobsha(raw)!=sha:raise RuntimeError("Previously source-verified EEA/WQ source CSV Git blob changed "+path)
    data=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"),newline="")))
    if len(data)!=expected:raise RuntimeError("Original source study station cohort size changed "+path)
    return data
bio={}
for yr,(filename,sha) in EEA.items():
    data=loadcsv(filename,sha,{2018:2634,2022:2475}[yr])
    bio[yr]={r["site_key"]:r for r in data}
    if len(bio[yr])!=len(data):
        raise RuntimeError("EEA original biological station duplicated in a year "+str(yr))
wq=loadcsv(WQ_DERIVED,WQ_DERIVED_BLOB,1291)
wq={r["WQ.ID"]:r for r in wq}
if len(wq)!=1291:raise RuntimeError("Original 2022 WQ station IDs duplicated")
with open(os.environ["ORIGINAL_454_H7_BIOTIC_AND_RECORDING_SOURCE"],encoding="utf-8-sig",newline="") as f:
    macro=list(csv.DictReader(f))
if len(macro)!=454 or len({(r["HYBAS_ID_L07"],r["PERIOD"]) for r in macro})!=454:
    raise RuntimeError("Original source H7 study frame changed")
macromap={(x["HYBAS_ID_L07"],x["PERIOD"]):x for x in macro}
if len(macromap)!=454:raise RuntimeError("Fixed original source basin period duplicated")
# Download exact original v5 author-certified archive, not estimated
# 2018 TP from previous 2022 data or from station coverage metadata.
url=f"https://ndownloader.figshare.com/files/{ZIP_ID}"
req=urllib.request.Request(url,headers={"User-Agent":"gesui-source-verified-2018-2022-TP-I2M2-real-longitudinal/1","Accept":"application/zip,application/octet-stream,*/*"})
data=bytearray()
with urllib.request.urlopen(req,timeout=100) as resp:
    while True:
        part=resp.read(131072)
        if not part:break
        data.extend(part)
        if len(data)>27_000_000:
            raise RuntimeError("Original Peters measured TP source archive unexpectedly large")
raw=bytes(data)
if len(raw)!=TP_ZIP_SIZE or hashlib.md5(raw).hexdigest()!=TP_ZIP_MD5 or hashlib.sha256(raw).hexdigest()!=TP_ZIP_SHA:
    raise RuntimeError("Original official 2018/2022 actually measured TP source archive exact MD5/SHA/bytes changed")
arc=zipfile.ZipFile(io.BytesIO(raw))
names=set(arc.namelist())
station=[]
invalid=collections.Counter()
for item in sites:
    wq_id=item["original_WQ_ID"]
    member="TP.2_-_Harmonised/"+wq_id+".csv"
    if member not in names:raise RuntimeError("Source-only eligible WQ TP datewise original file disappeared "+wq_id)
    rawzip=arc.read(member)
    rr=csv.DictReader(io.StringIO(rawzip.decode("utf-8-sig"),newline=""))
    if not {"date","value"}.issubset(rr.fieldnames or ()):
        raise RuntimeError("Genuine author harmonized TP.2 measurements fields changed")
    obs=collections.defaultdict(list)
    for q in rr:
        try:
            d=datetime.date.fromisoformat(str(q["date"]).strip()[:10])
            v=float(q["value"])
        except (TypeError,ValueError):
            continue
        if d.year in (2018,2022) and math.isfinite(v) and v>0:
            obs[d.year].append((d,v))
    for year in (2018,2022):
        raw_dates=obs[year]
        if len(raw_dates)<4 or len({d.month for d,v in raw_dates})<3 or len({d for d,v in raw_dates})<4:
            raise RuntimeError("Original source-frozen observed 2018/22 TP >=4 dates/3 months eligibility contradicted by original numeric zip")
    tp2018=float(statistics.median(v for _,v in obs[2018]))
    tp2022=float(statistics.median(v for _,v in obs[2022]))
    # Recompute and compare the original 2022 source TP median independently.
    published2022=float(wq[wq_id]["real_2022_station_median_mgL"])
    if not math.isclose(tp2022,published2022,rel_tol=2e-7,abs_tol=1e-7):
        raise RuntimeError(f"2018/2022 actual observed TP source reconstruction 2022 differs for {wq_id}: {tp2022}/{published2022}")
    r18=bio[2018].get(item["original_EEA2018_site_key"])
    r22=bio[2022].get(item["EEA2022_site_key"])
    if r18 is None or r22 is None:
        raise RuntimeError("Source-frozen 2018 or 2022 EEA physical monitoring station disappeared")
    if r18["site_key"]!=r22["site_key"]:
        raise RuntimeError("Original physical EEA key not identical across 2018–22; no hidden site replacements")
    if r18["countryCode"]!="FR" or r22["countryCode"]!="FR" or r18["procedureClassificationSystem"]!="I2M2" or r22["procedureClassificationSystem"]!="I2M2":
        raise RuntimeError("Original French EEA source BQE I2M2 fixed method changed")
    if r18["HYRIV_ID"]!=r22["HYRIV_ID"]:
        raise RuntimeError("Original 2018/22 EEA physical station HYRIV river ID changed")
    eqr18=float(r18["numeric_EQR_2018"]);eqr22=float(r22["normalised_EQR"])
    if not math.isfinite(eqr18) or not math.isfinite(eqr22) or eqr18<0 or eqr22<0 or eqr18>1.2 or eqr22>1.2:
        raise RuntimeError("Original 2018/2022 EEA I2M2 numeric scales incomparable or invalid")
    q18=float(r18["DIS_AV_CMS"]);q22=float(r22["DIS_AV_CMS"])
    if not math.isfinite(q18) or not math.isfinite(q22) or q18<=0 or q22<=0:
        raise RuntimeError("Original EEA site modeled flow unavailable")
    n18=len(obs[2018]);n22=len(obs[2022])
    h=item["HYBAS_ID_L07"]
    station.append({
      "EEA_identical_native_site_key_2018_2022":r22["site_key"],
      "Peters_identical_WQ_ID_2018_2022":wq_id,
      "HYBAS_ID_L07":h,"HYBAS_ID_L06":parents[h],"source_H6_original_geo_fold":fold_by[h],
      "real_original_2018_measured_TP_median_mg_L":tp2018,
      "real_original_2022_measured_TP_median_mg_L":tp2022,
      "real_original_delta_log_TP_2022_minus_2018":math.log(tp2022)-math.log(tp2018),
      "real_original_2018_EEA_I2M2_numeric_EQR":eqr18,
      "real_original_2022_EEA_I2M2_normalized_EQR":eqr22,
      "real_original_delta_I2M2_EQR_2022_minus_2018":eqr22-eqr18,
      "original_WQ_valid_2018_TP_observation_count":n18,
      "original_WQ_valid_2022_TP_observation_count":n22,
      "original_WQ_distinct_positive_2018_TP_months":len({d.month for d,v in obs[2018]}),
      "original_WQ_distinct_positive_2022_TP_months":len({d.month for d,v in obs[2022]}),
      "original_EEA2018_site_modelled_DIS_AV_CMS":q18,
      "original_EEA2022_site_modelled_DIS_AV_CMS":q22,
      "original_same_site_same_HYRIV_same_BQE_method_original_years":True
    })
assert len(station)==263
h7ids={r["HYBAS_ID_L07"] for r in station}
if len(h7ids)!=73 or len({r["HYBAS_ID_L06"] for r in station})!=24:
    raise RuntimeError("Original 263 sites 73 H7/24 H6 sampled data changed")
if len({r["EEA_identical_native_site_key_2018_2022"] for r in station})!=263:
    raise RuntimeError("Original site pairs reused across H7 or source dates")
by=collections.defaultdict(list)
for q in station:by[q["HYBAS_ID_L07"]].append(q)
def avg(a,key):
    return float(np.mean([float(r[key]) for r in a]))
def med(a,key):
    return float(np.median([float(r[key]) for r in a]))
rows=[]
for h in sorted(by):
    ss=by[h]
    A=atlas[h]
    if not all(A.get(k) is not None for k in ("crp_pc_use","urb_pc_use","dis_m3_pyr")):
        raise RuntimeError("Original actual HydroATLAS landuse/flow missing source")
    any_macro=macromap.get((h,"2020-2022"))
    if not any_macro:raise RuntimeError("Original source H7 with 2022 EEA I2M2 does not exist")
    both_gbif=((h,"2017-2019") in macromap and (h,"2020-2022") in macromap)
    M18=macromap.get((h,"2017-2019"))
    M22=macromap[(h,"2020-2022")]
    original_geo=dict(
       original_source_2020_22_H7_GBIF_supported_both_periods=both_gbif,
       mean_original_site_modelled_RIVER_Q_log1p=math.log1p((avg(ss,"original_EEA2018_site_modelled_DIS_AV_CMS")+avg(ss,"original_EEA2022_site_modelled_DIS_AV_CMS"))/2),
       log_original_H7_local_area_km2=float(M22["log_local_H7_area_km2"]),
       source_original_upstream_cropland_pct=float(A["crp_pc_use"]),
       source_original_upstream_urban_pct=float(A["urb_pc_use"]),
       log_original_H7_outlet_natural_DIS_m3s=math.log1p(float(A["dis_m3_pyr"])),
       log1p_original_same_site_2018_22_pair_n=math.log1p(len(ss)),
       mean_original_delta_log1p_TP_positive_sample_opportunities=float(np.mean([
          math.log1p(q["original_WQ_valid_2022_TP_observation_count"])-
          math.log1p(q["original_WQ_valid_2018_TP_observation_count"]) for q in ss
       ]))
    )
    r={
     "HYBAS_ID_L07":h,"HYBAS_ID_L06":parents[h],"original_H6_source_fold":fold_by[h],
     "original_same_physical_EEA_I2M2_and_WQ_TP_station_n":len(ss),
     "mean_site_delta_log_original_measured_TP_2022_minus_2018":avg(ss,"real_original_delta_log_TP_2022_minus_2018"),
     "mean_site_delta_original_EQR_2022_minus_2018":avg(ss,"real_original_delta_I2M2_EQR_2022_minus_2018"),
     "mean_original_2018_logTP":float(np.mean([math.log(q["real_original_2018_measured_TP_median_mg_L"]) for q in ss])),
     "mean_original_2022_logTP":float(np.mean([math.log(q["real_original_2022_measured_TP_median_mg_L"]) for q in ss])),
     "mean_original_2018_EQR":avg(ss,"real_original_2018_EEA_I2M2_numeric_EQR"),
     "mean_original_2022_EQR":avg(ss,"real_original_2022_EEA_I2M2_normalized_EQR"),
     **original_geo
    }
    if both_gbif:
        r.update({
         "delta_GBIF_original_plant_S50":float(M22["plant_S50_GBIF_speciesKey"])-float(M18["plant_S50_GBIF_speciesKey"]),
         "delta_GBIF_original_Lycaenidae_logodds":float(M22["lycaenidae_fivefamily_index_share_logodds"])-float(M18["lycaenidae_fivefamily_index_share_logodds"]),
         "delta_GBIF_plant_record_opportunity_log1p":float(M22["log1p_GBIF_Plantae_source_records"])-float(M18["log1p_GBIF_Plantae_source_records"]),
         "delta_GBIF_Animalia_record_opportunity_log1p":float(M22["log1p_GBIF_Animalia_source_records"])-float(M18["log1p_GBIF_Animalia_source_records"]),
         "delta_GBIF_fivefamily_record_opportunity_log1p":float(M22["log1p_GBIF_5family_butterfly_source_records"])-float(M18["log1p_GBIF_5family_butterfly_source_records"])
        })
    else:
        for key in ("delta_GBIF_original_plant_S50","delta_GBIF_original_Lycaenidae_logodds",
                    "delta_GBIF_plant_record_opportunity_log1p","delta_GBIF_Animalia_record_opportunity_log1p",
                    "delta_GBIF_fivefamily_record_opportunity_log1p"):r[key]=None
    rows.append(r)
assert len(rows)==73
H6,inv=np.unique([r["HYBAS_ID_L06"] for r in rows],return_inverse=True)
if len(H6)!=24:raise RuntimeError("Original H6 24 source group count changed")
F=np.array([r["original_H6_source_fold"] for r in rows],dtype=int)
if set(F)!=set(range(5)) or min(np.bincount(F,minlength=5))<10:
    raise RuntimeError("Original source-frozen H6 spatial holdout group support insufficient")
# One H7 mean response represents original equal physical site pairs once.
base_features=[
 "log_original_H7_local_area_km2",
 "source_original_upstream_cropland_pct",
 "source_original_upstream_urban_pct",
 "log_original_H7_outlet_natural_DIS_m3s",
 "mean_original_site_modelled_RIVER_Q_log1p",
 "log1p_original_same_site_2018_22_pair_n",
 "mean_original_delta_log1p_TP_positive_sample_opportunities"]
X0=np.array([[float(r[k]) for k in base_features] for r in rows],dtype=float)
Y=np.array([r["mean_site_delta_original_EQR_2022_minus_2018"] for r in rows])
T=np.array([r["mean_site_delta_log_original_measured_TP_2022_minus_2018"] for r in rows])
X1=np.column_stack([X0,T])
Xbase_2018=np.column_stack([X0,[r["mean_original_2018_logTP"] for r in rows]])
Xbase_2018_plusT=np.column_stack([Xbase_2018,T])
def cvpred(X,y,group):
    pred=np.empty(len(y))
    for k in sorted(set(group)):
        test=group==k;train=~test
        if not np.any(train) or not np.any(test):raise RuntimeError("Original 2018–2022 source group fold missing train/test")
        A=X[train];B=X[test]
        mu=A.mean(axis=0);sd=A.std(axis=0)
        sd=np.where(sd>1e-12,sd,1)
        r=Ridge(alpha=10.,fit_intercept=True,solver="svd")
        r.fit((A-mu)/sd,y[train])
        pred[test]=r.predict((B-mu)/sd)
    if not np.all(np.isfinite(pred)):raise RuntimeError("Original H6 geo holdout predicted source change nonfinite")
    return (y-pred)**2
rng=np.random.default_rng(20261113)
iid=rng.integers(0,len(H6),(2000,len(H6)))
def gain(e0,e1,subgroup_inv=inv):
    a=float(np.mean(e0));b=float(np.mean(e1))
    x=np.bincount(subgroup_inv,weights=e0,minlength=len(H6))
    y=np.bincount(subgroup_inv,weights=e1,minlength=len(H6))
    boot0=np.sum(x[iid],axis=1);boot1=np.sum(y[iid],axis=1)
    d=100*(boot0-boot1)/boot0
    return {"original_grouped_predictive_MSE_baseline":a,
            "original_grouped_predictive_MSE_with_delta_TP":b,
            "signed_MSE_improvement_pct":100*(a-b)/a,
            "H6_group_2000_bootstrap_descriptive95_pct":np.percentile(d,[2.5,97.5]).tolist(),
            "H6_group_bootstrap_positive_fraction":float(np.mean(d>0))}
def ols_beta(X,y):
    n=len(y);A=np.column_stack([np.ones(n),X])
    mu=A[:,1:].mean(axis=0);sd=A[:,1:].std(axis=0);sd=np.where(sd>1e-12,sd,1)
    A[:,1:]=(A[:,1:]-mu)/sd
    coef=np.linalg.lstsq(A,y,rcond=None)[0]
    res=y-A@coef
    br=np.linalg.pinv(A.T@A)
    meat=np.zeros((A.shape[1],A.shape[1]))
    for i in range(len(H6)):
        sel=inv==i; v=A[sel].T@res[sel];meat+=np.outer(v,v)
    rk=np.linalg.matrix_rank(A)
    var=br@meat@br*len(H6)/(len(H6)-1)*(n-1)/(n-rk)
    se=float(np.sqrt(max(var[-1,-1],0)))
    return {"beta_delta_EQR_per_1SD_delta_logTP":float(coef[-1]),
      "H6_parent_cluster_descriptive95": [float(coef[-1]-1.96*se),float(coef[-1]+1.96*se)],
      "original_delta_logTP_1SD":float(sd[-1]),
      "rank":int(rk)}
e0=cvpred(X0,Y,F);e1=cvpred(X1,Y,F)
elo0=cvpred(X0,Y,np.array([r["HYBAS_ID_L06"] for r in rows]))
elo1=cvpred(X1,Y,np.array([r["HYBAS_ID_L06"] for r in rows]))
eb0=cvpred(Xbase_2018,Y,F);eb1=cvpred(Xbase_2018_plusT,Y,F)
spr=spearmanr(T,Y)
# Site-level descriptive: not 263 independent river basins
site_delta_T=np.array([r["real_original_delta_log_TP_2022_minus_2018"] for r in station])
site_delta_E=np.array([r["real_original_delta_I2M2_EQR_2022_minus_2018"] for r in station])
site_spr=spearmanr(site_delta_T,site_delta_E)
within={}
for h,ss in by.items():
    if len(ss)<2:continue
    x=np.array([r["real_original_delta_log_TP_2022_minus_2018"] for r in ss])
    y=np.array([r["real_original_delta_I2M2_EQR_2022_minus_2018"] for r in ss])
    x=x-x.mean();y=y-y.mean()
    within[h]={"sxy":float(x@y),"sxx":float(x@x),"n":len(ss)}
xs=sum(r["sxy"] for r in within.values());xx=sum(r["sxx"] for r in within.values())
within_slope=xs/xx if xx>0 else None
rng2=np.random.default_rng(20261113)
bs=[]
h7keys=sorted(by)
for k in range(2000):
    hids=rng2.choice(h7keys,size=len(h7keys),replace=True)
    xy=sum(within[h]["sxy"] for h in hids if h in within)
    xx=sum(within[h]["sxx"] for h in hids if h in within)
    if xx>0:bs.append(xy/xx)
# Secondary original observational GBIF source index, not field services
terrestrial={}
good=[r for r in rows if r["original_source_2020_22_H7_GBIF_supported_both_periods"]]
for label,outcome,extra in [
 ("GBIF_Plant_speciesKey_50rarefied_S50_change","delta_GBIF_original_plant_S50",
    ["delta_GBIF_plant_record_opportunity_log1p","delta_GBIF_Animalia_record_opportunity_log1p"]),
 ("GBIF_Lycaenidae_5family_index_logodds_change","delta_GBIF_original_Lycaenidae_logodds",
    ["delta_GBIF_fivefamily_record_opportunity_log1p","delta_GBIF_Animalia_record_opportunity_log1p"])
]:
    q=good
    if len(q)<55 or len({r["HYBAS_ID_L06"] for r in q})<15:
        terrestrial[label]={"status":"SOURCE_STOP_TOO_FEW_TWO_YEAR_GBIF_RECORD_OPPORTUNITIES",
                            "n_original_same_H7_two_periods":len(q)}
        continue
    parentlist=np.array([r["HYBAS_ID_L06"] for r in q]);folds=np.array([r["original_H6_source_fold"] for r in q])
    XX=np.array([[float(z[k]) for k in base_features+extra] for z in q])
    TT=np.array([r["mean_site_delta_log_original_measured_TP_2022_minus_2018"] for r in q])
    YY=np.array([float(r[outcome]) for r in q])
    basee=cvpred(XX,YY,folds);adde=cvpred(np.column_stack([XX,TT]),YY,folds)
    baseLO=cvpred(XX,YY,parentlist);addLO=cvpred(np.column_stack([XX,TT]),YY,parentlist)
    corr=spearmanr(TT,YY)
    # Report 24 H6 parent block bootstrap only on source common nH7, no
    # repeated H7 pseudo-independent site-level responses.
    subset_H6,subinv=np.unique(parentlist,return_inverse=True)
    bb=np.random.default_rng(20261113).integers(0,len(subset_H6),(2000,len(subset_H6)))
    a=np.bincount(subinv,weights=basee,minlength=len(subset_H6))[bb].sum(axis=1)
    b=np.bincount(subinv,weights=adde,minlength=len(subset_H6))[bb].sum(axis=1)
    dv=100*(a-b)/a
    terrestrial[label]={
      "original_H7_both_periods_and_real_biological_TP_change":len(q),
      "original_source_hydrologic_H6_groups":len(subset_H6),
      "source_spearman_delta_logTP_vs_index_change":float(corr.statistic),
      "source_5fold_H6_group_cv_MSE_delta_TP_gain_pct":100*float(basee.mean()-adde.mean())/float(basee.mean()),
      "source_H6_2000_bootstrap_descriptive95_gain_pct":np.percentile(dv,[2.5,97.5]).tolist(),
      "source_24_parent_leaveout_MSE_delta_TP_gain_pct":100*float(baseLO.mean()-addLO.mean())/float(baseLO.mean()),
      "observed_GBIF_index_only_not_flowers_abundance_or_pollination":True}
summary={
 "schema":"GESUI_FR_TRUE_263_SAME_SITE_I2M2_MEASURED_TP_2018_22_DELTA_73H7_24H6_3BIOTA_v1",
 "original_pre_2018_22_measured_TP_eligibility_source_SHA256":SOURCE_TP_ELIGIBLE_SHA,
 "source_pre_outcome_native_H6_fold_SHA256":H6_SOURCE_FOLD_SHA,
 "source_publisher_measured_TP_zip_SHA256":TP_ZIP_SHA,
 "source_original_genuine_EEA_GitBlob_2018_2022":{yr:sha for yr,(name,sha) in EEA.items()},
 "source_original_fixed_263_identical_TP_and_EEA_I2M2_sites":len(station),
 "source_73_H7_equal_site_mean_independent_units":len(rows),
 "source_24_real_H6_parent_geographic_blocks":len(H6),
 "source_5fold_H6_sample_sizes":np.bincount(F,minlength=5).tolist(),
 "original_TP_and_I2M2_2018_2022_source_invalid_data_rejected":dict(invalid),
 "original_real_measured_TP_2018_median_range_mg_L":[min(q["real_original_2018_measured_TP_median_mg_L"] for q in station),max(q["real_original_2018_measured_TP_median_mg_L"] for q in station)],
 "original_real_measured_TP_2022_median_range_mg_L":[min(q["real_original_2022_measured_TP_median_mg_L"] for q in station),max(q["real_original_2022_measured_TP_median_mg_L"] for q in station)],
 "original_real_EQR_2018_range":[min(q["real_original_2018_EEA_I2M2_numeric_EQR"] for q in station),max(q["real_original_2018_EEA_I2M2_numeric_EQR"] for q in station)],
 "original_real_EQR_2022_range":[min(q["real_original_2022_EEA_I2M2_normalized_EQR"] for q in station),max(q["real_original_2022_EEA_I2M2_normalized_EQR"] for q in station)],
 "source_original_2022_Peters_reported_medians_recomputed_at_every_WQ_from_exact_original_dates":True,
 "primary_same_73_H7_basin_delta_TP_vs_delta_EQR_Spearman":{"rho":float(spr.statistic),"descriptive_p_not_independent_H6":float(spr.pvalue)},
 "secondary_site_raw_delta_Spearman_not263_independent_basins":{"rho":float(site_spr.statistic),"descriptive_p_not_independent_sites":float(site_spr.pvalue)},
 "within_H7_site_demeaned_delta_logTP_delta_EQR": {
    "eligible_H7_with_2plus_original_repeated_physical_sites":len(within),
    "eligible_EEA_original_physical_sites":sum(z["n"] for z in within.values()),
    "slope_delta_EQR_per_delta_logTP_within_H7":within_slope,
    "original_H7_2000bootstrap_descriptive95_slope":np.percentile(bs,[2.5,97.5]).tolist() if bs else None
 },
 "response_year_deltas":{
    "mean_siteH7_delta_logTP":float(T.mean()),"median_siteH7_delta_logTP":float(np.median(T)),
    "mean_siteH7_delta_EEA_I2M2_EQR":float(Y.mean()),"median_siteH7_delta_EEA_I2M2_EQR":float(np.median(Y))
 },
 "independent_geographic_real_H6_parent_out_of_region":{
    "original_fixed_H6_5fold_delta_TP_added_EQR_prediction":gain(e0,e1),
    "posthoc_leave_1_of_24_H6_parent_out_delta_TP_added_EQR_prediction":gain(elo0,elo1),
    "original_fixed_H6_5fold_2018_baseline_logTP_then_added_deltaTP_sensitivity":gain(eb0,eb1),
    "country_FR_and_EEA_benthic_method_I2M2_constant":True,
    "source_original_static_landuse_and_riverflow_covariates":base_features,
    "H6_clustered_adjusted_descriptive_delta_TP_delta_EQR_slope":ols_beta(X1,Y)
 },
 "original_same_H7_two_period_GBIF_index_change_vs_real_TP_change":terrestrial,
 "cautions":[
   "Observational time changes over 2018–22 among same river field sites are not causal TP interventions and can covary with seasonal weather, agriculture, local habitat and source monitoring differences.",
   "Although EEA I2M2 method labels agree, original 2018 numeric_EQR_2018 and 2022 normalised_EQR may not be completely identical reference/harmonization across collection years.",
   "The 2018 and 2022 TP are within-year sample medians across different dates, not simultaneous with EEA benthic sampling or controlled water-flow regimes.",
   "Whole H7 GBIF indexed Plant speciesKey-S50 and Lycaenidae relative five-family record observations are not sampled at the same physical river station or directly measured ecosystem services.",
   "The WWTP reported actual 2022 upstream P source is positive in only 4 of original 273 physical stations; this source gate FAILED, and zero recorded loads are NOT true untreated controls.",
   "H6 geographically blocked fold n24 all France single I2M2 study; no independent cross-country or continent transport claims.",
   "The 263 source and 73 H7 were selected outcome-blind for the new 2018→2022 changes but the research question is posthoc to known 2022 within-France TP-EQR source result."
 ]}
with (ROOT/"FR_REAL_2018_2022_I2M2_SAME_STATION_TP_DELTA_73H7_24H6_3BIOTA_RESULTS.json").open("w") as f:
    json.dump(summary,f,indent=2,sort_keys=True)
with (ROOT/"FR263_SAME_PHYSICAL_STATION_2018_2022_REAL_TP_AND_EEA_EQR.csv").open("w",newline="") as f:
    writer=csv.DictWriter(f,fieldnames=list(station[0]))
    writer.writeheader();writer.writerows(station)
with (ROOT/"FR73_HYBAS7_NATIVE_2018_2022_REAL_TP_I2M2_EQR_AND_GBIF_INDEX_DELTA.csv").open("w",newline="") as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]))
    writer.writeheader();writer.writerows(rows)
print("FR2018_2022_TRUE_TP_I2M2_DELTA_MACRO_START")
print(json.dumps({k:v for k,v in summary.items() if k!="cautions"},sort_keys=True),flush=True)
print("FR2018_2022_TRUE_TP_I2M2_DELTA_MACRO_END")
