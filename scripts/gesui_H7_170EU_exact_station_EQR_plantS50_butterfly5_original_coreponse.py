#!/usr/bin/env python3
"""SOURCE-FROZEN observational macro co-response: 2018-22 EEA EQI/EQR,
GBIF plant speciesKey S50 and Lycaenidae share of five butterfly families.

Original EEA site/method pairs, actual 234 H7 native basin source, new S50
and five-family source candidates were frozen BEFORE this ecological comparison.
Unit of inference = original H7 river basin, NOT 1047 repeated sites.
All results POSTHOC, observational, NOT a wastewater-mediated service chain.
"""
import collections,csv,hashlib,json,math,os,pathlib
import numpy as np
from scipy.stats import spearmanr
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_H7_170_EU_CO_RESPONSE"
ROOT.mkdir(parents=True,exist_ok=True)
CANDIDATES_SHA="07fa9744cd4e9574587a7b25c4ed4fa7f3183e377cac61c28723a4f24a86c875"
PAIRS_SHA="373a3a96ea790a7592cc4d10b0e4f9244fbdddab2c1af4901118792674874930"
BLOBS={
 "2018":("results/eea2018_independent_temporal_bqe_20261108/benthic_invertebrates_2018_WWTP_source_real_EQR.csv",
         "0f30bfbc3dac13f6a182c8a910e9899c04778322"),
 "2022":("results/eea2022_biology_bqe_20261108/BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv",
         "325f199e985d9862ffee5b0edb205c61a7be3833")
}
def require_file(path,expected):
    f=pathlib.Path(path)
    if hashlib.sha256(f.read_bytes()).hexdigest()!=expected:raise RuntimeError("Source snapshot SHA invalid "+str(f))
    return json.loads(f.read_text())
def gitsha(raw):
    return hashlib.sha1(("blob "+str(len(raw))+"\0").encode()+raw).hexdigest()
candidate_source=require_file(os.environ["EU234_FROZEN_CANDIDATES"],CANDIDATES_SHA)
same_source=require_file(os.environ["EU234_FROZEN_STATIONS"],PAIRS_SHA)
source_rows=candidate_source["original_source_eligible_basin_period_candidates"]
pairs=same_source["source_station_pairs"]
if len(source_rows)!=468 or len(pairs)!=1047:
    raise RuntimeError("Original frozen 234×2/1047 native biological source sample changed")
m={(r["HYBAS_ID_L07"],r["PERIOD"]):r for r in source_rows}
if len(m)!=468 or len({x[0] for x in m})!=234:
    raise RuntimeError("Original watershed/time unit duplicate")
idx=collections.defaultdict(list)
for r in pairs:
    idx[r["HYBAS_ID_L07"]].append(r)
if len(idx)!=174:raise RuntimeError("Original 174 unique method-comparable EEA river basins were altered")
eqr={}
for year,(path,digest) in BLOBS.items():
    raw=pathlib.Path(path).read_bytes()
    if gitsha(raw)!=digest:raise RuntimeError("Original EQR yearly source altered "+year)
    tab=list(csv.DictReader(raw.decode("utf-8-sig").splitlines()))
    expected_rows={"2018":2634,"2022":2475}[year]
    if len(tab)!=expected_rows:raise RuntimeError("Expected original biological source 2634/2475 differs")
    E={r["site_key"]:r for r in tab}
    if len(E)!=len(tab):raise RuntimeError("EEA original source monitoring site repeated within year")
    eqr[year]=E
def fnum(s):
    try:
        f=float(s)
        return f if math.isfinite(f) else None
    except (ValueError,TypeError):return None
out=[];excluded=collections.Counter()
for hid in sorted(idx):
    x=m[(hid,"2017-2019")]
    y=m[(hid,"2020-2022")]
    if not (x["original_5family_plant_S50_and_butterfly20_source_eligible"] and y["original_5family_plant_S50_and_butterfly20_source_eligible"]):
        excluded["new_original_GBIF_source_not_supported"]+=1
        continue
    # Preserve honest source missingness: predeclared fixed model requires
    # original indexed plant AND Animalia background in BOTH fixed periods.
    # Failed GBIF API source requests are NA, never biological zero, and no
    # source-absent row may be imputed to rescue an ecological coefficient.
    opportunity_keys=("original_GBIF_Animalia_record_opportunity",
                      "original_GBIF_Plantae_record_opportunity")
    if any(r.get(k) is None for r in (x,y) for k in opportunity_keys):
        excluded["missing_original_background_record_opportunity_source"]+=1
        continue
    country=collections.Counter(s["country"] for s in idx[hid])
    nation=sorted(country, key=lambda k:(-country[k],k))[0]
    vals=[]
    for s in idx[hid]:
        native=s["original_site_key"]
        u=eqr["2018"].get(native);v=eqr["2022"].get(native)
        if not u or not v or u["countryCode"]!=s["country"] or v["countryCode"]!=s["country"]:
            raise RuntimeError("Same original native station country/source mismatch")
        if u["procedureClassificationSystem"]!=s["method"] or v["procedureClassificationSystem"]!=s["method"]:
            raise RuntimeError("Same original native station method/source mismatch")
        if u["observedPropertyDeterminandBiologyEQRCode"]!=s["observed_bqe_code"] or v["observedPropertyDeterminandBiologyEQRCode"]!=s["observed_bqe_code"]:
            raise RuntimeError("Original EQR biological assessment indicator changed")
        p=fnum(u["numeric_EQR_2018"]);q=fnum(v["normalised_EQR"])
        if p is None or q is None:
            excluded["non_numeric_original_station_both_year"]+=1;continue
        if abs(p)>5 or abs(q)>5:raise RuntimeError("EQR abnormal source unit, not safe to mix")
        vals.append((p,q,s["country"],s["method"],native))
    if not vals:
        excluded["no_paired_original_finite_EQR"]+=1;continue
    def num(row,key):
        v=row.get(key)
        if v is None:raise RuntimeError("Source-record exposure field missing "+key)
        f=float(v)
        if not math.isfinite(f):raise RuntimeError("Nonfinite source variable "+key)
        return f
    for row in [x,y]:
        nbut=num(row,"original_source_5family_butterfly_records")
        nlyc=num(row,"original_source_lycaenidae_records")
        if nbut<20 or nlyc<0 or nlyc>nbut:raise RuntimeError("GBIF five family denominator inconsistent")
    # Haldane 0.5 pseudo only to handle a zero in a family count; fixed on whole cohort
    def logit(row):
        m=num(row,"original_source_5family_butterfly_records")
        l=num(row,"original_source_lycaenidae_records")
        return math.log((l+0.5)/(m-l+0.5))
    def dlog(key):
        return math.log1p(num(y,key))-math.log1p(num(x,key))
    static_flow=num(x,"registered_static_WWTP_outlet_modelled_flow_m3perday")
    area=num(x,"source_local_area_km2");uparea=num(x,"source_UP_AREA_KM2")
    if area<=0 or uparea<=0:raise RuntimeError("Invalid original watershed area")
    if static_flow!=num(y,"registered_static_WWTP_outlet_modelled_flow_m3perday"):
        raise RuntimeError("Static original modeled HydroWASTE discharge secretly varies between year windows")
    dS=num(y,"original_source_plant_S50")-num(x,"original_source_plant_S50")
    dLy=logit(y)-logit(x)
    dE=float(np.mean([b-a for a,b,*_ in vals]))
    out.append({
      "HYBAS_ID_L07":hid,"country":nation,"source_country_codes":";".join(sorted(country)),
      "n_original_same_site_both_year_EQR":len(vals),
      "n_original_method_labels":len({v[3] for v in vals}),
      "mean_EQR2018_same_station":float(np.mean([v[0] for v in vals])),
      "mean_EQR2022_same_station":float(np.mean([v[1] for v in vals])),
      "delta_EQR_2022_minus_2018":dE,
      "delta_plant_speciesKey_rarefied_S50":dS,
      "delta_lycaenidae_fivefamily_index_log_odds":dLy,
      "delta_lycaenidae_fivefamily_index_share":num(y,"original_source_lycaenidae_share_of_five_families")-num(x,"original_source_lycaenidae_share_of_five_families"),
      "delta_log1p_GBIF_Animalia_opportunity":dlog("original_GBIF_Animalia_record_opportunity"),
      "delta_log1p_GBIF_Plantae_opportunity":dlog("original_GBIF_Plantae_record_opportunity"),
      "delta_log1p_GBIF_fivefamily_opportunity":dlog("original_source_5family_butterfly_records"),
      "log_local_subbasin_area_km2":math.log(area),
      "log1p_modelled_WWTP_static_outflow_per_upstream_area":math.log1p(max(0.,static_flow)/(uparea+1.)),
      "log1p_original_EEA_station_count":math.log1p(len(vals)),
      "WWTP_static_same_both_periods_not_temporal_causal_dose":True
    })
if len(out)<120:raise RuntimeError("Predeclared 120 basin ecology gate failed after numeric outcome availability")
counts=collections.Counter(z["country"] for z in out)
if len(counts)<=5:raise RuntimeError("Predeclared distinct country source support <6")
if not all(math.isfinite(z["delta_EQR_2022_minus_2018"]) for z in out):
    raise RuntimeError("Missing original biological change")

# Exactly 1 row / source H7 basin, not 1047 pseudo-replicated river EQR comparisons
N=len(out)
if len({z["HYBAS_ID_L07"] for z in out})!=N:
    raise RuntimeError("Source one-basin-one-response grain violated")
Y=np.array([z["delta_EQR_2022_minus_2018"] for z in out])
p=np.array([z["delta_plant_speciesKey_rarefied_S50"] for z in out])
l=np.array([z["delta_lycaenidae_fivefamily_index_log_odds"] for z in out])
base_keys=["delta_log1p_GBIF_Animalia_opportunity","delta_log1p_GBIF_Plantae_opportunity",
 "delta_log1p_GBIF_fivefamily_opportunity","log_local_subbasin_area_km2",
 "log1p_modelled_WWTP_static_outflow_per_upstream_area","log1p_original_EEA_station_count"]
country=sorted(counts)
def build(keys):
    X=np.column_stack([
      np.ones(N),
      np.array([[z[k] for k in keys] for z in out],dtype=float),
      np.array([[float(z["country"]==c) for c in country[1:]] for z in out],dtype=float)
    ])
    return X
designs={
 "B_recording_country_area_static_flow":base_keys,
 "B_plus_plant_S50_change":base_keys+["delta_plant_speciesKey_rarefied_S50"],
 "B_plus_lycaenidae_share_change":base_keys+["delta_lycaenidae_fivefamily_index_log_odds"],
 "B_plus_plant_and_lycaenidae":base_keys+["delta_plant_speciesKey_rarefied_S50","delta_lycaenidae_fivefamily_index_log_odds"]
}
X={k:build(v) for k,v in designs.items()}
def fold(hid):return int(hashlib.sha256(hid.encode()).hexdigest()[:8],16)%5
folds=np.array([fold(x["HYBAS_ID_L07"]) for x in out])
if np.min(np.bincount(folds,minlength=5))<15:
    raise RuntimeError("Basin-disjoint cross-validation fold support inadequate")
def predict(XX):
    v=np.empty(N)
    for fold in range(5):
        test=folds==fold;train=~test
        Z=XX.copy()
        mu=Z[train,1:].mean(axis=0)
        sd=Z[train,1:].std(axis=0)
        sd=np.where(sd>1e-12,sd,1.)
        Z[:,1:]=(Z[:,1:]-mu)/sd
        A=Z[train];T=Z[test]
        penalty=np.eye(A.shape[1])*10;penalty[0,0]=0
        b=np.linalg.solve(A.T@A+penalty,A.T@Y[train])
        v[test]=T@b
    return v
mse={}
errs={}
for key,XX in X.items():
    pp=predict(XX)
    errs[key]=(Y-pp)**2
    mse[key]=float(np.mean(errs[key]))
baseline="B_recording_country_area_static_flow"
comparisons={}
for name in list(designs)[1:]:
    gain=100*(mse[baseline]-mse[name])/mse[baseline]
    rng=np.random.default_rng(20261112)
    ix=rng.integers(0,N,(2000,N))
    d=100*(errs[baseline][ix].mean(axis=1)-errs[name][ix].mean(axis=1))/errs[baseline][ix].mean(axis=1)
    comparisons[name]={"out_of_basin_5fold_MSE_gain_pct":float(gain),
       "bootstrap_original_basin_2000_descriptive95_gain_pct":np.percentile(d,[2.5,97.5]).tolist(),
       "bootstrap_positive_gain_fraction":float(np.mean(d>0))}
def rho(a,b):
    r=spearmanr(a,b)
    return float(r.statistic),float(r.pvalue)
corrs={}
for name,a,b in [
 ("delta_EQR_vs_delta_S50",Y,p),
 ("delta_EQR_vs_delta_lycaenidae_logodds",Y,l),
 ("delta_S50_vs_delta_lycaenidae_logodds",p,l)
]:
    r,pvalue=rho(a,b)
    corrs[name]={"spearman_r":r,"descriptive_asymptotic_p_nonindependent_rivers":pvalue}
def hc3(betaX,weights=None):
    X0=betaX.copy()
    if weights is None:weights=np.ones(len(Y))
    sw=np.sqrt(weights)
    A=X0*sw[:,None]
    z=Y*sw
    b=np.linalg.lstsq(A,z,rcond=None)[0]
    res=z-A@b
    XtX=np.linalg.pinv(A.T@A)
    leverage=np.sum(A*(A@XtX),axis=1)
    leverage=np.clip(leverage,0,.999999)
    G=A.T@((res/(1-leverage))[:,None]**2*A)
    V=XtX@G@XtX
    se=np.sqrt(np.maximum(V.diagonal(),0))
    # Coefficients raw units, also standardize by fixed sample SD of original predictor.
    found={}
    for var in ("delta_plant_speciesKey_rarefied_S50","delta_lycaenidae_fivefamily_index_log_odds"):
        if var not in designs["B_plus_plant_and_lycaenidae"]:continue
        pos=1+designs["B_plus_plant_and_lycaenidae"].index(var)
        rawbeta=float(b[pos]);ci=[float(rawbeta-1.96*se[pos]),float(rawbeta+1.96*se[pos])]
        v=np.array([a[var] for a in out],dtype=float)
        found[var]={"raw_beta":rawbeta,"HC3_descriptive95":ci,
            "beta_per_1SD_predictor":rawbeta*float(v.std(ddof=0)),
            "HC3_95_per_1SD_predictor":[t*float(v.std(ddof=0)) for t in ci]}
    return found
full=X["B_plus_plant_and_lycaenidae"]
fit=hc3(full)
# State unit = one H7 basin; original 7 countries heavily dominated by France.
weights=np.array([1/counts[z["country"]] for z in out],dtype=float)
weights/=weights.mean()
balanced=hc3(full,weights)
loo={}
for countryomit in counts:
    subset=[z for z in out if z["country"]!=countryomit]
    if len(subset)<25 or len({z["country"] for z in subset})<2:
        loo[countryomit]={"n":len(subset),"status":"insufficient_country_support"}
        continue
    rawY=np.array([z["delta_EQR_2022_minus_2018"] for z in subset])
    xx=full[np.array([z["country"]!=countryomit for z in out])]
    est=np.linalg.lstsq(xx,rawY,rcond=None)[0]
    loo[countryomit]={"n":len(subset),"country_left_out":countryomit,
         "coefs_unstandardized":{
            var:float(est[1+designs["B_plus_plant_and_lycaenidae"].index(var)])
            for var in ("delta_plant_speciesKey_rarefied_S50","delta_lycaenidae_fivefamily_index_log_odds")
         }}
# High-dimensional fixed country dummies with many tiny country groups create
# instability. These are diagnosis, not external prospective validation.
mode_fr=[z for z in out if z["country"]=="FR"]
non_fr=[z for z in out if z["country"]!="FR"]
subgroups={
 "FR": {"n":len(mode_fr),"delta_EQR_median":float(np.median([z["delta_EQR_2022_minus_2018"] for z in mode_fr])) if mode_fr else None},
 "non_FR":{"n":len(non_fr),"adequate_for_ecological_inference":len(non_fr)>=25}}
if len(mode_fr)>=20:
    subgroups["FR"]["unadjusted_spearman_delta_EQR_vs_S50"]=rho(
      [z["delta_EQR_2022_minus_2018"] for z in mode_fr],
      [z["delta_plant_speciesKey_rarefied_S50"] for z in mode_fr])[0]
    subgroups["FR"]["unadjusted_spearman_delta_EQR_vs_lycae"]=rho(
      [z["delta_EQR_2022_minus_2018"] for z in mode_fr],
      [z["delta_lycaenidae_fivefamily_index_log_odds"] for z in mode_fr])[0]
report={
 "schema":"GESUI_H7_170_EU_ORIGINAL_FIXED_SOURCE_EQR_PLANT_S50_BUTTERFLY5_COHANGE_v1",
 "exploratory_posthoc":True,
 "original_frozen_source_candidates_SHA256":CANDIDATES_SHA,
 "original_frozen_same_physical_station_method_SHA256":PAIRS_SHA,
 "EEA_original_GitBlob_Sha_by_year":{k:v[1] for k,v in BLOBS.items()},
 "sample":{"eligible_original_H7_basins":N,"native_2018_2022_stations":sum(x["n_original_same_site_both_year_EQR"] for x in out),
   "countries":dict(counts),"n_distinct_countries":len(counts),
   "fivefold_hash_basin_counts":np.bincount(folds,minlength=5).tolist(),
   "dropped_original_source_reasons":dict(excluded),
   "delta_EQR_median":float(np.median(Y)),
   "delta_EQR_mean":float(np.mean(Y)),
   "delta_plant_S50_median":float(np.median(p)),
   "delta_lycaenidae_logodds_median":float(np.median(l)),
   "n_in_FR":len(mode_fr),"n_outside_FR":len(non_fr)},
 "original_predictor_covariates":designs,
 "H7_disjoint_hash_5fold_ridge10_outcome_EQR_change_MSE":mse,
 "heldout_additional_predictor_gain":comparisons,
 "unadjusted_within_basin_co_response_Spearman":corrs,
 "all_basin_country_FEs_HC3_original_partial_slopes":fit,
 "all_country_equal_total_weight_HC3_sensitivity":balanced,
 "leave_one_country_out_recomputed_raw_OLS_sensitivity":loo,
 "France_vs_other_original_study_country_grain":subgroups,
 "not_actual_measured_ecosystem_services":True,
 "recording_adjustment_does_NOT_remove all GBIF_bias":True,
 "upstream_modelled_static_WWTP_pressure_not_time_change_or_station_receptor":True,
 "one_7country_posthoc_source_selection_not_pan_EU_evidence":True,
 "not_statistically_independent_country_or_larger_basin_blocks":True,
 "no_spatially_disjoint_regional_fold_claim":True,
 "limitations":[
    "Only exact original same EEA benthic stations with unchanged observed BQE method and country; euro monitoring representation heavily French.",
    "Original EEA 2018 numeric_EQR and 2022 normalized_EQR are distinct yearly source columns; procedural method matching reduces but does not eliminate harmonization differences.",
    "GBIF Plantae speciesKey 50-rarefied index richness is an opportunistic observation model, not standardized field flora.",
    "Five-family Lycaenidae share measures records, not true butterfly density, visit or pollination service.",
    "Basin hashed folds prevent same H7 basin in train/test but do not block geographically neighboring H7 basins; spatial transfer remains unvalidated.",
    "2015-era static HydroWASTE to BASIN OUTLET is not an observed 2018/22 WWTP effective dose at an interior EEA monitoring station.",
    "One basin contributes one delta; original 1047 station pairs are not independent 1047 catchments.",
    "No causal wastewater→river benthos→flora→butterfly mediation can be fit without shared actual nutrients, hosts or pollination observations."
 ]}
with (ROOT/"H7_170EU_EQR_PLANT_S50_LYCAENIDAE_ORIGINAL_TRUE_COHANGE_RESULTS.json").open("w") as f:
    json.dump(report,f,indent=2,sort_keys=True)
with (ROOT/"H7_170EU_ORIGINAL_BASIN_PAIRED_ECOLOGICAL_DATA.csv").open("w",newline="") as f:
    cols=list(out[0])
    writer=csv.DictWriter(f,fieldnames=cols)
    writer.writeheader();writer.writerows(out)
print("H7_170_EU_MACRO_CO_RESPONSE_RESULT_START")
brief={k:v for k,v in report.items() if k not in ("limitations","leave_one_country_out_recomputed_raw_OLS_sensitivity")}
brief["leave_one_country_out_recomputed_raw_OLS_sensitivity"]=loo
print(json.dumps(brief,sort_keys=True))
print("H7_170_EU_MACRO_CO_RESPONSE_RESULT_END")
