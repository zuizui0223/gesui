#!/usr/bin/env python3
"""FR178 POSTHOC 2-scene fixed-calendar-month S2 mosaic near–far NDVI × original measured river TP.

Entirely exploratory: after original single-scene 116/178 source STOP. Frozen
2-scene mosaics are an explicitly posthoc technical sensitivity (not a rescue
of the preregistered main null and not independent confirmation).
Minimum 120/40/15 and spatial bands remain unchanged. Source eligibility
before this numeric fit. STOP on a failed original source gate, never fit 14
pre-screened sites or change a radius/mask after viewing effects.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,statistics,zipfile,datetime
from pathlib import Path
from collections import defaultdict,Counter
import numpy as np
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge

PIX_SCHEMA="GESUI_FR178_S2_FROZEN_2SCENE_SAME_MONTH_CLOUD_MOSAIC_SOURCE_QA_POSTHOC_V1"
SOURCE_ID_SHA="76ba45169e58901a2cd82a18244d022773074963a6ac983b527f06302a4aedb8"
TP_ZIP_SHA="26f2af11f4972987ee6f7c3f497ad902820f145851f71919d4f45b589e1d809e"
TP_ZIP_MD5="b596dbe9833c49426147a99f10f3fe31"
TP_ZIP_SIZE=25003722
ORIGINAL_263_NUMERIC_SHA="1824b0d9b84d35a460b4d1d0513659f93a16c59d7f5666a2d8707be1c6da9942"
ORIGINAL_73_COVAR_SHA="f2ff4e56f245e1aa8f84f5ac099d68ec1c7e27a786c1067439db388b186dd90b"
MIN_SITES=120
MIN_H7=40
MIN_H6=15
NUM_BOOT=2000
SEED=20261117
SCHEMA="GESUI_FR178_2018_2022_TWO_SCENE_MOSAIC_REAL_TP_H6_EXPLORATORY_V1"

def sha(b):return hashlib.sha256(b).hexdigest()
def load_csv(path,expected):
    raw=path.read_bytes()
    if sha(raw)!=expected:raise RuntimeError("Original source numeric SHA mismatch: "+str(path))
    return list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))

def pixel_source(path):
    raw=[]
    for shard in range(4):
        p=path/f"FR178_2SCENE_MOSAIC_SOURCE_SHARD{shard}.json"
        q=json.loads(p.read_text())
        if q.get("schema")!=PIX_SCHEMA or q.get("original_source_shard")!=shard or q.get("original_source_178_pop")!=178:
            raise RuntimeError("Signed source true-riverline original shard mismatch")
        raw+=q["observations"]
    if len(raw)!=178 or len({x["site_key"] for x in raw})!=178:
        raise RuntimeError("Four exact original 178 source IDs missing or repeated")
    if sha("\n".join(sorted(x["site_key"] for x in raw)).encode())!=SOURCE_ID_SHA:
        raise RuntimeError("Original parent cohort not exactly preselected 178")
    eligible=[x for x in raw if x["source_two_years_real_mosaic_pixel_eligible"] and
              len(x["posthoc_composite_delta_NDVI_by_month"])>=2]
    h7=len({x["H7"] for x in eligible})
    h6=len({x["H6"] for x in eligible})
    qa={"original_178_input_source_site_n":178,
        "posthoc_2scene_fixed_original_TP_compatible_TWICE_NDVI_eligible_sites":len(eligible),
        "eligible_distinct_H7":h7,"eligible_distinct_H6":h6,
        "predeclared_minima":{"n_physical_sites":MIN_SITES,"n_native_H7":MIN_H7,"n_native_H6":MIN_H6},
        "all_178_original_site_IDs_authentic":True,
        "no_original_aquatic_concentration_or_EQR_opened_for_selecting_pixel_population":True,
        "posthoc_technical_mosaic_after_original_116_below_120_STOP":True,
        "original_116_site_single_scene_primary_source_gate_remains_FAILED":True}
    qa["gate_pass"]=len(eligible)>=MIN_SITES and h7>=MIN_H7 and h6>=MIN_H6
    return eligible,qa

def datewise_original_TP(archive,wq):
    name=f"TP.2_-_Harmonised/{wq}.csv"
    if name not in archive.namelist():raise RuntimeError("Original WQ site unavailable")
    r=csv.DictReader(io.StringIO(archive.read(name).decode("utf-8-sig")))
    q=defaultdict(list)
    for rec in r:
        try:
            d=datetime.date.fromisoformat(str(rec["date"]).strip()[:10])
            x=float(rec["value"])
        except (ValueError,TypeError,KeyError):continue
        if d.year in (2018,2022) and d.month in (6,7,8) and np.isfinite(x) and x>0:
            q[(d.year,d.month)].append(x)
    return q

def observations(eligible,a):
    E=load_csv(a.aquatic,ORIGINAL_263_NUMERIC_SHA)
    H=load_csv(a.h7,ORIGINAL_73_COVAR_SHA)
    E={r["EEA_identical_native_site_key_2018_2022"]:r for r in E}
    H={r["HYBAS_ID_L07"]:r for r in H}
    if len(E)!=263 or len(H)!=73:raise RuntimeError("Original 263 site or 73 H7 numeric source changed")
    tp_raw=a.original_TP.read_bytes()
    if len(tp_raw)!=TP_ZIP_SIZE or sha(tp_raw)!=TP_ZIP_SHA or hashlib.md5(tp_raw).hexdigest()!=TP_ZIP_MD5:
        raise RuntimeError("Original published TP.2 2018/2022 ZIP source not authenticated")
    archive=zipfile.ZipFile(io.BytesIO(tp_raw))
    rows=[]
    for x in eligible:
        e=E.get(x["site_key"]);h=H.get(x["H7"])
        if e is None or h is None or x["H6"]!=e["HYBAS_ID_L06"] or x["H7"]!=e["HYBAS_ID_L07"]:
            raise RuntimeError("Same-site original numeric source ID/H6/H7 mismatches fixed pixel site")
        bymonth=datewise_original_TP(archive,e["Peters_identical_WQ_ID_2018_2022"])
        pairs=[]
        for m in x["posthoc_composite_delta_NDVI_by_month"]:
            month=int(m["same_calendar_month"])
            r18=bymonth[(2018,month)];r22=bymonth[(2022,month)]
            if not r18 or not r22:
                raise RuntimeError("A pixel-qualified month had no truly measured TP in both years")
            d_tp=math.log(statistics.median(r22))-math.log(statistics.median(r18))
            pairs.append((d_tp,float(m["delta_near_minus_far_2022_minus_2018"]),
              float(m["ndvi_near_minus_far_2018"]),month))
        if len(pairs)<2:raise RuntimeError("Preselected original pixel site lost 2 paired months")
        rows.append({"site_key":x["site_key"],"H7":x["H7"],"H6":x["H6"],
          "source_original_same_calendar_months":[p[3] for p in pairs],
          "delta_ln_original_measured_river_TP":float(np.mean([p[0] for p in pairs])),
          "delta_true_riverline_near100_minus_far1000_NDVI":float(np.mean([p[1] for p in pairs])),
          "original_2018_true_riverline_near_far_NDVI":float(np.mean([p[2] for p in pairs])),
          "actual_delta_I2M2_benthic_EQR":float(e["real_original_delta_I2M2_EQR_2022_minus_2018"]),
          "h7_upstream_cropland":float(h["source_original_upstream_cropland_pct"]),
          "h7_upstream_urban":float(h["source_original_upstream_urban_pct"]),
          "h7_log_area":float(h["log_original_H7_local_area_km2"]),
          "h7_log_flow":float(h["mean_original_site_modelled_RIVER_Q_log1p"]),
          "source_original_H6_preassigned_fold":int(h["original_H6_source_fold"])})
    return rows

def aggregate_to_h7(rows):
    by=defaultdict(list)
    for x in rows:by[x["H7"]].append(x)
    result=[]
    for group,v in sorted(by.items()):
        ref=v[0]
        if len({x["H6"] for x in v})!=1 or len({x["source_original_H6_preassigned_fold"] for x in v})!=1:
            raise RuntimeError("H7 observations incorrectly span H6 parent folds")
        out={"H7":group,"H6":ref["H6"],"source_original_H6_preassigned_fold":ref["source_original_H6_preassigned_fold"],
             "original_distinct_physical_sites_per_H7":len(v)}
        for k in ("delta_ln_original_measured_river_TP",
                  "delta_true_riverline_near100_minus_far1000_NDVI",
                  "original_2018_true_riverline_near_far_NDVI",
                  "actual_delta_I2M2_benthic_EQR",
                  "h7_upstream_cropland","h7_upstream_urban","h7_log_area","h7_log_flow"):
            out[k]=float(np.mean([x[k] for x in v]))
        result.append(out)
    return result

def bootstrap_spearman(z):
    x=np.array([v["delta_ln_original_measured_river_TP"] for v in z])
    y=np.array([v["delta_true_riverline_near100_minus_far1000_NDVI"] for v in z])
    parents=np.array([v["H6"] for v in z])
    keys=np.unique(parents)
    stat=float(spearmanr(x,y).statistic)
    rng=np.random.default_rng(SEED)
    draws=[]
    for _ in range(NUM_BOOT):
        picks=rng.choice(keys,size=len(keys),replace=True)
        ii=np.concatenate([np.flatnonzero(parents==p) for p in picks])
        if len(np.unique(x[ii]))<2 or len(np.unique(y[ii]))<2:continue
        val=float(spearmanr(x[ii],y[ii]).statistic)
        if np.isfinite(val):draws.append(val)
    return {"n_original_H7":len(z),"n_original_H6":len(keys),
            "rho":stat,"descriptive_parent_bootstrap95":list(np.quantile(draws,[.025,.975])) if draws else None,
            "this_is_exploratory_association_not_a_confirmatory_p_or_causal_estimate":True}

def H6_geocv(z):
    X=np.array([[r[k] for k in ("h7_upstream_cropland","h7_upstream_urban",
                "h7_log_area","h7_log_flow","original_2018_true_riverline_near_far_NDVI")]
                for r in z],float)
    y=np.array([r["delta_true_riverline_near100_minus_far1000_NDVI"] for r in z],float)
    tp=np.array([r["delta_ln_original_measured_river_TP"] for r in z],float)
    folds=np.array([r["source_original_H6_preassigned_fold"] for r in z],int)
    parents=np.array([r["H6"] for r in z])
    if len(set(folds))!=5 or any(len(set(folds[parents==p]))!=1 for p in set(parents)):
        raise RuntimeError("Original H6-group folds not preserved")
    pred=np.full((len(z),2),np.nan)
    for fold in sorted(set(folds)):
        train=np.where(folds!=fold)[0];test=np.where(folds==fold)[0]
        if len(test)<3 or len(train)<20:raise RuntimeError("Source fixed H6 out-of-region fold too sparse")
        for col,features in enumerate((X,np.column_stack((X,tp)))):
            scaler=StandardScaler()
            fit=scaler.fit_transform(features[train])
            model=Ridge(alpha=10)
            model.fit(fit,y[train])
            pred[test,col]=model.predict(scaler.transform(features[test]))
    if not np.isfinite(pred).all():raise RuntimeError("Original geographical holdout predicted missing result")
    errors=(y[:,None]-pred)**2
    mse=errors.mean(axis=0)
    k=sorted(set(parents))
    agg=np.array([[errors[parents==p,c].mean() for c in (0,1)] for p in k])
    gain=100*(1-agg[:,1].mean()/agg[:,0].mean())
    rng=np.random.default_rng(SEED)
    r=rng.integers(0,len(k),size=(NUM_BOOT,len(k)))
    gg=100*(1-agg[r,1].mean(axis=1)/agg[r,0].mean(axis=1))
    return {"original_true_H6_holdout_folds":5,"original_Ridge_alpha":10,
      "source_site_H7_equal_original_model_n":len(z),
      "site_MSE_base_and_addTP":mse.tolist(),
      "site_weighted_MSE_gain_pct":float(100*(1-mse[1]/mse[0])),
      "H6_equal_MSE_gain_pct":float(gain),
      "fixed_parent_error_descriptive95_pct":np.quantile(gg,[.025,.975]).tolist(),
      "not_a_wastewater_or_actual_nutrient_transmission_effect":True}

def run(a):
    a.out.mkdir(parents=True,exist_ok=True)
    good,qa=pixel_source(a.pixel_dir)
    if not qa["gate_pass"]:
        output={"schema":SCHEMA,"status":"STOP_POSTHOC_MOSAIC_SOURCE_UNDERSUPPORTED_NO_NUMERIC_ECOLOGY",
                "source_gate":qa,"no_tp_or_EQR_magnitudes_read":True}
    else:
        rows=observations(good,a)
        h7=aggregate_to_h7(rows)
        outcome=bootstrap_spearman(h7)
        prediction=H6_geocv(h7)
        output={"schema":SCHEMA,"status":"POSTHOC_EXPLORATORY_MOSAIC_RIVER_TP_NDVI_FITTED_NOT_ORIGINAL_MAIN_MODEL",
          "source_gate":qa,"original_river_source_chemistry_method":"physical same WQ ID and matched summer real TP medians",
          "original_h7_spearman":outcome,"original_h6_grouped_prediction":prediction,
          "n_original_biological_monitor_sites":len(rows),
          "original_single_scene_source_gate_failed_116_below_120":True,
          "this_composite_is_after_known_original_source_failure":True,
          "original_h7_site_equal_values":h7,"original_native_site_month_matched_values":rows,
          "no_actual_municipal_source_specific_water_to_land_NP_transport":True,
          "no_confirmatory_or_independent_reproduction_claim":True}
    (a.out/"FR178_POSTHOC_TWO_SCENE_MOSAIC_REAL_TP_NDVI_RESULTS.json").write_text(
        json.dumps(output,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("FR178_POSTHOC_2SCENE_MOSAIC_NUTRIENT_TERRESTRIAL="+json.dumps({
      "status":output["status"],"source":qa,
      "if_go_rho":output.get("original_h7_spearman",{}).get("rho")},ensure_ascii=False),flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--pixel-dir",type=Path,required=True)
    p.add_argument("--original-TP",type=Path)
    p.add_argument("--aquatic",type=Path)
    p.add_argument("--h7",type=Path)
    p.add_argument("--out",type=Path,required=True)
    run(p.parse_args())
if __name__=="__main__":main()
