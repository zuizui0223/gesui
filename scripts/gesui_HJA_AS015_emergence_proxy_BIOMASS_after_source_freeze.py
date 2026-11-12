#!/usr/bin/env python3
"""Independent empirical adult-emergence calibration of EPT proxy vs export amount.

Original HJA LTER AS015 v2 actual CSV, MD5/SHA verified, PRE-numeric stream /
period / taxon manifest verified from distinct earlier GH Actions job.

This DOES NOT join NRSA (different US sites / years), estimate sewage impact,
sample a new independent basin, or measure pollination service.
"""
import collections,csv,hashlib,io,itertools,json,math,os,pathlib,statistics,time,urllib.request
import numpy as np
from scipy.stats import rankdata

ORIGINAL_MD5="f3e6e400990da883ae0d341731431749"
ORIGINAL_SHA="c8ae00e91e7b210fbc86f5b08c6047ceaeade30bb2d410c3748f495fea0de4b2"
SOURCE_MANIFEST_SHA="51354ce5cc3ab3acc63aecc28dd93534438ab3ffaf09b86ab2b0994b51874923"
mfile=pathlib.Path(os.environ["AS015_FROZEN_SOURCE_MANIFEST"])
mraw=mfile.read_bytes()
if hashlib.sha256(mraw).hexdigest()!=SOURCE_MANIFEST_SHA:
    raise RuntimeError("Source-only AS015 cohort changed after pre-numeric gate")
manifest=json.loads(mraw)
url=manifest["original_url"]
if url!="https://andrewsforest.oregonstate.edu/fsdbdata/cache/AS01501_v1.csv":
    raise RuntimeError("Attempted to use source other than original official EML-linked table")
if manifest["original_MD5"]!=ORIGINAL_MD5 or manifest["original_SHA256"]!=ORIGINAL_SHA:
    raise RuntimeError("Original publication source hashes changed")
def download():
    err=None
    for i in range(4):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"Mozilla/5.0 (compatible; gesui-lter-source-pinned-analysis/1)",
                    "Accept":"text/csv,*/*"})
            with urllib.request.urlopen(req,timeout=70) as resp:raw=resp.read()
            if hashlib.md5(raw).hexdigest()!=ORIGINAL_MD5 or hashlib.sha256(raw).hexdigest()!=ORIGINAL_SHA:
                raise RuntimeError("Original CC-BY AS015 body SHA/MD5 changed — stop")
            return raw
        except Exception as e:
            err=e;time.sleep(1+2*i)
    raise RuntimeError("Unable to retrieve original AS015 verified bytes "+str(err))
raw=download()
rdr=csv.DictReader(io.StringIO(raw.decode("utf-8-sig"),newline=""))
rows=list(rdr)
if len(rows)!=1564:raise RuntimeError("AS015 v2 1564 rows count changed")
key=lambda r:(r["SITE"].strip(),r["TIME_PERIOD"].strip(),r["SAMPLE_DATE"].strip(),r["TAXA_ORDER"].strip())
keys=set(key(r) for r in rows)
if len(keys)!=len(rows) or keys!=set(tuple(a) for a in manifest["source_complete_keys"]):
    raise RuntimeError("Original row identity differs from SHA-locked source-only manifest")
orders=("Ephemeroptera","Plecoptera","Trichoptera","Diptera")
groups=collections.defaultdict(dict)
regimes={}
for row in rows:
    site,period,date,order=key(row)
    if site in regimes and regimes[site]!=row["STREAM_REGIME"]:
        raise RuntimeError("Contradictory stream hydrological regime")
    regimes[site]=row["STREAM_REGIME"]
    abu=float(row["EMERG_ABUNDANCE"])
    mass=float(row["EMERG_BIOMASS"])
    if not math.isfinite(abu) or not math.isfinite(mass) or abu<0 or mass<0:
        raise RuntimeError("Invalid original quantified adult emergence amount")
    groups[(site,period,date)][order]=(abu,mass)
if set(regimes.values())!={"I","P"} or len(regimes)!=17 or len(groups)!=391:
    raise RuntimeError("Frozen original study field units changed")
for k,g in groups.items():
    if set(g)!=set(orders):raise RuntimeError("Original AS015 complete four-order event lost")
events=[]
for (site,period,date),g in sorted(groups.items()):
    abu=sum(v[0] for v in g.values())
    biomass=sum(v[1] for v in g.values())
    E=sum(g[o][0] for o in orders[:3])
    M=sum(g[o][1] for o in orders[:3])
    diptera_ab=g["Diptera"][0];diptera_mass=g["Diptera"][1]
    events.append({
      "site":site,"period":period,"date":date,"regime":regimes[site],
      "total_adult_abundance_per_m2_interval":abu,
      "total_adult_biomass_g_m2_interval":biomass,
      "EPT_adult_abundance_per_m2_interval":E,
      "EPT_adult_biomass_g_m2_interval":M,
      "Diptera_Chironomidae_abundance_per_m2_interval":diptera_ab,
      "Diptera_Chironomidae_biomass_g_m2_interval":diptera_mass,
      "EPT_adult_abundance_fraction":E/abu if abu>0 else None,
      "EPT_adult_biomass_fraction":M/biomass if biomass>0 else None
    })
sby=collections.defaultdict(list)
for e in events:sby[e["site"]].append(e)
def mean(x):return sum(x)/len(x)
def med(x):return statistics.median(x)
stats=[]
for site,ev in sorted(sby.items()):
    if len(ev)!=23 or len(set(r["period"] for r in ev))!=23:
        raise RuntimeError("Site has different source sampling opportunities; stop")
    abu=[r["total_adult_abundance_per_m2_interval"] for r in ev]
    bio=[r["total_adult_biomass_g_m2_interval"] for r in ev]
    shares=[r["EPT_adult_abundance_fraction"] for r in ev if r["EPT_adult_abundance_fraction"] is not None]
    bshares=[r["EPT_adult_biomass_fraction"] for r in ev if r["EPT_adult_biomass_fraction"] is not None]
    E=[r["EPT_adult_abundance_per_m2_interval"] for r in ev]
    M=[r["EPT_adult_biomass_g_m2_interval"] for r in ev]
    D=[r["Diptera_Chironomidae_abundance_per_m2_interval"] for r in ev]
    Db=[r["Diptera_Chironomidae_biomass_g_m2_interval"] for r in ev]
    stats.append({
     "site":site,"regime":regimes[site],"n_interval":len(ev),
     "intervals_with_positive_adult_abundance":len(shares),
     "intervals_with_positive_adult_biomass":len(bshares),
     "mean_interval_total_adult_abundance_per_m2":mean(abu),
     "mean_interval_total_adult_biomass_g_per_m2":mean(bio),
     "mean_interval_EPT_adult_abundance_fraction":mean(shares) if shares else None,
     "mean_interval_EPT_adult_biomass_fraction":mean(bshares) if bshares else None,
     "pooled_abundance_weighted_EPT_fraction":sum(E)/sum(abu) if sum(abu)>0 else None,
     "pooled_biomass_weighted_EPT_fraction":sum(M)/sum(bio) if sum(bio)>0 else None,
     "mean_interval_EPT_adult_abundance_per_m2":mean(E),
     "mean_interval_EPT_adult_biomass_g_per_m2":mean(M),
     "mean_interval_Diptera_adult_abundance_per_m2":mean(D),
     "mean_interval_Diptera_adult_biomass_g_per_m2":mean(Db),
     "biomass_g_per_recorded_individual":sum(bio)/sum(abu) if sum(abu)>0 else None
    })
if any(z["mean_interval_EPT_adult_abundance_fraction"] is None for z in stats):
    raise RuntimeError("Some stream has no observed positive adult emergence event")
site_array=lambda k:np.array([z[k] for z in stats],dtype=float)
def spearman(a,b):
    if len(a)<3 or np.std(a)==0 or np.std(b)==0:return float("nan")
    return float(np.corrcoef(rankdata(a,method="average"),rankdata(b,method="average"))[0,1])
x=site_array("mean_interval_EPT_adult_abundance_fraction")
biomass=site_array("mean_interval_total_adult_biomass_g_per_m2")
abundance=site_array("mean_interval_total_adult_abundance_per_m2")
regime=np.array([z["regime"] for z in stats])
rng=np.random.default_rng(20261112)
def perm_spearman(a,b,n=20000):
    obs=spearman(a,b)
    vals=np.array([spearman(a,rng.permutation(b)) for i in range(n)])
    return {"rho":obs,"MC_draws":n,
       "two_sided_MC_reference_fraction":float((1+np.count_nonzero(np.abs(vals)>=abs(obs)-1e-14))/(n+1))}
def boot_spearman(a,b,n=10000):
    indices=rng.integers(0,len(a),size=(n,len(a)))
    draws=[spearman(a[ii],b[ii]) for ii in indices]
    valid=np.asarray([v for v in draws if np.isfinite(v)],dtype=float)
    return {"n_draws_requested":n,"n_valid":len(valid),
      "descriptive_percentile95":[float(v) for v in np.quantile(valid,[0.025,0.975])] if len(valid) else None}
primary=perm_spearman(x,biomass)
primary["bootstrap_stream_resample"]=boot_spearman(x,biomass)
secondary=perm_spearman(x,abundance)
secondary["bootstrap_stream_resample"]=boot_spearman(x,abundance)
# Regime comparisons: exactly 17 site means, 9 intermittent and 8 perennial,
# exhaustive 17 choose 9 = 24310 label arrangements.
I=np.where(regime=="I")[0];P=np.where(regime=="P")[0]
if len(I)!=9 or len(P)!=8:raise RuntimeError("Original stream regimes do not match source metadata")
def exact_regime_difference(v):
    obs=float(np.mean(v[I])-np.mean(v[P]))
    allidx=np.arange(17)
    diffs=[]
    for combo in itertools.combinations(range(17),9):
        ids=np.array(combo);ids2=np.setdiff1d(allidx,ids,assume_unique=True)
        diffs.append(float(np.mean(v[ids])-np.mean(v[ids2])))
    dd=np.asarray(diffs)
    assert len(dd)==24310
    # Not a randomized-experiment p; merely conditional stream-label exchangeability.
    return {"mean_I":float(np.mean(v[I])),"mean_P":float(np.mean(v[P])),
            "mean_I_minus_P":obs,
            "exact_combinations":len(dd),
            "two_sided_exact_label_permutation_reference_fraction":float(np.mean(np.abs(dd)>=abs(obs)-1e-14))}
regime_contrasts={
 "EPT_abundance_fraction":exact_regime_difference(x),
 "total_emergence_biomass_g_m2_per_interval":exact_regime_difference(biomass),
 "total_emergence_abundance_individuals_m2_per_interval":exact_regime_difference(abundance)
}
# The ecological interpretation of a high/low median classification is
# secondary only; prevalence and direction cannot be called hypothesis tests.
x_med=med(x);bm_med=med(biomass)
discordant=[z for z in stats if
    ((z["mean_interval_EPT_adult_abundance_fraction"]>x_med) !=
     (z["mean_interval_total_adult_biomass_g_per_m2"]>bm_med))]
# No event pseudoreplication; within-site, report site-specific correlations
# of event-level adult share and biomass, not one spurious n=391 p-value.
within=[]
for site,es in sorted(sby.items()):
    positive=[a for a in es if a["EPT_adult_abundance_fraction"] is not None]
    rho=spearman(np.array([e["EPT_adult_abundance_fraction"] for e in positive]),
                 np.array([e["total_adult_biomass_g_m2_interval"] for e in positive])) if len(positive)>=5 else None
    within.append({"site":site,"n_positive_intervals":len(positive),"rho":rho})
positive_within=[z["rho"] for z in within if z["rho"] is not None and np.isfinite(z["rho"])]
result={
 "study":"HJA_AS015_2022_ORIGINAL_ADULT_EMERGENCE_PROPORTION_VS_BIOMASS_CALIBRATION",
 "source_original_HJA_CC_BY_4":True,
 "citation":"Latorre-B T & Arismendi I (2026), HJ Andrews Experimental Forest AS015 v2, DOI:10.6073/pasta/6e9b83e2547e3a82d92dab068476f1d7",
 "original_2022_CSV_SHA256":ORIGINAL_SHA,
 "original_2022_CSV_MD5":ORIGINAL_MD5,
 "original_source_manifest_SHA256":SOURCE_MANIFEST_SHA,
 "study_period":"2022-06-25 to 2022-09-01",
 "number_streams":len(stats),"streams_regime":{"I":len(I),"P":len(P)},
 "site_x_period_events":len(events),"source_taxon_rows":len(rows),
 "zero_total_abundance_site_periods":sum(e["total_adult_abundance_per_m2_interval"]==0 for e in events),
 "zero_total_biomass_site_periods":sum(e["total_adult_biomass_g_m2_interval"]==0 for e in events),
 "primary_stream_level_EPT_adult_share_vs_total_emergence_biomass":primary,
 "secondary_stream_level_EPT_share_vs_total_emergence_abundance":secondary,
 "secondary_stream_level_EPT_share_vs_mean_EPT_biomass":spearman(x,site_array("mean_interval_EPT_adult_biomass_g_per_m2")),
 "secondary_stream_level_EPT_share_vs_mean_chironomid_biomass":spearman(x,site_array("mean_interval_Diptera_adult_biomass_g_per_m2")),
 "secondary_stream_level_EPT_adult_share_vs_EPT_biomass_share":spearman(x,site_array("mean_interval_EPT_adult_biomass_fraction")),
 "regime_site_mean_differences":regime_contrasts,
 "median_EPT_share":x_med,"median_total_biomass_g_m2_per_interval":bm_med,
 "discordance_cross_median_sites":len(discordant),
 "within_stream_temporal_correlations":{
   "n_streams_with_positive_intervals_rho":len(positive_within),
   "median_rho":med(positive_within) if positive_within else None,
   "streams_with_positive_rho":sum(v>0 for v in positive_within),
   "streams_with_negative_rho":sum(v<0 for v in positive_within),
   "each_stream":within},
 "site_level_17_independent_stream_records":stats,
 "restrictions":[
   "This tests emergence-proxy suitability, NOT whether NRSA benthic EPT fraction statistically causes emergence at the same sites.",
   "The 17 emergence sites share one Oregon watershed and calendar summer; not globally independent catchments.",
   "AS015 collected only 4 insect orders including Chironomidae within Diptera, not all freshwater or terrestrial insects.",
   "Biomass was measured or estimated by size for four orders and is g dry mass m-2 per collection interval, not an independent continuous daily flux.",
   "Zero total emergent insect observations are not assigned an EPT fraction of zero.",
   "No WWTP N/P, true source-to-river nutrient dose, riparian plants, pollinators, predators, or realized ecosystem services measured in this 17-stream source.",
   "Monte Carlo and exhaustive permutation fractions rely on questionable between-stream exchangeability; report as descriptive only."
  ]
}
root=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_HJA_AS015_emergence_actual"
root.mkdir(parents=True,exist_ok=True)
(root/"AS015_EMERGENCE_EPT_SHARE_TO_BIOMASS_CALIBRATION_20261112.json").write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False))
brief={k:result[k] for k in [
 "number_streams","streams_regime","site_x_period_events","zero_total_abundance_site_periods",
 "primary_stream_level_EPT_adult_share_vs_total_emergence_biomass",
 "secondary_stream_level_EPT_share_vs_total_emergence_abundance",
 "secondary_stream_level_EPT_share_vs_mean_EPT_biomass",
 "secondary_stream_level_EPT_share_vs_mean_chironomid_biomass",
 "regime_site_mean_differences","discordance_cross_median_sites",
 "within_stream_temporal_correlations"]}
brief["within_stream_temporal_correlations"].pop("each_stream",None)
print("AS015_EMERGENCE_CALIBRATION_RESULT_START",flush=True)
print(json.dumps(brief,sort_keys=True,allow_nan=False),flush=True)
print("AS015_EMERGENCE_CALIBRATION_RESULT_END",flush=True)
