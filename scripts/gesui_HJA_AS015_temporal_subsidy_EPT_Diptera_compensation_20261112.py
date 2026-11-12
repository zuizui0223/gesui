#!/usr/bin/env python3
"""Posthoc 17-stream EPT vs chironomid seasonal covariance/pulse decomposition.

Original HJA Andrews LTER AS015 v2 original SHA and PRE-numeric 17x23x4
source site+period manifest required. Nothing inferred about WWTP or consumers.
"""
import collections,csv,hashlib,io,json,math,os,pathlib,statistics,time,urllib.request
import numpy as np
URL="https://andrewsforest.oregonstate.edu/fsdbdata/cache/AS01501_v1.csv"
SOURCE_SHA256="c8ae00e91e7b210fbc86f5b08c6047ceaeade30bb2d410c3748f495fea0de4b2"
SOURCE_MD5="f3e6e400990da883ae0d341731431749"
MANIFEST_SHA256="51354ce5cc3ab3acc63aecc28dd93534438ab3ffaf09b86ab2b0994b51874923"
src=pathlib.Path(os.environ["AS015_FROZEN_SOURCE_MANIFEST"]).read_bytes()
if hashlib.sha256(src).hexdigest()!=MANIFEST_SHA256:
    raise RuntimeError("Pre-outcome original AS015 stream/period/taxa manifest SHA changed")
source=json.loads(src)
if source["original_url"]!=URL or source["original_SHA256"]!=SOURCE_SHA256:
    raise RuntimeError("Original HJA source metadata does not match frozen protocol")
err=None
for n in range(4):
    try:
        req=urllib.request.Request(URL,headers={"User-Agent":"Mozilla/5.0 (compatible; gesui-AS015-seasonal-20261112/1)","Accept":"text/csv,*/*"})
        with urllib.request.urlopen(req,timeout=70) as resp:raw=resp.read()
        if hashlib.sha256(raw).hexdigest()!=SOURCE_SHA256 or hashlib.md5(raw).hexdigest()!=SOURCE_MD5:
            raise RuntimeError("Publisher HJA original source bytes changed; stop")
        break
    except Exception as e:
        err=e;time.sleep(1+2*n)
else:raise RuntimeError("Source unavailable "+str(err))
rows=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"),newline="")))
if len(rows)!=1564:raise RuntimeError("Source rows not matching original v2")
def key(r):return (r["SITE"].strip(),r["TIME_PERIOD"].strip(),r["SAMPLE_DATE"].strip(),r["TAXA_ORDER"].strip())
if set(key(r) for r in rows)!={tuple(k) for k in source["source_complete_keys"]}:
    raise RuntimeError("AS015 original complete event identities changed")
by=collections.defaultdict(dict)
regime={}
for r in rows:
    site,period,date,order=key(r)
    if site in regime and regime[site]!=r["STREAM_REGIME"]:raise RuntimeError("Stream regime changed")
    regime[site]=r["STREAM_REGIME"]
    a=float(r["EMERG_ABUNDANCE"]);b=float(r["EMERG_BIOMASS"])
    if not all(math.isfinite(x) and x>=0 for x in [a,b]):
        raise RuntimeError("Invalid emerging insect reported original value")
    z=by[(site,period,date)]
    if order in z:raise RuntimeError("Original order duplicates")
    z[order]={"n":a,"g":b}
assert len(by)==391 and len(regime)==17
orders={"Ephemeroptera","Plecoptera","Trichoptera","Diptera"}
streams=collections.defaultdict(list)
for (site,period,date),g in by.items():
    if set(g)!=orders:raise RuntimeError("Incomplete 4-group native event")
    E=sum(g[x]["g"] for x in ("Ephemeroptera","Plecoptera","Trichoptera"))
    D=g["Diptera"]["g"]
    EA=sum(g[x]["n"] for x in ("Ephemeroptera","Plecoptera","Trichoptera"))
    DA=g["Diptera"]["n"]
    streams[site].append({"period":period,"date":date,"E":E,"D":D,"T":E+D,
                         "EA":EA,"DA":DA,"TA":EA+DA})
outrows=[]
for site,rawsite in sorted(streams.items()):
    es=sorted(rawsite,key=lambda x:(x["date"],x["period"]))
    if len(es)!=23 or len({x["period"] for x in es})!=23:
        raise RuntimeError("Original 23 nonrepeating site periods not present")
    E=np.asarray([x["E"] for x in es])
    D=np.asarray([x["D"] for x in es])
    T=E+D
    vari=np.var(E,ddof=1);vard=np.var(D,ddof=1)
    vart=np.var(T,ddof=1)
    cov=float(np.cov(E,D,ddof=1)[0,1])
    denom=vari+vard
    comp=float(-2*cov/denom) if denom>0 else None
    stde=math.sqrt(vari);stdd=math.sqrt(vard)
    sync=float(vart/(stde+stdd)**2) if (stde+stdd)>0 else None
    total=float(np.sum(T))
    totalabu=float(sum(x["TA"] for x in es))
    epttot=float(np.sum(E))
    eptabu=float(sum(x["EA"] for x in es))
    if total<=0:raise RuntimeError("Site has all-zero biomass, undefined pulse fraction")
    cv=float(math.sqrt(vart)/np.mean(T))
    if not np.isfinite(cv):raise RuntimeError("CV nonfinite")
    shiftvars=[float(np.var(E+np.roll(D,j),ddof=1)) for j in range(23)]
    if abs(shiftvars[0]-vart)>1e-12:raise RuntimeError("Circular-shift 0 changed total biomass")
    streamrow={
        "site":site,"regime":regime[site],"n_periods":23,
        "zero_total_emergence_biomass_periods":int(np.sum(T==0)),
        "mean_total_emergent_biomass_g_m2_per_interval":float(np.mean(T)),
        "mean_EPT_emergent_biomass_g_m2_per_interval":float(np.mean(E)),
        "mean_Diptera_emergent_biomass_g_m2_per_interval":float(np.mean(D)),
        "EPT_adult_biomass_fraction_pooled":float(epttot/total),
        "EPT_adult_abundance_fraction_pooled":float(eptabu/totalabu) if totalabu>0 else None,
        "variance_total_biomass_per_interval":float(vart),
        "variance_EPT_biomass":float(vari),
        "variance_Diptera_biomass":float(vard),
        "within_site_EPT_Diptera_biomass_covariance":cov,
        "temporal_compensation_index_C_positive_means_opposing_peaks":comp,
        "temporal_synchrony_phi":sync,
        "total_biomass_CV_23_intervals":cv,
        "largest_one_interval_share_total_biomass":float(np.max(T)/total),
        "top5_intervals_share_total_biomass":float(np.sum(np.sort(T)[-5:])/total),
        "actual_vartot_percentile_among_23_circular_shifts":float(np.mean(np.array(shiftvars)<=shiftvars[0]+1e-14)),
        "actual_vartot_less_than_median_circular_shift":bool(shiftvars[0]<np.median(shiftvars)),
        "sampled_daily_flux_NOT_measured":True
    }
    outrows.append(streamrow)
if any(v["temporal_compensation_index_C_positive_means_opposing_peaks"] is None for v in outrows):
    raise RuntimeError("Some streams lack variation; no posthoc exclusion allowed")
def median(key,rows=outrows):
    return float(statistics.median(z[key] for z in rows))
def mean(key,rows=outrows):
    return float(statistics.mean(z[key] for z in rows))
def rangeinfo(key):
    a=[z[key] for z in outrows]
    return {"min":float(min(a)),"median":float(statistics.median(a)),"max":float(max(a))}
N=len(outrows)
C=np.asarray([z["temporal_compensation_index_C_positive_means_opposing_peaks"] for z in outrows])
A=np.asarray([z["top5_intervals_share_total_biomass"] for z in outrows])
R=np.asarray([z["total_biomass_CV_23_intervals"] for z in outrows])
rng=np.random.default_rng(20261112)
samples=rng.integers(0,N,size=(10000,N))
bootstrap_C=np.median(C[samples],axis=1)
bootstrap_positive=np.mean(C[samples]>0,axis=1)
def interval(t):return [float(x) for x in np.quantile(t,[.025,.975])]
regime_out={}
for k in ["I","P"]:
    sub=[z for z in outrows if z["regime"]==k]
    regime_out[k]={
      "stream_count":len(sub),
      "median_temporal_compensation_C":median("temporal_compensation_index_C_positive_means_opposing_peaks",sub),
      "mean_top5_pulse_fraction":mean("top5_intervals_share_total_biomass",sub),
      "mean_total_biomass_CV":mean("total_biomass_CV_23_intervals",sub),
      "n_positive_compensation":sum(z["temporal_compensation_index_C_positive_means_opposing_peaks"]>0 for z in sub),
      "n_actual_variance_lower_than_circular_shift_median":sum(z["actual_vartot_less_than_median_circular_shift"] for z in sub)
    }
out={
 "schema":"AS015_2022_adult_emergence_biomass_temporal_subsidy_continuity_posthoc_v1",
 "original_HJA_AS015_MD5":SOURCE_MD5,"source_sha256":SOURCE_SHA256,
 "site_period_order_manifest_SHA256":MANIFEST_SHA256,
 "stream_independent_n":N,"all_in_one_watershed":True,
 "n_stream_periods":len(by),"source_taxon_rows":len(rows),
 "site_stats":outrows,
 "group_stability":{
   "n_stream_C_positive":int(np.sum(C>0)),
   "n_stream_C_negative":int(np.sum(C<0)),
   "median_C":float(np.median(C)),
   "mean_C":float(np.mean(C)),
   "descriptive_bootstrap95_median_C":interval(bootstrap_C),
   "descriptive_bootstrap95_fraction_C_positive":interval(bootstrap_positive),
   "n_stream_actual_variance_less_than_circular_shift_median":sum(z["actual_vartot_less_than_median_circular_shift"] for z in outrows),
   "median_temporal_synchrony_phi":median("temporal_synchrony_phi"),
   "total_emergence_biomass_CV":rangeinfo("total_biomass_CV_23_intervals"),
   "top5_of23_intervals_total_biomass_fraction":rangeinfo("top5_intervals_share_total_biomass"),
   "largest_one_interval_biomass_fraction":rangeinfo("largest_one_interval_share_total_biomass"),
   "pooled_EPT_biomass_fraction_across_intervals":rangeinfo("EPT_adult_biomass_fraction_pooled"),
   "pooled_EPT_abundance_fraction_across_intervals":rangeinfo("EPT_adult_abundance_fraction_pooled"),
   "regimes":regime_out},
 "scientific_limits":[
   "This is a posthoc within-watershed descriptive study of 17 streams, not 391 independent sites.",
   "Only EPT and Diptera/Chironomidae adults measured, not complete terrestrial prey or pollinator assemblages.",
   "g/m² per collection interval, not directly per-day flux or predator consumption; periods typically 3–4 days.",
   "Temporal covariance of two component biomasses is mathematical and depends on sampling and phenology. It does not show beneficial ecosystem services.",
   "Circular shifts preserve each group's marginal series and autocorrelation but are not a randomization of habitat or season; no causal p is estimated.",
   "No river water nutrients, wastewater source, spider/bird diet or floral reproduction measured in AS015; no cross-paper row joins."
 ]}
outpath=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_AS015_seasonal_continuity"
outpath.mkdir(parents=True,exist_ok=True)
(outpath/"AS015_EPT_DIPTERA_EMERGENCE_TEMPORAL_COMPENSATION_SOURCE_ORIGINAL_20261112.json").write_text(json.dumps(out,indent=2,sort_keys=True))
print("AS015_TEMPORAL_SUBSIDY_RESULT_START")
print(json.dumps({"n":N,"site_periods":len(by),"summary":out["group_stability"]},sort_keys=True))
print("AS015_TEMPORAL_SUBSIDY_RESULT_END")
