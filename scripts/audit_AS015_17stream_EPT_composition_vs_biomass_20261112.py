#!/usr/bin/env python3
"""Independent 17-stream LTER AS015 sensitivity using previously verified aggregate.

No requery, outcome-aware model selection, or merging NRSA with Oregon.
Manual rank calculations (no scipy) challenge the original mean-interval
EPT-share interpretation with abundance-weighted alternative.
"""
import hashlib,json,math,os,pathlib,statistics
from collections import Counter

path=pathlib.Path(os.environ["AS015_VERIFIED_NUMERIC_RESULT"])
data=json.loads(path.read_text())
if data["original_2022_CSV_SHA256"]!="c8ae00e91e7b210fbc86f5b08c6047ceaeade30bb2d410c3748f495fea0de4b2" or \
   data["original_source_manifest_SHA256"]!="51354ce5cc3ab3acc63aecc28dd93534438ab3ffaf09b86ab2b0994b51874923":
    raise RuntimeError("Original source results were changed")
sites=data["site_level_17_independent_stream_records"]
assert len(sites)==17 and len({x["site"] for x in sites})==17
def rank(arr):
    indexed=sorted(enumerate(arr),key=lambda x:x[1])
    ranks=[0.]*len(arr)
    i=0
    while i<len(arr):
        j=i+1
        while j<len(arr) and indexed[j][1]==indexed[i][1]:j+=1
        avg=(i+1+j)/2
        for k in range(i,j):ranks[indexed[k][0]]=avg
        i=j
    return ranks
def corr(x,y):
    if len(x)<3:return None
    a=rank(x);b=rank(y)
    ma=sum(a)/len(a);mb=sum(b)/len(b)
    da=[v-ma for v in a];db=[v-mb for v in b]
    denom=math.sqrt(sum(v*v for v in da)*sum(v*v for v in db))
    return sum(u*v for u,v in zip(da,db))/denom if denom else None
def get(arr,key):return [z[key] for z in arr]
K="mean_interval_EPT_adult_abundance_fraction"
B="mean_interval_total_adult_biomass_g_per_m2"
A="mean_interval_total_adult_abundance_per_m2"
x=get(sites,K);m=get(sites,B)
rho=corr(x,m)
assert abs(rho-data["primary_stream_level_EPT_adult_share_vs_total_emergence_biomass"]["rho"])<1e-12
variants={
 "original_equal_interval_mean_share":rho,
 "pooled_count_weighted_EPT_abundance_share":corr(get(sites,"pooled_abundance_weighted_EPT_fraction"),m),
 "pooled_biomass_weighted_EPT_share":corr(get(sites,"pooled_biomass_weighted_EPT_fraction"),m),
 "mean_interval_EPT_biomass_share":corr(get(sites,"mean_interval_EPT_adult_biomass_fraction"),m),
 "original_EPT_abundance_share_vs_adult_abundance":corr(x,get(sites,A)),
 "original_EPT_abundance_share_vs_average_Chironomidae_biomass":corr(x,get(sites,"mean_interval_Diptera_adult_biomass_g_per_m2")),
}
loo={z["site"]:corr(get(sites[:i]+sites[i+1:],K),get(sites[:i]+sites[i+1:],B)) for i,z in enumerate(sites)}
regime_r={}
for reg in ("I","P"):
    sub=[z for z in sites if z["regime"]==reg]
    regime_r[reg]={"sites":len(sub),"within_regime_rank_rho":corr(get(sub,K),get(sub,B)),
        "mean_EPT_fraction":statistics.mean(get(sub,K)),
        "mean_biomass_gm2_per_event":statistics.mean(get(sub,B))}
counts=[z["intervals_with_positive_adult_abundance"] for z in sites]
supported=[z for z in sites if z["intervals_with_positive_adult_abundance"]>=20]
ratios={}
for var in [B,A,K]:
    vals=get(sites,var)
    ratios[var]={"min":min(vals),"median":statistics.median(vals),"max":max(vals),
      "mean":statistics.mean(vals),"standard_deviation_site":statistics.stdev(vals)}
ratio=max(m)/min(v for v in m if v>0)
median_x=statistics.median(x);median_m=statistics.median(m)
n_low_ept_high_biomass=sum(z[K]<median_x and z[B]>median_m for z in sites)
n_high_ept_low_biomass=sum(z[K]>median_x and z[B]<median_m for z in sites)
out={
 "source_dataset":"HJA AS015 2022 adult emergence, original v2; DOI 10.6073/pasta/6e9b83e2547e3a82d92dab068476f1d7",
 "source_report_original_HJA_sha256":data["original_2022_CSV_SHA256"],
 "site_n":17,"events":data["site_x_period_events"],"source_taxon_rows":data["source_taxon_rows"],
 "manual_pure_python_reproduces_original_Spearman":True,
 "EPT_share_biomass_rank_correlation_variants":variants,
 "leave_one_stream_out_Spearman":{
   "min":min(loo.values()),"max":max(loo.values()),
   "negative_count":sum(r<0 for r in loo.values()),"positive_count":sum(r>0 for r in loo.values()),
   "per_stream_omission":loo},
 "regime_subgroup_Spearman":regime_r,
 "stream_level_distributions":ratios,
 "biomass_max_min_positive_ratio":ratio,
 "stream_periods_with_positive_adults":{"min":min(counts),"max":max(counts),"median":statistics.median(counts)},
 "at_least_20_active_periods_sites":len(supported),
 "at_least_20_active_periods_spearman":corr(get(supported,K),get(supported,B)) if len(supported)>=3 else None,
 "discordant_counts":{
     "low_EPT_abundance_share_high_total_emergence_biomass":n_low_ept_high_biomass,
     "high_EPT_share_low_total_emergence_biomass":n_high_ept_low_biomass,
     "half_split_conditions":"Original site-equal 17-stream medians, posthoc descriptive, no site labels needed"},
 "limitations":[
  "17 stream records all in one watershed, not independent continents, not 391 replicated rivers",
  "Source collects only 4 orders (Diptera mostly Chironomidae), not all adult insect taxa",
  "Saturated high EPT fraction can obscure any real correlation; absence of positive correlation is NOT proof of ecological independence",
  "Within adult-emergence data we can compare relative EPT and total biomass; cannot validate that benthic EPT fraction in NRSA is adult emergence EPT fraction on same sites",
  "River regime and elevation not experimentally assigned; cross-regime differences not causal",
  "Biomass in g/m2 per sampled 3–4-day interval, not annual production or dry mass flux per day"
 ]}
outdir=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_AS015_independent"
outdir.mkdir(exist_ok=True)
(outdir/"HJA_AS015_17STREAM_EPT_PROXY_INDEPENDENT_QA_20261112.json").write_text(json.dumps(out,sort_keys=True,indent=2))
compact={k:v for k,v in out.items() if k not in ("limitations","leave_one_stream_out_Spearman")}
compact["leave_one_stream_out_Spearman"]={k:v for k,v in out["leave_one_stream_out_Spearman"].items() if k!="per_stream_omission"}
print("AS015_17STREAM_PROXY_INDEPENDENT_AUDIT_START")
print(json.dumps(compact,sort_keys=True))
print("AS015_17STREAM_PROXY_INDEPENDENT_AUDIT_END")
