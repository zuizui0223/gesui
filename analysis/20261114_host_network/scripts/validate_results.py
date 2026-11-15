"""Independent arithmetic and taxonomy audit of candidate host opportunities."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import pathlib,json,hashlib,collections
import numpy as np,pandas as pd
from scipy.stats import hypergeom
from sklearn.cluster import KMeans
R=pathlib.Path(__file__).resolve().parents[1];d=pd.read_csv(R/'output/H7_candidate_host_network_panel.csv',dtype={'h':str,'h6':str})
z=np.load(R/'output/candidate_network_matrices.npz');B,P,A,NB,NP=[z[k] for k in ['B','P','A','NB','NP']]
# Independent SciPy hypergeometric survival probabilities, not original product code.
bprob=np.array([hypergeom.sf(0,int(n),b,20) for n,b in zip(NB,B)])
pprob=np.array([hypergeom.sf(0,int(n),p,50) for n,p in zip(NP,P)])
known=A.sum(1)>0
# Sparse edge sums differ from matrix dot-product construction.
i,j=np.nonzero(A)
links=np.sum(bprob[:,i]*pprob[:,j],axis=1)
counts=np.column_stack([P[:,np.flatnonzero(A[b])].sum(1) for b in range(len(A))])
qp=np.array([hypergeom.sf(0,int(n),c,50) for n,c in zip(NP,counts)])
support=np.sum(bprob[:,known]*qp[:,known],axis=1)/bprob[:,known].sum(1)
assert np.max(abs(links-d.expected_links))<1e-8
assert np.max(abs(support-d.host_support))<1e-8
# Raw response content hashes were checked independently of the collector.
receipts=[json.loads(s) for s in (R/'source/source_receipts.jsonl').read_text().splitlines()]
for r in receipts:
 raw=(R/'source/responses'/r['cache_file']).read_bytes()
 assert hashlib.sha256(raw).hexdigest()==r['response_sha256']
# Same graph null, now ALSO constrain similar observed H6 geographic distributions.
# This is a labelled posthoc sensitivity: exactly six unsupervised clusters, no selection by result.
hs=sorted(set(d.h6));inc=np.array([(P[d.h6==h]>0).any(0) for h in hs]).T.astype(float)
norm=np.linalg.norm(inc,axis=1);nonzero=norm>0
labels=np.full(len(norm),-1);labels[nonzero]=KMeans(n_clusters=6,random_state=20261114,n_init=10).fit_predict(inc[nonzero]/norm[nonzero,None])
strata=np.array([s+'|geo'+str(k) for s,k in zip(z['strata'],labels)]);groups=[np.flatnonzero(strata==s) for s in sorted(set(strata))]
weights=1/d.groupby('h6').h6.transform('size').to_numpy();weights/=weights.sum()
rng=np.random.default_rng(20261114);ba=bprob@A;null=np.empty((499,len(d)))
for it in range(499):
 perm=np.arange(len(strata))
 for g in groups:perm[g]=rng.permutation(g)
 assert np.array_equal(strata[perm],strata)
 assert np.array_equal(A[:,perm].sum(1),A.sum(1))
 null[it]=np.sum(ba[:,perm]*pprob,axis=1)
na=null@weights;obs=float(weights@links)
qa={'source_response_sha_checks':len(receipts),'max_expected_link_difference_SciPy_sparse':float(np.max(abs(links-d.expected_links))),'max_host_support_difference_SciPy':float(np.max(abs(support-d.host_support))),'all_metric_rows_checked':len(d),'geography_constrained_null':{'fixed_clusters':6,'draws':499,'movable_host_keys':sum(len(g) for g in groups if len(g)>1),'host_keys_total':len(strata),'observed_H6_equal_expected_links':obs,'null_H6_equal_mean':float(na.mean()),'null_H6_equal95':np.quantile(na,[.025,.975]).tolist(),'relative_excess_percent':100*(obs/na.mean()-1),'null_at_least_observed_plus1':float((1+sum(na>=obs))/500),'posthoc_not_causal':True}}
(R/'output/independent_numerical_source_QA.json').write_text(json.dumps(qa,indent=2));print(json.dumps(qa,indent=2))
# Descriptive source contributions: no local feeding or wastewater slope inferred.
plants=json.loads((R/'source/plant_name_crosswalk.json').read_text());names={}
for p in plants:
 k=p['accepted_speciesKey']
 if k is None:continue
 n=p['match'].get('species') or p['match'].get('canonicalName')
 if k not in names or p['match'].get('status')=='ACCEPTED':names[k]=n
contrib=weights@((bprob@A)*pprob)
keys=z['pkeys'];out=pd.DataFrame({'plantKey':keys,'accepted_name':[names[k] for k in keys],'H6_equal_expected_link_contribution':contrib,'contribution_fraction':contrib/contrib.sum(),'catalogue_butterfly_degree':A.sum(0).astype(int),'recorded_basin_periods':(P>0).sum(0)})
out=out.sort_values('H6_equal_expected_link_contribution',ascending=False);out.to_csv(R/'output/host_plant_candidate_contributions.csv',index=False)
print(out.head(12).to_string(index=False))
