#!/usr/bin/env python3
"""Independent reviewer implementation retained for reproducible verification.

Reconstructs the original source, enumerates every fixed-count alternative,
and uses analytic Hessians/Newton updates. Does not import the production
source builder, likelihood kernel, or inference module. Writes only stdout.
"""
import os, itertools, collections, hashlib, json, pathlib,sys,time
import numpy as np, pandas as pd
from scipy import stats
from scipy.special import logsumexp
import argparse
parser=argparse.ArgumentParser(description='Independent raw-source enumeration/Newton audit; fits the frozen checks only.')
parser.add_argument('--source',type=pathlib.Path,required=True)
parser.add_argument('--output',type=pathlib.Path,required=True)
args=parser.parse_args()
cache=args.source; out=args.output
res=json.loads((out/'results.json').read_text()); p=res['primary']; names=list(p['coefficients']); controls=names[1:]; GATE=res['source']
assert all(hashlib.sha256((cache/row['file']).read_bytes()).hexdigest()==row['sha256'] for row in GATE['source']['source_files'])
tax=pd.read_csv(cache/'GESUI_303_HOST_NUTRIENT_SCORES.csv',dtype={'butterflyKey':str}).set_index('butterflyKey')
tax=tax.loc[(tax.taxon_status=='ACCEPTED') & ~tax.actual_used_alias_conflict & tax.host_nutrient_mean.notna()].copy();tax['gen']=tax.accepted_name.str.split().str[0]
assert tax.index.is_unique and tax['accepted_name'].is_unique and len(tax)==186 and len(set(tax.gen))==88
assert np.max(abs(tax.host_nutrient_mean-tax.score_sum/tax.n_scored_hosts))<1e-12
assert np.max(abs(tax.trait_coverage-tax.n_scored_hosts/tax.n_documented_host_concepts))<1e-12
plant=pd.read_csv(cache/'EU454_PLANT_11511_SOURCE_EIV_EXACT_CROSSWALK.csv',dtype={'speciesKey':str}).dropna(subset=['EIV_NUTRI_exact']).set_index('speciesKey')
assert len(plant)==4011 and set(plant.gbif_taxon_status)=={'ACCEPTED'} and set(plant.gbif_rank)=={'SPECIES'}
panel=pd.read_csv(cache/'GENUINE_TERRESTRIAL_ONLY_183_H7_40_H6_PAIRED_H7_MEANS.csv',dtype={'h':str,'h6':str}).set_index('h')
facets={}
for fname in [GATE['source']['source_files'][0]['file'],GATE['source']['source_files'][1]['file']]:
 for row in json.loads((cache/fname).read_text()):
  if str(row['HYBAS_ID_L07']) in panel.index:
   key=str(row['HYBAS_ID_L07']),row['period'],row['group']
   assert key not in facets; facets[key]=row['source_speciesKey_counts']
period=['2017-2019','2020-2022']
X=[]; occ=[]; env=[]; cid=[]; cl=[]
for h,q in panel.iterrows():
 for per in period:
  raw=facets[h,per,'Plantae']; sel=[k for k,v in raw.items() if v>0 and k in plant.index]; cwm=float(plant.loc[sel,'EIV_NUTRI_exact'].mean())
  sid=set()
  for f in tax.family.unique():
   for k,v in facets[h,per,f].items():
    if v>0 and k in tax.index:
     assert tax.loc[k,'family']==f;sid.add(k)
  env.append(cwm);X.append([cwm,*q[controls]]);occ.append(sid);cid.append((h,per));cl.append((str(q.h6),per))
X=np.array(X);X=(X-X.mean(axis=0))/X.std(axis=0); reg=np.array([c[0] for c in cl]);sid=np.array(tax.index.tolist()); vals=tax.host_nutrient_mean.to_dict();gens=tax.gen.to_dict();M=np.array([[k in d for k in sid] for d in occ])
assert abs(np.array(env).std()-p['source_scaling'][names[0]][1])<1e-12
sets={c:np.where(np.array([v==c for v in cl]))[0] for c in sorted(set(cl))}
def strata(tax):
 st=[]
 for g,t in tax.groupby('gen'):
  if t.host_nutrient_mean.nunique()<2:continue
  for a,b in itertools.combinations(sorted(t.index),2):
   ca,cb=np.where(sid==a)[0][0],np.where(sid==b)[0][0];sel=M[:,ca]^M[:,cb]
   for (region,per),ix in sets.items():
    row=ix[sel[ix]];obs=M[row,ca]
    if obs.size and obs.min()<obs.max():
     st.append((region,per,g,a,b,row,obs))
 return st
st=strata(tax); hi=strata(tax[tax.trait_coverage>=.8])
def sc(st):
 diff=[s for s in st if vals[s[3]]!=vals[s[4]]];rows=set([c for s in diff for c in s[5]]);return dict(strata=len(st),informative_nonzero=len(diff),rows=sum(len(s[5]) for s in st),rows_nonzero=sum(len(s[5]) for s in diff),pairs=len(set((s[3],s[4]) for s in st)),pairs_nonzero=len(set((s[3],s[4]) for s in diff)),species=len(set(k for s in st for k in s[3:5])),species_nonzero=len(set(k for s in diff for k in s[3:5])),genera=len(set(s[2] for s in diff)),h6=len(set(s[0] for s in diff)),h7=len(set(cid[c][0] for c in rows)),h7period=len(rows),all_h7period=len(set(c for s in st for c in s[5])), stratum_hist=dict(collections.Counter(len(s[5]) for s in st)),pair_strata_hist=dict(collections.Counter(collections.Counter((s[3],s[4]) for s in st).values())))
print('PRIMARY SUPPORT',sc(st),'COVERAGE80',sc(hi),flush=True)
weightcnt=collections.Counter((s[3],s[4]) for s in st);w=np.array([len(st)/len(weightcnt)/weightcnt[s[3],s[4]] for s in st]);clreg=np.array([s[0] for s in st]);allreg=np.unique(clreg);genera=np.array([s[2] for s in st]);pgen=np.unique(genera)
used=sorted({k for s in st for k in s[3:5]}); grouped={g:[k for k in used if gens[k]==g] for g in pgen}
za=np.array([used.index(s[3]) for s in st]);zb=np.array([used.index(s[4]) for s in st]);v0=np.array([vals[k] for k in used]);d0=v0[za]-v0[zb]
allalts=[]; ix=[]
for i,s in enumerate(st):
 row=s[5]; obs=s[6]; o=X[row[obs]].sum(axis=0)
 for ss in itertools.combinations(row,int(sum(obs))):
  allalts.append(X[list(ss)].sum(axis=0)-o);ix.append(i)
ix=np.array(ix);al=np.array(allalts);beg=np.r_[0,np.where(np.diff(ix))[0]+1]
def valgrad(theta,A,W,hess=False):
 eta=A@theta; vmax=np.maximum.reduceat(eta,beg);e=np.exp(eta-vmax[ix]);sumE=np.add.reduceat(e,beg);prob=e/sumE[ix];expA=np.add.reduceat(prob[:,None]*A,beg,axis=0)
 loss=np.dot(W,vmax+np.log(sumE));grad=W@expA
 if hess:
  diff=A-expA[ix]
  return loss,grad,diff.T@(diff*(prob*W[ix])[:,None]),-W[:,None]*expA
 return loss,grad
def refit(A,W,theta):
 theta=theta.copy()
 for i in range(100):
  loss,g,H,S=valgrad(theta,A,W,True)
  if abs(g).max()/sum(W)<4e-10:break
  inc=-np.linalg.solve(H,g); l=1
  while l>1e-7 and valgrad(theta+l*inc,A,W)[0]>loss+1e-11:l*=.5
  theta+=l*inc
 else:raise ValueError('Newton failed')
 return theta,H,S
t,H,S=refit(al*d0[ix,None],w,np.zeros(len(names)));eig=np.linalg.eigvalsh(H)
U=np.stack([S[clreg==c].sum(axis=0) for c in allreg]);Inv=np.linalg.inv(H); se=np.sqrt((len(U)/(len(U)-1)*Inv@U.T@U@Inv)[0,0]);print('PRIMARY EXACT beta,se,eigs,normalized solution difference',t[0],se,eig[[0,-1]],abs(t-np.array(p['fit']['params'])).max(),flush=True)
assert abs(t-np.array(p['fit']['params'])).max()<2e-6 and eig.min()>0
assert abs(se-p['H6_cluster_standard_error'])<2e-6
keep=[i for i,n in enumerate(names) if n not in ['nButterfly','hns_coverage','hns_host_completeness','hns_scored_species']]; t2,H2,S2=refit(al[:,keep]*d0[ix,None], w,np.zeros(len(keep)));print('SENS EXACT',t2[0],flush=True)
assert abs(t2[0]-res['sensitivity_without_four_butterfly_derived_controls']['beta'])<2e-6
def run_draws(kind,reps,refpath):
 prev=pd.read_csv(out/refpath); diff=[];mins=[];maxcoefs=[];rng=np.random.default_rng(20261122 if kind=='bootstrap' else 20261123);group_ix=np.searchsorted(allreg,clreg); zix={g:[used.index(k) for k in members] for g,members in grouped.items()}
 for j in range(reps):
  if kind=='bootstrap':
   mult=rng.multinomial(len(allreg),np.full(len(allreg),1/len(allreg)));W=w*mult[group_ix];D=d0
  elif kind=='remove':
   W=w*(genera!=pgen[j]);D=d0;assert prev.iloc[j].removed_butterfly_genus==pgen[j]
  else:
   v=v0.copy()
   for g in pgen:v[zix[g]]=rng.permutation(v0[zix[g]])
   W=w;D=v[za]-v[zb]
  b,I,S=refit(al*D[ix,None],W,np.zeros(len(names)))
  diff.append(b[0]-prev.iloc[j].beta);mins.append(np.linalg.eigvalsh(I)[0]);maxcoefs.append(abs(b).max())
 assert max(abs(np.array(diff)))<2e-6 and min(mins)>0
 print('INDEPENDENT',kind,reps,'max_beta_abs_difference',max(abs(np.array(diff))),'minimum_information_eigenvalue',min(mins),'max_absolute_nuisance_or_beta',max(maxcoefs),flush=True)
for a in [('bootstrap',399,'bootstrap_H6.csv'),('labels',999,'within_butterfly_genus_label_reference.csv'),('remove',29,'leave_one_butterfly_genus_out.csv')]:run_draws(*a)
