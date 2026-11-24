#!/usr/bin/env python3
"""Expanded Clarke support audit; fixed prior cohorts, plus separate expanded-support sensitivity."""
from pathlib import Path
import sys,hashlib,json
P=Path(sys.argv[1]);prior=P/'gesui_strict_reference_falsifier.py'
exec(compile(prior.read_text().split('rows=[];foldrows=[];')[0],str(prior),'exec'))
# Restore variable P in case older code defaulted to a different directory.
P=Path(sys.argv[1]); rr=json.loads((P/'clarke_expanded_mapping_receipt.json').read_text());assert rr['status']=='PASS_SOURCE_MAPPING'
for file,key in [('clarke_edges.csv','original_edge_sha256'),('clarke_expanded_edges.csv','expanded_edge_sha256'),('clarke_expanded_plant_keys.csv','plant_keys_sha256')]:assert hashlib.sha256((P/file).read_bytes()).hexdigest()==rr[key]
oldedge=pd.read_csv(P/'clarke_edges.csv',dtype={'butterflyKey':str,'plantKey':str});newedge=pd.read_csv(P/'clarke_expanded_edges.csv',dtype={'butterflyKey':str,'plantKey':str});nt=pd.read_csv(P/'clarke_expanded_plant_keys.csv',dtype={'plantKey':str});assert len(nt)==rr['new_plants']
newpk=nt.plantKey.tolist(); assert not set(newpk)&set(pk)
# Preserve original column order for exact KMeans replay. Add new columns at end.
upk=pk+newpk; upidx={k:i for i,k in enumerate(upk)};families=np.array([pf[k] for k in pk]+nt.family.tolist());nc=len(upk);nn=len(newpk)
extra=np.zeros((len(d),nn),int);ei={k:i for i,k in enumerate(newpk)}
for ii,x in enumerate(d.itertuples()):
 for k,v in fac[x.h,x.period,'Plantae']['source_speciesKey_counts'].items():
  if k in ei:extra[ii,ei[k]]=v
assert (extra.sum(0)>0).all()
UPcount=np.column_stack((Pcount,extra));UP=incl(np.array(NP)[:,None],UPcount,50);np.testing.assert_allclose(UP[:,:990],PP,atol=1e-13,rtol=0)
oldU=np.pad(A,((0,0),(0,nn))); accepted=L&np.array([bmeta[k]['status']=='ACCEPTED' for k in bk]);oldmapped=L&(A.sum(1)>0);denom=PB[:,oldmapped].sum(1)
cat=np.array([f'{c}::{m}' for c,m in zip(country,method)])
def graph(ed,cutoff):
 g=np.zeros_like(oldU)
 for x in ed[ed.status<=cutoff].itertuples():g[bidx[x.butterflyKey],upidx[x.plantKey]]=1
 return g
def strats(tr,n):
 ix=a[tr];hs=sorted(set(d.h6.iloc[ix]));assert len(hs)>=15
 C=UPcount[:,:n];freq=(C[ix]>0).sum(0);bins=np.digitize(freq,[0,1,5,20,100,468],right=True)
 z=np.array([(C[ix[d.h6.iloc[ix].to_numpy()==h]]>0).any(0) for h in hs]).T.astype(float);norm=np.linalg.norm(z,axis=1);lab=np.full(n,-1);good=norm>0
 lab[good]=KMeans(n_clusters=6,random_state=20261114,n_init=10).fit_predict(z[good]/norm[good,None])
 ss=np.array([f+'|'+str(v)+'|'+str(g) for f,v,g in zip(families[:n],bins,lab)]);return [np.flatnonzero(ss==s) for s in sorted(set(ss))]
def decomp(g,b0,b1,st,pa,pc):
 N0=np.zeros_like(g)
 for j in st:N0[:,j]=g[:,j].mean(1)[:,None]
 np.testing.assert_allclose(g.sum(1),N0.sum(1),atol=1e-10,rtol=0)
 M=g-N0;u=b0@M;v=b1@M;dp=pc-pa;dc=v-u;early=(u*pa).sum(1);late=(v*pc).sum(1);p=(u*dp).sum(1);c=(dc*pa).sum(1);i=(dc*dp).sum(1)
 np.testing.assert_allclose(late-early,p+c+i,atol=1e-12,rtol=0)
 return early,{'plant_symmetric':p+.5*i,'consumer_symmetric':c+.5*i,'plant_early_consumer_fixed':p}
rows=[];losses=[];contrasts=[];gates=[];statics=[]
for cutoff in [3,2]:
 restricted=graph(oldedge,cutoff);expanded=graph(newedge,cutoff);assert (expanded>=restricted).all()
 for scope,maskgraph in [('fixed_prior',restricted),('expanded_support',expanded)]:
  both=accepted&(A.sum(1)>0)&(maskgraph.sum(1)>0);b0,b1=PB[a][:,both],PB[b][:,both];frac=np.divide(PB[:,both].sum(1),denom,out=np.zeros(len(d)),where=denom>0);cov=(frac[a]>=.8)&(frac[b]>=.8)
  if scope=='fixed_prior':assert both.sum()==(46 if cutoff==3 else 45)
  graphs={'old990':oldU[both,:990],'clarke990':restricted[both,:990],'old_union':oldU[both],'clarke_restricted_union':restricted[both],'clarke_expanded_union':expanded[both]}
  for name,g in graphs.items():statics.append(dict(cutoff=cutoff,scope=scope,graph=name,butterflies=int(both.sum()),plants=int((g.sum(0)>0).sum()),edges=int(g.sum()),universe=g.shape[1]))
  for tag,sm in [('quality',np.asarray(quality)),('all',np.ones(len(q),bool))]:
   ids=np.flatnonzero(sm&cov);ng=len(set(H6[ids]));gate=dict(cutoff=cutoff,scope=scope,sample=tag,common_butterflies=int(both.sum()),H7=len(ids),H6=ng,median_early_coverage=float(np.median(frac[a][sm])),median_late_coverage=float(np.median(frac[b][sm])))
   gate['status']='FIT' if both.sum()>=30 and min((graphs['old_union'].sum(0)>0).sum(),(graphs['clarke_expanded_union'].sum(0)>0).sum())>=60 and len(ids)>=40 and ng>=20 else 'STOP_SUPPORT';gates.append(gate)
   if scope=='fixed_prior':assert len(ids)==(159 if tag=='quality' else 200)
   if gate['status']!='FIT':continue
   store=collections.defaultdict(list)
   for f in sorted(set(fold[ids])):
    te=ids[fold[ids]==f];tr=ids[fold[ids]!=f];assert not set(H6[tr])&set(H6[te]); st0=strats(tr,990);stu=strats(tr,nc)
    for gn,g in graphs.items():
     n=g.shape[1];st=st0 if n==990 else stu;pa=UP[a,:n];pc=UP[b,:n];early,ys=decomp(g,b0,b1,st,pa,pc);X=np.column_stack((X0,early));X1=np.column_stack((X,E0))
     for ep,y in ys.items():
      p0,e0,_=fit(tr,te,y,X,cat);p1,e1,coef=fit(tr,te,y,X1,cat)
      for h6 in sorted(set(H6[te])):
       zz=H6[te]==h6;l0=float(np.mean((y[te][zz]-p0[zz])**2));l1=float(np.mean((y[te][zz]-p1[zz])**2));store[gn,ep].append((h6,l0,l1));losses.append(dict(cutoff=cutoff,scope=scope,sample=tag,graph=gn,endpoint=ep,h6=h6,mse0=l0,mse1=l1,solver_diff=max(e0,e1)))
   boot={};idx=np.random.default_rng(20261124).integers(0,ng,(5000,ng))
   for (gn,ep),v in store.items():
    df=pd.DataFrame(v,columns=['h6','mse0','mse1']).sort_values('h6');assert len(df)==ng and not df.h6.duplicated().any();m=df[['mse0','mse1']].to_numpy();x0,x1=m.mean(0)
    if x0<=1e-20:rows.append(dict(cutoff=cutoff,scope=scope,sample=tag,graph=gn,endpoint=ep,H7=len(ids),H6=ng,status='STOP_CONSTANT'));continue
    v0=m[idx].mean(1);bs=100*(v0[:,0]-v0[:,1])/v0[:,0];gain=100*(x0-x1)/x0;lo,hi=np.quantile(bs,[.025,.975]);boot[gn,ep]=(gain,bs);rows.append(dict(cutoff=cutoff,scope=scope,sample=tag,graph=gn,endpoint=ep,H7=len(ids),H6=ng,gain=gain,lo=lo,hi=hi,status='FIT'))
   for x,y in [('clarke_expanded_union','clarke_restricted_union'),('clarke_expanded_union','old_union')]:
    ep='plant_symmetric';g1,b1s=boot[x,ep];g0,b0s=boot[y,ep];diff=b1s-b0s;lo,hi=np.quantile(diff,[.025,.975]);l2,h2=np.quantile(diff,[.0125,.9875]);contrasts.append(dict(cutoff=cutoff,scope=scope,sample=tag,A=x,B=y,difference=g1-g0,lo=lo,hi=hi,bonf2_lo=l2,bonf2_hi=h2))
for name,v in [('gate',gates),('structure',statics),('results',rows),('losses',losses),('contrasts',contrasts)]:pd.DataFrame(v).to_csv(P/f'clarke_expanded_{name}.csv',index=False)
# Verify old 990 and restricted 990 against archived results, not only approximate agreement.
p0=pd.read_csv(P/'archived_clarke_results.csv')
r=pd.DataFrame(rows);cmp=r[(r.scope=='fixed_prior')&r.graph.isin(['old990','clarke990'])].copy();cmp['catalogue']=cmp.graph.map({'old990':'old','clarke990':'clarke'});cmp=cmp.rename(columns={'cutoff':'max_status'}).merge(p0,on=['max_status','sample','catalogue','endpoint'],suffixes=('_new','_archived'),validate='one_to_one');assert len(cmp)==24
for k in ['gain','lo','hi']:np.testing.assert_allclose(cmp[k+'_new'],cmp[k+'_archived'],rtol=0,atol=1e-8)
print('REPLAY 24 prior numerical rows PASS',flush=True);print(pd.DataFrame(gates).to_string(index=False));print(pd.DataFrame(statics).to_string(index=False));print(r[r.endpoint=='plant_symmetric'].round(3).to_string(index=False));print(pd.DataFrame(contrasts).round(3).to_string(index=False))
(P/'clarke_expanded_fit_receipt.json').write_text(json.dumps(dict(status='FIT_COMPLETE',replayed_prior_990_rows=24,fit_rows=len(rows),loss_rows=len(losses),maximum_solver_diff=max(x['solver_diff'] for x in losses),source_union_keys=nc,source_mapping_receipt=rr,notes='Fold-specific targets, fixed error bootstrap; same observations. Distinguish graph and null-universe changes.'),indent=2)+'\n')
