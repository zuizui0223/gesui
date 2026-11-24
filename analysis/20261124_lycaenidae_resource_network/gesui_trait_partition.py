#!/usr/bin/env python3
"""Exploratory source-locked plant trait projection under training-only early null.
Fixed before trait model results. No observed resource biomass or local management measurements.
"""
from pathlib import Path
import sys
P=Path(sys.argv[1])
prior=P/'gesui_strict_reference_falsifier.py'
exec(compile(prior.read_text().split('rows=[];foldrows=[];')[0],str(prior),'exec'))
import requests,io
# Byte hashes from verified source-only run 45381165806, artifact 14915055021.
sources={
'life':('https://files.ibot.cas.cz/cevs/downloads/floraveg/Life_form.xlsx','cdbe95457bfb0a3154f1d676f73b965f42122d996a02b7c2362ff5f2f27df14a'),
'disturbance':('https://zenodo.org/records/7116802/files/disturbance_indicator_values.csv?download=1','7596c9354018816d408e051d59437b3303666d14ec1cecd170f38b9f7206a04d')}
raws={}
for k,(url,sha) in sources.items():
 r=requests.get(url,timeout=(15,60),headers={'User-Agent':'gesui academic trait analysis'});r.raise_for_status()
 assert hashlib.sha256(r.content).hexdigest()==sha,k
 raws[k]=r.content
life=pd.read_excel(io.BytesIO(raws['life']))
dist=pd.read_csv(io.BytesIO(raws['disturbance']))
assert len(life)==16402 and len(dist)==6382
names=defaultdict(set)
for x in cross:
 k=x['accepted_speciesKey'];m=x['match']
 if k in pidx and m.get('species'):names[k].add(' '.join(m['species'].casefold().split()))
hostnames=[next(iter(names[k])) if len(names[k])==1 else None for k in pk]
families=np.array([pf[k] for k in pk])
def trait_matrix(table,key,cols,binary=False):
 assert key in table and all(c in table for c in cols)
 table=table.copy();table['name']=table[key].astype(str).str.casefold().str.replace(r'\s+',' ',regex=True).str.strip()
 table[cols]=table[cols].apply(pd.to_numeric,errors='coerce')
 lookup={}
 for name,g in table.groupby('name'):
  vals=g[cols].dropna().drop_duplicates()
  if len(vals)!=1:continue
  v=vals.iloc[0].to_numpy(float)
  if not np.isfinite(v).all():continue
  if binary and (not np.isin(v,[0,1]).all() or v.sum()==0):continue
  lookup[name]=v
 return np.array([lookup.get(k,np.full(len(cols),np.nan)) for k in hostnames])
specs={
 'life_form':(life,'FloraVeg.Taxon',['Phanerophyte','Chamaephyte','Hemicryptophyte','Geophyte','Hydrophyte','Therophyte','Epiphyte'],True),
 'management_affinity':(dist,'species',['Mowing.Frequency','Grazing.Pressure','Soil.Disturbance'],False),
 'disturbance_affinity':(dist,'species',['Disturbance.Severity','Disturbance.Frequency'],False)}
rows=[];gates=[];lossrows=[];contrasts=[]
for domain,(table,key,cols,binary) in specs.items():
 T=trait_matrix(table,key,cols,binary);known=np.isfinite(T).all(1);groups=[];valid=np.zeros(len(pk),bool)
 for fam in sorted(set(families)):
  idx=np.flatnonzero((families==fam)&known)
  if len(idx)<max(5,len(cols)+2):continue
  z=T[idx]-T[idx].mean(0);rank=np.linalg.matrix_rank(z)
  if rank<1 or len(idx)<rank+2:continue
  proj=z@np.linalg.pinv(z)
  np.testing.assert_allclose(proj@proj,proj,atol=1e-8,rtol=0)
  groups.append((idx,proj));valid[idx]=True
 den0=PA.sum(1);den1=PC.sum(1)
 cov0=np.divide(PA[:,valid].sum(1),den0,out=np.zeros(len(q)),where=den0>0)
 cov1=np.divide(PC[:,valid].sum(1),den1,out=np.zeros(len(q)),where=den1>0)
 for tag,mask0 in [('quality',np.asarray(quality)),('all',np.ones(len(q),bool))]:
  # Support gate uses source/trait coverage only. Retain all STOPs and do not relax thresholds.
  mask=mask0&(cov0>=.5)&(cov1>=.5)
  ids=np.flatnonzero(mask);ng=len(set(H6[ids]));gate=dict(domain=domain,sample=tag,exact_scored_keys=int(known.sum()),eligible_keys=int(valid.sum()),eligible_families=len(groups),median_coverage_early=float(np.median(cov0[mask0])),median_coverage_late=float(np.median(cov1[mask0])),H7=len(ids),H6=ng,H7_at_80=int((mask0&(cov0>=.8)&(cov1>=.8)).sum()))
  gate['status']='FIT' if known.sum()>=100 and len(ids)>=40 and ng>=20 else 'STOP_SUPPORT';gates.append(gate)
  if gate['status']!='FIT':continue
  losses=defaultdict(list);errmax=0
  for f in sorted(set(fold[ids])):
   te=ids[fold[ids]==f];tr=ids[fold[ids]!=f];assert not set(H6[tr])&set(H6[te]);assert len(set(H6[tr]))>=15
   M,_=reference(tr);C=.5*(BA+BB)@M
   full=((PC-PA)*C).sum(1);matched=((PC[:,valid]-PA[:,valid])*C[:,valid]).sum(1)
   identity=np.zeros(len(q));trait=identity.copy()
   for idx,proj in groups:
    x=PA[:,idx];y=PC[:,idx];c=C[:,idx];m0=x.sum(1);m1=y.sum(1)
    u0=np.divide(x,m0[:,None],out=np.zeros_like(x),where=m0[:,None]>1e-14);u1=np.divide(y,m1[:,None],out=np.zeros_like(y),where=m1[:,None]>1e-14)
    u0=np.where(m0[:,None]>1e-14,u0,u1);u1=np.where(m1[:,None]>1e-14,u1,u0)
    change=.5*(m0+m1)[:,None]*(u1-u0)
    identity+=(change*c).sum(1);trait+=(change*(c@proj)).sum(1)
   ys={'full_plant_delta':full,'matched_delta':matched,'unmatched_delta':full-matched,'matched_quantity':matched-identity,'trait_identity':trait,'residual_identity':identity-trait}
   np.testing.assert_allclose(ys['matched_quantity']+ys['trait_identity']+ys['residual_identity']+ys['unmatched_delta'],full,rtol=0,atol=1e-12)
   early=(BA@M*PA).sum(1);X=np.column_stack((X0,early));X1=np.column_stack((X,E0));cat=np.array([f'{c}::{m}' for c,m in zip(country,method)])
   for name,y in ys.items():
    p0,e0,_=fit(tr,te,y,X,cat);p1,e1,c=fit(tr,te,y,X1,cat);errmax=max(errmax,e0,e1)
    for g in sorted(set(H6[te])):
     z=H6[te]==g;l0=float(np.mean((y[te][z]-p0[z])**2));l1=float(np.mean((y[te][z]-p1[z])**2));losses[name].append((g,l0,l1));lossrows.append(dict(domain=domain,sample=tag,endpoint=name,h6=g,mse0=l0,mse1=l1))
  boot={}
  for name,ls in losses.items():
   df=pd.DataFrame(ls,columns=['h6','mse0','mse1']).sort_values('h6');assert not df.h6.duplicated().any()
   m=df[['mse0','mse1']].to_numpy();l0,l1=m.mean(0);idx=np.random.default_rng(20261124).integers(0,len(m),(5000,len(m)));v=m[idx].mean(1)
   if l0<=1e-20:rows.append(dict(domain=domain,sample=tag,endpoint=name,H7=len(ids),H6=ng,status='STOP_CONSTANT'));continue
   bs=100*(v[:,0]-v[:,1])/v[:,0];boot[name]=bs;lo,hi=np.quantile(bs,[.025,.975]);rows.append(dict(domain=domain,sample=tag,endpoint=name,H7=len(ids),H6=ng,status='FIT',gain=100*(l0-l1)/l0,lo=lo,hi=hi,independent_prediction_max_diff=errmax))
  if 'trait_identity' in boot and 'residual_identity' in boot:
   bs=boot['trait_identity']-boot['residual_identity'];lo,hi=np.quantile(bs,[.025/3,1-.025/3]);contrasts.append(dict(domain=domain,sample=tag,contrast='trait_minus_residual',bonf3_lo=lo,bonf3_hi=hi))
out=P/'trait_partition_output';out.mkdir(exist_ok=True)
for f,rr in [('source_coverage',gates),('results',rows),('losses',lossrows),('contrasts',contrasts)]:
 pd.DataFrame(rr).to_csv(out/(f+'.csv'),index=False)
print('TRAIT SOURCE GATES');print(pd.DataFrame(gates).to_string(index=False))
print('TRAIT PROJECTION RESULTS');print(pd.DataFrame(rows).to_string(index=False))
(out/'source_receipt.json').write_text(json.dumps({'sources':sources,'source_gate_artifact':14915055021,'protocol':'source matched 50% mass coverage both periods, min40 H7/20 H6, min100 scored hosts; training-only null; fixed 3 trait domains; all source failures retained','notes':'trait projection is algebraic, not causal; source EIV/disturbance are ecological optima; no direct management, N-fixation or actual feeding measurement'},indent=2))
