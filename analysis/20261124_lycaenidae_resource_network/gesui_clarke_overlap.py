#!/usr/bin/env python3
"""Post-expansion fixed-cohort edge overlap diagnostic; no actual feeding inference."""
from pathlib import Path
import sys
P=Path(sys.argv[1]);prior=P/'gesui_clarke_expanded_compare.py'
exec(compile(prior.read_text().split('rows=[];losses=[];contrasts=[];gates=[];statics=[]')[0],str(prior),'exec'))
rows=[];losses=[];contrasts=[];structures=[];maxid=0.
for cutoff in [3,2]:
 restricted=graph(oldedge,cutoff);expanded=graph(newedge,cutoff);both=accepted&(A.sum(1)>0)&(restricted.sum(1)>0);assert both.sum()==(46 if cutoff==3 else 45)
 ao=oldU[both];ac=expanded[both];inter=ao*ac;gr={'shared':inter,'old_only':ao-inter,'clarke_only':ac-inter};assert (ao==gr['shared']+gr['old_only']).all() and (ac==gr['shared']+gr['clarke_only']).all()
 for k,g in gr.items():structures.append(dict(cutoff=cutoff,component=k,edges=int(g.sum()),plants=int((g.sum(0)>0).sum()),butterflies=int((g.sum(1)>0).sum())))
 frac=np.divide(PB[:,both].sum(1),denom,out=np.zeros(len(d)),where=denom>0);cov=(frac[a]>=.8)&(frac[b]>=.8);b0,b1=PB[a][:,both],PB[b][:,both];pa,pc=UP[a],UP[b]
 for tag,sm in [('quality',np.asarray(quality)),('all',np.ones(len(q),bool))]:
  ids=np.flatnonzero(sm&cov);ng=len(set(H6[ids]));assert len(ids)==(159 if tag=='quality' else 200);store=collections.defaultdict(list)
  for f in sorted(set(fold[ids])):
   te=ids[fold[ids]==f];tr=ids[fold[ids]!=f];assert not set(H6[tr])&set(H6[te]);st=strats(tr,nc);cs={k:decomp(g,b0,b1,st,pa,pc) for k,g in gr.items()}
   for n,full in [('old_only',ao),('clarke_only',ac)]:
    e,y=decomp(full,b0,b1,st,pa,pc);maxid=max(maxid,float(np.max(np.abs(e-cs['shared'][0]-cs[n][0]))));np.testing.assert_allclose(e,cs['shared'][0]+cs[n][0],atol=1e-12,rtol=0)
    for ep in y:
     maxid=max(maxid,float(np.max(np.abs(y[ep]-cs['shared'][1][ep]-cs[n][1][ep]))));np.testing.assert_allclose(y[ep],cs['shared'][1][ep]+cs[n][1][ep],atol=1e-12,rtol=0)
   for gn,(early,ys) in cs.items():
    y=ys['plant_symmetric'];X=np.column_stack((X0,early));X1=np.column_stack((X,E0));p0,e0,_=fit(tr,te,y,X,cat);p1,e1,_=fit(tr,te,y,X1,cat)
    for h6 in sorted(set(H6[te])):
     zz=H6[te]==h6;l0=float(np.mean((y[te][zz]-p0[zz])**2));l1=float(np.mean((y[te][zz]-p1[zz])**2));store[gn].append((h6,l0,l1));losses.append(dict(cutoff=cutoff,sample=tag,component=gn,h6=h6,mse0=l0,mse1=l1,solver_diff=max(e0,e1)))
  idx=np.random.default_rng(20261124).integers(0,ng,(5000,ng));boot={}
  for gn,v in store.items():
   df=pd.DataFrame(v,columns=['h6','mse0','mse1']).sort_values('h6');assert len(df)==ng and not df.h6.duplicated().any();m=df[['mse0','mse1']].to_numpy();l0,l1=m.mean(0);assert l0>1e-20;gain=100*(l0-l1)/l0;vv=m[idx].mean(1);bs=100*(vv[:,0]-vv[:,1])/vv[:,0];lo,hi=np.quantile(bs,[.025,.975]);boot[gn]=gain,bs;rows.append(dict(cutoff=cutoff,sample=tag,component=gn,H7=len(ids),H6=ng,gain=gain,lo=lo,hi=hi))
  for gn in ['old_only','clarke_only']:
   v=boot[gn][1]-boot['shared'][1];lo,hi=np.quantile(v,[.0125,.9875]);contrasts.append(dict(cutoff=cutoff,sample=tag,A=gn,B='shared',difference=boot[gn][0]-boot['shared'][0],bonf2_lo=lo,bonf2_hi=hi))
for name,v in [('structure',structures),('results',rows),('losses',losses),('contrasts',contrasts)]:pd.DataFrame(v).to_csv(P/f'clarke_overlap_{name}.csv',index=False)
print(pd.DataFrame(structures).to_string(index=False));print(pd.DataFrame(rows).round(3).to_string(index=False));print(pd.DataFrame(contrasts).round(3).to_string(index=False));print('maximum decomposition error',maxid,'maximum independent solver difference',max(x['solver_diff'] for x in losses))
(P/'clarke_overlap_receipt.json').write_text(json.dumps({'fits':len(rows),'loss_rows':len(losses),'max_decomposition_error':maxid,'max_solver_difference':max(x['solver_diff'] for x in losses)},indent=2)+'\n')
