#!/usr/bin/env python3
"""Posthoc exact plant resource decomposition, not nutrient causality or abundance."""
from pathlib import Path
import sys
P=Path(sys.argv[1] if len(sys.argv)>1 else '/mnt/data')
prior=P/'gesui_fab_resource_consumer.py'
# Rebuild signed source and old model design without running old model loop.
exec(compile(prior.read_text().split('base=q[K].to_numpy(float);ee=q[E].to_numpy(float)')[0],str(prior),'exec'))
Ez='gesui_eiv_source.zip';z=zipfile.ZipFile(P/Ez)
assert hashlib.sha256((P/Ez).read_bytes()).hexdigest()=='1d3051d4c527421bec3ec3ba917422a9162afb81cc7f38db0d8130287174a14c'
m='EU454_PLANT_11511_SOURCE_EIV_EXACT_CROSSWALK.csv';assert hashlib.sha256(z.read(m)).hexdigest()=='ab45f36c0dad047442bb958bcef3135b45014bafec8a70811bcafb6109a02679'
eiv=pd.read_csv(z.open(m),dtype={'speciesKey':str});assert len(eiv)==11511 and eiv.speciesKey.nunique()==11511 and eiv.EIV_NUTRI_exact.notna().sum()==4011
nut=eiv.set_index('speciesKey').EIV_NUTRI_exact
s=np.array([nut.get(str(k),np.nan) for k in pk]);scored=np.isfinite(s);fam=np.array([pf[k] for k in pk]);C=(.5*(PB[a][:,L]+PB[b][:,L]))@Dmat[L];p0=PP[a];p1=PP[b];dp=p1-p0
actual=np.sum(C*dp,axis=1);np.testing.assert_allclose(actual,q.L_plant_delta,atol=1e-12,rtol=0)
meta={'signed_eiv_zip_sha256':hashlib.sha256((P/Ez).read_bytes()).hexdigest(),'source_host_keys':len(pk),'scored_host_keys':int(scored.sum()),'fab_keys':int(fab.sum()),'fab_scored_keys':int((fab&scored).sum()),'decomposition':{}}
# Fixed group definitions. Empty group's mean leverage is set to the other period's mean; zero mass contributes zero.
def partition(label,groups):
    assert np.all(np.sum(groups,axis=0)==1)
    qty=np.zeros(len(q));comp=qty.copy()
    for g in groups:
        x=p0[:,g];y=p1[:,g];c=C[:,g];m0=x.sum(1);m1=y.sum(1);t0=(x*c).sum(1);t1=(y*c).sum(1)
        u0=np.divide(t0,m0,out=np.zeros_like(t0),where=m0>1e-14);u1=np.divide(t1,m1,out=np.zeros_like(t1),where=m1>1e-14)
        u0=np.where(m0>1e-14,u0,u1);u1=np.where(m1>1e-14,u1,u0)
        qg=(m1-m0)*.5*(u0+u1);cg=(u1-u0)*.5*(m0+m1)
        np.testing.assert_allclose(qg+cg,t1-t0,atol=1e-12,rtol=0);qty+=qg;comp+=cg
    np.testing.assert_allclose(qty+comp,actual,atol=1e-12,rtol=0)
    q[label+'_quantity']=qty;q[label+'_identity']=comp
    meta['decomposition'][label]={'groups':len(groups),'max_error':float(np.max(np.abs(qty+comp-actual))),'quantity_sd':float(qty.std()),'identity_sd':float(comp.std())}
partition('all_hosts',[np.ones(len(pk),bool)])
partition('fab_nonfab',[fab,~fab])
partition('plant_family',[fam==f for f in sorted(set(fam))])
partition('family_coverage',[(fam==f)&sc for f in sorted(set(fam)) for sc in [scored,~scored] if np.any((fam==f)&sc)])
# Linear projection of the per-host matching leverage within each scored plant family; >=5 scored keys, >=2 distinct NUTRI values.
# This is an algebraic descriptive split, NOT a fitted biological effect of nitrogen.
trait=np.zeros(len(q));eligible=[];trait_detail=[]
for f in sorted(set(fam)):
    g=(fam==f)&scored
    if g.sum()<5 or len(set(s[g]))<2:continue
    eligible.append(f);x=p0[:,g];y=p1[:,g];m0=x.sum(1);m1=y.sum(1);ss=s[g]-s[g].mean();slope=C[:,g]@ss/(ss@ss)
    v0=np.divide(x@ss,m0,out=np.zeros_like(m0),where=m0>1e-14);v1=np.divide(y@ss,m1,out=np.zeros_like(m1),where=m1>1e-14)
    v0=np.where(m0>1e-14,v0,v1);v1=np.where(m1>1e-14,v1,v0)
    contribution=.5*(m0+m1)*(v1-v0)*slope;trait+=contribution
    trait_detail.append({'family':f,'scored_keys':int(g.sum()),'mean_mass':float((m0+m1).mean()/2),'sd_contribution':float(contribution.std())})
q['nutri_identity']=trait;q['other_identity']=q.family_coverage_identity-trait
np.testing.assert_allclose(q.family_coverage_quantity+q.nutri_identity+q.other_identity,actual,atol=1e-12,rtol=0)
meta['nutri_projection']={'eligible_families':eligible,'eligible_host_keys':int(np.sum(scored&np.isin(fam,eligible))),'family_details':trait_detail}
def marginal(g,label):
    x=p0[:,g];y=p1[:,g];m0=x.sum(1);m1=y.sum(1);q[label+'_mass0']=m0;q[label+'_mass1']=m1;q[label+'_mass_delta']=m1-m0
    ss=s[g];good=np.isfinite(ss);x=x[:,good];y=y[:,good];ss=ss[good];den0=x.sum(1);den1=y.sum(1)
    q[label+'_scorecov0']=np.divide(den0,m0,out=np.full_like(m0,np.nan),where=m0>1e-14);q[label+'_scorecov1']=np.divide(den1,m1,out=np.full_like(m1,np.nan),where=m1>1e-14)
    c0=np.divide(x@ss,den0,out=np.full_like(den0,np.nan),where=den0>1e-14);c1=np.divide(y@ss,den1,out=np.full_like(den1,np.nan),where=den1>1e-14)
    q[label+'_nutri0']=c0;q[label+'_nutri1']=c1;q[label+'_nutri_delta']=c1-c0
for label,g in [('host',np.ones(len(pk),bool)),('fab',fab),('nonfab',~fab)]:marginal(g,label)
rawtrait=pd.read_csv(z.open('EU454_H7PERIOD_NUTRI_TRAIT_MEANS_SOURCEONLY.csv'),dtype={'h':str});assert len(rawtrait)==451
for t0,t1 in [('Plant_N_incidence','flora_nutri'),('Plant_N_record','flora_nutri_record')]:
    g=rawtrait.pivot(index='h',columns='period',values=t0);q[t1+'0']=q.h.map(g['2017-2019']);q[t1+'1']=q.h.map(g['2020-2022']);q[t1+'_delta']=q[t1+'1']-q[t1+'0']
meta['coverage']={k:{'median_early':float(q.loc[quality,k+'_scorecov0'].median()),'median_late':float(q.loc[quality,k+'_scorecov1'].median()),'both_80_H7':int((quality&(q[k+'_scorecov0']>=.8)&(q[k+'_scorecov1']>=.8)).sum())} for k in ['host','fab','nonfab']}
base=q[K].to_numpy(float);ee=q[E].to_numpy(float);B0=np.column_stack((base,q.L_early.to_numpy()))
checks=['L_plant_delta','all_hosts_quantity','all_hosts_identity','fab_nonfab_quantity','fab_nonfab_identity','plant_family_quantity','plant_family_identity','family_coverage_quantity','family_coverage_identity','nutri_identity','other_identity']
allmask=np.ones(len(q),bool)
for tag,mask in [('quality',quality),('all',allmask)]:
    for k in checks:evaluate(tag+'__'+k,mask,q[k].to_numpy(),B0,ee)
    for k in ['host','fab','nonfab']:
        for v in ['mass','nutri']:
            m=mask if v=='mass' else mask&(q[k+'_scorecov0']>=.8)&(q[k+'_scorecov1']>=.8)
            out=evaluate(tag+'__'+k+'_'+v+'_delta',m,q[k+'_'+v+'_delta'].to_numpy(),np.column_stack((base,q[k+'_'+v+'0'].to_numpy())),ee)
            if out.get('status')!='FIT':results.append(out)
    for k in ['flora_nutri','flora_nutri_record']:
        out=evaluate(tag+'__'+k+'_delta',mask,q[k+'_delta'].to_numpy(),np.column_stack((base,q[k+'0'].to_numpy())),ee)
        if out.get('status')!='FIT':results.append(out)
for k in ['L_plant_delta','fab_nonfab_quantity','fab_nonfab_identity','nutri_identity','other_identity']:
    evaluate('quality_W__'+k,quality,q[k].to_numpy(),np.column_stack((q[[k for k in K if k!='W']].to_numpy(float),q.L_early.to_numpy())),q.W.to_numpy())
    for pred in ['mean_original_2018_logTP','mean_original_2018_EQR']:
        evaluate('FR73_'+pred+'__'+k,frmask,q[k].to_numpy(),B0,q[pred].to_numpy(),'loho')
comparisons=[]
for tag in ['quality','all']:
    for x,y in [('all_hosts_identity','all_hosts_quantity'),('fab_nonfab_identity','fab_nonfab_quantity'),('plant_family_identity','plant_family_quantity'),('nutri_identity','other_identity')]:
        g1,l1,b1=lossstore[tag+'__'+x];g2,l2,b2=lossstore[tag+'__'+y];assert g1==g2
        gain=lambda l:100*(l[:,0].mean()-l[:,1].mean())/l[:,0].mean()
        dd=b1-b2;lo,hi=np.quantile(dd,[.025,.975]);lo2,hi2=np.quantile(dd,[.00625,.99375]);comparisons.append(dict(sample=tag,A=x,B=y,diff=gain(l1)-gain(l2),lo=lo,hi=hi,bonf4_lo=lo2,bonf4_hi=hi2))
meta['model_count']=len(results);meta['max_independent_prediction_difference']=max(v['independent_prediction_diff'] for v in results if v.get('status')=='FIT');meta['prior_exact_plant_term_match']=float(np.max(np.abs(actual-q.L_plant_delta)))
pd.DataFrame(results).to_csv(P/'gesui_plant_resource_partition_results.csv',index=False);pd.DataFrame(comparisons).to_csv(P/'gesui_plant_resource_partition_contrasts.csv',index=False)
cols=['h','h6','country','L_plant_delta']+[k for k in q if k.startswith(('all_hosts_','fab_nonfab_','plant_family_','family_coverage_','nutri_identity','other_identity','host_','fab_','nonfab_','flora_nutri')) and k not in ['h','h6']];q[cols].to_csv(P/'gesui_plant_resource_partition_panel.csv',index=False)
(P/'gesui_plant_resource_partition_receipt.json').write_text(json.dumps(meta,indent=2))
print(json.dumps({k:v for k,v in meta.items() if k!='nutri_projection'},indent=2));print(pd.DataFrame(results)[['name','n','H6','gain','lo','hi','coef_min','coef_max','status']].round(3).to_string(index=False));print(pd.DataFrame(comparisons).round(3).to_string(index=False))
