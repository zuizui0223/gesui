#!/usr/bin/env python3
"""Exploratory independent 93-site within-butterfly-genus relative abundance test.

Requires source-locked Guariento (2023) Springer ESM3 and original gesui HNS.
Pooled site-by-genus conditional multinomial composition, species fixed effects;
relative abundance association, not sewage impacts, temporal tracking, or causal N.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,sys,time
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
from scipy.stats import norm,t
from tyrol93_source_gate import (fetch_xlsx,table,xlsx_sheets,normal,count_source,
                                  EXPECTED_HNS_SHA256,EXPECTED_XLSX_SHA256)


def get_data(blob,hns_raw,source_summary):
    sheets=xlsx_sheets(blob)
    _,all_sites=table(sheets['Community']); _,vrows=table(sheets['Variables'])
    variable={normal(r['site_code']):r for r in vrows}
    hs=list(csv.DictReader(io.StringIO(hns_raw.decode('utf-8-sig'))))
    scores={normal(r['accepted_name']):float(r['host_nutrient_mean']) for r in hs
       if r['taxon_status']=='ACCEPTED' and r.get('actual_used_alias_conflict','').lower() not in ('true','1','yes')
       and r['host_nutrient_mean']}
    support=source_summary['genus_support']
    genus_species={r['genus']:r['candidate_names'] for r in support}
    species_counts={}
    for site in all_sites:
        k=normal(site['site_code'])
        species_counts[k]={normal(s):float(count) for s,count in site.items()
                           if s not in ('site_code','Land use type')}
    siteids=sorted(variable,key=lambda x:int(x))
    df=[]
    for site in siteids:
        r=variable[site]
        df.append({'site':site,'N':float(r['Nutrients (Elenberg plant indicator)']),
          'light':float(r['Light (Elenberg plant indicator)']),
          'temp':float(r['T_Mean (site mean annual temperature °C)']),
          'elevation':float(r['Elevation (m a.s.l.)']),
          'landuse':normal(r['Land use type']),'year':normal(r['Year of sampling'])})
    for row in df:
        for k in ('N','light','temp','elevation'):
            if not math.isfinite(row[k]):raise ValueError('MISSING_COVARIATE '+k)
    landuses=sorted(set(r['landuse'] for r in df));years=sorted(set(r['year'] for r in df))
    scales={k:(float(np.mean([r[k] for r in df])),float(np.std([r[k] for r in df],ddof=0)))
      for k in ('N','light','temp','elevation')}
    if any(sd<1e-9 for mu,sd in scales.values()):raise ValueError('CONSTANT_NATIVE_COVARIATE')
    X=np.array([[((r[k]-scales[k][0])/scales[k][1]) for k in ('N','light','temp','elevation')]
          +[float(r['landuse']==l) for l in landuses[1:]]
          +[float(r['year']==y) for y in years[1:]] for r in df])
    names=['N','light','temp','elevation']+['landuse:'+l for l in landuses[1:]]+['year:'+y for y in years[1:]]
    obs=[];total=0;unique_species=set()
    for i,k in enumerate(siteids):
        for gen,sp in sorted(genus_species.items()):
            counts=np.array([species_counts[k][s] for s in sp],float)
            n=counts.sum()
            if n<=0:continue
            obs.append({'site_index':i,'site_id':k,'genus':gen,'species':sp,'y':counts,
                         'hns':np.array([scores[s] for s in sp],float),'n':n})
            total+=n
            unique_species.update(s for s,c in zip(sp,counts) if c>0)
    return {'sites':df,'site_ids':siteids,'X':X,'xnames':names,'scales':scales,
            'landuses':landuses,'years':years,'genus_species':genus_species,'obs':obs,
            'total_eligible_count':total,'eligible_species_recorded':len(unique_species)}


def design(data, names, site_keep=None, excluded_genus=None):
    site_keep=set(range(len(data['site_ids']))) if site_keep is None else set(site_keep)
    observations=[g for g in data['obs'] if g['site_index'] in site_keep and g['genus']!=excluded_genus]
    valid_genera=sorted(g for g in data['genus_species'] if g!=excluded_genus)
    all_sp={g:data['genus_species'][g] for g in valid_genera}
    alphas={};offset=0
    for gen in valid_genera:
        for name in all_sp[gen][1:]:
            alphas[(gen,name)]=offset;offset+=1
    cov_indices=[data['xnames'].index(name) for name in names]
    npar=offset+len(names)
    constructed=[]
    for g in observations:
        hns_c=g['hns']-g['hns'].mean()
        Z=np.zeros((len(g['species']),npar),dtype=float)
        for j,sp in enumerate(g['species']):
            if (g['genus'],sp) in alphas:Z[j,alphas[(g['genus'],sp)]]=1.
        Z[:,offset:]=np.outer(hns_c,data['X'][g['site_index'],cov_indices])
        constructed.append({'Z':Z,'y':g['y'],'n':g['n'],'site_index':g['site_index'],
                            'genus':g['genus'],'species':g['species']})
    return {'groups':constructed,'p':npar,'names':names,'nalpha':offset,
            'sites':len(set(g['site_index'] for g in constructed)),
            'genera':valid_genera,'alphas':alphas}


def evaluate(theta,groups,weights=None,return_hessian=False,return_site_scores=False):
    npar=len(theta);grad=np.zeros(npar);nll=0.;H=np.zeros((npar,npar)) if return_hessian else None
    scores=defaultdict(lambda:np.zeros(npar)) if return_site_scores else None
    for j,g in enumerate(groups):
        Z=g['Z'];y=g['y'];n=g['n']; w=1. if weights is None else weights[j]
        if w<=0:continue
        logits=Z@theta; lse=logsumexp(logits);probs=np.exp(logits-lse)
        mu=probs@Z
        nll+=w*(n*lse-np.dot(y,logits))
        contrib=w*(n*mu-y@Z)
        grad+=contrib
        if return_site_scores:scores[g['site_index']]-=contrib
        if return_hessian:
            H+=w*n*((Z.T*probs)@Z-np.outer(mu,mu))
    if return_hessian and return_site_scores:return nll,grad,H,scores
    if return_hessian:return nll,grad,H
    if return_site_scores:return nll,grad,scores
    return nll,grad


def fit(d,weights=None,initial=None):
    p=d['p']; init=np.zeros(p) if initial is None else np.asarray(initial,float)
    def fn(theta):return evaluate(theta,d['groups'],weights=weights)
    res=minimize(fn,init,jac=True,method='BFGS',options={'maxiter':700,'gtol':1e-6})
    loss,g,H,sc=evaluate(res.x,d['groups'],weights=weights,return_hessian=True,return_site_scores=True)
    eig=np.linalg.eigvalsh(H)
    site_scores=np.stack(list(sc.values()))
    nsites=len(sc)
    converged=(np.isfinite(res.x).all() and eig.min()>1e-8 and
      np.max(np.abs(g))/max(1.,sum(x['n'] for x in d['groups']))<1e-7)
    output={'status':'FIT' if converged else 'NUMERICAL_STOP','success_optimizer':bool(res.success),
      'optimizer_message':str(res.message),'maximum_grad_per_count':float(max(abs(g))/max(1,sum(x['n'] for x in d['groups']))),
      'npar':p,'nll':float(loss),'eigen_min':float(eig.min()),'eigen_max':float(eig.max()),
      'information_condition':float(eig.max()/eig.min()) if eig.min()>0 else None,
      'sampled_sites':d['sites'],'sampled_genera':len(d['genera']),
      'params':res.x.tolist(),'site_count_for_sandwich':nsites}
    if converged and nsites>1:
        inv=np.linalg.inv(H)
        robust=(nsites/(nsites-1))*inv@(site_scores.T@site_scores)@inv
        output['robust_se']=[math.sqrt(max(0.,robust[i,i])) for i in range(p)]
        output['H_site_cluster_score_sum']=[float(v) for v in site_scores.sum(axis=0)]
    return output


def model_result(data,names,site_keep=None,excluded_genus=None):
    d=design(data,names,site_keep,excluded_genus)
    result=fit(d)
    result['features']=names
    if 'N' in names and result['status']=='FIT':
        i=d['nalpha']+names.index('N')
        b=float(result['params'][i]);se=float(result['robust_se'][i]);q=float(t.ppf(.975,result['site_count_for_sandwich']-1))
        result.update(beta_N_per_1SDplantN_1HNSpoint=b,odds_ratio=math.exp(b),
          cluster_interval_beta=[b-q*se,b+q*se],
          cluster_interval_OR=[math.exp(b-q*se),math.exp(b+q*se)],
          nominal_site_cluster_t_p=float(2*t.sf(abs(b/se),result['site_count_for_sandwich']-1)) if se>0 else None,
          N_coefficient_index=i)
    return result


def heldout_score(model,d,holdouts,training):
    # Pooled site groups; parameterization uses the same globally frozen genus pool
    # and fixed means/SD from all 93 source sites. Check all species supported.
    names=model['features'];train=design(d,names,site_keep=training)
    fitobj=fit(train)
    if fitobj['status']!='FIT':return None,'TRAIN_NUMERIC_STOP'
    if train['p']!=model['npar']:
        return None,'SPECIES_NOT_SUPPORTED_IN_TRAIN'
    test=design(d,names,site_keep=holdouts)
    if test['p']!=train['p']:
        return None,'HOLDOUT_INVALID_MODEL_DIMENSIONS'
    nll=evaluate(np.array(fitobj['params']),test['groups'])[0]
    return -nll,None


def main():
    p=argparse.ArgumentParser();p.add_argument('--hns',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.parent.mkdir(parents=True,exist_ok=True)
    report={'status':'SOURCE_STOP'}
    try:
        blob,url,_=fetch_xlsx();hns=a.hns.read_bytes()
        if hashlib.sha256(hns).hexdigest()!=EXPECTED_HNS_SHA256:raise RuntimeError('HNS_SHA_MISMATCH')
        s=count_source(blob,hns)
        if not s['gate']['eligible_to_proceed_exploratory_model']:
            report={'status':'SOURCE_STOP_BY_FROZEN_GATE','source':s}
        else:
            d=get_data(blob,hns,s)
            # fixed rank evaluation for nutrient beyond confounders in native sites
            z=d['X']; ix=d['xnames'].index('N')
            nc=np.delete(z,ix,axis=1)
            n_res=z[:,ix]-np.c_[np.ones(len(z)),nc]@np.linalg.lstsq(np.c_[np.ones(len(z)),nc],z[:,ix],rcond=None)[0]
            resid_fraction=float(np.mean(n_res**2)/np.mean((z[:,ix]-z[:,ix].mean())**2))
            base=['light','temp','elevation']+[x for x in d['xnames'] if x.startswith('landuse:')]
            primary=['N']+base
            year=['N']+base+[x for x in d['xnames'] if x.startswith('year:')]
            if resid_fraction<0.05:report={'status':'IDENTIFICATION_STOP_NUTRIENT_TOO_COLLINEAR',
              'source':s,'N_variance_fraction_after_landuse_light_temp_elevation_year':resid_fraction}
            else:
                specs={'unadjusted_N':['N'],'adjusted_N':primary,'adjusted_without_N':base,
                   'no_landuse_N':['N','light','temp','elevation'],'with_year_N':year}
                fits={key:model_result(d,names) for key,names in specs.items()}
                print('MODELS '+json.dumps({k:{'status':v['status'],'OR':v.get('odds_ratio'),
                      'interval':v.get('cluster_interval_OR')} for k,v in fits.items()}),flush=True)
                cv=[]; seeds=[20261122,20261123,20261124]
                for seed in seeds:
                    rng=np.random.default_rng(seed)
                    landuse=np.array([r['landuse'] for r in d['sites']]);folds=np.zeros(len(landuse),dtype=int)
                    for g in d['landuses']:
                        sites=np.flatnonzero(landuse==g);rng.shuffle(sites)
                        folds[sites]=np.arange(len(sites))%5
                    for fold in range(5):
                        test=np.flatnonzero(folds==fold);train=np.flatnonzero(folds!=fold)
                        values={};errors={}
                        for label in ('adjusted_without_N','adjusted_N'):
                            ll,err=heldout_score(fits[label],d,test,train)
                            values[label]=ll;errors[label]=err
                        gain=(values['adjusted_N']-values['adjusted_without_N']) if all(v is not None for v in values.values()) else None
                        cv.append({'seed':seed,'fold':fold,'holdout_sites':len(test),
                            'heldout_N_added_loglik_gain':gain,'status':'FIT' if gain is not None else str(errors)})
                outcomes=[x['heldout_N_added_loglik_gain'] for x in cv if x['heldout_N_added_loglik_gain'] is not None]
                # Posthoc genus influence diagnostic; no independent genus-specific hypothesis testing.
                influence=[]
                for excluded in sorted(d['genus_species']):
                    item={'excluded_genus':excluded}
                    for label,names in [('adjusted_N',primary),('without_landuse_N',['N','light','temp','elevation'])]:
                        rr=model_result(d,names,excluded_genus=excluded)
                        item[label]={'status':rr['status'],'odds_ratio':rr.get('odds_ratio'),
                          'beta':rr.get('beta_N_per_1SDplantN_1HNSpoint'),
                          'cluster_interval_OR':rr.get('cluster_interval_OR'),
                          'p_descriptive':rr.get('nominal_site_cluster_t_p')}
                    influence.append(item)
                report={'status':'COMPLETE_EXPLORATORY_INDEPENDENT_SITE_NUTRIENT_ASSOCIATION',
                    'source':{k:v for k,v in s.items() if k!='genus_support'},
                    'source_supported_genera':s['genus_support'],
                    'design':{'site_n':len(d['sites']),'site_genera_combinations_observed':len(d['obs']),
                        'total_individual_count_in_10_eligible_genera':d['total_eligible_count'],
                        'landuses':d['landuses'],'years':d['years'],'trait_scaled_per_1_HNS_point':True,
                        'plantN_scaled_per_1SD':d['scales']['N'][1],
                        'N_fraction_residual_var_after_controls_and_year':resid_fraction,
                        'site_year_is_cross_sectional':True,
                        'fixed_source_species_pool_used_at_all_sites':True},
                    'fits':fits,'posthoc_leave_one_genus_out':influence,'heldout_site_stratified_landuse_fivefold':cv,
                    'cv_summary':{'total_folds':len(cv),'valid_folds':len(outcomes),
                       'positive_folds':sum(v>0 for v in outcomes),
                       'total_N_added_loglik_gain':float(sum(outcomes)) if outcomes else None,
                       'failure_folds':len(cv)-len(outcomes)},
                    'no_causal_sewage_or_population_dynamics_claim':True}
    except Exception as e:report={'status':'SOURCE_OR_NUMERIC_STOP','exception':str(e)[:700]}
    a.output.write_text(json.dumps(report,indent=2,ensure_ascii=False,
        default=lambda obj: obj.item() if isinstance(obj,np.generic) else (_ for _ in ()).throw(TypeError(type(obj).__name__)))+'\n')
    print('FINAL '+json.dumps({'status':report['status'],'cv':report.get('cv_summary'),
      'primary':{k:report.get('fits',{}).get('adjusted_N',{}).get(k) for k in ('status','beta_N_per_1SDplantN_1HNSpoint','odds_ratio','cluster_interval_OR')}},ensure_ascii=False,
      default=lambda obj: obj.item() if isinstance(obj,np.generic) else (_ for _ in ()).throw(TypeError(type(obj).__name__))),flush=True)
    if 'STOP' in report['status'] or (report.get('fits',{}).get('adjusted_N',{}).get('status')=='NUMERICAL_STOP'):sys.exit(4)

if __name__=='__main__':main()
