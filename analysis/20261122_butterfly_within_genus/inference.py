"""Fixed exploratory conditional-pair inference; no unreported model search."""
from __future__ import annotations
from collections import defaultdict
import json
import math
import time

import numpy as np
import pandas as pd
from scipy import linalg, optimize, stats

from conditional import conditional_loglike


def fit_conditional(X,y,lengths,weights,start=None):
    p=X.shape[-1]
    total=float(np.sum(weights))
    if total<=0:
        return {'converged':False,'reason':'zero total weight','params':[0.0]*p}
    def objective(theta):
        loss,grad=conditional_loglike(theta,X,y,lengths,weights)
        return loss/total,grad/total
    initial=np.zeros(p) if start is None else np.asarray(start,float)
    result=optimize.minimize(objective,initial,jac=True,method='BFGS',
        options={'gtol':1e-7,'maxiter':600})
    loss,gradient=objective(result.x)
    maxgrad=float(np.max(np.abs(gradient)))
    # A precision-loss status with a tiny actual score is recorded, not hidden.
    converged=bool(np.isfinite(result.x).all() and np.isfinite(loss) and maxgrad<=1e-6)
    return {'converged':converged,'params':result.x.tolist(),
        'negative_loglik':float(loss*total),'optimizer_status':int(result.status),
        'optimizer_message':str(result.message),'iterations':int(result.nit),
        'max_normalized_gradient':maxgrad}


def observed_information(theta,X,y,lengths,weights,epsilon=1e-5):
    theta=np.asarray(theta,float)
    out=np.zeros((len(theta),len(theta)))
    for j in range(len(theta)):
        step=epsilon*max(1.,abs(theta[j]))
        a=theta.copy();b=theta.copy();a[j]+=step;b[j]-=step
        ga=conditional_loglike(a,X,y,lengths,weights)[1]
        gb=conditional_loglike(b,X,y,lengths,weights)[1]
        out[:,j]=(ga-gb)/(2*step)
    return (out+out.T)/2


def cluster_covariance(information,scores,clusters):
    grouped=defaultdict(lambda:np.zeros(scores.shape[1]))
    for c,s in zip(clusters,scores):
        grouped[str(c)]+=s
    if len(grouped)<2:
        raise ValueError('cluster covariance requires at least2 H6')
    meat=sum(np.outer(s,s) for s in grouped.values())
    inv=np.linalg.inv(information)
    return (len(grouped)/(len(grouped)-1))*inv@meat@inv


def permute_trait_labels(original,genus,rng):
    grouped=defaultdict(list)
    for species in sorted(original):
        grouped[genus[species]].append(species)
    out={}
    for g in sorted(grouped):
        members=grouped[g]
        values=rng.permutation([original[s] for s in members])
        out.update({s:float(v) for s,v in zip(members,values)})
    return out


def prepare_design(strata,cells,panel,controls):
    names=['plant_n_incidence']+list(controls)
    scaling={names[0]:[float(cells[names[0]].mean()),float(cells[names[0]].std(ddof=0))]}
    for k in controls:
        scaling[k]=[float(panel[k].mean()),float(panel[k].std(ddof=0))]
    if any(v[1]<=1e-12 for v in scaling.values()):
        raise RuntimeError('constant original source covariate')
    n=len(strata);mx=max(len(r['h']) for r in strata)
    base=np.zeros((n,mx,len(names)))
    y=np.zeros((n,mx),dtype=int);lengths=np.zeros(n,dtype=int)
    for i,r in enumerate(strata):
        count=len(r['h']);lengths[i]=count;y[i,:count]=r['y']
        for j,k in enumerate(names):
            raw=r[k] if j==0 else r['controls'][k]
            mu,sd=scaling[k]
            values=(np.asarray(raw,float)-mu)/sd
            base[i,:count,j]=values-values.mean()
    difference=np.array([r['trait_difference'] for r in strata])
    weights=np.array([r['weight'] for r in strata])
    z=(base*difference[:,None,None]*np.sqrt(weights)[:,None,None]).reshape(-1,len(names))
    # Only exactly redundant nuisance directions may be removed; protect X.
    nuisance=z[:,1:]
    if nuisance.shape[1]:
        _,R,pivot=linalg.qr(nuisance,mode='economic',pivoting=True)
        diag=np.abs(np.diag(R));tol=max(nuisance.shape)*np.finfo(float).eps*(diag.max() if len(diag) else 1.)
        rank=int(np.sum(diag>tol))
        keep=sorted((np.asarray(pivot[:rank])+1).tolist())
    else:
        keep=[]
    residual=z[:,0].copy()
    if keep:
        residual-=z[:,keep]@np.linalg.lstsq(z[:,keep],z[:,0],rcond=None)[0]
    residual_norm=float(np.linalg.norm(residual))
    if residual_norm<=1e-10*max(1.,float(np.linalg.norm(z[:,0]))):
        raise RuntimeError('identification STOP: focal plant N adds no rank beyond nuisance')
    chosen=[0]+keep
    return dict(base=base[:,:,chosen],y=y,lengths=lengths,difference=difference,weights=weights,
        clusters=np.array([r['h6'] for r in strata]),genera=np.array([r['genus'] for r in strata]),
        species_a=[r['species_a'] for r in strata],species_b=[r['species_b'] for r in strata],
        names=[names[j] for j in chosen],scaling=scaling,
        dropped_nuisance=[k for j,k in enumerate(names) if j not in chosen],
        design_rank=len(chosen),focal_residual_norm=residual_norm)


def model_arrays(data,difference=None,weights=None):
    difference=data['difference'] if difference is None else np.asarray(difference,float)
    return data['base']*difference[:,None,None],data['y'],data['lengths'],data['weights'] if weights is None else weights


def report_fit(data,fit=None):
    X,y,lengths,weights=model_arrays(data)
    fit=fit_conditional(X,y,lengths,weights) if fit is None else fit
    if not fit['converged']:
        return dict(status='NUMERICAL_STOP',fit=fit)
    theta=np.asarray(fit['params'])
    info=observed_information(theta,X,y,lengths,weights)
    info2=observed_information(theta,X,y,lengths,weights,epsilon=2e-5)
    eig=np.linalg.eigvalsh(info)
    if eig.min()<=0 or not np.isfinite(info).all():
        return dict(status='INFORMATION_STOP',fit=fit,information_eigenvalues=eig.tolist())
    _,_,scores=conditional_loglike(theta,X,y,lengths,weights,return_scores=True)
    # Zero-HNS-difference strata contribute zero score and no information.
    active=np.abs(data['difference'])>1e-12
    cov=cluster_covariance(info,scores[active],data['clusters'][active])
    groups=len(np.unique(data['clusters'][active]))
    se=math.sqrt(max(0.,float(cov[0,0])))
    q=float(stats.t.ppf(.975,groups-1));beta=float(theta[0]);ci=[beta-q*se,beta+q*se]
    return dict(status='FIT',fit=fit,beta=beta,odds_ratio=math.exp(beta),
        H6_cluster_standard_error=se,H6_clusters=groups,cluster_df=groups-1,
        beta95_H6_cluster_t=ci,odds_ratio95_H6_cluster_t=[math.exp(x) for x in ci],
        nominal_two_sided_H6_cluster_t_p=float(2*stats.t.sf(abs(beta/se),groups-1)) if se else None,
        coefficients=dict(zip(data['names'],fit['params'])),
        information_condition_number=float(eig.max()/eig.min()),
        information_eigenvalues=eig.tolist(),
        information_epsilon_check_max_relative=float(np.max(np.abs(info-info2))/max(1.,np.max(np.abs(info)))),
        nuisance_removed_only_if_exactly_redundant=data['dropped_nuisance'],
        source_scaling=data['scaling'],fixed_observed_species_inference=True,
        not_true_occupancy_abundance_or_wastewater_causality=True)


def _save(path,value):
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n')


def run_all(cells,counts,scores,panel,strata,gate,output,controls,build_pair_strata,source_summary):
    if not gate['primary']['support_gate_pass']:
        _save(output/'results.json',dict(status='SOURCE_STOP',gate=gate));return
    data=prepare_design(strata,cells,panel,controls)
    primary=report_fit(data)
    _save(output/'primary.json',primary)
    print('PRIMARY='+json.dumps(primary,ensure_ascii=False,allow_nan=False),flush=True)
    if primary['status']!='FIT':
        _save(output/'results.json',dict(status=primary['status'],primary=primary));return
    started=time.monotonic()
    X,y,lengths,weights=model_arrays(data)
    start=primary['fit']['params']
    regions=sorted(np.unique(data['clusters']))
    group_index={g:i for i,g in enumerate(regions)}
    index=np.array([group_index[g] for g in data['clusters']])
    bootstrap=[];rng=np.random.default_rng(20261122)
    for i in range(399):
        multiplicity=rng.multinomial(len(regions),np.repeat(1/len(regions),len(regions)))
        w=weights*multiplicity[index]
        f=fit_conditional(X,y,lengths,w,start=start)
        bootstrap.append(dict(replicate=i+1,converged=f['converged'],beta=f['params'][0] if f['converged'] else None,
                              iterations=f.get('iterations'),max_gradient=f.get('max_normalized_gradient'),status=f.get('optimizer_status')))
        if (i+1)%50==0:
            print('BOOTSTRAP '+str(i+1)+'/399 elapsed='+str(round(time.monotonic()-started,1)),flush=True)
    pd.DataFrame(bootstrap).to_csv(output/'bootstrap_H6.csv',index=False)
    good=[r['beta'] for r in bootstrap if r['converged']]
    bsummary=dict(draws=399,seed=20261122,converged=len(good),failed=399-len(good),
                  fixed_source_pair_weights=True,model_refitted_per_draw=True)
    if len(good)>=math.ceil(.95*399):
        bsummary.update(beta95_percentile=np.quantile(good,[.025,.975]).tolist(),
                        odds_ratio95_percentile=np.exp(np.quantile(good,[.025,.975])).tolist())
    else:
        bsummary['status']='BOOTSTRAP_UNRESOLVED'
    _save(output/'bootstrap_summary.json',bsummary)
    print('BOOTSTRAP_SUMMARY='+json.dumps(bsummary),flush=True)

    used=sorted(set(data['species_a'])|set(data['species_b']))
    score=scores.set_index('butterflyKey')
    original={k:float(score.loc[k,'host_nutrient_mean']) for k in used}
    genera={k:str(score.loc[k,'genus']) for k in used}
    references=[];rng=np.random.default_rng(20261123)
    for i in range(999):
        shuffled=permute_trait_labels(original,genera,rng)
        d=np.array([shuffled[a]-shuffled[b] for a,b in zip(data['species_a'],data['species_b'])])
        permX=data['base']*d[:,None,None]
        f=fit_conditional(permX,y,lengths,weights)
        references.append(dict(replicate=i+1,converged=f['converged'],beta=f['params'][0] if f['converged'] else None,
                               iterations=f.get('iterations'),max_gradient=f.get('max_normalized_gradient'),status=f.get('optimizer_status')))
        if (i+1)%100==0:
            print('TRAIT_REFERENCE '+str(i+1)+'/999 elapsed='+str(round(time.monotonic()-started,1)),flush=True)
    pd.DataFrame(references).to_csv(output/'within_butterfly_genus_label_reference.csv',index=False)
    goodref=[r['beta'] for r in references if r['converged']]
    rsummary=dict(draws=999,seed=20261123,converged=len(goodref),failed=999-len(goodref),
        descriptive_full_within_genus_exchangeability_not_exact_beta_only_null=True,
        labels_global_once_per_species=True,eligible_strata_and_tied_pairs_fixed=True,
        participating_species=len(used),participating_genera=len(set(genera.values())))
    if len(goodref)==999:
        extreme=int(np.sum(np.abs(goodref)>=abs(primary['beta'])-1e-12))
        rsummary.update(two_sided_absolute_beta_reference_p=(1+extreme)/1000,extreme_draws=extreme,
                        reference_beta_quantiles=np.quantile(goodref,[.025,.5,.975]).tolist())
    else:
        rsummary['status']='REFERENCE_INCOMPLETE_NO_FULL_REFERENCE_P'
    _save(output/'trait_reference_summary.json',rsummary)
    print('TRAIT_REFERENCE_SUMMARY='+json.dumps(rsummary),flush=True)

    influence=[]
    for genus in sorted(np.unique(data['genera'])):
        w=weights*(data['genera']!=genus)
        f=fit_conditional(X,y,lengths,w,start=start)
        influence.append(dict(removed_butterfly_genus=genus,converged=f['converged'],
                              beta=f['params'][0] if f['converged'] else None))
    pd.DataFrame(influence).to_csv(output/'leave_one_butterfly_genus_out.csv',index=False)
    valid=[r['beta'] for r in influence if r['converged']]
    isummary=dict(genera=len(influence),converged=len(valid),failed=len(influence)-len(valid),
                  positive=sum(b>0 for b in valid),negative=sum(b<0 for b in valid),
                  beta_range=[min(valid),max(valid)] if valid else None)
    remaining=[k for k in controls if k not in ['nButterfly','hns_coverage','hns_host_completeness','hns_scored_species']]
    ds=prepare_design(strata,cells,panel,remaining)
    no_butterfly_controls=report_fit(ds)
    high=scores[scores.trait_coverage>=.8]
    hs=build_pair_strata(cells,counts,high)
    hg=source_summary(hs)
    high_result=dict(status='SOURCE_STOP',source=hg)
    if hg['support_gate_pass']:
        high_result=dict(source=hg,**report_fit(prepare_design(hs,cells,panel,controls)))
    summary=dict(status='COMPLETE_EXPLORATORY',specification_commit=gate['specification_commit'],
        source=gate,primary=primary,H6_bootstrap=bsummary,within_butterfly_genus_trait_reference=rsummary,
        leave_one_butterfly_genus_out=isummary,
        sensitivity_without_four_butterfly_derived_controls=no_butterfly_controls,
        sensitivity_host_coverage80=high_result,
        total_inference_seconds=float(time.monotonic()-started))
    _save(output/'results.json',summary)
    print('FINAL_SUMMARY='+json.dumps({k:v for k,v in summary.items() if k not in ['source','primary']},ensure_ascii=False),flush=True)
    return summary
