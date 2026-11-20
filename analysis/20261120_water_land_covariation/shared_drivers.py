"""Order-averaged covariance allocation and held-parent-basin robustness."""
from pathlib import Path
import importlib.util
import itertools
import json
import math

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('covariation_source', HERE / 'analyse.py')
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
OUT = base.ROOT / 'results/20261120_water_land_covariation/shared_drivers'
BLOCKS = {'urban': ['urban'], 'crop': ['crop'], 'hydrology': ['area', 'riverflow', 'naturalout'], 'W': ['W']}
PAIRS = [(0, 2), (0, 3), (1, 2), (1, 3), (2, 3), (0, 1)]


def allocations(cov):
    result = []
    for k in range(4):
        contribution = np.zeros_like(cov[0])
        for mask in range(16):
            if mask & (1 << k):
                continue
            n = mask.bit_count()
            factor = math.factorial(n)*math.factorial(3-n)/math.factorial(4)
            contribution += factor*(cov[mask]-cov[mask | (1 << k)])
        result.append(contribution)
    result = np.array(result)
    np.testing.assert_allclose(result.sum(0), cov[0]-cov[15], atol=1e-10)
    return result


def main():
    OUT.mkdir(exist_ok=True)
    d, source = base.load(base.ROOT.parent)
    frame, _, outcomes, _ = base.paired(d, False)
    group, labels = pd.factorize(frame.h6, sort=True)
    G = len(labels)
    w = 1/np.bincount(group)[group]
    v = frame[outcomes].to_numpy(float)
    center = np.average(v, weights=w, axis=0)
    scale = np.sqrt(np.average((v-center)**2, weights=w, axis=0))
    v = (v-center)/scale
    columns = {mask: base.OBS + [c for k, cs in enumerate(BLOCKS.values()) if mask & (1 << k) for c in cs]
               for mask in range(16)}
    designs = {mask: base.matrix(frame, cs)[0] for mask, cs in columns.items()}

    def fit(weights):
        result = [base.model(designs[mask], v, weights) for mask in range(16)]
        if any(r is None for r in result):
            return None
        cov = np.array([r[0] for r in result])
        return cov, allocations(cov), np.array([r[1] for r in result]), result[-1][2]

    point = fit(w)
    assert point is not None
    original = pd.read_csv(OUT.parent/'correlations.csv')
    original = original[(original.cohort=='all_butterflies') & (original['mode']=='spatial')]
    for name,mask in [('observation',0),('observation_W',8),('landscape',7),('landscape_W',15)]:
        for row in original[original.stage==name].itertuples():
            np.testing.assert_allclose(point[2][mask,outcomes.index(row.left),outcomes.index(row.right)],row.r,atol=1e-12)
    rng = np.random.default_rng(20261121)
    boot = []
    for _ in range(3000):
        mult = np.bincount(rng.integers(G, size=G), minlength=G)
        value = fit(w*mult[group])
        if value is not None:
            boot.append(value)
    assert len(boot) >= 2850
    bc, ba = [np.array([b[k] for b in boot]) for k in [0, 1]]
    # Independent ordering algorithm checks allocation weights and masks.
    check = np.zeros_like(point[1])
    for order in itertools.permutations(range(4)):
        mask = 0
        for k in order:
            next_mask = mask | (1 << k)
            check[k] += (point[0][mask]-point[0][next_mask])/24
            mask = next_mask
    np.testing.assert_allclose(check, point[1], atol=1e-12)
    rows = []
    for k, block in enumerate(BLOCKS):
        for j, l in PAIRS:
            specs = [('shapley', point[1][k,j,l], ba[:,k,j,l]),
                     ('standalone', point[0][0,j,l]-point[0][1<<k,j,l], bc[:,0,j,l]-bc[:,1<<k,j,l]),
                     ('unique_conditional', point[0][15^(1<<k),j,l]-point[0][15,j,l], bc[:,15^(1<<k),j,l]-bc[:,15,j,l])]
            for kind, estimate, values in specs:
                lo, hi = np.quantile(values, [.025,.975])
                rows.append(dict(block=block,kind=kind,left=outcomes[j],right=outcomes[l],estimate=estimate,low=lo,high=hi))
    pd.DataFrame(rows).to_csv(OUT/'covariance_allocation.csv', index=False)
    rows=[]
    coefficients=np.array([b[3] for b in boot])
    for k,term in enumerate(['intercept']+columns[15]):
        if term not in base.LAND+['W']:
            continue
        for j,y in enumerate(outcomes):
            lo,hi=np.quantile(coefficients[:,k,j],[.025,.975])
            rows.append(dict(term=term,outcome=y,beta=point[3][k,j],low=lo,high=hi))
    pd.DataFrame(rows).to_csv(OUT/'full_model_slopes.csv',index=False)

    # W single-block minus W unique is its shared/redundant covariance allocation.
    rows = []
    for j,l in PAIRS:
        values = bc[:,0,j,l]-bc[:,8,j,l]-bc[:,7,j,l]+bc[:,15,j,l]
        lo,hi = np.quantile(values,[.025,.975])
        rows.append(dict(left=outcomes[j],right=outcomes[l],
                         estimate=point[0][0,j,l]-point[0][8,j,l]-point[0][7,j,l]+point[0][15,j,l],low=lo,high=hi))
    pd.DataFrame(rows).to_csv(OUT/'W_landscape_shared.csv', index=False)

    rows = []
    stages = {'raw':None, 'observation':0, 'landscape':7, 'landscape_W':15}
    for g in range(G):
        weights = w*(group != g)
        fitted = fit(weights)
        assert fitted is not None
        for stage, mask in stages.items():
            corr = base.moments(v,weights)[1] if mask is None else fitted[2][mask]
            for j,l in PAIRS:
                rows.append(dict(omitted_H6=labels[g],stage=stage,left=outcomes[j],right=outcomes[l],r=corr[j,l]))
    pd.DataFrame(rows).to_csv(OUT/'leave_one_H6_correlations.csv',index=False)

    # Training-only scaling; fixed penalty and no cross-validation tuning.
    predictions = {}
    for stage, mask in [('observation',0),('landscape',7),('landscape_W',15)]:
        pred = np.full_like(v,np.nan)
        for g in range(G):
            train, test = group != g, group == g
            assert set(group[train]).isdisjoint(group[test])
            scaler = StandardScaler().fit(frame.loc[train,columns[mask]])
            weights = w[train]/w[train].mean()
            # Original response units, so full-data outcome scaling cannot leak.
            model = Ridge(alpha=10).fit(scaler.transform(frame.loc[train,columns[mask]]),
                                        frame.loc[train,outcomes],sample_weight=weights)
            pred[test] = model.predict(scaler.transform(frame.loc[test,columns[mask]]))
        assert np.isfinite(pred).all()
        err = (pred-frame[outcomes].to_numpy())**2
        predictions[stage] = np.array([err[group==g].mean(0) for g in range(G)])
    draws = np.random.default_rng(20261122).integers(G,size=(5000,G))
    rows, errors = [], []
    for stage, err in predictions.items():
        for g in range(G):
            for j,y in enumerate(outcomes):
                errors.append(dict(stage=stage,H6=labels[g],outcome=y,mse=err[g,j]))
    for old,new in [('observation','landscape'),('landscape','landscape_W')]:
        a,b = predictions[old],predictions[new]
        values = 100*(1-b[draws].mean(1)/a[draws].mean(1))
        for j,y in enumerate(outcomes):
            lo,hi=np.quantile(values[:,j],[.025,.975])
            rows.append(dict(comparison=new+'_over_'+old,outcome=y,
                             improvement_pct=100*(1-b[:,j].mean()/a[:,j].mean()),low=lo,high=hi))
    pd.DataFrame(rows).to_csv(OUT/'heldout_prediction.csv',index=False)
    pd.DataFrame(errors).to_csv(OUT/'heldout_H6_errors.csv',index=False)

    W=frame[['W']].to_numpy()
    m0=base.model(designs[0],W,w)
    m1=base.model(designs[7],W,w)
    overlap=1-m1[0][0,0]/m0[0][0,0]
    wm=[]
    rng=np.random.default_rng(20261123)
    for _ in range(3000):
        mult=np.bincount(rng.integers(G,size=G),minlength=G)
        weights=w*mult[group]
        a,b=base.model(designs[0],W,weights),base.model(designs[7],W,weights)
        if a is not None and b is not None:
            wm.append(1-b[0][0,0]/a[0][0,0])
    corr=base.moments(frame[base.LAND+['W']].to_numpy(float),w)[1]
    pd.DataFrame(corr,index=base.LAND+['W'],columns=base.LAND+['W']).to_csv(OUT/'candidate_correlations.csv')
    receipt={'n_H7':len(frame),'n_H6':G,'bootstrap_valid':len(boot),'requested':3000,
             'blocks':BLOCKS,'source':source,'scripts':{p.name:base.sha(p) for p in [Path(__file__),HERE/'FOLLOWUP_PROTOCOL.md']},
             'W_partial_R2_landscape_after_observation':overlap,'W_partial_R2_marginal95':np.quantile(wm,[.025,.975]).tolist(),
             'checks':['Shapley efficiency every bootstrap','Independent 24-order allocation check','Original saved correlations matched',
                       'Full-rank subsets','24 disjoint heldout H6 folds','Finite predictions'],
             'limits':['Posthoc spatial follow-up; temporal uncertainty unchanged.',
                       'Shapley is an allocation convention, not causal attribution.',
                       'Marginal CIs, not corrected across candidates.',
                       'Prediction intervals resample fixed H6 errors, not full model uncertainty.',
                       'Partial R2 in sample, not predictive accuracy; climate unmeasured.']}
    (OUT/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['source','scripts','limits','checks']}))


if __name__=='__main__':main()
