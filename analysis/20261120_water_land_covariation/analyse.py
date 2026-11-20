"""Matched catchment covariation; no causal pathway or annual synchrony inference."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'analysis/20261115_hns_family_taxon'))
from analyse import load_reconstruction, FAMS

LAND = ['crop', 'urban', 'area', 'riverflow', 'naturalout']
OBS = ['nPlant', 'nButterfly', 'NUTRI_coverage_all_Plantae_speciesKeys',
       'primary_record_coverage_raw', 'water_support']
Y = ['aquatic_stress', 'TP', 'plant', 'butterfly', 'lyca']
STAGES = ['raw', 'observation', 'observation_W', 'landscape', 'landscape_W']


def sha(p):
    data = p.read_bytes()
    if p.suffix in {'.py', '.md'}:
        data = data.replace(b'\r\n', b'\n')
    return hashlib.sha256(data).hexdigest()


def load(work):
    panel = work / 'original/H7_454_HNS_panel.csv'
    scores = work / 'original/GESUI_303_HOST_NUTRIENT_SCORES.csv'
    d, taxa, B, score, known, fam, raw, _, provenance = load_reconstruction(
        panel, scores, work / 'facets')
    m = known & (fam == 'Lycaenidae')
    counts = B[:, m]
    incidence = counts > 0
    ok = ((counts.sum(1) >= 10) & (incidence.sum(1) >= 2)
          & (counts.sum(1) >= .7 * raw[:, FAMS.index('Lycaenidae')]))
    d['lyca'] = np.where(ok, incidence @ score[m] / np.maximum(incidence.sum(1), 1), np.nan)
    p = next((work / 'fr73').rglob('FR73*.csv'))
    assert sha(p) == 'f2ff4e56f245e1aa8f84f5ac099d68ec1c7e27a786c1067439db388b186dd90b'
    fr = pd.read_csv(p, dtype={'HYBAS_ID_L07': str, 'HYBAS_ID_L06': str})
    d = d.merge(fr, left_on=['h', 'h6'], right_on=['HYBAS_ID_L07', 'HYBAS_ID_L06'], validate='many_to_one')
    inputs = [panel, scores, p]
    for pat in ['*ELIGIBILITY*.csv', '*MEANS*.csv']:
        p = next((work / 'plant').glob(pat))
        d = d.merge(pd.read_csv(p, dtype={'h': str}), on=['h', 'period'],
                    suffixes=('', '_source'), validate='one_to_one')
        inputs.append(p)
    assert len(d) == 146 and d.full_source_trait_eligible.all()
    assert not d.duplicated(['h', 'period']).any()
    assert d.groupby('h').h6.nunique().max() == 1
    late = d.period == '2020-2022'
    d['TP'] = np.where(late, d.mean_original_2022_logTP, d.mean_original_2018_logTP)
    d['aquatic_stress'] = -np.where(late, d.mean_original_2022_EQR, d.mean_original_2018_EQR)
    d['plant'] = d.Plant_N_incidence
    d['butterfly'] = d.primary_incidence_mean
    d['water_support'] = d.log1p_original_same_site_2018_22_pair_n
    assert d.lyca.notna().sum() == 139
    inputs.extend([Path(__file__), Path(__file__).with_name('PROTOCOL.md')])
    inputs.extend(sorted((ROOT / 'analysis/20261115_hns_family_taxon').glob('*.py')))
    inputs.extend([ROOT / 'analysis/20261115_hns_family_decomposition/analyse.py',
                   ROOT / 'analysis/20261115_hns/scripts/score_core.py'])
    provenance['inputs'] = [{'path': str(p.relative_to(work)) if p.is_relative_to(work) else str(p),
                              'sha256': sha(p),
                              'hash_mode': 'LF-normalized UTF-8' if p.suffix in {'.py', '.md'} else 'original bytes'} for p in inputs]
    return d, provenance


def paired(d, lyca):
    early = d[d.period == '2017-2019'].set_index('h').sort_index()
    late = d[d.period == '2020-2022'].set_index('h').sort_index()
    assert early.index.equals(late.index)
    if lyca:
        keep = early.lyca.notna() & late.lyca.notna()
        early, late = early.loc[keep], late.loc[keep]
    outcomes = Y if lyca else Y[:-1]
    cols = list(dict.fromkeys(outcomes + OBS + LAND + ['W']))
    mean = (early[cols] + late[cols]) / 2
    delta = late[cols] - early[cols]
    static = {c: float((late[c] - early[c]).abs().max()) for c in LAND + ['W']}
    assert static['W'] < 1e-12
    for c in LAND + ['W']:
        delta[c] = mean[c]
    delta['water_support'] = early.mean_original_delta_log1p_TP_positive_sample_opportunities
    np.testing.assert_allclose(delta.TP, early.mean_site_delta_log_original_measured_TP_2022_minus_2018, atol=1e-12)
    np.testing.assert_allclose(delta.aquatic_stress, -early.mean_site_delta_original_EQR_2022_minus_2018, atol=1e-12)
    for frame in (mean, delta):
        frame['h6'] = early.h6
        assert np.isfinite(frame[cols].to_numpy()).all()
    return mean, delta, outcomes, static


def matrix(frame, columns):
    active = [c for c in columns if frame[c].std(ddof=0) > 1e-10]
    v = frame[active].to_numpy(float)
    if active:
        v = (v - v.mean(0)) / v.std(0)
    return np.column_stack([np.ones(len(frame)), v]), active


def moments(v, w):
    center = v - np.average(v, weights=w, axis=0)
    cov = (center * w[:, None]).T @ center / w.sum()
    corr = cov / np.sqrt(np.outer(np.diag(cov), np.diag(cov)))
    return cov, corr


def model(X, v, w):
    sw = np.sqrt(w)
    b, _, rank, _ = np.linalg.lstsq(X * sw[:, None], v * sw[:, None], rcond=None)
    if rank != X.shape[1]:
        return None
    resid = v - X @ b
    cov, corr = moments(resid, w)
    return cov, corr, b


def run(frame, outcomes, draws, seed):
    groups, labels = pd.factorize(frame.h6, sort=True)
    sizes = np.bincount(groups)
    G = len(labels)
    w = 1 / sizes[groups]
    v = frame[outcomes].to_numpy(float)
    center = np.average(v, weights=w, axis=0)
    scale = np.sqrt(np.average((v-center)**2, weights=w, axis=0))
    v = (v-center) / scale
    specs = [[], OBS, OBS + ['W'], OBS + LAND, OBS + LAND + ['W']]
    designs = [matrix(frame, cols) for cols in specs]

    def fit(w):
        fitted = [model(X, v, w) for X, active in designs]
        if any(m is None for m in fitted):
            return None
        cov = np.array([m[0] for m in fitted])
        corr = np.array([m[1] for m in fitted])
        beta = np.array([fitted[k][2][-1] for k in [2, 4]])
        # Frisch-Waugh-Lovell identity: covariance removed by W equals
        # beta_j * beta_k * conditional variance(W) on these same weights.
        for before, after in [(1, 2), (3, 4)]:
            W = designs[after][0][:, -1:]
            wr = model(designs[before][0], W, w)
            explained = np.outer(fitted[after][2][-1], fitted[after][2][-1]) * wr[0][0, 0]
            np.testing.assert_allclose(cov[before] - cov[after], explained, atol=1e-9)
        return cov, corr, beta

    point = fit(w)
    assert point is not None
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(draws):
        mult = np.bincount(rng.integers(G, size=G), minlength=G)
        result = fit(w * mult[groups])
        if result is not None:
            boots.append(result)
    assert len(boots) >= .95 * draws, (len(boots), draws)
    boot = [np.array([b[i] for b in boots]) for i in range(3)]
    meta = {'n_H7': len(frame), 'n_H6': G, 'valid_bootstraps': len(boots),
            'requested_bootstraps': draws, 'controls': dict(zip(STAGES, [x[1] for x in designs])),
            'outcome_scale': dict(zip(outcomes, scale)), 'H6_sizes': sizes.tolist()}
    return point, boot, meta


def synthetic_check():
    # Exact orthogonal noise makes the known shared-driver covariance recoverable.
    q, _ = np.linalg.qr(np.column_stack([np.ones(200), np.random.default_rng(7).normal(size=(200, 3))]))
    W, a, b = q[:, 1], q[:, 2], q[:, 3]
    v = np.column_stack([2*W+a, 3*W+b])
    w = np.ones(200)
    raw = model(np.ones((200, 1)), v, w)
    adj = model(np.column_stack([np.ones(200), W]), v, w)
    assert raw[1][0, 1] > .8 and abs(adj[1][0, 1]) < 1e-12
    np.testing.assert_allclose(raw[0][0, 1] - adj[0][0, 1], 6*np.var(W), atol=1e-12)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--work', type=Path, default=ROOT.parent)
    p.add_argument('--out', type=Path, default=ROOT / 'results/20261120_water_land_covariation')
    p.add_argument('--bootstraps', type=int, default=3000)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    synthetic_check()
    d, provenance = load(a.work.resolve())
    rows, contributions, land_contributions, slopes, meta, primary = [], [], [], [], {}, {}
    for lyca in [False, True]:
        spatial, temporal, outcomes, static = paired(d, lyca)
        cohort = 'all_butterflies' if not lyca else 'lyca_supported'
        for mode, frame in [('spatial', spatial), ('temporal', temporal)]:
            key = cohort + '_' + mode
            point, boot, info = run(frame, outcomes, a.bootstraps, 20261120)
            meta[key] = dict(info, static_max_period_difference=static)
            print(key, json.dumps(info), flush=True)
            pairs = [(0, 2), (0, 3), (1, 2), (1, 3), (0, 1), (2, 3)]
            if lyca:
                pairs += [(0, 4), (1, 4), (2, 4)]
            for stage, name in enumerate(STAGES):
                for j, k in pairs:
                    lo, hi = np.quantile(boot[1][:, stage, j, k], [.025, .975])
                    row = dict(cohort=cohort, mode=mode, stage=name, left=outcomes[j], right=outcomes[k],
                               r=point[1][stage, j, k], low=lo, high=hi, covariance=point[0][stage, j, k],
                               n_H7=len(frame), n_H6=info['n_H6'])
                    rows.append(row)
                    if not lyca and (j, k) in [(0, 2), (0, 3)]:
                        primary[(mode, name, k)] = (row, boot[1][:, stage, j, k])
            for ix, (before, after) in enumerate([(1, 2), (3, 4)]):
                for j, k in pairs:
                    for metric, m in [('covariance_removed', 0), ('correlation_removed', 1)]:
                        values = boot[m][:, before, j, k] - boot[m][:, after, j, k]
                        lo, hi = np.quantile(values, [.025, .975])
                        contributions.append(dict(cohort=cohort, mode=mode, baseline=STAGES[before],
                                                  left=outcomes[j], right=outcomes[k], metric=metric,
                                                  estimate=point[m][before, j, k]-point[m][after, j, k], low=lo, high=hi))
                for j, outcome in enumerate(outcomes):
                    lo, hi = np.quantile(boot[2][:, ix, j], [.025, .975])
                    slopes.append(dict(cohort=cohort, mode=mode, stage=STAGES[after], outcome=outcome,
                                       beta=point[2][ix, j], low=lo, high=hi))
            for j, k in pairs:
                values = boot[0][:, 1, j, k] - boot[0][:, 3, j, k]
                lo, hi = np.quantile(values, [.025, .975])
                land_contributions.append(dict(cohort=cohort, mode=mode, left=outcomes[j], right=outcomes[k],
                                               estimate=point[0][1, j, k]-point[0][3, j, k], low=lo, high=hi))
    # Approximate simultaneous max-|bootstrap error| intervals per specification:
    # the same resampling seed draws the same H6 groups in the two primary modes.
    for stage in STAGES:
        selected = [primary[(mode, stage, k)] for mode in ['spatial', 'temporal'] for k in [2, 3]]
        assert len({len(b) for row, b in selected}) == 1
        error = np.column_stack([b-row['r'] for row, b in selected])
        critical = np.quantile(np.max(np.abs(error), axis=1), .95)
        for row, b in selected:
            row['simultaneous_low'] = max(-1., row['r']-critical)
            row['simultaneous_high'] = min(1., row['r']+critical)
    pd.DataFrame(rows).to_csv(a.out / 'correlations.csv', index=False)
    pd.DataFrame(contributions).to_csv(a.out / 'W_covariance_contribution.csv', index=False)
    pd.DataFrame(land_contributions).to_csv(a.out / 'landscape_covariance_contribution.csv', index=False)
    pd.DataFrame(slopes).to_csv(a.out / 'W_response_slopes.csv', index=False)
    receipt = {'source': provenance, 'cohorts': meta, 'python': platform.python_version(),
               'numpy': np.__version__, 'pandas': pd.__version__, 'seed': 20261120,
               'checks': ['frozen source hashes and uniqueness', 'same-site water difference identities',
                          'static W verified', 'full-rank models', 'shared-driver synthetic check',
                          'FWL covariance decomposition on every valid bootstrap'],
               'limitations': ['Exploratory, not preregistered.', 'W static: no change-in-wastewater estimate.',
                               'No climate adjustment or causal attribution.', 'Two windows only, no long-term synchrony.',
                               'H6 equal weighted, marginal CIs unless simultaneous columns present.',
                               'Spatial and temporal results reuse the same sites; not independent replication.']}
    (a.out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
