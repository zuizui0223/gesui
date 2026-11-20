"""Independent statsmodels check of every saved point correlation."""
from pathlib import Path
import importlib.util
import json
import numpy as np
import pandas as pd
import statsmodels.api as sm

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('water_land_analysis', HERE / 'analyse.py')
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
OUT = mod.ROOT / 'results/20261120_water_land_covariation'
d, provenance = mod.load(mod.ROOT.parent)
saved = pd.read_csv(OUT / 'correlations.csv')
errors = []
for lyca in [False, True]:
    spatial, temporal, outcomes, static = mod.paired(d, lyca)
    cohort = 'lyca_supported' if lyca else 'all_butterflies'
    for mode, frame in [('spatial', spatial), ('temporal', temporal)]:
        weights = 1 / frame.h6.map(frame.h6.value_counts()).to_numpy(float)
        controls = [[], mod.OBS, mod.OBS + ['W'], mod.OBS + mod.LAND, mod.OBS + mod.LAND + ['W']]
        for stage, columns in zip(mod.STAGES, controls):
            # Build the model independently with unstandardized predictors.
            active = [c for c in columns if frame[c].std(ddof=0) > 1e-10]
            X = sm.add_constant(frame[active], has_constant='add').to_numpy(float)
            residuals = np.column_stack([sm.WLS(frame[c], X, weights=weights).fit().resid for c in outcomes])
            cov = np.cov(residuals, rowvar=False, aweights=weights, ddof=0)
            corr = cov / np.sqrt(np.outer(cov.diagonal(), cov.diagonal()))
            matched = saved[(saved.cohort == cohort) & (saved['mode'] == mode) & (saved.stage == stage)]
            for row in matched.itertuples():
                actual = corr[outcomes.index(row.left), outcomes.index(row.right)]
                errors.append(abs(actual-row.r))
assert len(errors) == len(saved) == 150
assert max(errors) < 1e-9
assert not saved.duplicated(['cohort', 'mode', 'stage', 'left', 'right']).any()
receipt = json.loads((OUT / 'receipt.json').read_text())
assert all(m['valid_bootstraps'] == 3000 for m in receipt['cohorts'].values())
for item in receipt['source']['inputs']:
    p = mod.ROOT.parent / item['path']
    assert mod.sha(p) == item['sha256'], p
result = {'point_correlations_verified': len(errors), 'max_absolute_error': max(errors),
          'statsmodels_version': sm.__version__, 'input_hashes_verified': True,
          'validation_script_sha256': mod.sha(Path(__file__))}
(OUT / 'validation.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
print(json.dumps(result))
