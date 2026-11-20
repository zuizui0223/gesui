"""Replay the unchanged historical France model and retain joint uncertainty.

Requires the original external workspace inputs. Does not download, impute,
change a cohort/model, or pool unrelated studies.
"""
import argparse
import contextlib
import hashlib
import io
import json
import os
from pathlib import Path
import runpy
import shutil
import sys
import platform

import numpy as np
import pandas as pd
import scipy
import sklearn

REPO = Path(__file__).resolve().parents[2]
DOC = REPO / 'docs/manuscript_riparian_interactions'

def record(path, root):
    return {'path': str(path.relative_to(root)).replace('\\','/'),
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'bytes': path.stat().st_size}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--scratch', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    workspace, scratch, output = [p.resolve() for p in (args.workspace,args.scratch,args.output)]
    if (workspace/'work/gesui').resolve() != REPO:
        raise ValueError('Historical script expects this checkout at workspace/work/gesui')
    scratch.mkdir(parents=True, exist_ok=True); output.mkdir(parents=True, exist_ok=True)
    original = Path(__file__).with_name('fr_interface_paths_original.py')
    replay = scratch / original.name
    shutil.copyfile(original, replay)
    stdout = io.StringIO(); previous = Path.cwd()
    try:
        os.chdir(workspace)
        with contextlib.redirect_stdout(stdout):
            state = runpy.run_path(str(replay),run_name='__main__')
    finally:
        os.chdir(previous)
    expected = pd.read_csv(DOC/'evidence/fr_interface_paths.csv',index_col='path')
    actual = pd.read_csv(scratch/'fr_interface_paths.csv',index_col='path')
    if not expected.index.equals(actual.index) or not expected.columns.equals(actual.columns):
        raise ValueError('Historical output schema changed')
    np.testing.assert_allclose(actual.to_numpy(),expected.to_numpy(),rtol=1e-10,atol=1e-12)
    draws = pd.DataFrame(state['boot'])[expected.index]
    assert len(draws)==3000 and np.isfinite(draws.to_numpy()).all()
    covariance = draws.cov()
    assert np.linalg.eigvalsh(covariance.to_numpy()).min() > -1e-10
    covariance.to_csv(output/'FR139_JOINT_BOOTSTRAP_COVARIANCE_20261120.csv',index_label='path')
    draws.to_csv(output/'FR139_JOINT_BOOTSTRAP_DRAWS_20261120.csv',index_label='draw')
    (output/'FR139_REPLAY_STDOUT_20261120.txt').write_text(stdout.getvalue(),encoding='utf-8')
    inputs = [workspace/'work/original/H7_454_HNS_panel.csv',workspace/'work/original/GESUI_303_HOST_NUTRIENT_SCORES.csv']
    inputs += sorted((workspace/'work/fr73').rglob('FR73*.csv'))
    inputs += sorted((workspace/'work/plant').glob('*ELIGIBILITY*.csv')) + sorted((workspace/'work/plant').glob('*MEANS*.csv'))
    inputs += sorted(p for p in (workspace/'work/facets').rglob('*') if p.suffix in ('.zip','.json'))
    modules = sorted({Path(m.__file__).resolve() for m in tuple(sys.modules.values())
                      if getattr(m,'__file__',None) and Path(m.__file__).resolve().is_relative_to(REPO)
                      and Path(m.__file__).suffix=='.py'})
    receipt = {'cohort': {'rows':len(state['z']),'H7':int(state['z'].h.nunique()),'H6':state['G']},
               'replicates':len(draws),'historical_paths_verified':len(actual),
               'maximum_absolute_numeric_difference':float(np.abs(actual.to_numpy()-expected.to_numpy()).max()),
               'original_script':record(original,REPO),'inputs':[record(p,workspace) for p in inputs],
               'loaded_repository_modules':[record(p,REPO) for p in modules],
               'versions':dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,sklearn=sklearn.__version__),
               'interpretation':'Joint H6 bootstrap of the same historical observational equations; not causal mediation, new independent evidence, or cross-study covariance.',
               'input_access':'Original external inputs required; this replay does not constitute fresh source acquisition.'}
    (output/'FR139_REPLAY_RECEIPT_20261120.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:receipt[k] for k in ('cohort','replicates','historical_paths_verified','maximum_absolute_numeric_difference')}))

if __name__ == '__main__':
    main()
