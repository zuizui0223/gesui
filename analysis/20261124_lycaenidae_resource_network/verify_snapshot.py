#!/usr/bin/env python3
"""Verify archive bytes and optional replay of two preserved analysis scripts."""
from pathlib import Path
import argparse, hashlib, json, shutil, subprocess, sys
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
M = json.loads((HERE / 'SOURCE_MANIFEST.json').read_text())

def digest(path):
    b = path.read_bytes()
    return len(b), hashlib.sha256(b).hexdigest(), hashlib.sha1(f'blob {len(b)}\0'.encode()+b).hexdigest()

def verify_snapshot():
    for x in M['files']:
        assert digest(HERE/x['name']) == (x['bytes'],x['sha256'],x['git_blob_sha1']), x['name']
    print(f"SNAPSHOT: {len(M['files'])}/{len(M['files'])} exact bytes SHA256 and Git blob SHA1 PASS",flush=True)

def compare_csv(original, replay):
    a=pd.read_csv(original,dtype={'h':str,'h6':str,'butterfly_key':str},keep_default_na=False)
    b=pd.read_csv(replay,dtype={'h':str,'h6':str,'butterfly_key':str},keep_default_na=False)
    assert a.shape==b.shape and list(a)==list(b), original.name
    maxdiff=0.
    for k in a:
        x,y=a[k],b[k]
        # Identity and boolean columns require exact equality; numpy booleans cannot be subtracted.
        if k in {'h','h6','butterfly_key','butterflyKey','plantKey','speciesKey','HYBAS_ID_L07','HYBAS_ID_L06'} or pd.api.types.is_bool_dtype(x) or pd.api.types.is_bool_dtype(y):
            assert (x.astype(str)==y.astype(str)).all(), (original.name,k,'identity or boolean mismatch')
            continue
        # Numeric columns compare with full row alignment; explicit missingness preserved.
        n1,n2=pd.to_numeric(x,errors='coerce'),pd.to_numeric(y,errors='coerce')
        if (n1.notna() | (x=='')).all() and (n2.notna() | (y=='')).all():
            assert (n1.isna()==n2.isna()).all(), (original.name,k,'missingness')
            good=n1.notna()
            if good.any():
                np.testing.assert_allclose(n1[good],n2[good],atol=1e-8,rtol=1e-8,err_msg=f'{original.name}:{k}')
                maxdiff=max(maxdiff,float(np.max(np.abs(n1[good]-n2[good]))))
        else:
            assert (x.astype(str)==y.astype(str)).all(), (original.name,k)
    print(f'REPLAY {original.name}: {len(a)} rows, max absolute numeric difference {maxdiff:.3g}',flush=True)

def replay(work):
    work.mkdir(parents=True,exist_ok=True)
    for x in M['sources']:
        assert hashlib.sha256((work/x['name']).read_bytes()).hexdigest()==x['sha256'],x['name']
    for name in ['gesui_fab_resource_consumer.py','gesui_plant_resource_partition.py','gesui_rebuilt_network.csv']:
        if (HERE/name).resolve() != (work/name).resolve(): shutil.copyfile(HERE/name,work/name)
    for name in ['gesui_fab_resource_consumer.py','gesui_plant_resource_partition.py']:
        subprocess.run([sys.executable,str(work/name),str(work)],check=True)
    for name in ['gesui_fab_consumer_group_summary.csv','gesui_fab_consumer_model_results.csv','gesui_fab_paired_contrasts.csv','gesui_fab_temporal_decomposition.csv','gesui_lycaenidae_catalogue_groups.csv','gesui_plant_resource_partition_results.csv','gesui_plant_resource_partition_contrasts.csv','gesui_plant_resource_partition_panel.csv']:
        compare_csv(HERE/name,work/name)
    print('REPLAY PASS: both analysis scripts and 8 derived CSVs; no independent field samples added',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--work',type=Path);a=p.parse_args()
    verify_snapshot()
    if a.work:
        assert a.work.resolve()!=HERE, 'Replay must not overwrite archived snapshot'
        replay(a.work)
