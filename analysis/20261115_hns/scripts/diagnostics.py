"""Temporal protocol and labelled post-result geography/missingness diagnostics."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import hashlib,json,zipfile
from urllib.parse import urlparse,parse_qs
import numpy as np,pandas as pd
from shapely import wkt
from shapely.ops import transform
from pyproj import Transformer
from model_core import model,FAMS
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'output'

def main():
    d=pd.read_csv(OUT/'H7_454_HNS_panel.csv',dtype={'h':str,'h6':str},keep_default_na=False)
    with zipfile.ZipFile(ROOT/'inputs/gesui_h7_234_complete_source429_repaired.zip') as z:raw=json.loads(z.read('GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json'))
    polys={};checks=0
    for row in raw:
        if row['group']!='Plantae' or row['period']!='2017-2019':continue
        h=str(row['HYBAS_ID_L07']);polys[h]={}
        for rec in row['source_request_receipts']:
            q=parse_qs(urlparse(rec['url']).query)['geometry'][0];sh=hashlib.sha256(q.encode()).hexdigest();assert sh==rec['component_SHA'];polys[h][sh]=q;checks+=1
    fwd=Transformer.from_crs('EPSG:4326','EPSG:6933',always_xy=True).transform;inv=Transformer.from_crs('EPSG:6933','EPSG:4326',always_xy=True).transform
    rows=[]
    for h,qs in polys.items():
        gs=[transform(fwd,wkt.loads(q)) for q in qs.values()];aa=np.array([abs(g.area) for g in gs]);assert aa.sum()>0
        xx=float(np.average([g.centroid.x for g in gs],weights=aa));yy=float(np.average([g.centroid.y for g in gs],weights=aa));lon,lat=inv(xx,yy)
        rows.append({'h':h,'longitude':lon,'latitude':lat,'n_polygon_parts':len(gs),'n_invalid_query_parts':sum(not g.is_valid for g in gs)})
    geo=pd.DataFrame(rows);geo.to_csv(OUT/'source_polygon_centroids.csv',index=False)
    d=d.merge(geo,on='h',validate='many_to_one',sort=False);assert len(d)==454
    d['lon2']=d.longitude**2;d['lat2']=d.latitude**2;d['lonlat']=d.longitude*d.latitude
    extras=['longitude','latitude','lon2','lat2','lonlat'];rs={}
    for kind in ['record','incidence']:
        rs['geography_'+kind]=model(d,'primary_'+kind+'_mean',kind=kind,extras=extras,save_prediction=OUT/f'geography_{kind}_predictions.csv')
    rs['geography_source80']=model(d[d.minimum_species_coverage>=.8],'primary_record_mean',extras=extras,save_prediction=OUT/'geography_source80_predictions.csv')
    # Original paired periods; static source cannot stand for exposure change.
    p1=d[d.period=='2017-2019'].set_index('h');p2=d[d.period=='2020-2022'].set_index('h');shared=sorted(set(p1.index)&set(p2.index));assert len(shared)==222
    a=p1.loc[shared].reset_index();b=p2.loc[shared].reset_index();np.testing.assert_allclose(a.W,b.W,atol=1e-15);assert (a.h6==b.h6).all()
    pair=a.copy();pair['period']='2017-2019_to_2020-2022'
    paired_desc={}
    for kind in ['record','incidence']:
        y='primary_'+kind+'_mean';pair['delta_'+y]=b[y].to_numpy(float)-a[y].to_numpy(float);pair['initial_HNS']=a[y].to_numpy(float)
        changes=['nAnimal','nPlant','nButterfly','primary_record_coverage_raw','primary_species_coverage','primary_host_completeness']+[f'primary_{"rf" if kind=="record" else "if"}_{f}' for f in FAMS[:-1]]
        de=[]
        for c in changes:pair['delta_'+c]=b[c].to_numpy(float)-a[c].to_numpy(float);de.append('delta_'+c)
        rs['temporal_'+kind]=model(pair,'delta_'+y,kind=kind,extras=['initial_HNS']+de,save_prediction=OUT/f'temporal_{kind}_predictions.csv')
        yy=pair['delta_'+y];paired_desc[kind]={'H7_pairs':len(pair),'H6':int(pair.h6.nunique()),'median_change':float(yy.median()),'mean_change':float(yy.mean()),'positive_changes':int(sum(yy>0)),'negative_changes':int(sum(yy<0)),'source_W_unchanged_in_all_pairs':True}
    pair.to_csv(OUT/'paired_222_H7_HNS_changes.csv',index=False)
    # Fixed-design adverse missing-host bounds; do not equate all-low/all-high to extrema.
    z=np.load(OUT/'trait_matrices.npz');s=z['scores_primary'];known=np.isfinite(s);lo=z['scores_missing_hosts_low'];hi=z['scores_missing_hosts_high'];sd=pd.read_csv(OUT/'species_slope_decomposition.csv',dtype={'butterflyKey':str})
    bounds={}
    for kind,g in sd.groupby('kind',sort=False):
        g=g.set_index('butterflyKey').loc[z['keys']];c=g.conditional_W_record_composition_coefficient.to_numpy()
        low=float(np.sum(c[known]*np.where(c[known]>=0,lo[known],hi[known])));high=float(np.sum(c[known]*np.where(c[known]>=0,hi[known],lo[known])))
        bounds[kind]={'fixed_design_partial_beta_min':low,'fixed_design_partial_beta_max':high,'includes_zero':low<=0<=high,'not_sampling_CI':True,'no_bounds_for_unlisted_hosts_or_unscored_species':True}
    summary=[]
    for name,r in rs.items():
        summary.append({'analysis':name,'n':r['n_rows'],'H7':r['H7'],'H6':r['n_H6'],'H6_gain_pct':r['H6']['gain_pct'],'H6_low':r['H6']['fixed_error_bootstrap95'][0],'H6_high':r['H6']['fixed_error_bootstrap95'][1],'W_beta':r['beta']['W_beta_per1SD'],'beta_low':r['beta']['W_beta_cluster95'][0],'beta_high':r['beta']['W_beta_cluster95'][1],'country_gain':r['country']['gain_pct'],'country_low':r['country']['fixed_error_bootstrap95'][0],'country_high':r['country']['fixed_error_bootstrap95'][1]})
    (OUT/'diagnostic_results.json').write_text(json.dumps({'models':rs,'paired_change':paired_desc,'missing_host_bounds':bounds,'source_geometry_SHA_checks':checks},indent=2,allow_nan=False));pd.DataFrame(summary).to_csv(OUT/'diagnostic_summary.csv',index=False)
    print(pd.DataFrame(summary).round(5).to_string(index=False));print(json.dumps({'paired':paired_desc,'missing_host_bounds':bounds},indent=2))
if __name__=='__main__':main()
