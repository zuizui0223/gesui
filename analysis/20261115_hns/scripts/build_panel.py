"""Rebuild same H7 panel and fixed host-trait summaries from frozen originals."""
from __future__ import annotations
import json,zipfile,hashlib,collections,math,os
from pathlib import Path
import numpy as np,pandas as pd
from score_core import aggregate_traits
ROOT=Path(__file__).resolve().parents[1]
DATA=Path(os.environ.get('GESUI_DATA',str(ROOT/'inputs')))
FAMS=['Lycaenidae','Nymphalidae','Pieridae','Papilionidae','Hesperiidae']

def zj(name,member):
    with zipfile.ZipFile(DATA/name) as z:return json.loads(z.read(member))

def build():
    out=ROOT/'output';out.mkdir(exist_ok=True)
    manifest=json.loads((DATA/'SOURCE_MANIFEST.json').read_text())
    for name in ['gesui_303_host_nutrient_scores.csv','gesui_H7_posthoc_S20_source.csv','gesui_h7_234_complete_source429_repaired.zip','gesui_h7_454_original.zip','gesui_h7_h6_parent_original.zip','gesui_h7_atlas_original.zip']:
        assert hashlib.sha256((DATA/name).read_bytes()).hexdigest()==manifest[name]['sha256'],name
    raw=zj('gesui_h7_234_complete_source429_repaired.zip','GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json')
    assert len(raw)==2808 and all(x['source_status']=='ok' for x in raw)
    rr={(str(x['HYBAS_ID_L07']),x['period'],x['group']):x for x in raw};assert len(rr)==len(raw)
    with zipfile.ZipFile(DATA/'gesui_h7_454_original.zip') as z:
        d=pd.read_csv(z.open('EU234_H7_ORIGINAL_STUDY_NATIVE_BASINPERIOD_EQR_PLANT_BUTTERFLY_SOURCE_COHORT.csv'),dtype={'HYBAS_ID_L07':str},keep_default_na=False)
    d=d.rename(columns={'HYBAS_ID_L07':'h','PERIOD':'period','original_majority_BQE_method':'method','W_static_upstream_MODEL_flow_per_UpArea':'W','log_local_H7_area_km2':'area','log1p_original_EQR_site_n':'nEEA','log1p_GBIF_Animalia_source_records':'nAnimal','log1p_GBIF_Plantae_source_records':'nPlant','log1p_GBIF_5family_butterfly_source_records':'nButterfly','log1p_original_EEA_site_DIS_AV_CMS':'riverflow'})
    h6=zj('gesui_h7_h6_parent_original.zip','HYBAS7_232_TO_OFFICIAL_HYBAS6_SOURCE_PARENT_EXPLICIT_REPAIR_v2.json')['original_454_H7_study_basin_geographic_parent_mapping']
    mp={str(x['HYBAS_ID_L07']):str(x['official_parent_HYBAS_ID_L06']) for x in h6}
    atlas=zj('gesui_h7_atlas_original.zip','HYDROATLAS_V1_L07_1219_H7_BASIN_RAW_LANDUSE_ATTRIBUTES_SOURCE_ONLY.json')['original_1219_H7_covariate_by_real_ID']
    d['h6']=d.h.map(mp);d['crop']=[atlas[h]['crp_pc_use'] for h in d.h];d['urban']=[atlas[h]['urb_pc_use'] for h in d.h];d['naturalout']=[math.log1p(float(atlas[h]['dis_m3_pyr'])) for h in d.h]
    prev=pd.read_csv(DATA/'gesui_H7_posthoc_S20_source.csv',dtype={'h':str,'h6':str},keep_default_na=False)
    for col in ['h','period','country','method','h6']:assert (d[col]==prev[col]).all(),col
    for col in ['W','area','nEEA','nAnimal','nPlant','nButterfly','riverflow','crop','urban','naturalout']:np.testing.assert_allclose(d[col],prev[col],atol=1e-12,rtol=1e-12)
    d['butterfly_S20']=prev.butterfly_S20;d['plant_S50']=prev.plant_S50_new
    d['method']=d.method.replace('','__UNKNOWN__')
    country={h:sorted(collections.Counter(g.country),key=lambda c:(-sum(g.country==c),c))[0] for h,g in d.groupby('h')}
    d['country_holdout']=d.h.map(country)
    t=pd.read_csv(DATA/'gesui_303_host_nutrient_scores.csv',dtype={'butterflyKey':str})
    assert len(t)==303 and t.butterflyKey.is_unique
    keys=t.butterflyKey.tolist();ki={k:i for i,k in enumerate(keys)}
    B=np.zeros((len(d),len(t)),dtype=np.int64);raw_n=[];sourcecov=[]
    for r,x in enumerate(d.itertuples()):
        total=0;seen=set()
        for f in FAMS:
            q=rr[(x.h,x.period,f)];total+=q['source_query_occurrence_count']
            for k,n in q['source_speciesKey_counts'].items():
                assert k in ki and k not in seen,(x.h,k)
                assert t.loc[ki[k],'family']==f
                B[r,ki[k]]=int(n);seen.add(k)
        raw_n.append(total)
        # A missing coverage in a genuine zero-record family is not taxonomic quality >=80%.
        sourcecov.append(min(rr[(x.h,x.period,g)]['source_facet_coverage'] or 0. for g in FAMS+['Plantae']))
    raw_n=np.array(raw_n);d['minimum_species_coverage']=sourcecov
    np.testing.assert_array_equal(B.sum(1),prev.butterfly_N_source)
    np.testing.assert_array_equal(raw_n,prev.butterfly_N_newQuery)
    # Independent prior retained source matrix, keyed rather than positional.
    old=np.load(ROOT/'prior/gesui_host_network_20261114/output/candidate_network_matrices.npz')
    oldix={str(k):i for i,k in enumerate(old['bkeys'])}
    np.testing.assert_array_equal(B,old['B'][:,[oldix[k] for k in keys]])
    family=t.family.to_numpy();coverage=t.trait_coverage.to_numpy(float)
    accepted=(t.taxon_status=='ACCEPTED') & ~t.actual_used_alias_conflict.astype(bool)
    score_maps={
      'primary':np.where(accepted,t.host_nutrient_mean,np.nan),
      'host80':np.where(accepted & (t.trait_coverage>=.8),t.host_nutrient_mean,np.nan),
      'strict_names':np.where(accepted,t.strict_accepted_name_only_mean,np.nan),
      'no_wikipedia':np.where(accepted,t.no_wikipedia_only_mean,np.nan),
      'drop_Aglais_io':np.where(accepted & (t.accepted_name!='Aglais io'),t.host_nutrient_mean,np.nan),
      'missing_hosts_low':np.where(accepted,t.documented_host_missing_lower,np.nan),
      'missing_hosts_high':np.where(accepted,t.documented_host_missing_upper,np.nan)}
    # Bounds must retain the same scored set (not supply traits to previously unscored species).
    for v in ['missing_hosts_low','missing_hosts_high']:score_maps[v][~np.isfinite(score_maps['primary'])]=np.nan
    families_ref=None;variant_qa={};derived={}
    for name,s in score_maps.items():
        q=aggregate_traits(B,s,family,coverage,raw_n)
        for col in ['record_mean','incidence_mean','record_coverage_raw','record_coverage_identified','species_coverage','scored_records','scored_species','host_completeness','between_family','within_family']:
            derived[f'{name}_{col}']=q[col]
        for j,f in enumerate(q['families']):
            derived[f'{name}_rf_{f}']=q['record_family_shares'][:,j]
            derived[f'{name}_if_{f}']=q['incidence_family_shares'][:,j]
        variant_qa[name]={'n_scored_keys':int(np.isfinite(s).sum()),'total_raw_record_score_coverage':float(B[:,np.isfinite(s)].sum()/raw_n.sum()),'median_raw_record_score_coverage':float(np.nanmedian(q['record_coverage_raw'])),'n_finite_rows':int(np.isfinite(q['record_mean']).sum()),'family_reference_means':dict(zip(q['families'].tolist(),[float(v) if np.isfinite(v) else None for v in q['family_reference_means']]))}
        if name=='primary':families_ref=variant_qa[name]['family_reference_means']
    d=pd.concat([d,pd.DataFrame(derived,index=d.index)],axis=1)
    # All species kept in roster; only positive observations in compact long table.
    rows=[]
    for i,k in enumerate(keys):
        for r in np.flatnonzero(B[:,i]):
            rows.append({'h':d.h.iloc[r],'period':d.period.iloc[r],'h6':d.h6.iloc[r],'butterflyKey':k,'accepted_name':t.accepted_name.iloc[i],'family':family[i],'n_records':int(B[r,i]),'primary_HNS':score_maps['primary'][i],'host_completeness':coverage[i]})
    long=pd.DataFrame(rows);long.to_csv(out/'species_records_with_HNS.csv',index=False)
    d.to_csv(out/'H7_454_HNS_panel.csv',index=False)
    np.savez_compressed(out/'trait_matrices.npz',B=B,raw_totals=raw_n,keys=np.array(keys),family=family.astype(str),host_coverage=coverage,**{f'scores_{k}':v for k,v in score_maps.items()})
    t.to_csv(out/'species_score_roster_303.csv',index=False)
    qa={'n_original_rows':len(d),'H7':int(d.h.nunique()),'H6':int(d.h6.nunique()),'n_country_labels':int(d.country_holdout.nunique()),'n_source_facets':len(raw),'all_original_303_species_retained':len(t),'positive_species_record_rows':len(long),'all_source_species_records':int(B.sum()),'all_raw_family_records':int(raw_n.sum()),'historical_family_count_drift_rows':int(sum(raw_n!=prev.butterfly_N_old)),'max_historical_count_drift':int(max(abs(raw_n-prev.butterfly_N_old))),'source_coverage80_rows':int(sum(d.minimum_species_coverage>=.8)),'blank_original_method_rows':int(sum(d.method=='__UNKNOWN__')),'score_variants':variant_qa,'traits_are_not_measured_leaf_nitrogen_or_actual_sewage_dose':True}
    (out/'SOURCE_QA.json').write_text(json.dumps(qa,indent=2));print(json.dumps(qa,indent=2),flush=True)
    return d,t,score_maps,B
if __name__=='__main__':build()
