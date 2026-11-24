#!/usr/bin/env python3
"""Diagnostic of exact plant score overlap, not a causal mediation test."""
from pathlib import Path
import sys
P=Path(sys.argv[1]);p=P/'gesui_filter_compare.py'
exec(compile(p.read_text().split('for cutoff in [3,2]:')[0],str(p),'exec'))
new=pd.read_csv(P/'filter_clarke_edges.csv',dtype={'plantKey':str,'butterflyKey':str});assert hashlib.sha256((P/'filter_clarke_edges.csv').read_bytes()).hexdigest()==r['edge_sha256']
blocked=set(pk)|set(new.plantKey);taxnames=tx.set_index('speciesKey').gbif_canonical_name.to_dict();genera={str(taxnames[k]).split()[0].casefold() for k in blocked if k in taxnames and isinstance(taxnames[k],str)};unlisted_species=set(S)-blocked;unlisted_genera={k for k in S if str(taxnames[k]).split()[0].casefold() not in genera};assert unlisted_genera<=unlisted_species
for lab,keys in [('exclude_host_species',unlisted_species),('exclude_host_genera',unlisted_genera)]:
 vals=[];ns=[]
 for x in d.itertuples():
  c=fac[x.h,x.period,'Plantae']['source_speciesKey_counts'];s=[S[k] for k in c if k in keys];vals.append(np.mean(s) if s else np.nan);ns.append(len(s))
 d[lab]=vals;d[lab+'_n']=ns
# reconstruct the exact same eligibility and scores from previous panel; compare all original covariates with verified source.
panel=pd.read_csv(P/'filter_panel.csv',dtype={'h':str,'h6':str});panel=panel[panel.cutoff==3].copy();assert len(panel)==414
q0=d.merge(panel[['h','period','old_HNS','clarke_HNS','common_record_cov','common_n','plant_N']],on=['h','period'],suffixes=('','_archived'),validate='one_to_one');np.testing.assert_allclose(q0.plant_N,q0.plant_N_archived,atol=1e-12,rtol=0);q0['common_log_n']=np.log1p(q0.common_n)
gate=[]
for cohort,df in [('EU',q0),('FR73',q0[q0.frh6.notna()].copy())]:
 if cohort=='FR73':df['category']='FR::I2M2'
 for test in ['exclude_host_species','exclude_host_genera']:
  good=df[df[test+'_n']>=30].groupby('h').period.nunique();keep=set(good[good==2].index);s=df[df.h.isin(keep)].copy();g={'cohort':cohort,'test':test,'original_H7':df.h.nunique(),'H7':s.h.nunique(),'H6':s.h6.nunique(),'min_remaining':int(s[test+'_n'].min()) if len(s) else 0,'median_remaining':float(s[test+'_n'].median()),'eligible_plant_keys':len(unlisted_species if test=='exclude_host_species' else unlisted_genera)};g['status']='FIT' if g['H7']>=.8*g['original_H7'] and g['H6']>=20 else 'STOP_SUPPORT';gate.append(g)
  if g['status']!='FIT':continue
  bcols=base+['W']+(['TP','log_sites'] if cohort=='FR73' else [])
  for out in ['old_HNS','clarke_HNS']:
   for pred in ['plant_N',test]:evaluate(cohort,3,test,s,out,pred,bcols)
for lab,rr in [('gate',gate),('results',result),('losses',losses)]:pd.DataFrame(rr).to_csv(P/f'filter_nonhost_{lab}.csv',index=False)
print('GATES',pd.DataFrame(gate).to_string(index=False));print('RESULTS',pd.DataFrame(result)[['cohort','mode','predictor','outcome','H7','H6','gain','lo','hi']].round(4).to_string(index=False))
