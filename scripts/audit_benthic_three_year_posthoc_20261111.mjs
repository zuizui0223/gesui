// Read-only, POSTHOC audit; uses only source-pinned archived files from gesui PR #56.
// From repository root: node scripts/audit_benthic_three_year_posthoc_20261111.mjs
// Node.js >= 18; no external packages; emits JSON to stdout.
import fs from 'node:fs';
import crypto from 'node:crypto';
const sources = {
  region2016: ['results/eea2016_thirdyear_source_20261111/2016_EEA_RIVER_BENTHIC_EQR_JRC_GREEN_MODEL_POINTNP_NUTS2_HELDOUT_RESULTS.json', '108ac711022921af8cfb2b5477d046f386d43f17'],
  site2018: ['results/eea2018_independent_temporal_bqe_20261108/benthic_invertebrates_2018_WWTP_source_real_EQR.csv', '0f30bfbc3dac13f6a182c8a910e9899c04778322'],
  site2022: ['results/eea2022_biology_bqe_20261108/BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv', '325f199e985d9862ffee5b0edb205c61a7be3833']
};
function sourceText([path, expected]) {
  const bytes = fs.readFileSync(path);
  const sha = crypto.createHash('sha1').update('blob ' + bytes.length + '\0').update(bytes).digest('hex');
  if (sha !== expected) throw Error('Pinned source SHA mismatch: ' + path + ' ' + sha);
  return bytes.toString('utf8');
}
function csv(text) {
  const lines = text.trim().split(/\r?\n/), columns = lines.shift().split(',');
  return lines.map(line => {
    const fields = line.split(',');
    if (fields.length !== columns.length) throw Error('CSV fields changed');
    return Object.fromEntries(columns.map((h,i) => [h,fields[i]]));
  });
}
const mean = a => a.reduce((sum,x) => sum+x,0)/a.length;
function covarianceWithin(blocks, predictor) {
  const sums = {xx:0,yy:0,xy:0};
  for (const g of blocks) {
    const xm=mean(g.map(z=>z[predictor]));
    const ym=mean(g.map(z=>z.mean_EQR_original_numeric));
    for (const z of g) {
      const dx=z[predictor]-xm,dy=z.mean_EQR_original_numeric-ym;
      sums.xx+=dx*dx; sums.yy+=dy*dy; sums.xy+=dx*dy;
    }
  }
  return sums;
}
const r = s => s.xy/Math.sqrt(s.xx*s.yy);
function rng(seed) {
  let s=seed>>>0;
  return () => {s=(Math.imul(s,1664525)+1013904223)>>>0;return s/4294967296;};
}
const orig16=JSON.parse(sourceText(sources.region2016));
const regions=orig16.results.benthic_invertebrates.original_NUTS2_disjoint_holdout_Ridge10.source_JRC_green_region_sampled_EEA_EQR_summary;
const groups={};
for(const z of regions)(groups[z.country]??=[]).push(z);
const countries=Object.keys(groups);
if(regions.length!==73 || countries.length!==12) throw Error('2016 region cohort changed');
const random=rng(20261111);
const stats16={};
for(const [label,key] of Object.entries({
  N:'mean_logPointN_modelled2014to18_density',
  P:'mean_logPointP_modelled2014to18_density'
})) {
  const original=r(Object.values(groups).reduce((a,g)=>a.concat(g),[]).length ? covarianceWithin(Object.values(groups),key) : null);
  const omit=Object.fromEntries(countries.map(c=>[c,r(covarianceWithin(countries.filter(x=>x!==c).map(x=>groups[x]),key))]));
  const individual=Object.fromEntries(countries.filter(c=>groups[c].length>=3).map(c=>[c,r(covarianceWithin([groups[c]],key))]));
  const bootstrap=[];
  for(let i=0;i<10000;i++) {
    const sampled=Array.from({length:countries.length},()=>groups[countries[Math.floor(random()*countries.length)]]);
    bootstrap.push(r(covarianceWithin(sampled,key)));
  }
  bootstrap.sort((a,b)=>a-b);
  const percentile=q=>bootstrap[Math.floor((bootstrap.length-1)*q)];
  stats16[label]={
    regional_country_centered_pearson:original,
    leave_one_country_out_min:Math.min(...Object.values(omit)),
    leave_one_country_out_max:Math.max(...Object.values(omit)),
    leave_one_country_out:omit,
    n_countries_with_at_least_3_regions:Object.keys(individual).length,
    n_negative_countries:Object.values(individual).filter(z=>z<0).length,
    mean_of_country_correlations:mean(Object.values(individual)),
    per_country_correlation:individual,
    posthoc_country_block_bootstrap_10000_descriptive95:[percentile(.025),percentile(.975)]
  };
}
const data18=csv(sourceText(sources.site2018));
const data22=csv(sourceText(sources.site2022));
const index22=new Map(data22.map(z=>[z.site_key,z]));
const matched=[];
for (const a of data18) {
  const b=index22.get(a.site_key);
  if (!b) continue;
  const n18=Number(a['2018_measured_upstream_N_t_per_year'])>0;
  const n22=Number(b.upstream_2022_EEA_N_measured_t_per_yr)>0;
  if(n18!==(Number(a['2018_measured_upstream_P_t_per_year'])>0) ||
     n22!==(Number(b.upstream_2022_EEA_P_measured_t_per_yr)>0))
    throw Error('N and P source-presence mismatch');
  matched.push({
    status:String(Number(n18))+String(Number(n22)),
    delta:Number(b.normalised_EQR)-Number(a.numeric_EQR_2018),
    methodStable:a.procedureClassificationSystem===b.procedureClassificationSystem
  });
}
if(data18.length!==2634 || data22.length!==2475 || matched.length!==1167)
  throw Error('Expected original site cohort changed');
function table(z) {
  return Object.fromEntries(['00','01','10','11'].map(status=>{
    const a=z.filter(t=>t.status===status);
    return [status,{sites:a.length,mean_EQR_change:a.length?mean(a.map(t=>t.delta)):null}];
  }));
}
const report={
  schema:'gesui_posthoc_benthic_three_year_20161111_v1',
  not_causal:true, not_independent_third_measured_annual_WWTP_year:true,
  source_git_blob_shas:Object.fromEntries(Object.entries(sources).map(([k,v])=>[k,{path:v[0],git_sha:v[1]}])),
  year2016:{regions:regions.length,countries:countries.length,region_mean_EQR_vs_JRC2014_2018_modelled_point_NP:stats16},
  matched2018_2022:{site_count_2018:data18.length,site_count_2022:data22.length,shared_sites:matched.length,
    full:table(matched),method_stable_count:matched.filter(x=>x.methodStable).length,
    method_stable:table(matched.filter(x=>x.methodStable)),
    method_changed_count:matched.filter(x=>!x.methodStable).length},
  limits:[
    'Region-mean correlations use NUTS2 as unit; within-country centering removes country means, not land-use or method confounding.',
    '2016 JRC GREEN modeled 2014-2018 regional domestic+industrial loads are not actual 2016 measured wastewater N/P.',
    '2018/2022 detected WWTP report-present changes are not validated treatment changes; only 7 gain and 4 loss sites.',
    'Country bootstrap is posthoc descriptive and does not establish independent spatial or causal units.',
    'All data preexist and no new ecological response was collected.'
  ]
};
process.stdout.write(JSON.stringify(report,null,2)+'\n');
