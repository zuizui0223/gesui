#!/usr/bin/env node
/**
 * Reproduce original 454 H7-period WLS within-family host-N slope randomization.
 * Dependencies: Node >=18 only. Explicitly NOT a wastewater causal/field test.
 *
 * This is EXACT algebraic coefficient recomputation: original saved species
 * slope factors project each source taxon's proportion through the fixed
 * WLS residualized exposure. Family-reference HNS in PR #95 was the
 * UNWEIGHTED mean over scored taxa; permuting labels within a family
 * preserves this reference exactly. The original design and all records,
 * strata, covariates, factor scores' missingness and family memberships
 * remain untouched. If source factor CSV or reference means change, FAIL.
 */
"use strict";
const fs=require("node:fs");
const path=require("node:path");
const crypto=require("node:crypto");
const root=path.resolve(__dirname,"../..");
const srcRel="results/20261115_hns_family_taxon/SPECIES_W_SLOPE_ATTRIBUTION.csv";
const savedRel="results/20261115_hns_trait_shuffle/PERMUTATION_RESULTS.json";
const expectedBlob="02b0ade550c7cb2a7689757519a7fcf5c1f3703a";
const source=fs.readFileSync(path.join(root,srcRel));
const gitBlob=crypto.createHash("sha1").update(Buffer.from("blob "+source.length+"\0")).update(source).digest("hex");
if(gitBlob!==expectedBlob)throw Error("Original SHA-1 Git blob mismatch: "+gitBlob);
const [header,...lines]=source.toString("utf8").trim().split(/\r?\n/).map(s=>s.split(","));
const ix=Object.fromEntries(header.map((v,i)=>[v,i]));
for(const key of ["family","host_N_fixed_score","within_family_reference_HNS","W_species_record_composition_slope_factor","W_beta_component_algebraic"])if(!(key in ix))throw Error("Missing "+key);
const rows=lines.map(a=>({family:a[ix.family],trait:Number(a[ix.host_N_fixed_score]),ref:Number(a[ix.within_family_reference_HNS]),slope:Number(a[ix.W_species_record_composition_slope_factor]),reported:Number(a[ix.W_beta_component_algebraic])}));
const familyNames=[...new Set(rows.map(r=>r.family))].sort();
const gr=Object.fromEntries(familyNames.map(f=>[f,rows.filter(r=>r.family===f)]));
if(rows.length!==186 || familyNames.length!==5)throw Error("Scored source species/families changed");
for(const f of familyNames){
  const a=gr[f],mean=a.reduce((s,r)=>s+r.trait,0)/a.length;
  if(a.some(r=>Math.abs(r.ref-mean)>1e-10))throw Error("Reference mean not invariant within "+f);
  if(a.some(r=>Math.abs(r.reported-r.slope*(r.trait-r.ref))>1e-11))throw Error("Source original attribution mismatch "+f);
}
const obsFam=Object.fromEntries(familyNames.map(f=>[f,gr[f].reduce((s,r)=>s+r.slope*(r.trait-r.ref),0)]));
const obs=familyNames.reduce((s,f)=>s+obsFam[f],0);
if(Math.abs(obs-0.12708367495390768)>1e-9)throw Error("Original source W beta not reproduced");
let seed=20261115>>>0;
function rand(){seed^=seed<<13;seed^=seed>>>17;seed^=seed<<5;return(seed>>>0)/4294967296;}
const dist=[],byFam=Object.fromEntries(familyNames.map(f=>[f,[]]));
for(let b=0;b<9999;b++){
  let all=0;
  for(const f of familyNames){
    const a=gr[f],s=a.map(r=>r.trait);
    for(let k=s.length-1;k>0;k--){const j=Math.floor(rand()*(k+1));[s[k],s[j]]=[s[j],s[k]];}
    let beta=0;
    for(let k=0;k<a.length;k++)beta+=a[k].slope*(s[k]-a[k].ref);
    byFam[f].push(beta);all+=beta;
  }
  dist.push(all);
}
function quant(a,p){const s=[...a].sort((x,y)=>x-y),i=(s.length-1)*p,l=Math.floor(i),u=Math.ceil(i);return s[l]+(s[u]-s[l])*(i-l);}
function stat(a,v){const exceed=a.filter(x=>Math.abs(x)>=Math.abs(v)-1e-12).length;const mu=a.reduce((s,x)=>s+x,0)/a.length;const sd=Math.sqrt(a.reduce((s,x)=>s+(x-mu)**2,0)/a.length);return{observed_beta:v,extreme_permutations:exceed,p_two_sided:(exceed+1)/10000,permutation_null_mean:mu,permutation_null_sd:sd,permutation_null_quantile_025:quant(a,.025),permutation_null_median:quant(a,.5),permutation_null_quantile_975:quant(a,.975)};}
const resultFam=Object.fromEntries(familyNames.map(f=>[f,{n_species:gr[f].length,...stat(byFam[f],obsFam[f])}]));
const ranked=familyNames.map(f=>({family:f,p:resultFam[f].p_two_sided})).sort((a,b)=>a.p-b.p);let q=1;
for(let j=ranked.length-1;j>=0;j--){q=Math.min(q,ranked[j].p*ranked.length/(j+1));resultFam[ranked[j].family].q_BH_five_families=q;}
const result={
  schema:"gesui-20261115-hns-within-family-label-permutation-v1",
  generated_on:"2026-11-15",
  run_status:"independently_recomputed_from_frozen_git_species_slope_csv",
  original_source_git_blob_sha1:gitBlob,
  source_csv_path:srcRel,
  source_species_n:186,
  seed:20261115,rng:"xorshift32",iterations:9999,
  original_beta_reconstruction_error:obs-0.12708367495390768,
  primary:stat(dist,obs),
  families:resultFam,
  interpretation:"Permutation concerns trait assignment conditional on unchanged original residualized WLS design and species records. It does not test true wastewater delivery, phylogenetic exchangeability, opportunistic GBIF bias, prediction transfer or causal effect. Family q values BH adjusted across five posthoc secondary comparisons.",
  plant_nutrient_index_not_fitted:true
};
if(process.argv.includes("--verify")){
  const expected=JSON.parse(fs.readFileSync(path.join(root,savedRel),"utf8"));
  if(JSON.stringify(expected)!==JSON.stringify(result))throw Error("Reproduction mismatch; original permutation results changed");
  process.stdout.write("SHA-preserved 186 species, 9999 shuffles, and stored result verified\n");
}else if(process.argv.includes("--write")){
  fs.mkdirSync(path.dirname(path.join(root,savedRel)),{recursive:true});
  fs.writeFileSync(path.join(root,savedRel),JSON.stringify(result,null,2)+"\n");
}else{
  process.stdout.write(JSON.stringify(result,null,2)+"\n");
}
