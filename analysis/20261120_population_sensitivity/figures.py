"""Three manuscript review figures; keep original and final sensitivity distinct."""
from pathlib import Path
import hashlib,json,shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
RESULT=ROOT/'results/20261120_population_sensitivity'
OUT=ROOT/'docs/manuscript_riparian_interactions/figures/population_final'
OUT.mkdir(parents=True,exist_ok=True)
pop=json.loads((RESULT/'result.json').read_text())
original=json.loads((ROOT/'results/20261120_eu451_W_confirmation/result.json').read_text())
rows=[('Original 451 rows\n232 H7 / 90 H6',original,'#929292'),
      ('Same-sample landscape controls\n437 rows / 225 H7 / 88 H6',pop['models']['same_sample_original'],'#197e95'),
      ('Same sample + upstream population\n437 rows / 225 H7 / 88 H6',pop['models']['population_adjusted'],'#a04a6c')]
fig,axes=plt.subplots(1,2,figsize=(12,5),sharey=True)
for i,(label,r,color) in enumerate(rows):
    y=2-i
    a=r['association']; b=r['prediction']
    for ax,point,ci in [(axes[0],a['adjusted_beta_per_1SD_focal'],a['H6_cluster_descriptive95']),
                        (axes[1],b['H6_equal_added_MSE_gain_pct'],b['fixed_out_of_fold_error_H6_resampling95_pct'])]:
        ax.plot(ci,[y,y],color=color,lw=2)
        ax.plot(point,y,'o',color=color,ms=7)
        ax.annotate(f'{point:.4f}' if ax==axes[0] else f'{point:.2f}%',(point,y),xytext=(0,13),textcoords='offset points',ha='center',fontsize=10)
axes[0].set_yticks([2,1,0],[r[0] for r in rows],fontsize=9)
axes[0].set_xlabel('W coefficient: HNS units per W SD\nH6-cluster 95% interval')
axes[1].set_xlabel('Heldout MSE improvement from adding W (%)\nFixed-error H6 bootstrap 95% interval')
axes[0].set_title('Conditional association')
axes[1].set_title('Additional predictive information')
for ax in axes:
    ax.axvline(0,color='#555555',ls='--',lw=.8)
    ax.set_ylim(-.5,2.6);ax.spines[['top','right']].set_visible(False);ax.grid(axis='x',alpha=.15)
fig.suptitle('Figure 1 | Butterfly host affinity: final population sensitivity',fontsize=14,y=.96)
fig.text(.025,.035,'GHS-POP 2020 upstream population: a separate product, with potentially shared census inputs.\nSeven H7 excluded by frozen upstream-area consistency rule. Static population; no causal wastewater effect identified.',fontsize=9)
fig.tight_layout(rect=[0,.13,1,.91])
for ext in ['png','pdf']:fig.savefig(OUT/f'FIGURE1_butterfly_population.{ext}',dpi=180)
plt.close(fig)

# Figure 2 retains the already verified French comparison without refitting.
for ext in ['png','pdf']:
    shutil.copyfile(ROOT/f'results/20261120_water_land_covariation/water_land_covariation.{ext}',OUT/f'FIGURE2_water_land_correspondence.{ext}')

fig,axes=plt.subplots(1,2,figsize=(12,5),gridspec_kw={'width_ratios':[1.4,1]})
uk=[('Three water-quality measures\nadded after effluent fraction',10.56,5.06,14.29),
    ('Effluent fraction\nadded after water quality',-.519,-.761,-.253)]
for i,(label,x,lo,hi) in enumerate(uk):
    axes[0].plot([lo,hi],[1-i,1-i],color='#197e95',lw=2)
    axes[0].plot(x,1-i,'o',color='#197e95',ms=7)
    axes[0].annotate(f'{x:.2f}%',(x,1-i),xytext=(0,13),textcoords='offset points',ha='center')
axes[0].set_yticks([1,0],[r[0] for r in uk],fontsize=9)
axes[0].set_title('England | 591 sites, 2010–2018\n2,624 site-years / 7 regional holdouts',fontsize=11)
for y,x in [(1,16.30),(0,19.27)]:
    axes[1].plot(x,y,'D',color='#d18535',ms=7)
    axes[1].annotate(f'{x:.2f}%',(x,y),xytext=(0,13),textcoords='offset points',ha='center')
axes[1].set_yticks([1,0],['River TP','Benthic EQR'],fontsize=10)
axes[1].set_title('Austria | original 68 sites, 2022\nReported P / flow as added predictor',fontsize=11)
for ax in axes:
    ax.axvline(0,color='#555555',ls='--',lw=.8)
    ax.set_xlabel('Prediction MSE improvement (%)')
    ax.set_ylim(-.5,1.6);ax.spines[['top','right']].set_visible(False);ax.grid(axis='x',alpha=.15)
axes[0].set_xlim(-3,18);axes[1].set_xlim(-2,25)
fig.suptitle('Figure 3 | Aquatic evidence in its own samples',fontsize=14,y=.97)
fig.text(.025,.035,'England intervals: reported 95% descriptive intervals. Austria: point estimates only; intervals unavailable in the cited summary.\nDifferent samples, outcomes and models: no pooled effect or scale comparison. This does not validate exposure measurement for butterflies.',fontsize=9)
fig.tight_layout(rect=[0,.14,1,.9],w_pad=3)
for ext in ['png','pdf']:fig.savefig(OUT/f'FIGURE3_aquatic_context.{ext}',dpi=180)
plt.close(fig)
source={'figure1':{'population_result':'results/20261120_population_sensitivity/result.json',
                   'original_result':'results/20261120_eu451_W_confirmation/result.json'},
        'figure2':{'copied_from':'results/20261120_water_land_covariation/water_land_covariation.png'},
        'figure3':{'source_table':'docs/manuscript_riparian_interactions/PATH_EVIDENCE_SCALES_20261120.md',
                   'UK_fixed_source':'https://github.com/zuizui0223/gesui/blob/76b0dd0368a14589ef496c1dd8cfc09809b02717/docs/RESULTS_20261120_CHEMPOP_WQ_EDF_AFTER_NP.md',
                   'Austria_source':'https://github.com/zuizui0223/gesui/pull/118',
                   'UK_MSE_gain':[10.56,-.519],'UK_intervals':[[5.06,14.29],[-.761,-.253]],
                   'Austria_MSE_gain':[16.30,19.27],'Austria_intervals':None}}
source['outputs']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.iterdir() if p.suffix in ['.png','.pdf']}
(OUT/'receipt.json').write_text(json.dumps(source,indent=2)+'\n',encoding='utf-8')
