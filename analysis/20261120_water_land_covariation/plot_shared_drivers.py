"""Show every candidate block and heldout response, without significance filtering."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

OUT=Path(__file__).resolve().parents[2]/'results/20261120_water_land_covariation/shared_drivers'
d=pd.read_csv(OUT/'covariance_allocation.csv')
cv=pd.read_csv(OUT/'heldout_prediction.csv')
fig,axes=plt.subplots(1,3,figsize=(14,5.5),gridspec_kw={'width_ratios':[1,1,1.3]})
for ax,right,title in zip(axes[:2],['plant','butterfly'],['Aquatic state / plants','Aquatic state / butterflies']):
    for i,(block,label,color) in enumerate([('urban','Urban','#555e68'),('crop','Cropland','#47905a'),('hydrology','Basin / flow','#238ca3'),('W','Wastewater W','#8b5294')]):
        r=d[(d.kind=='shapley')&(d.left=='aquatic_stress')&(d.right==right)&(d.block==block)].iloc[0]
        y=3-i
        ax.plot([r.low,r.high],[y,y],color=color,lw=2)
        ax.plot(r.estimate,y,'o',color=color,ms=7)
    ax.set_yticks([3,2,1,0],['Urban','Cropland','Basin / flow','Wastewater W'])
    ax.set_xlim(-.06,.23)
    ax.set_ylim(-.55,3.55)
    ax.set_xlabel('Allocated covariance reduction\n(order-averaged, standardized)')
    ax.set_title(title,pad=15)
for stage,offset,color,label in [('landscape_over_observation',.12,'#238ca3','Add landscape'),('landscape_W_over_landscape',-.12,'#8b5294','Then add W')]:
    for i,(outcome,name) in enumerate([('aquatic_stress','Aquatic state'),('TP','River TP'),('plant','Plants'),('butterfly','Butterflies')]):
        r=cv[(cv.comparison==stage)&(cv.outcome==outcome)].iloc[0]
        y=3-i+offset
        axes[2].plot([r.low,r.high],[y,y],color=color,lw=2)
        axes[2].plot(r.improvement_pct,y,'o',color=color,ms=6,label=label if i==0 else None)
axes[2].set_yticks([3,2,1,0],['Aquatic state','River TP','Plants','Butterflies'])
axes[2].set_xlim(-25,70)
axes[2].set_ylim(-.55,3.55)
axes[2].set_xlabel('Heldout H6 MSE improvement (%)\nPositive = better prediction')
axes[2].set_title('Leave-one-parent-basin-out prediction',pad=15)
axes[2].legend(loc='lower right',fontsize=8,frameon=False)
for ax in axes:
    ax.axvline(0,color='#555555',ls='--',lw=.8)
    ax.spines[['top','right']].set_visible(False)
    ax.grid(axis='x',alpha=.15)
fig.suptitle('Shared landscape signal, uncertain additional wastewater information',fontsize=15,y=.97)
fig.text(.025,.035,'France: 73 catchments / 24 H6. Spatial two-window means; not temporal synchrony. EQR is sign-reversed.\nLeft: marginal H6-bootstrap intervals, descriptive allocation, not causal effects. Right: intervals resample fixed heldout H6 errors.',fontsize=9)
fig.tight_layout(rect=[0,.12,1,.91],w_pad=2)
fig.savefig(OUT/'shared_drivers.png',dpi=180)
fig.savefig(OUT/'shared_drivers.pdf')
