"""Plot all specified cross-realm pairs, including uncertain results."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'results/20261120_water_land_covariation'
d = pd.read_csv(OUT / 'correlations.csv')
d = d[d.cohort == 'all_butterflies']
stages = ['raw', 'observation', 'landscape', 'landscape_W']
labels = ['Unadjusted', 'Observation controls', '+ Landscape', '+ Wastewater proxy W']
colors = ['#909090', '#007f86', '#dc8b25', '#653a92']
pairs = [('aquatic_stress', 'plant'), ('aquatic_stress', 'butterfly'), ('TP', 'plant'), ('TP', 'butterfly')]
names = ['Reversed EQR / plant NUTRI', 'Reversed EQR / butterfly HNS', 'River logTP / plant NUTRI', 'River logTP / butterfly HNS']
fig, axes = plt.subplots(1, 2, figsize=(13, 7), sharey=True)
for ax, mode, title in zip(axes, ['spatial', 'temporal'], ['Spatial correspondence', 'Paired changes between two windows']):
    for s, (stage, label, color) in enumerate(zip(stages, labels, colors)):
        for i, (left, right) in enumerate(pairs):
            z = d[(d['mode'] == mode) & (d.stage == stage) & (d.left == left) & (d.right == right)].iloc[0]
            y = (3-i)*1.4 + (1.5-s)*.2
            ax.plot([z.low, z.high], [y, y], color=color, lw=1.8)
            ax.plot(z.r, y, 'o', color=color, ms=5, label=label if i == 0 else None)
    ax.axvline(0, color='#555555', lw=.9, ls='--')
    ax.set_xlim(-.55, .9)
    ax.set_xticks([-.5, -.25, 0, .25, .5, .75])
    ax.set_xlabel('Residual correlation (marginal 95% H6 bootstrap interval)')
    ax.set_title(title, fontsize=13, pad=18)
    ax.grid(axis='x', alpha=.15)
    ax.spines[['top', 'right']].set_visible(False)
axes[0].set_yticks(np.arange(4)*1.4, names[::-1], fontsize=10)
axes[1].legend(loc='upper right', fontsize=9, frameon=False)
fig.suptitle('Matched aquatic and terrestrial indicators: France, 73 catchments / 24 parent basins', fontsize=14, y=.97)
fig.text(.02, .025, 'Water: identical sites in 2018 and 2022. Plants and butterflies: 2017–2019 and 2020–2022.\nW is static geography, not a measured change in sewage input. Intervals here are marginal; simultaneous primary intervals are in correlations.csv.', fontsize=9)
fig.tight_layout(rect=[0, .09, 1, .93])
fig.savefig(OUT / 'water_land_covariation.png', dpi=180)
fig.savefig(OUT / 'water_land_covariation.pdf')
