#!/usr/bin/env python3
"""Render frozen results; this script fits no model and adds no analysis."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def render(output: Path) -> None:
    result = json.loads((output / 'results.json').read_text())
    primary = result['primary']
    sensitivity = result['sensitivity_without_four_butterfly_derived_controls']
    boot = result['H6_bootstrap']
    reference = result['within_butterfly_genus_trait_reference']
    draws = pd.read_csv(output / 'within_butterfly_genus_label_reference.csv')
    assert result['status'] == 'COMPLETE_EXPLORATORY'
    assert draws.converged.all() and len(draws) == 999

    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 10,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.titleweight': 'bold', 'svg.fonttype': 'none',
        'svg.hashsalt': 'gesui-within-genus-89c349c',
    })
    fig, (left, right) = plt.subplots(1, 2, figsize=(12.6, 6.6),
                                    gridspec_kw={'width_ratios': [1.02, 1]})
    fig.subplots_adjust(left=.245, right=.97, bottom=.30, top=.77, wspace=.43)
    fig.suptitle('Within-genus butterfly recording association',
                 x=.04, ha='left', y=.96, fontsize=18, fontweight='bold')
    fig.text(.04, .885,
             'Exploratory analysis: positive primary estimate, sensitivity to controls and trait labels',
             fontsize=11, color='#444444')

    entries = [
        ('Primary: 13 controls\nH6-cluster t interval', primary['odds_ratio'],
         primary['odds_ratio95_H6_cluster_t'], '#17647c'),
        ('Same primary estimate\nH6-bootstrap percentile interval', primary['odds_ratio'],
         boot['odds_ratio95_percentile'], '#609fb0'),
        ('Nine controls\nFour butterfly-derived controls omitted', sensitivity['odds_ratio'],
         sensitivity['odds_ratio95_H6_cluster_t'], '#b45d27'),
    ]
    for y, (label, estimate, interval, color) in zip([2, 1, 0], entries):
        lo, hi = interval
        left.errorbar(estimate, y, xerr=[[estimate-lo], [hi-estimate]],
                      fmt='o', color=color, capsize=5, markersize=7, lw=2)
        left.text(estimate, y+.20, f'{estimate:.2f} [{lo:.2f}, {hi:.2f}]',
                  color=color, ha='center', fontsize=10)
    left.axvline(1, color='#777777', ls='--', lw=1)
    left.set_xscale('log')
    left.set_xlim(.9, 2.4)
    left.set_xticks([1, 1.25, 1.5, 2])
    left.set_xticklabels(['1', '1.25', '1.5', '2'])
    left.set_yticks([2, 1, 0])
    left.set_yticklabels([r[0] for r in entries], fontsize=9)
    left.set_ylim(-.45, 2.55)
    left.tick_params(axis='y', length=0, pad=10)
    left.set_xlabel('Conditional relative recording odds ratio\n+1 plant NUTRI SD; +1 HNS point between species')
    left.set_title('A  Effect and 95% intervals', loc='left', pad=18, fontsize=11)

    beta = primary['beta']
    values = draws.beta.to_numpy()
    heights, _, _ = right.hist(values, bins=np.linspace(min(values), max(values), 31),
                              color='#b9c5ce', edgecolor='white')
    right.set_ylim(0, max(heights)*1.55)
    right.axvline(beta, color='#17647c', lw=2.2)
    right.axvline(-abs(beta), color='#17647c', ls=':', lw=1.2)
    right.text(.98, .95,
               f'Observed beta = {beta:.3f}\n144 / 999 at least as extreme\nPlus-one reference fraction = {reference["two_sided_absolute_beta_reference_p"]:.3f}',
               transform=right.transAxes, ha='right', va='top', fontsize=9,
               bbox={'facecolor':'white', 'edgecolor':'none', 'alpha':.92, 'pad':4})
    right.set_xlabel('Fitted beta after within-butterfly-genus HNS shuffling')
    right.set_ylabel('Number of label assignments')
    right.set_title('B  Descriptive trait-label reference', loc='left', pad=18, fontsize=11)

    fig.text(.04, .09,
             '34 H6 catchments; 29 genera. All resampling fits converged. The bootstrap uses the same primary estimate.\n'
             'The trait-label reference is not an exact adjusted beta-only null test. High host-score coverage analysis stopped for insufficient support.\n'
             'Plant NUTRI SD = 0.3841 score units. Recorded identity contrasts do not measure abundance or temporal population change.',
             fontsize=9, linespacing=1.55, color='#444444')
    fig.text(.04, .025, 'Source: results/20261122_butterfly_within_genus/results.json | specification: 89c349c | 2026-11-22',
             fontsize=8, color='#666666')
    fig.savefig(output / 'within_genus_summary.png', dpi=190, facecolor='white')
    svg_path = output / 'within_genus_summary.svg'
    fig.savefig(svg_path, facecolor='white', metadata={'Date': None})
    svg_path.write_text('\n'.join(line.rstrip() for line in svg_path.read_text().splitlines())+'\n')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    render(parser.parse_args().output)
