"""Figures for chapter 05 (thermal spin-up). Run from this folder: python make_ch05.py [preview_dir]

The drift-vs-brute-force figure (ch05_drift_vs_bruteforce.svg) is made by
notebooks/05_spinup_drift_correction.ipynb.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt

OUT = Path(__file__).parent
PREVIEW = Path(sys.argv[1]) if len(sys.argv) > 1 else None
plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none'})

SURFACE, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']   # categorical slots 1-4, fixed order
MARKERS = ['o', 's', '^', 'D']

# Model minus glenglat, final-year profiles, same parameters in every run (score_runs.py).
STEPS = ['before\nPhase 2', '1961-90\nclimate', 'snow, rain,\nfull melt', 'no melt heat\ninto ice',
         'firn mask\nfrom balance', 'lapse rate\nrule', 'burial from\nbalance', 'relative\ndepth',
         'advection\nin spin-up']
BIAS = {  # zone -> bias per step [C]
    'Monte Rosa 3000-4200 m': [4.47, 4.46, 4.47, 4.47, 4.47, 4.32, 2.35, 2.26, 2.25],
    'Colle Gnifetti > 4200 m': [2.97, 2.77, 2.83, 2.83, 2.83, 2.08, 0.74, 0.74, 0.68],
    'Grenz tongue < 3000 m': [1.64, 1.59, 1.80, 1.40, 1.54, 1.53, 1.43, 0.40, 0.50],
    'Tete Rousse': [1.02, 1.01, 1.12, 0.96, 1.06, 1.05, 1.03, 1.03, 1.03],
}


def save(fig, name):
    fig.savefig(OUT / f'{name}.svg', facecolor=SURFACE)
    if PREVIEW:
        fig.savefig(PREVIEW / f'{name}.png', dpi=100, facecolor=SURFACE)
    plt.close(fig)


def fig_phase2_steps():
    fig, ax = plt.subplots(figsize=(10, 4.8), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    x = range(len(STEPS))
    for (zone, vals), col, mk in zip(BIAS.items(), SERIES, MARKERS):
        ax.plot(x, vals, color=col, lw=2, marker=mk, ms=7, mec=SURFACE, mew=1.5, label=zone)
        ax.text(len(STEPS) - 0.8, vals[-1], zone, color=INK, va='center', fontsize=9)
    ax.axhline(0, color=INK2, lw=1)
    ax.set_xticks(list(x))
    ax.set_xticklabels(STEPS, fontsize=8, color=INK2)
    ax.set_xlim(-0.3, len(STEPS) + 1.9)
    ax.set_ylabel('bias, model - borehole [°C]', color=INK)
    ax.tick_params(colors=INK2)
    ax.grid(axis='y', color=GRID, lw=0.8)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(GRID)
    ax.legend(loc='upper right', fontsize=8, frameon=False, labelcolor=INK)
    ax.set_title('How each Phase 2 step moved the fit to the boreholes (same parameters throughout)',
                 fontsize=10, color=INK, loc='left')
    fig.tight_layout()
    save(fig, 'ch05_phase2_steps')


if __name__ == '__main__':
    fig_phase2_steps()
    print('written:', sorted(p.name for p in OUT.glob('ch05_phase2*.svg')))
