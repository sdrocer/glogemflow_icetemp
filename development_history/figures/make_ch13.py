"""Figures for chapter 13 (the air from above: ERA5 free air for the thermal model).
Run from this folder: python make_ch13.py [preview_dir]

ch13_air_products: product minus measured air by month, mean of the 7 stations above 2800 m, the months 2010-2016
  that every product covers (week3/climate_products/tables/summary_above_2800m.csv, window 'common all'): GloGEM's
  air today, ERA5 free air at the station height, ICON-DREAM-EU and COSMO-REA6 (cell closest in height).
ch13_mechanism: ERA5's 2 m air minus its own free air at the grid cell's height, by month, at the nine grid cells
  of the station test, 1995-2022 (week3/climate_products/tables/mechanism_2m_minus_freeair_wide.csv).
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path(__file__).parent
PREVIEW = Path(sys.argv[1]) if len(sys.argv) > 1 else None
W3 = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3')
plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none'})
SURFACE, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
BLUE, ORANGE, GREEN, LIGHT = '#2a78d6', '#eb6834', '#2f9e6b', '#9cc3ef'
MONTHS = list('JFMAMJJASOND')


def save(fig, name):
    fig.savefig(OUT / f'{name}.svg', facecolor=SURFACE)
    if PREVIEW:
        fig.savefig(PREVIEW / f'{name}.png', dpi=100, facecolor=SURFACE)
    plt.close(fig)


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.tick_params(colors=INK2)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(INK2)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def air_products():
    S = pd.read_csv(W3 / 'climate_products/tables/summary_above_2800m.csv')
    S = S[S.window == 'common all']
    rows = [('GloGEM today', 'primary node', INK2, '--', 'GloGEM today'),
            ('COSMO-REA6', 'closest-height 3x3', ORANGE, '-', 'COSMO-REA6'),
            ('ICON-DREAM-EU', 'closest-height 3x3', GREEN, '-', 'ICON-DREAM'),
            ('ERA5 free air', 'node (each station)', BLUE, '-', 'ERA5 free air')]
    fig, ax = plt.subplots(figsize=(7.2, 3.8), facecolor=SURFACE)
    style(ax)
    x = np.arange(12)
    ax.axhline(0, color=INK, lw=0.8)
    for prod, var, col, ls, lab in rows:
        r = S[(S['product'] == prod) & (S.variant == var)].iloc[0]
        y = np.array([r[f'm{m:02d}'] for m in range(1, 13)])
        ax.plot(x, y, color=col, ls=ls, lw=2, marker='o', ms=4, label=f'{lab} (year {r.annual:+.1f} °C)')
        ax.text(11.25, y[-1], lab, color=INK2, fontsize=9, va='center')
    ax.set_xticks(x)
    ax.set_xticklabels(MONTHS)
    ax.set_xlim(-0.4, 13.3)
    ax.set_ylabel('product minus measured air (°C)')
    ax.legend(loc='lower center', frameon=False, fontsize=8.5, ncol=2)
    fig.tight_layout()
    save(fig, 'ch13_air_products')


def mechanism():
    M = pd.read_csv(W3 / 'climate_products/tables/mechanism_2m_minus_freeair_wide.csv')
    fig, ax = plt.subplots(figsize=(7.2, 3.4), facecolor=SURFACE)
    style(ax)
    x = np.arange(12)
    ax.axhline(0, color=INK, lw=0.8)
    Y = M[[str(m) for m in range(1, 13)]].to_numpy()
    for y in Y:
        ax.plot(x, y, color=LIGHT, lw=1)
    ax.plot(x, Y.mean(axis=0), color=BLUE, lw=2.2, marker='o', ms=4,
            label=f'mean of the {len(Y)} grid cells (cell heights {M.hclim.min():.0f}–{M.hclim.max():.0f} m)')
    ax.plot([], [], color=LIGHT, lw=1, label='each grid cell')
    ax.set_xticks(x)
    ax.set_xticklabels(MONTHS)
    ax.set_ylabel('2 m air minus free air (°C)')
    ax.legend(loc='lower center', frameon=False, fontsize=8.5)
    fig.tight_layout()
    save(fig, 'ch13_mechanism')


if __name__ == '__main__':
    air_products()
    mechanism()
