"""Figures for chapter 08 (the thermal forcing and the tests a forecast needs).
Run from this folder: python make_ch08.py [preview_dir]

ch08_toff_by_elevation: T_off of the 4073 RGI7 CentralEurope glaciers in the mass-balance calibration the
thermal runs use (mhuss r7spec_global_results, monthly, ERA5), against the glacier's median elevation (GloGEM
inventory thick_centraleurope.dat); "both bounds": c_prec at its lower bound 1.3 and DDF at its upper bound 4.5,
so T_off had to close the budget.
ch08_massif_split: in-sample M1 (outline-equal RMSE below 15 m, outlines with at least 3 points) of the fair
run set by f_ref, for H and FRH, on the Mont Blanc and Monte Rosa massifs (week2/fair_test/points.csv, scored
with forcing_fair_test.m1).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUT = Path(__file__).parent
PREVIEW = Path(sys.argv[1]) if len(sys.argv) > 1 else None
REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'scripts' / 'icetemp_analysis'))
CAL = Path('/scratch_net/iceberg_second/mhuss/r7spec_global_results/monthly/CentralEurope/calibration/'
           'calibrate_m1_cID9_centraleurope_final_era5.dat')
THICK = Path('/itet-stor/jabeer/glogem/geometricdata/rgiv7/files/thick_centraleurope.dat')
POINTS = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week2/fair_test/points.csv')
plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none'})

SURFACE, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
AFTER, BEFORE = '#2a78d6', '#eb6834'   # FRH (the central forcing) / H, as the corrected-first order of ch06-07
NEUTRAL = '#bdbbb5'


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
        ax.spines[s].set_color(GRID)


def inventory():
    rows = []
    for ln in open(THICK).readlines()[1:]:
        p = ln.replace(';', ' ').split()
        rows.append((p[0], float(p[3]), float(p[9])))
    return pd.DataFrame(rows, columns=['ID', 'area', 'hmed']).set_index('ID')


def fig_toff_by_elevation():
    c = pd.read_csv(CAL, sep=r'\s+', dtype={'ID': str}).set_index('ID')
    t = c.join(inventory(), how='inner')
    both = (t.Cprec <= 1.3 + 0.005) & (t.DDFsnow >= 4.5 - 0.005)
    fig, ax = plt.subplots(figsize=(9.5, 5.2), facecolor=SURFACE)
    style(ax)
    ax.grid(color=GRID, lw=0.6)
    ax.scatter(t.hmed[~both], t.T_off[~both], s=5, color=NEUTRAL, alpha=0.6, lw=0, label='other glaciers', zorder=2)
    ax.scatter(t.hmed[both], t.T_off[both], s=6, color=BEFORE, alpha=0.6, lw=0,
               label='c_prec and DDF at their bounds: T_off closes the budget', zorder=3)
    bins = np.arange(1900, 4100, 100)
    mid, med = [], []
    for lo in bins:
        s = t.T_off[(t.hmed >= lo) & (t.hmed < lo + 100)]
        if len(s) >= 10:
            mid.append(lo + 50)
            med.append(s.median())
    ax.plot(mid, med, color=INK, lw=2, label='median per 100 m', zorder=4)
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xlabel('median glacier elevation [m a.s.l.]', color=INK)
    ax.set_ylabel('T_off [°C]', color=INK)
    ax.legend(fontsize=8, frameon=False, labelcolor=INK, loc='upper left', markerscale=3)
    fig.text(0.01, 0.005, f'{len(t)} RGI7 CentralEurope glaciers, mass-balance calibration of the thermal runs '
             '(ERA5 monthly, geodetic targets 2000–2020).', color=INK2, fontsize=8)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, 'ch08_toff_by_elevation')
    return t, both


def fig_massif_split():
    import forcing_fair_test as ft
    P = pd.read_csv(POINTS, dtype={'rgi': str})
    frefs = ft.FREFS
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.2), sharey=True, facecolor=SURFACE)
    for ax, m in zip(axes, ('Mont Blanc', 'Monte Rosa')):
        style(ax)
        ax.grid(color=GRID, lw=0.6)
        sub = P[P.massif == m]
        for v, col in (('H', BEFORE), ('FRH', AFTER)):
            y = [ft.m1(sub, ft.col(v, fr)) for fr in frefs]
            k = int(np.nanargmin(y))
            ax.plot(frefs, y, color=col, lw=2, marker='o', ms=5, mec=SURFACE, mew=1, label=v)
            ax.plot([frefs[k]], [y[k]], marker='o', ms=12, mfc='none', mec=col, mew=1.5)
            ax.annotate(f'best {frefs[k]:.2f}', (frefs[k], y[k]), xytext=(0, -16), textcoords='offset points',
                        ha='center', color=col, fontsize=8)
        ax.set_title(f'{m} ({sub.rgi.nunique()} outlines, {len(sub)} points)', loc='left', color=INK, fontsize=10)
        ax.set_xlabel('f_ref', color=INK)
    axes[0].set_ylabel('M1 in-sample [°C]', color=INK)
    axes[0].legend(fontsize=9, frameon=False, labelcolor=INK)
    fig.text(0.01, 0.005, 'Without T_off (FRH) the two massifs want opposite f_ref; with T_off (H) both sit at '
             '0.25. One Alps-wide f_ref cannot serve both.', color=INK2, fontsize=8)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save(fig, 'ch08_massif_split')


if __name__ == '__main__':
    t, both = fig_toff_by_elevation()
    hi = t.hmed > 3400
    print(f'{len(t)} glaciers; above 3400 m: {int(hi.sum())}, both bounds {100 * both[hi].mean():.1f} %, '
          f'median T_off {t.T_off[hi].median():.2f}')
    fig_massif_split()
