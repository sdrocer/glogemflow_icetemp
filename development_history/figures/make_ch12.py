"""Figures for chapter 12 (sun and shade, burial, air: three tests and the orientation build).
Run from this folder: python make_ch12.py [preview_dir]

ch12_shade_slope: the shade error inside sites (within-site slope of the band bias on summer sunshine, OLS with site
  fixed effects; 90 % interval from the wild cluster bootstrap; all without the three entities Stage 1 left out of
  the within-site measures) for FRHS, Stage 1 (C, W, S), with the measured
  burial (B, B+C) and the build (three parts), against the rule-B line (half of the FRHS slope).
ch12_burial: the thermal model's burial rate (b_ref) against the measured long-term accumulation at the 41 borehole
  bands with a measurement (week3/burial_tables/burial.csv).
ch12_air_bias: the model's air temperature minus the station's, by month, at the high stations of Monte Rosa, Mont
  Blanc and the Bernese Alps / Engadin (week3/airtemp_check/tables/bias_by_station.csv, full overlap).
ch12_parts: the part factors F (central, dT 0) by elevation on the Gorner complex (01225) and Fieschergletscher
  (02611), with each part's share of the band (week3/orient_parts/parts.csv).
ch12_results: the error on left-out massifs overall, in the zone and above 3800 m for FRHS, Stage 1, the measured
  burial, burial with orientation, and the build.
"""
import json
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
BLUE, ORANGE, NEUTRAL, LIGHT, GREEN = '#2a78d6', '#eb6834', '#9a9893', '#9cc3ef', '#2f9e6b'
PART_COL = {'G': NEUTRAL, 'N': BLUE, 'S': ORANGE}


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


def load(path):
    return json.load(open(path))['variants']


def shade_slope():
    s1, rc = load(W3 / 'orient_stage1/orient_stage1.json'), load(W3 / 'orient_recheck/orient_stage1.json')
    cp = load(W3 / 'orient_parts_score/parts_score.json')
    rows = [('FRHS', s1['FRHS']['within']), ('Stage 1, central', s1['C']['within']),
            ('Stage 1, weaker', s1['W']['within']), ('Stage 1, stronger', s1['S']['within']),
            ('measured burial', rc['FRHS']['within']), ('burial + orientation', rc['C']['within']),
            ('the build (three parts)', cp['CP']['within (Stage-1 exclusions)'])]
    fig, ax = plt.subplots(figsize=(7.2, 3.8), facecolor=SURFACE)
    style(ax)
    y = np.arange(len(rows))[::-1]
    for yi, (lab, w) in zip(y, rows):
        lo, hi = w['slope_90']
        col = NEUTRAL if lab in ('FRHS', 'measured burial') else BLUE
        ax.plot([lo, hi], [yi, yi], color=col, lw=2, alpha=0.5)
        ax.plot(w['slope'], yi, 'o', color=col, ms=7)
    half = 0.5 * s1['FRHS']['within']['slope']
    ax.axvline(half, color=ORANGE, lw=1.2, ls='--')
    ax.axvline(0, color=INK2, lw=0.8)
    ax.text(half, len(rows) - 0.4, ' rule B: at least half of FRHS', color=ORANGE, fontsize=9, va='bottom')
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows])
    ax.set_xlabel('shade error inside sites (°C per +0.1 of summer sunshine)')
    ax.set_ylim(-0.7, len(rows) - 0.1)
    fig.tight_layout()
    save(fig, 'ch12_shade_slope')


def burial():
    B = pd.read_csv(W3 / 'burial_tables/burial.csv', dtype={'rgi': str})
    massif = {'00773': 'Mont Blanc', '00774': 'Mont Blanc', '01195': 'Monte Rosa', '01205': 'Monte Rosa',
              '01225': 'Monte Rosa', '01230': 'Monte Rosa'}
    B['massif'] = B.rgi.map(massif).fillna('other')
    fig, ax = plt.subplots(figsize=(6.4, 5.0), facecolor=SURFACE)
    style(ax)
    for m, col in (('Mont Blanc', ORANGE), ('Monte Rosa', BLUE), ('other', NEUTRAL)):
        d = B[B.massif == m]
        ax.errorbar(d.acc_central, d.model_bref, xerr=[d.acc_central - d.acc_low, d.acc_high - d.acc_central],
                    fmt='o', color=col, ms=5, lw=0.8, alpha=0.85, label=m)
    ax.plot([0, 6.5], [0, 6.5], color=INK2, lw=0.8)
    ax.text(3.0, 3.25, 'model = measured', color=INK2, fontsize=9, rotation=28, ha='center')
    ax.set_xlim(0, 3.6)
    ax.set_ylim(0, 6.5)
    ax.set_xlabel('measured long-term accumulation (m w.e. per year)')
    ax.set_ylabel("the thermal model's burial rate (m w.e. per year)")
    ax.legend(frameon=False, loc='upper right')
    fig.tight_layout()
    save(fig, 'ch12_burial')


def air_bias():
    A = pd.read_csv(W3 / 'airtemp_check/tables/bias_by_station.csv', dtype={'node': str})
    A = A[A.window == 'full overlap'].drop_duplicates('station')
    keep = {'Capanna Margherita': BLUE, 'Plateau Rosa': BLUE, 'Gornergrat': BLUE, 'Aiguille du Midi': ORANGE,
            'Jungfraujoch': NEUTRAL, 'Corvatsch': NEUTRAL}
    ls = {'Capanna Margherita': '-', 'Plateau Rosa': '--', 'Gornergrat': ':', 'Aiguille du Midi': '-',
          'Jungfraujoch': '-', 'Corvatsch': '--'}
    fig, ax = plt.subplots(figsize=(7.0, 3.8), facecolor=SURFACE)
    style(ax)
    months = np.arange(1, 13)
    for st, col in keep.items():
        r = A[A.station == st]
        if r.empty:
            continue
        r = r.iloc[0]
        ax.plot(months, [r[f'm{m:02d}'] for m in months], color=col, ls=ls[st], lw=1.8,
                label=f'{st} ({int(r.z)} m): year {r.annual:+.1f} °C')
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xticks(months)
    ax.set_xticklabels(list('JFMAMJJASOND'))
    ax.set_ylabel('model air minus measured (°C)')
    ax.legend(frameon=False, fontsize=8, loc='lower center', ncol=2)
    fig.tight_layout()
    save(fig, 'ch12_air_bias')


def parts():
    R = pd.read_csv(W3 / 'orient_parts/parts.csv', dtype={'rgi': str})
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 4.4), sharey=False, facecolor=SURFACE)
    for ax, (gid, name) in zip(axes, (('01225', 'Gorner complex (01225)'), ('02611', 'Fieschergletscher (02611)'))):
        style(ax)
        d = R[R.rgi == gid]
        for p, lab in (('G', 'gentle'), ('N', 'north-facing'), ('S', 'south-facing')):
            q = d[(d.part == p) & (d.flag.isna())].sort_values('band')
            ax.scatter(q['F+0'], q.band, s=8 + 60 * q.share, color=PART_COL[p], alpha=0.75, label=lab)
        ax.axvline(1, color=INK2, lw=0.8)
        ax.set_title(name, fontsize=10, color=INK)
        ax.set_xlabel('orientation factor F (dot size: share of the band)')
    axes[0].set_ylabel('band elevation (m)')
    axes[0].legend(frameon=False, loc='lower left', fontsize=9)
    fig.tight_layout()
    save(fig, 'ch12_parts')


def results():
    s1, rc = load(W3 / 'orient_stage1/orient_stage1.json'), load(W3 / 'orient_recheck/orient_stage1.json')
    cp = load(W3 / 'orient_parts_score/parts_score.json')
    rows = [('FRHS', s1['FRHS']['lomo']), ('Stage 1', s1['C']['lomo']), ('measured burial', rc['FRHS']['lomo']),
            ('burial + orientation', rc['C']['lomo']), ('the build', cp['CP']['lomo'])]
    keys = (('M1', 'all bands'), ('zone_M1', 'zone 3000-3800 m'), ('above_M1', 'above 3800 m'))
    fig, ax = plt.subplots(figsize=(7.4, 3.6), facecolor=SURFACE)
    style(ax)
    x = np.arange(len(keys))
    w = 0.16
    cols = [NEUTRAL, LIGHT, '#c9c7c2', BLUE, GREEN]
    for i, ((lab, L), col) in enumerate(zip(rows, cols)):
        vals = [L[k] for k, _ in keys]
        ax.bar(x + (i - 2) * w, vals, w, color=col, label=lab)
        for xi, v in zip(x, vals):
            ax.text(xi + (i - 2) * w, v + 0.03, f'{v:.2f}', ha='center', fontsize=7, color=INK2)
    ax.set_xticks(x)
    ax.set_xticklabels([k[1] for k in keys])
    ax.set_ylabel('error on left-out massifs (°C)')
    ax.legend(frameon=False, fontsize=8, ncol=3, loc='upper left')
    ax.set_ylim(0, 3.6)
    fig.tight_layout()
    save(fig, 'ch12_results')


if __name__ == '__main__':
    burial()
    air_bias()
    parts()
    if (W3 / 'orient_parts_score/parts_score.json').exists():
        shade_slope()
        results()
