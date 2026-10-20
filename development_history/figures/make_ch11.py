"""Figures for chapter 11 (the missing heat).
Run from this folder: python make_ch11.py [preview_dir]

Both from stage 2 (tag s2, scored with the depth bins pooled per borehole and year; firn/stage2_rule.py, outputs in
week3/stage2/hole_year_{bands,points}.csv), each variant at the f_ref chosen leaving one massif out.
ch11_s_per_band: per band, the RMSE change against FRH (readings below 15 m) by band elevation, for S with the model's
own firn (FRHS) and S with the right firn (FRHLOS; observed firn only on the 9 GLAMOS glaciers).
ch11_profiles: measured and modelled temperatures at depth, band by band, on the small glaciers of the zone and the
Grenz tongue and Tête Rousse 3135 m (its frozen bed).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUT = Path(__file__).parent
PREVIEW = Path(sys.argv[1]) if len(sys.argv) > 1 else None
S2 = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3/stage2')
plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none'})
SURFACE, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
BLUE, ORANGE, NEUTRAL, LIGHT = '#2a78d6', '#eb6834', '#9a9893', '#9cc3ef'
SMALL = ('Hohsaas', 'Corvatsch', 'Alphubel South', 'Sex Rouge', 'Tortin')
ZONE = (3005, 3795)


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


def fig_s_per_band():
    t = pd.read_csv(S2 / 'hole_year_bands.csv', dtype={'rgi': str})
    small = t.band.str.contains('|'.join(SMALL))
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), facecolor=SURFACE, sharey=True)
    for ax, (v, title) in zip(axes, (('FRHS', 'S with the model\'s own firn (FRHS)'),
                                     ('FRHLOS', 'S with the right firn (FRHLOS)'))):
        style(ax)
        ax.axvspan(*ZONE, color=GRID, alpha=0.5, lw=0)
        ax.axhspan(-0.2, 0.2, color=GRID, alpha=0.9, lw=0)
        ax.axhline(0, color=INK2, lw=0.8)
        d = t[f'rmse_{v}'] - t.rmse_FRH
        for sel, col, lab in ((~small, NEUTRAL, 'other bands'), (small, ORANGE, 'small zone glaciers: ' + ', '.join(SMALL))):
            ax.scatter(t.elev[sel], d[sel], s=12 + 6 * np.sqrt(t.n[sel]), color=col, edgecolor=SURFACE, lw=0.8,
                       label=lab, zorder=3)
        ax.set_title(title, loc='left', color=INK, fontsize=10)
        ax.set_xlabel('band elevation [m]', color=INK)
    axes[0].set_ylabel('RMSE change against FRH [°C]\n← better        worse →', color=INK)
    axes[0].legend(fontsize=8, frameon=False, labelcolor=INK, loc='upper left')
    fig.text(0.01, 0.005, 'Each dot one band (site × 10 m band), readings below 15 m, f_ref chosen leaving one massif out; '
             'dot size: number of points. Grey strip: within 0.2 °C (thermistor precision);\nshaded column: the zone '
             '3000–3800 m. The right firn exists only on the 9 GLAMOS glaciers; elsewhere FRHLOS is the ledger on the '
             'model\'s balance with S.', color=INK2, fontsize=8)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    save(fig, 'ch11_s_per_band')


def fig_profiles():
    p = pd.read_csv(S2 / 'hole_year_points.csv', dtype={'rgi': str})
    keep = p.entity.str.contains('|'.join(SMALL)) | p.entity.isin(['Grenzgletscher@2505', 'Grenzgletscher@2605', 'Glacier de Tête Rousse@3135'])
    p = p[keep]
    order = sorted(p.entity.unique(), key=lambda e: (e.split('@')[0], int(e.split('@')[1].split('~')[0])))
    nc = 6
    nr = int(np.ceil(len(order) / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(13, 3.0 * nr), facecolor=SURFACE)
    axes = np.atleast_2d(axes)
    variants = (('FRH', 'FRH (today)', NEUTRAL, 'o'), ('FRHS', 'S, model firn', ORANGE, 's'),
                ('FRHLO', 'right firn', LIGHT, 'v'), ('FRHLOS', 'right firn + S', BLUE, 'D'))
    for ax, e in zip(axes.flat, order):
        style(ax)
        d = p[p.entity == e].sort_values('depth')
        ax.axvline(0, color=GRID, lw=0.8)
        for v, lab, col, mk in variants:
            ax.plot(d[v], d.depth, ls='none', marker=mk, ms=8 if v == 'FRH' else 5, color=col, mec=SURFACE, mew=0.6,
                    label=lab)
        ax.plot(d.obs, d.depth, ls='none', marker='o', ms=6.5, mfc='none', mec=INK, mew=1.4, label='measured')
        ax.invert_yaxis()
        name, band = e.split('@')
        ax.set_title(f'{name.replace("Vadret dal ", "").replace("gletscher", "").replace("Glacier du ", "").replace("Glacier de ", "")} '
                     f'{band.split("~")[0]} m', loc='left', color=INK, fontsize=9)
        ax.tick_params(labelsize=8)
    for ax in axes.flat[len(order):]:
        ax.set_visible(False)
    for ax in axes[:, 0]:
        ax.set_ylabel('depth [m]', color=INK)
    for c in range(nc):
        vis = [ax for ax in axes[:, c] if ax.get_visible()]
        if vis:
            vis[-1].set_xlabel('temperature [°C]', color=INK)
    h, l = axes.flat[0].get_legend_handles_labels()
    fig.legend(h, l, fontsize=8.5, frameon=False, labelcolor=INK, loc='lower center', ncol=5, bbox_to_anchor=(0.5, 0.01))
    fig.text(0.01, 0.005, 'Readings below 15 m, one point per borehole, year and 1 m depth bin; each variant at its '
             'out-of-sample f_ref. Where markers coincide, the larger grey one (FRH) lies underneath.', color=INK2, fontsize=8)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    save(fig, 'ch11_profiles')


if __name__ == '__main__':
    fig_s_per_band()
    fig_profiles()
