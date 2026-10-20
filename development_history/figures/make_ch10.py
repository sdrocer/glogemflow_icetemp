"""Figures for chapter 10 (the firn ledger).
Run from this folder: python make_ch10.py [preview_dir]

ch10_check5: the firn ledger against the firn model (Beer et al. 2026) in the borehole bands of the four glaciers
with GLAMOS grids over the whole run period (scripts/icetemp_analysis/firn/ledger_vs_firnmodel.py, outputs in
week2/firnmodel/check5/): firn [m w.e.] by year, from the firn model (GLAMOS cells, its own rules), from the ledger's
rules on the GLAMOS band balance (started from the steady state of the first ten observed years), and from the
ledger GloGEM keeps (its own band balance; FRHL run of 14 Oct).
(The site-mean figure of the perfect-firn test was removed on 20 Oct: comparisons are per band at depth; the depth
profiles per band, the right firn included, are in chapter 11, figures/make_ch11.py.)
"""
import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

OUT = Path(__file__).parent
PREVIEW = Path(sys.argv[1]) if len(sys.argv) > 1 else None
CHECK5 = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week2/firnmodel/check5')
plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none'})

SURFACE, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
BLUE, ORANGE = '#2a78d6', '#eb6834'
BH = {'hohsaas': ('Hohsaas', '3275–3345 m'), 'corvatsch': ('Corvatsch', '3295–3345 m'),
      'alphubel': ('Alphubel', '3755–3825 m'), 'aletsch': ('Gr. Aletsch (Jungfraufirn)', '3415–3635 m')}
VARIANTS = (('A_firn_model', 'firn model (GLAMOS cells)', INK, '-'),
            ("C'_band_steady_start", 'ledger rules, GLAMOS balance', BLUE, '-'),
            ('D_glogem_ledger', "ledger rules, GloGEM's balance", ORANGE, '-'))


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


def fig_check5():
    t = pd.read_csv(CHECK5 / 'borehole_band_firn_by_year.csv')
    t = t[t.year.between(1972, 2025)]
    fig, axes = plt.subplots(1, 4, figsize=(12, 4.3), facecolor=SURFACE, sharex=True)
    for ax, (g, (name, band)) in zip(axes, BH.items()):
        style(ax)
        ax.grid(axis='y', color=GRID, lw=0.6)
        for v, lab, col, ls in VARIANTS:
            s = t[(t.glacier == g) & (t.variant == v)].sort_values('year')
            ax.plot(s.year, s.firn_we, color=col, lw=2 if v != 'A_firn_model' else 1.6, ls=ls, label=lab)
        ax.set_title(f'{name}\n{band}', loc='left', color=INK, fontsize=9.5)
        ax.set_ylim(bottom=0)
    axes[0].set_ylabel('firn in the borehole bands [m w.e.]', color=INK)
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, fontsize=8.5, frameon=False, labelcolor=INK, loc='lower center', ncol=3, bbox_to_anchor=(0.5, 0.075))
    fig.text(0.01, 0.005, 'Band mean of the firn left after each hydrological year. The two ledger curves differ only '
             'in their input. The firn model and the blue ledger differ in their rules and their start\n(the firn '
             'model starts from no firn in 1970, the ledger from the steady state of the first ten observed years).',
             color=INK2, fontsize=8)
    fig.tight_layout(rect=(0, 0.13, 1, 1))
    save(fig, 'ch10_check5')


if __name__ == '__main__':
    fig_check5()
