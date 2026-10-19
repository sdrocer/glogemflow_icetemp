"""Figures for chapter 10 (the firn ledger).
Run from this folder: python make_ch10.py [preview_dir]

ch10_check5: the firn ledger against the firn model (Beer et al. 2026) in the borehole bands of the four glaciers
with GLAMOS grids over the whole run period (scripts/icetemp_analysis/firn/ledger_vs_firnmodel.py, outputs in
week2/firnmodel/check5/): firn [m w.e.] by year, from the firn model (GLAMOS cells, its own rules), from the ledger's
rules on the GLAMOS band balance (started from the steady state of the first ten observed years), and from the
ledger GloGEM keeps (its own band balance; FRHL run of 14 Oct).
ch10_perfect_firn: the perfect-firn test of the batch of 19 Oct (b19, f_ref 0.85), site by site on the 9 GLAMOS
glaciers: bias (model minus measured) of FRH, FRHL, FRHLo (ledger fed the observed balance) and FRHLO (also the
observed firn mask); scripts/icetemp_analysis/firn/perfect_firn_sites.py, week3/batch_rule/perfect_firn_sites.csv.
"""
import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

OUT = Path(__file__).parent
PREVIEW = Path(sys.argv[1]) if len(sys.argv) > 1 else None
CHECK5 = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week2/firnmodel/check5')
BATCH = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3/batch_rule')
NEUTRAL = '#bdbbb5'
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


def fig_perfect_firn():
    t = pd.read_csv(BATCH / 'perfect_firn_sites.csv', dtype={'rgi': str}).sort_values('bias_FRH')
    variants = (('FRH', 'FRH (today)', NEUTRAL, 'o'), ('FRHL', 'ledger, GloGEM balance', ORANGE, 's'),
                ('FRHLo', 'ledger, observed balance', BLUE, '^'),
                ('FRHLO', 'ledger and firn mask, observed', INK, 'D'))
    fig, ax = plt.subplots(figsize=(9.5, 5.2), facecolor=SURFACE)
    style(ax)
    ax.grid(axis='x', color=GRID, lw=0.6)
    ax.axvline(0, color=INK2, lw=1)
    for k, (v, lab, col, mk) in enumerate(variants):
        y = [i + (k - 1.5) * 0.17 for i in range(len(t))]
        ax.plot(t[f'bias_{v}'], y, ls='none', marker=mk, ms=7, color=col, mec=SURFACE, mew=1.2, label=lab, zorder=3)
    ax.set_yticks(range(len(t)), [f"{r.site} ({r.band_lo:.0f}{'' if r.band_hi == r.band_lo else f'–{r.band_hi:.0f}'} m)"
                                  for r in t.itertuples()])
    ax.set_xlabel('bias at the borehole points, model minus measured [°C]', color=INK)
    ax.legend(fontsize=8, frameon=False, labelcolor=INK, loc='lower right', bbox_to_anchor=(1.0, 0.06))
    ax.text(0.01, 0.01, '← too cold', transform=ax.transAxes, color=INK2, fontsize=8, va='bottom')
    ax.text(0.99, 0.01, 'too warm →', transform=ax.transAxes, color=INK2, fontsize=8, va='bottom', ha='right')
    fig.text(0.01, 0.005, 'Batch b19, f_ref 0.85. Sites of the 9 glaciers with GLAMOS grids; Fiescher, Gorner and Morteratsch '
             'are observed only from 2005-2008,\nso their observed firn barely acts. Giving the small zone glaciers their '
             'observed firn (none) makes them too cold.', color=INK2, fontsize=8)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    save(fig, 'ch10_perfect_firn')


if __name__ == '__main__':
    fig_check5()
    fig_perfect_firn()
