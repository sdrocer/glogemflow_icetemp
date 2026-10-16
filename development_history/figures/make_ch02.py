"""Figures for chapter 02 (parameter redesign).
Run from this folder: python make_ch02.py [preview_dir]

ch02_glacier_optima: each glacier's own best (perm_frac, dT_scale), fitted alone against campaign 8's
emulator on a 41 x 41 grid over its training range, all depths (notebooks/archive/
06_parameter_identifiability.ipynb, Section 2; icetemp ca37122, 56d0e08). Per glacier from
data/figures/per_glacier_fit.csv, Grenzgletscher's bands from data/figures/per_entity_fit.csv, the
design box from data/bayescal/centraleurope_covfix/design.csv. The one shared set (0.119, 0.052) is
the notebook's printed best over all 695 observations; the notebook does not store it.
ch02_conductivity: thermal conductivity against density, before (linear ice/air mix, kice 2.33) and
after GloGEM 1aef7f6 (Calonne 2011 x f_ins, capped at Yen's ice value, Yen 1981 ice law at and above
900 kg m-3), with Sturm (1997) and Yen's (1981) snow relation inside their published ranges; and the
snow/ice diffusivity ratio k/C over k_ice/C_ice that the heat equation uses, with C the linear
ice/air mix (cice 1.89e6, cair 1297 J m-3 K-1) and ice taken at 900 kg m-3 and 0 °C, as in
settings.pro. Constants from settings.pro and the 1aef7f6 commit message.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

OUT = Path(__file__).parent
PREVIEW = Path(sys.argv[1]) if len(sys.argv) > 1 else None
REPO = Path(__file__).resolve().parents[2]
FITS = REPO / 'data' / 'figures'
plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none'})

SURFACE, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
AFTER, BEFORE = '#2a78d6', '#eb6834'   # categorical slots 1-2, fixed order; corrected state first
NEUTRAL = '#bdbbb5'
AFTER_LIGHT = '#cde2fb'                 # blue ramp step 100

SHARED = (0.119, 0.052)   # notebook 06, cell 6: one shared set, best over all 695 observations

# settings.pro (GloGEM GloGEMflow, from 1aef7f6)
KICE_OLD, KAIR, CICE, CAIR = 2.33, 0.001, 1890000., 1297.
K_CA = (0.024, -1.23e-4, 2.5e-6)
K_ICE_A, K_ICE_B, RHO_ICE = 9.828, -5.7e-3, 900.

SHORT = {'Glacier des Bossons': 'Bossons', 'Taconnaz Glacier': 'Taconnaz',
         'Vedretta Alta dell’Ortles': 'Ortles (Alta)', 'Vadret dal Corvatsch': 'Corvatsch',
         'Glacier du Sex Rouge': 'Sex Rouge'}


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


def fig_glacier_optima():
    g = pd.read_csv(FITS / 'per_glacier_fit.csv')
    e = pd.read_csv(FITS / 'per_entity_fit.csv')
    grenz = e[e.base == 'Grenzgletscher']
    d = pd.read_csv(REPO / 'data' / 'bayescal' / 'centraleurope_covfix' / 'design.csv')
    lo, hi = d[['perm_frac', 'dT_scale']].min().values, d[['perm_frac', 'dT_scale']].max().values
    size = lambda n: 7 * np.asarray(n, float)

    fig, ax = plt.subplots(figsize=(9.2, 6.2), facecolor=SURFACE)
    style(ax)
    ax.add_patch(Rectangle(lo, *(hi - lo), fc='none', ec=INK2, lw=1, ls=(0, (4, 3)), zorder=1))
    ax.text(hi[0], hi[1] + 0.06, 'design box: range of campaign 8\'s 150 training runs', ha='right',
            va='bottom', color=INK2, fontsize=8.5)
    # Grenzgletscher's bands (jittered slightly: several share a grid cell)
    rng = np.random.default_rng(1)
    jx, jy = rng.normal(0, 0.006, len(grenz)), rng.normal(0, 0.03, len(grenz))
    ax.scatter(grenz.pf_opt + jx, grenz.ds_opt + jy, s=size(grenz.n), facecolor='none', edgecolor=INK,
               lw=0.8, zorder=4)
    ax.scatter(g.pf_opt, g.ds_opt, s=size(g.n), color=AFTER, alpha=0.75, edgecolor=SURFACE, lw=1.8, zorder=3)
    ax.scatter(*SHARED, s=260, marker='*', color=BEFORE, edgecolor=SURFACE, lw=1.5, zorder=5, clip_on=False)
    ax.annotate('one shared set for all 11\n(costs a median +0.67 °C)', SHARED, xytext=(0.2, 0.55),
                color=INK, fontsize=9, va='center', arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    offs = {'Grenzgletscher': (-6, -30), 'Glacier des Bossons': (-14, 14), 'Lysgletscher': (-14, 14),
            'Hintereisferner': (12, 10), 'Taconnaz Glacier': (8, -14), 'Gornergletscher': (12, -4),
            'Vedretta Alta dell’Ortles': (10, -14), 'Vadret dal Corvatsch': (10, -16),
            'Fieschergletscher': (10, 6), 'Breithornplateau': (10, 0), 'Glacier du Sex Rouge': (-12, 0)}
    for _, r in g.iterrows():
        dx, dy = offs[r.glacier]
        ax.annotate(f'{SHORT.get(r.glacier, r.glacier)} ({r.n})', (r.pf_opt, r.ds_opt), xytext=(dx, dy),
                    textcoords='offset points', ha='left' if dx > 0 or r.glacier == 'Grenzgletscher' else 'right', va='center', color=INK,
                    fontsize=8.5)
    ax.annotate('open circles: Grenzgletscher\'s 21 bands,\neach fitted alone', (0.99, 0.68),
                xytext=(0.70, 1.65), color=INK2, fontsize=8.5, va='center', ha='left',
                arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    ax.set_xlim(-0.04, 1.06)
    ax.set_ylim(-0.6, 5.35)
    ax.set_xlabel('best perm_frac (scales the percolation depth) [–]', color=INK)
    ax.set_ylabel('best dT_scale (scales the decision-tree surface offset) [–]', color=INK)
    ax.set_title('Each glacier fitted alone wants its own parameters: the optima spread over 92 % of '
                 'both ranges', fontsize=10, color=INK, loc='left')
    handles = [plt.Line2D([], [], marker='o', ls='', ms=np.sqrt(7 * 40), color=AFTER, mec=SURFACE,
                          label='one glacier; area ∝ number of observations (in brackets)'),
               plt.Line2D([], [], marker='o', ls='', ms=np.sqrt(7 * 15), mfc='none', mec=INK,
                          label='one Grenzgletscher band')]
    ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.0, 0.93), fontsize=8.5,
              frameon=False, labelcolor=INK, labelspacing=1.2)
    fig.text(0.01, 0.012, 'Fitted against campaign 8\'s emulator, all depths, 41 × 41 grid '
             '(notebooks/archive/06_parameter_identifiability.ipynb; data/figures/per_glacier_fit.csv,\n'
             'per_entity_fit.csv). Glaciers with at least 8 observations; Grenzgletscher bands with at '
             'least 5, jittered slightly where they share a grid cell.', color=INK2, fontsize=8)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.94, bottom=0.15)
    save(fig, 'ch02_glacier_optima')


def calonne(rho):
    return K_CA[0] + K_CA[1] * rho + K_CA[2] * rho ** 2


def k_ice(t_c):
    return K_ICE_A * np.exp(K_ICE_B * (t_c + 273.15))


def heat_cap(rho):
    return (1 - rho / 1000.) * CAIR + rho / 1000. * CICE


def fig_conductivity():
    rho = np.linspace(100, 900, 400)
    kcap = k_ice(0.0)
    old = (1 - rho / 1000.) * KAIR + rho / 1000. * KICE_OLD
    new = {f: np.minimum(f * calonne(rho), kcap) for f in (0.6, 1.0, 1.4)}
    sturm_r = np.linspace(156, 600, 200)
    sturm = 0.138 - 1.01 * sturm_r / 1000 + 3.233 * (sturm_r / 1000) ** 2
    yen_r = np.linspace(100, 600, 200)
    yen = 2.22362 * (yen_r / 1000) ** 1.885

    fig, (ax, axd) = plt.subplots(1, 2, figsize=(10.5, 4.9), facecolor=SURFACE,
                                  gridspec_kw=dict(width_ratios=[1.15, 1], wspace=0.28))
    for a in (ax, axd):
        style(a)
        a.grid(color=GRID, lw=0.8)
        a.set_axisbelow(True)

    ax.fill_between(rho, new[0.6], new[1.4], color=AFTER_LIGHT, lw=0, zorder=1)
    ax.plot(rho, old, color=BEFORE, lw=2, zorder=3)
    ax.plot(rho, new[1.0], color=AFTER, lw=2, zorder=3)
    for f in (0.6, 1.4):
        ax.plot(rho, new[f], color=AFTER, lw=0.8, zorder=2)
    ax.plot(sturm_r, sturm, color=INK2, lw=1.2, ls=(0, (5, 3)), zorder=2)
    ax.plot(yen_r, yen, color=INK2, lw=1.2, ls=(0, (1, 2)), zorder=2)
    ax.axhline(kcap, color=INK, lw=0.8, ls=':', zorder=1)
    ax.scatter([RHO_ICE], [kcap], s=45, color=INK, edgecolor=SURFACE, lw=1.2, zorder=4)
    ax.text(105, kcap + 0.05, f'cap: ice, Yen (1981), {kcap:.2f} at 0 °C (2.32 at −20 °C)', color=INK,
            fontsize=8.5, va='bottom')
    ax.text(470, 1.42, 'old: linear ice/air mix', color=INK, fontsize=9, rotation=0, ha='left')
    ax.text(700, 0.66, 'Calonne (2011) × f_ins\nshaded: f_ins 0.6 – 1.4', color=INK, fontsize=9, ha='left')
    ax.annotate('Sturm (1997)', (600, sturm[-1]), xytext=(640, 0.3), color=INK2, fontsize=8.5,
                va='center', arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    ax.annotate('Yen (1981) snow', (600, yen[-1]), xytext=(440, 1.05), color=INK2, fontsize=8.5,
                va='center', ha='right', arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    # values at 250 kg m-3
    r0 = 250.
    ax.axvline(r0, color=INK2, lw=0.8, ls=(0, (2, 2)), zorder=1)
    k_old0 = (1 - r0 / 1000.) * KAIR + r0 / 1000. * KICE_OLD
    ax.annotate(f'{k_old0:.3f}', (r0, k_old0), xytext=(-6, 4), textcoords='offset points', ha='right',
                color=INK, fontsize=8.5)
    ax.text(108, 1.85, f'at 250 kg m⁻³:\nCalonne {calonne(r0):.3f}\nYen {2.22362 * 0.25 ** 1.885:.3f}\n'
            f'Sturm {0.138 - 1.01 * 0.25 + 3.233 * 0.0625:.3f}', color=INK, fontsize=8.5, va='top')
    ax.set_xlim(100, 930)
    ax.set_ylim(0, 2.45)
    ax.set_xlabel('density ρ [kg m⁻³]', color=INK)
    ax.set_ylabel('thermal conductivity k [W m⁻¹ K⁻¹]', color=INK)
    ax.set_title('a  Conductivity: the old mix was ~4× too high in snow', fontsize=10, color=INK, loc='left')

    kap_ice_new = kcap / heat_cap(RHO_ICE)
    kap_ice_old = ((1 - RHO_ICE / 1000.) * KAIR + RHO_ICE / 1000. * KICE_OLD) / heat_cap(RHO_ICE)
    ratio_old = old / heat_cap(rho) / kap_ice_old
    ratio = {f: new[f] / heat_cap(rho) / kap_ice_new for f in new}
    axd.fill_between(rho, ratio[0.6], ratio[1.4], color=AFTER_LIGHT, lw=0, zorder=1)
    axd.plot(rho, ratio_old, color=BEFORE, lw=2, zorder=3)
    axd.plot(rho, ratio[1.0], color=AFTER, lw=2, zorder=3)
    for f in (0.6, 1.4):
        axd.plot(rho, ratio[f], color=AFTER, lw=0.8, zorder=2)
    axd.axvline(r0, color=INK2, lw=0.8, ls=(0, (2, 2)), zorder=1)
    i0 = np.argmin(np.abs(rho - r0))
    axd.text(r0 + 10, ratio_old[i0] - 0.03, f'old: {ratio_old[i0]:.2f}× everywhere (density cancels)',
             color=INK, fontsize=8.5, va='top')
    axd.scatter([r0] * 3, [ratio[f][i0] for f in (0.6, 1.0, 1.4)], s=22, color=AFTER, edgecolor=SURFACE,
                lw=1, zorder=4)
    axd.text(262, 0.88, 'at 250 kg m⁻³:\n' + '\n'.join(f'f_ins {f:.1f}: {ratio[f][i0]:.2f}×' for f in (1.4, 1.0, 0.6)),
             color=INK, fontsize=8.5, va='top')
    axd.text(560, 0.42, 'after 1aef7f6\n(Calonne × f_ins)', color=INK, fontsize=9, ha='left')
    axd.set_xlim(100, 930)
    axd.set_ylim(0, 1.15)
    axd.set_xlabel('density ρ [kg m⁻³]', color=INK)
    axd.set_ylabel('diffusivity k/C, relative to ice [–]', color=INK)
    axd.set_title('b  How fast heat diffuses, relative to ice', fontsize=10, color=INK, loc='left')
    fig.text(0.01, 0.012, 'Constants from settings.pro and GloGEM commit 1aef7f6. C = linear ice/air mix '
             '(the same before and after); ice at 900 kg m⁻³ and 0 °C. Sturm and Yen snow drawn only\n'
             'inside their published ranges (156–600 and 80–600 kg m⁻³). The old ice conductivity was the '
             'fixed kice = 2.33 W m⁻¹ K⁻¹.', color=INK2, fontsize=8)
    fig.subplots_adjust(left=0.07, right=0.99, top=0.92, bottom=0.2)
    save(fig, 'ch02_conductivity')


if __name__ == '__main__':
    fig_glacier_optima()
    fig_conductivity()
    print('written:', sorted(p.name for p in OUT.glob('ch02_*.svg')))
