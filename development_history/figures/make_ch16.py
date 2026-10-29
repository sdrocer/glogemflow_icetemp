"""Figures for chapter 16 (Step 5, the calibration).
Run from this folder: python make_ch16.py [preview_dir]

All numbers come from data/bayescal/step5_parts5 (the overnight run of 28 Oct on the rules commit 3467851):
ch16_parameters: left, the 60 design runs, the stage-2 box, the stage-1 (40 runs) and final (60 runs) results as 90 %
  regions, the seven left-out-massif results, the grid choice and the starting values; right, a zoom on the result
  with the variants and the result without each decisive site (fit/, fit40/, design.csv, benchmark/folds.json,
  fit/variants.json, fit/site_checks.csv).
ch16_split: the error at depth of each massif against f_ref, from the real grid runs at f_ins 1.017
  (benchmark/benchmark_points.csv).
ch16_benchmark: the benchmark, each massif predicted from the other six: error at depth overall, in the zone and above
  3800 m, and the polythermal bands right (benchmark/benchmark.json).
ch16_bands: measured against predicted, per band at depth, each massif from its left-out fit (real runs;
  benchmark/benchmark_points.csv).
ch16_coverage: the share of the left-out readings inside their 90 % range for four parameter spreads
  (benchmark/coverage_spread.json).
"""
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Ellipse, Rectangle

OUT = Path(__file__).parent
PREVIEW = Path(sys.argv[1]) if len(sys.argv) > 1 else None
R = OUT.parents[1] / 'data' / 'bayescal' / 'step5_parts5'
plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none'})
SURFACE, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
BLUE, ORANGE, NEUTRAL, LIGHT, GREEN = '#2a78d6', '#eb6834', '#9a9893', '#9cc3ef', '#2f9e6b'
MASSIF_COL = {'Mont Blanc': ORANGE, 'Monte Rosa': BLUE}
CHI2_90 = 4.605   # chi-square with 2 degrees of freedom at 90 %
GRID_FREF = (0.03, 0.10, 0.25, 0.40, 0.55, 0.70, 0.85, 1.00)
START, GRID_CHOICE = (0.55, 1.0), (0.85, 1.017)
SHORT = {'Valais west and Diablerets': 'Valais west', 'Mischabel-Saas': 'Mischabel-Saas'}


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


def region(ax, g, **kw):
    """The 90 % region of a normal with the grid posterior's mean, sd and correlation."""
    sd, r = np.asarray(g['sd']), np.asarray(g['corr'])[0, 1]
    cov = np.array([[sd[0] ** 2, r * sd[0] * sd[1]], [r * sd[0] * sd[1], sd[1] ** 2]])
    val, vec = np.linalg.eigh(cov)
    ang = np.degrees(np.arctan2(vec[1, 1], vec[0, 1]))
    w, h = 2 * np.sqrt(CHI2_90 * val[::-1])
    ax.add_patch(Ellipse(g['mean'][:2], w, h, angle=ang, **kw))


def m1(df, c):
    g = ((df[c] - df.obs) ** 2).groupby(df.rgi)
    mse = g.mean()[g.size() >= 3]
    return float(np.sqrt(mse.mean()))


def points():
    return pd.read_csv(R / 'benchmark' / 'benchmark_points.csv', dtype={'rgi': str, 'hole': str})


def parameters():
    F, F40 = json.load(open(R / 'fit' / 'fit.json')), json.load(open(R / 'fit40' / 'fit.json'))
    D = pd.read_csv(R / 'design.csv')
    folds = json.load(open(R / 'benchmark' / 'folds.json'))['folds']
    V = json.load(open(R / 'fit' / 'variants.json'))
    S = pd.read_csv(R / 'fit' / 'site_checks.csv')
    fig, axs = plt.subplots(1, 2, figsize=(10.0, 4.9), facecolor=SURFACE, gridspec_kw=dict(width_ratios=[1.15, 1]))

    ax = axs[0]
    style(ax)
    ax.scatter(D.refreeze_frac[:40], D.insul_scale[:40], s=22, facecolor='none', edgecolor=NEUTRAL, lw=1,
               label='first 40 design runs')
    ax.scatter(D.refreeze_frac[40:], D.insul_scale[40:], s=22, color=NEUTRAL, edgecolor=SURFACE, lw=0.8,
               label='20 runs near the stage-1 result')
    c, h = (0.9747, 1.2865), (0.4478, 0.2865)
    lo, hi = max(c[0] - h[0], 0.03), min(c[0] + h[0], 1.0)
    ax.add_patch(Rectangle((lo, c[1] - h[1]), hi - lo, 2 * h[1], fill=False, ec=INK2, ls='--', lw=1))
    ax.text(lo + 0.01, c[1] + h[1] + 0.03, 'stage-2 box', color=INK2, fontsize=8)
    region(ax, F40['grid'], fill=False, ec=INK2, lw=1.2, ls=':')
    region(ax, F['grid'], fc=LIGHT, ec=BLUE, lw=1.5)
    for f in folds:
        col = MASSIF_COL.get(f['massif'], INK2)
        ax.scatter(*f['theta'][:2], s=60, marker='o', color=col, edgecolor=SURFACE, lw=1.5, zorder=4)
    mr = next(f for f in folds if f['massif'] == 'Monte Rosa')
    ax.annotate('result without\nMonte Rosa', mr['theta'][:2], (0.22, 0.48), color=INK, fontsize=8,
                arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    ax.annotate('the other six left-out\nresults (Mont Blanc orange)', (0.925, 1.21), (0.60, 1.82), color=INK,
                fontsize=8, arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    ax.scatter(*GRID_CHOICE, s=60, marker='s', color=INK, edgecolor=SURFACE, lw=1.5, zorder=5)
    ax.annotate('grid choice', GRID_CHOICE, (0.62, 0.78), color=INK, fontsize=8,
                arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    ax.scatter(*START, s=70, marker='^', color=INK, edgecolor=SURFACE, lw=1.5, zorder=5)
    ax.annotate('starting values', START, (0.30, 1.10), color=INK, fontsize=8,
                arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    ax.text(0.985, 1.66, 'stage 1\n(40 runs)', color=INK2, fontsize=8, ha='right')
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0.3, 2.0)
    ax.set_xlabel('f_ref (share of the meltwater that refreezes and heats)')
    ax.set_ylabel('f_ins (scale of the snow\'s insulation)')
    ax.set_title('the runs and where the fits landed', fontsize=10, color=INK)
    ax.legend(loc='upper left', frameon=False, fontsize=8)

    ax = axs[1]
    style(ax)
    region(ax, F40['grid'], fill=False, ec=INK2, lw=1.2, ls=':')
    region(ax, F['grid'], fc=LIGHT, ec=BLUE, lw=1.5)
    ax.text(0.958, 1.45, 'stage 1', color=INK2, fontsize=8)
    ax.text(0.903, 1.30, 'result\n(60 runs)', color=BLUE, fontsize=8, ha='right')
    for f in folds:
        if f['massif'] != 'Monte Rosa':
            ax.scatter(*f['theta'][:2], s=40, color=MASSIF_COL.get(f['massif'], INK2), edgecolor=SURFACE, lw=1.2,
                       zorder=4)
    labels = dict(deep='deep readings only', info_cap='no site above 20 %', eq_pub='equilibrium and published',
                  no_1999='without the 1999 bands')
    offs = dict(deep=(-0.075, -0.02), info_cap=(-0.06, -0.06), eq_pub=(0.012, 0.0), no_1999=(-0.075, -0.07))
    for k, lab in labels.items():
        m = V[k]['mean']
        ax.scatter(*m, s=60, marker='D', facecolor=SURFACE, edgecolor=GREEN, lw=1.6, zorder=5)
        ax.annotate(lab, m, (m[0] + offs[k][0], m[1] + offs[k][1]), color=INK, fontsize=8, ha='left',
                    arrowprops=dict(arrowstyle='-', color=INK2, lw=0.6))
    names = {'Breithornplateau 3810-3810': 'without Breithornplateau',
             'Glacier des Bossons 3930-4304': 'without Bossons 3930–4304 m',
             'Grenzgletscher 3898-4485': 'without the Colle Gnifetti saddles'}
    offs = {'Breithornplateau 3810-3810': (-0.06, 0.12), 'Glacier des Bossons 3930-4304': (-0.02, -0.13),
            'Grenzgletscher 3898-4485': (-0.005, -0.06)}
    for site, lab in names.items():
        r = S[S.site == site].iloc[0]
        p = (r.without_refreeze_frac, r.without_insul_scale)
        ax.scatter(*p, s=60, marker='X', color=INK, edgecolor=SURFACE, lw=1.2, zorder=5)
        ax.annotate(lab, p, (p[0] + offs[site][0], p[1] + offs[site][1]), color=INK, fontsize=8,
                    arrowprops=dict(arrowstyle='-', color=INK2, lw=0.6))
    ax.scatter(*GRID_CHOICE, s=60, marker='s', color=INK, edgecolor=SURFACE, lw=1.5, zorder=5)
    ax.text(GRID_CHOICE[0] + 0.006, GRID_CHOICE[1] - 0.035, 'grid choice', color=INK, fontsize=8)
    ax.scatter([], [], s=60, marker='D', facecolor=SURFACE, edgecolor=GREEN, lw=1.6, label='variants (rule 5)')
    ax.scatter([], [], s=60, marker='X', color=INK, label='without a decisive site (check 3)')
    ax.scatter([], [], s=40, color=INK2, label='left-out-massif results')
    ax.set_xlim(0.78, 1.0)
    ax.set_ylim(0.7, 1.6)
    ax.set_xlabel('f_ref')
    ax.set_ylabel('f_ins')
    ax.set_title('zoom: what moves the result', fontsize=10, color=INK)
    ax.legend(loc='upper left', frameon=False, fontsize=8)
    fig.tight_layout()
    save(fig, 'ch16_parameters')


def split():
    P = points()
    fig, ax = plt.subplots(figsize=(7.2, 4.6), facecolor=SURFACE)
    style(ax)
    x = np.array(GRID_FREF)
    order = sorted(P.massif.unique(), key=lambda m: m in MASSIF_COL)
    ends = []
    for m in order:
        d = P[P.massif == m]
        y = np.array([m1(d, f'G@{fr:.2f}') for fr in GRID_FREF])
        col = MASSIF_COL.get(m, NEUTRAL)
        lw = 2.2 if m in MASSIF_COL else 1.2
        ax.plot(x, y, color=col, lw=lw, marker='o', ms=4, zorder=3 if m in MASSIF_COL else 2)
        k = int(np.argmin(y))
        ax.scatter(x[k], y[k], s=70, color=col, edgecolor=SURFACE, lw=1.5, zorder=4)
        ends.append([y[-1], SHORT.get(m, m), INK if m in MASSIF_COL else INK2, 'normal'])
    y = np.array([m1(P, f'G@{fr:.2f}') for fr in GRID_FREF])
    ax.plot(x, y, color=INK, lw=2.6, marker='o', ms=4, zorder=5)
    k = int(np.argmin(y))
    ax.scatter(x[k], y[k], s=80, color=INK, edgecolor=SURFACE, lw=1.5, zorder=6)
    ends.append([y[-1], 'all 7 massifs', INK, 'bold'])
    ends.sort(key=lambda e: e[0])
    for i in range(1, len(ends)):   # labels at least 0.3 °C apart, pushed upwards
        ends[i][0] = max(ends[i][0], ends[i - 1][0] + 0.3)
    for yl, lab, col, wt in ends:
        ax.text(1.015, yl, lab, color=col, fontsize=8, va='center', weight=wt)
    ax.annotate('Mont Blanc\'s best: 0.25', (0.25, m1(P[P.massif == 'Mont Blanc'], 'G@0.25')), (0.27, 0.95),
                color=INK, fontsize=8, arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    ax.annotate('Monte Rosa\'s best: 1.0', (1.0, m1(P[P.massif == 'Monte Rosa'], 'G@1.00')), (0.60, 1.05),
                color=INK, fontsize=8, arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    ax.set_xlim(0, 1.2)
    ax.set_xticks(x)
    ax.set_xticklabels([f'{v:.2f}' for v in x])
    ax.set_ylim(0, 8)
    ax.set_xlabel('f_ref (real GloGEM runs, f_ins 1.017)')
    ax.set_ylabel('error at depth (°C, below 15 m)')
    ax.set_title('each massif\'s error against f_ref: the dot marks its best value', fontsize=10, color=INK)
    fig.tight_layout()
    save(fig, 'ch16_split')


def benchmark():
    B = json.load(open(R / 'benchmark' / 'benchmark.json'))['methods']
    rows = [('step5', 'Step 5', BLUE), ('grid', 'grid choice', INK), ('start', 'starting values', GREEN),
            ('nearest5', '5 nearest bands', NEUTRAL), ('elev_line', 'line in elevation', NEUTRAL),
            ('air_line', 'line in the free air', NEUTRAL)]
    panels = [('M1', 'error at depth, all (°C)', (0, 5)), ('zone_M1', 'in the zone 3000–3800 m (°C)', (0, 5)),
              ('above_M1', 'above 3800 m (°C)', (0, 5)), ('poly', 'polythermal bands right (of 18)', (0, 18))]
    fig, axs = plt.subplots(1, 4, figsize=(11.0, 3.4), facecolor=SURFACE, sharey=True)
    yy = np.arange(len(rows))[::-1]
    for ax, (key, title, lim) in zip(axs, panels):
        style(ax)
        ax.grid(axis='y', visible=False)
        for y, (m, lab, col) in zip(yy, rows):
            v = int(B[m]['classes']['counts']['polythermal'].split()[0]) if key == 'poly' else B[m][key]
            ax.plot([lim[0], v], [y, y], color=GRID, lw=1, zorder=1)
            ax.scatter(v, y, s=70, color=col, edgecolor=SURFACE, lw=1.5, zorder=3)
            ax.text(v + (lim[1] - lim[0]) * 0.04, y, f'{v:d}' if key == 'poly' else f'{v:.2f}', color=INK,
                    fontsize=8, va='center')
        ax.set_xlim(*lim)
        ax.set_title(title, fontsize=9, color=INK)
    axs[0].set_yticks(yy)
    axs[0].set_yticklabels([r[1] for r in rows])
    fig.tight_layout()
    save(fig, 'ch16_benchmark')


def bands():
    P = points()
    fig, ax = plt.subplots(figsize=(6.4, 6.0), facecolor=SURFACE)
    style(ax)
    t = np.array([-17, 1])
    ax.fill_between(t, t - 1.7, t + 1.7, color=GRID, alpha=0.7, lw=0, zorder=0)
    ax.plot(t, t, color=INK2, lw=1, zorder=1)
    for m in sorted(P.massif.unique(), key=lambda m: m in MASSIF_COL):
        d = P[P.massif == m]
        col = MASSIF_COL.get(m, NEUTRAL)
        ax.scatter(d.obs, d.step5, s=16, color=col, edgecolor=SURFACE, lw=0.5, zorder=3 if m in MASSIF_COL else 2,
                   label=m if m in MASSIF_COL else None)
    ax.scatter([], [], s=16, color=NEUTRAL, label='the other five massifs')
    notes = {'Taconnaz Glacier 3680-3830': ('Taconnaz 3680–3830 m: 6.2 °C too warm\n(the model at or near 0 °C)',
                                            (-16.5, -1.6)),
             'Glacier des Bossons 3930-4304': ('Bossons 3930–4304 m:\n3.1 °C too warm', (-16.5, -5.6)),
             'Fieschergletscher 3889-3893': ('Fiescher:\n2.7 °C too warm', (-11.0, -3.4)),
             'Breithornplateau 3810-3810': ('Breithornplateau:\n4.0 °C too cold', (-2.6, -10.5))}
    for site, (txt, xy) in notes.items():
        d = P[P.site == site]
        k = ((d.obs - d.obs.median()) ** 2 + (d.step5 - d.step5.median()) ** 2).idxmin()   # a real point near the middle
        ax.annotate(txt, (d.obs[k], d.step5[k]), xy, color=INK, fontsize=8,
                    arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    ax.text(-16.8, -16.3, '±1.7 °C: the best a band model can do (ch. 11)', color=INK2, fontsize=8)
    ax.set_xlim(-17, 1)
    ax.set_ylim(-17, 1)
    ax.set_aspect('equal')
    ax.set_xlabel('measured (°C)')
    ax.set_ylabel('Step 5, predicted from the other six massifs (°C)')
    ax.set_title('each band at each depth below 15 m, each hole and year', fontsize=10, color=INK)
    ax.legend(loc='lower right', frameon=False, fontsize=8)
    fig.tight_layout()
    save(fig, 'ch16_bands')


def coverage():
    C = json.load(open(R / 'benchmark' / 'coverage_spread.json'))
    rows = [('rule', 'each left-out fit\'s own (±0.02 / ±0.06):\nrule 9 as coded'),
            ('fold_sd', 'the sd of the 7 left-out results\n(±0.20 / ±0.20)'),
            ('map', 'the map\'s spread of rule 10\n(±0.20 / ±0.37)'),
            ('single_massif', 'one massif\'s own spread,\n(n − 1) × sd (±1.2 / ±1.2)')]
    fig, ax = plt.subplots(figsize=(8.6, 3.6), facecolor=SURFACE)
    style(ax)
    ax.grid(axis='y', visible=False)
    ax.axvspan(80, 97, color=GRID, alpha=0.8, lw=0, zorder=0)
    ax.axvline(90, color=INK2, lw=1, ls='--', zorder=1)
    ax.text(80.5, len(rows) - 0.45, 'rule 9 asks for 80–97 %', color=INK2, fontsize=8)
    yy = np.arange(len(rows))[::-1]
    for y, (k, lab) in zip(yy, rows):
        v = 100 * C['overall'][k]
        ok = 80 <= v <= 97
        ax.barh(y, v, height=0.5, color=BLUE if ok else NEUTRAL, zorder=2)
        s = C['side'][k]
        ax.text(v + 1.0, y, f'{v:.0f} % {"inside the window" if ok else "fails"}\n(model too warm '
                f'{100 * s["below"]:.0f} %, too cold {100 * s["above"]:.0f} %)', color=INK, fontsize=8, va='center')
    ax.set_yticks(yy)
    ax.set_yticklabels([r[1] for r in rows], fontsize=8)
    ax.set_xlim(0, 140)
    ax.set_xticks([0, 20, 40, 60, 80, 90, 100])
    ax.set_xlabel('left-out readings inside their 90 % range (%)')
    ax.set_title('spread given to f_ref / f_ins', fontsize=10, color=INK, loc='left')
    fig.tight_layout()
    save(fig, 'ch16_coverage')


if __name__ == '__main__':
    parameters()
    split()
    benchmark()
    bands()
    coverage()
