"""Figures for chapter 01 (the starting point).
Run from this folder: python make_ch01.py [preview_dir]

ch01_cart: the depth-2 regression tree (CART) for the firn warming offset dT_firn, and its training
vs 5-fold cross-validated r², with the one-variable linear fit on T_amplitude for comparison.
Recomputed as in notebooks/archive/02_glenglat_ML.ipynb (cells 3, 10, 11; scikit-learn, random_state
42) from data/glenglat_profiles_derived.csv. That file now holds one profile more than when the
notebook ran (209 instead of 208 after the filter): borehole 678, profile 3. Without it every number
the notebook printed comes back exactly (training r² 0.707, CV r² 0.104 ± 0.797, linear CV r²
0.332 ± 0.040, leaves 0.54 / 4.40 / 7.26 / 10.34 °C on 31 / 82 / 59 / 36 profiles), so it is left
out here; the script checks this. The rounded thresholds are the ones GloGEM used (commit 035c46e).
ch01_paramspace: balanced score of the 25 x 25 (perm_frac, dT_scale) sweep through campaign 7's
emulator (data/figures/sweep25.csv, scripts/figures/sweep25.py), the settings where recall and
specificity both exceed 50 %, campaign 8's posterior samples (centraleurope_covfix) and campaign 7's
posterior mean (centraleurope_elevsplit); the book version of figs/paramspace_target_vs_posterior.png
(scripts/figures/fig_paramspace.py).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch

OUT = Path(__file__).parent
PREVIEW = Path(sys.argv[1]) if len(sys.argv) > 1 else None
REPO = Path(__file__).resolve().parents[2]
RUNS = REPO / 'data' / 'bayescal'
plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none'})

SURFACE, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
AFTER, BEFORE = '#2a78d6', '#eb6834'   # categorical slots 1-2, fixed order; corrected state first
NEUTRAL = '#bdbbb5'

FEATURES = ['T_amplitude', 'T_maat', 'PDD', 'elevation', 'latitude']
LATER_PROFILE = (678, 3)   # (borehole_id, profile_id) added to the csv after the notebook ran


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


def cart_data():
    """The notebook's 208 filtered profiles: firn/equilibrium zone, shallow T at or below -0.5 °C."""
    p = pd.read_csv(REPO / 'data' / 'glenglat_profiles_derived.csv')
    p = p.dropna(subset=['T_amplitude', 'dT_firn'])
    p = p[(p.T_shallow <= -0.5) & (p.mass_balance_area != 'ablation')]
    later = (p.borehole_id == LATER_PROFILE[0]) & (p.profile_id == LATER_PROFILE[1])
    d = p[~later].dropna(subset=FEATURES)
    assert len(d) == 208, len(d)
    return d


def fig_cart():
    from sklearn.tree import DecisionTreeRegressor
    from sklearn.linear_model import LinearRegression
    from sklearn.model_selection import cross_val_score

    d = cart_data()
    X, y = d[FEATURES].values, d.dT_firn.values
    tree = DecisionTreeRegressor(max_depth=2, random_state=42).fit(X, y)
    r2_train = tree.score(X, y)
    cv_tree = cross_val_score(DecisionTreeRegressor(max_depth=2, random_state=42), X, y, cv=5,
                              scoring='r2')
    lin = LinearRegression()
    cv_lin = cross_val_score(lin, d[['T_amplitude']].values, y, cv=5, scoring='r2')
    r2_lin_train = lin.fit(d[['T_amplitude']].values, y).score(d[['T_amplitude']].values, y)
    assert round(r2_train, 3) == 0.707 and round(cv_tree.mean(), 3) == 0.104
    assert round(cv_tree.std(), 3) == 0.797 and round(cv_lin.mean(), 3) == 0.332

    t = tree.tree_
    # node ids: 0 root (T_amplitude), 1 / 4 elevation splits, leaves 2, 3, 5, 6
    assert t.feature[0] == FEATURES.index('T_amplitude')
    assert t.feature[1] == t.feature[4] == FEATURES.index('elevation')
    leaf = {i: (t.value[i, 0, 0], t.n_node_samples[i]) for i in (2, 3, 5, 6)}

    fig = plt.figure(figsize=(10.5, 5.0), facecolor=SURFACE)
    ax = fig.add_axes([0.0, 0.12, 0.6, 0.8])
    ax.set_facecolor(SURFACE)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')

    def box(x, y, text, w=3.5, h=1.45, fc=SURFACE, ec=INK2, bold=False, size=9):
        ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle='round,pad=0.05,rounding_size=0.18',
                                    fc=fc, ec=ec, lw=1))
        ax.text(x, y, text, ha='center', va='center', color=INK, fontsize=size,
                fontweight='bold' if bold else 'normal', linespacing=1.35)

    def edge(x0, y0, x1, y1, label, side):
        ax.plot([x0, x1], [y0, y1], color=INK2, lw=1, zorder=0)
        xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
        ax.text(xm + (-0.15 if side == 'l' else 0.15), ym, label, ha='right' if side == 'l' else 'left',
                va='center', color=INK2, fontsize=8.5)

    box(5, 8.7, 'seasonal amplitude T_amplitude\n≤ 20 °C ?   (fitted 20.3 °C)', w=4.6)
    box(2.5, 5.6, 'elevation ≤ 4300 m ?\n(fitted 4296 m)', w=3.2)
    box(7.5, 5.6, 'elevation ≤ 1500 m ?\n(fitted 1492 m)', w=3.2)
    edge(4.0, 7.97, 2.9, 6.33, 'yes: "maritime"', 'l')
    edge(6.0, 7.97, 7.1, 6.33, 'no: "continental"', 'r')
    leaves = [(1.25, 2, 'yes'), (3.75, 3, 'no'), (6.25, 5, 'yes'), (8.75, 6, 'no')]
    for x, i, ans in leaves:
        v, n = leaf[i]
        parent = 2.5 if i in (2, 3) else 7.5
        ax.plot([parent + (-0.6 if ans == 'yes' else 0.6), x], [4.87, 3.15], color=INK2, lw=1, zorder=0)
        ax.text((parent + x) / 2 + (-0.35 if ans == 'yes' else 0.35), 4.05, ans,
                ha='right' if ans == 'yes' else 'left', va='center', color=INK2, fontsize=8.5)
        box(x, 2.35, f'+{v:.2f} °C\n{n} profiles', w=2.2, h=1.5, fc='#eef4fb', ec=AFTER, bold=False,
            size=9.5)
    ax.text(5, 0.85, 'leaf = offset added to the air temperature at the firn surface (ΔT_firn)',
            ha='center', va='center', color=INK2, fontsize=8.5)
    ax.text(0.15, 9.95, 'a  The tree GloGEM used (fitted to 208 glenglat profiles)', ha='left',
            va='top', color=INK, fontsize=10)

    # b: training vs cross-validated r²
    axr = fig.add_axes([0.66, 0.2, 0.32, 0.68])
    style(axr)
    axr.grid(axis='x', color=GRID, lw=0.8)
    axr.set_axisbelow(True)
    rows = [(1.0, 'depth-2 tree', cv_tree, r2_train), (0.0, 'linear fit on\nT_amplitude', cv_lin, r2_lin_train)]
    for yy, _, cv, tr in rows:
        axr.scatter(cv, np.full(5, yy) + np.linspace(-0.12, 0.12, 5), s=26, color=NEUTRAL,
                    edgecolor=SURFACE, lw=1, zorder=3)
        axr.errorbar(cv.mean(), yy - 0.3, xerr=cv.std(), fmt='o', color=AFTER, ms=7, mec=SURFACE,
                     mew=1.5, ecolor=AFTER, elinewidth=1.5, capsize=3, zorder=4)
        axr.scatter(tr, yy + 0.3, marker='D', s=45, color=BEFORE, edgecolor=SURFACE, lw=1.2, zorder=4)
    axr.text(r2_train, 1.47, f'training {r2_train:.3f}', ha='center', va='bottom', color=INK, fontsize=8.5)
    axr.text(cv_tree.mean(), 0.6, f'cross-validated {cv_tree.mean():.3f} ± {cv_tree.std():.3f}',
             ha='center', va='top', color=INK, fontsize=8.5)
    worst = cv_tree.min()
    axr.annotate('one test fold:\nr² = ' + f'{worst:.2f}'.replace('-', '−'), xy=(worst, 1.0 - 0.06),
                 xytext=(worst + 0.12, 1.32), color=INK, fontsize=8.5, va='center',
                 arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    axr.text(cv_lin.mean(), -0.4, f'{cv_lin.mean():.3f} ± {cv_lin.std():.3f}', ha='center', va='top',
             color=INK, fontsize=8.5)
    axr.axvline(0, color=INK2, lw=0.8)
    axr.set_yticks([r[0] for r in rows])
    axr.set_yticklabels([r[1] for r in rows], color=INK)
    axr.set_ylim(-1.3, 1.75)
    axr.set_xlim(-1.65, 1.0)
    axr.set_xlabel('r² (explained share of the variance) [–]', color=INK)
    handles = [plt.Line2D([], [], marker='D', ls='', ms=7, color=BEFORE, mec=SURFACE, label='on the training data'),
               plt.Line2D([], [], marker='o', ls='', ms=5, color=NEUTRAL, mec=SURFACE, label='each of the 5 test folds'),
               plt.Line2D([], [], marker='o', ls='-', ms=7, color=AFTER, mec=SURFACE, label='5-fold mean ± sd')]
    axr.legend(handles=handles, loc='lower left', fontsize=8, frameon=False, labelcolor=INK)
    axr.set_title('b  Skill on profiles it was not fitted to', fontsize=10, color=INK, loc='left')
    fig.text(0.01, 0.015, 'Recomputed as in notebooks/archive/02_glenglat_ML.ipynb from '
             'data/glenglat_profiles_derived.csv (firn and equilibrium zone, shallow temperature ≤ −0.5 °C; '
             '208 profiles).\nTree depth 2, 5-fold cross-validation (scikit-learn, random_state 42). '
             'Thresholds as implemented in GloGEM 035c46e.', color=INK2, fontsize=8)
    save(fig, 'ch01_cart')


def fig_paramspace():
    sw = pd.read_csv(REPO / 'data' / 'figures' / 'sweep25.csv')
    pf, dt = np.sort(sw.perm_frac.unique()), np.sort(sw.dT_scale.unique())
    piv = lambda v: sw.pivot_table(index='dT_scale', columns='perm_frac', values=v).values
    bal = piv('balanced')
    meets = (piv('recall') > 0.5) & (piv('specificity') > 0.5)
    best = sw[(sw.recall > 0.5) & (sw.specificity > 0.5)].nlargest(1, 'balanced').iloc[0]
    p8 = pd.read_csv(RUNS / 'centraleurope_covfix' / 'posterior_samples.csv')
    p7 = pd.read_csv(RUNS / 'centraleurope_elevsplit' / 'posterior_samples.csv').mean()
    cloud = p8.sample(n=2000, random_state=0)

    fig, ax = plt.subplots(figsize=(9.2, 6.2), facecolor=SURFACE)
    style(ax)
    cmap = LinearSegmentedColormap.from_list('grey', ['#fcfcfb', '#d9d8d3', '#a3a19b', '#6e6d68'])
    im = ax.pcolormesh(pf, dt, bal, cmap=cmap, shading='gouraud', vmin=0.5, vmax=0.75, rasterized=True)
    plt.rcParams['hatch.color'] = INK
    plt.rcParams['hatch.linewidth'] = 0.6
    ax.contourf(pf, dt, meets.astype(float), levels=[0.5, 1.5], colors='none', hatches=['////'], zorder=2)
    ax.contour(pf, dt, meets.astype(float), levels=[0.5], colors=[INK], linewidths=1.3, zorder=2.5)
    ax.scatter(cloud.perm_frac, cloud.dT_scale, s=9, color=AFTER, alpha=0.35, edgecolor='none', zorder=3,
               rasterized=True)
    ax.scatter([p8.perm_frac.mean()], [p8.dT_scale.mean()], s=120, color=AFTER, edgecolor=SURFACE, lw=2,
               zorder=5)
    ax.scatter([p7.perm_frac], [p7.dT_scale], s=120, marker='D', color=BEFORE, edgecolor=SURFACE, lw=2,
               zorder=5, clip_on=False)
    ax.scatter([best.perm_frac], [best.dT_scale], s=90, marker='*', color=INK, edgecolor=SURFACE, lw=1,
               zorder=5, clip_on=False)
    lab = dict(color=INK, fontsize=9, bbox=dict(boxstyle='round,pad=0.25', fc=SURFACE, ec='none', alpha=0.9))
    ax.annotate(f'campaign 8 (honest likelihood)\nposterior {p8.perm_frac.mean():.3f} ± {p8.perm_frac.std():.3f}, '
                f'{p8.dT_scale.mean():.3f} ± {p8.dT_scale.std():.3f}', (p8.perm_frac.mean(), p8.dT_scale.mean()),
                xytext=(0.36, 0.55), va='center', arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8), **lab)
    ax.annotate(f'campaign 7 (over-confident)\npinned at the top wall: {p7.perm_frac:.3f}, {p7.dT_scale:.2f}',
                (p7.perm_frac, p7.dT_scale), xytext=(0.5, 4.25), va='center', ha='left',
                arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8), **lab)
    ax.annotate(f'best setting meeting both:\nrecall {best.recall:.1%}, specificity {best.specificity:.1%}',
                (best.perm_frac, best.dT_scale), xytext=(0.06, 4.25), va='center', ha='left',
                arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8), **lab)
    ax.text(0.62, 2.55, f'hatched: recall and specificity\nboth above 50 % ({int(meets.sum())} of {meets.size} settings)',
            **lab)
    ax.set_xlim(pf.min(), pf.max())
    ax.set_ylim(dt.min(), dt.max())
    ax.set_xlabel('perm_frac (scales the percolation depth) [–]', color=INK)
    ax.set_ylabel('dT_scale (scales the decision-tree surface offset) [–]', color=INK)
    ax.set_title('The model can meet the classification bar; the temperature fit does not choose it',
                 fontsize=10, color=INK, loc='left')
    cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.04)
    cb.set_label('balanced score: (recall + specificity) / 2 [–]', color=INK)
    cb.ax.tick_params(colors=INK2)
    cb.outline.set_visible(False)
    fig.text(0.01, 0.012, 'Background: 25 × 25 sweep scored on the 107 entities through campaign 7\'s emulator '
             '(sweep25.csv, scripts/figures/sweep25.py).\nBlue dots: 2000 of campaign 8\'s posterior samples '
             '(centraleurope_covfix); diamond: campaign 7\'s posterior mean (centraleurope_elevsplit).',
             color=INK2, fontsize=8)
    fig.subplots_adjust(left=0.09, right=0.93, top=0.94, bottom=0.17)
    save(fig, 'ch01_paramspace')


if __name__ == '__main__':
    fig_cart()
    fig_paramspace()
    print('written:', sorted(p.name for p in OUT.glob('ch01_*.svg')))
