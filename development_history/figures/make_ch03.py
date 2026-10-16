"""Figures for chapter 03 (why the parameters look alike).
Run from this folder: python make_ch03.py [preview_dir]

ch03_depth_shapes: the temperature change per 1 sd of each parameter against depth, on three
Grenzgletscher entities of the 48-run sweep grenz_tuning3 (tongue, 2505 m band; mid, 4005 m;
Colle Gnifetti, 4455 m). As in notebooks/02_parameter_sensitivity.ipynb section C: per depth, the
modelled temperatures of all runs are regressed on the standardised parameters; depths unresolved
in any run are dropped. Computed here from data/bayescal/grenz_tuning3 (design.csv and the
training_runs/*.output.pkl files). f_ins is drawn with its sign flipped. The same code on all 24
entities, with the observation table of icetemp 43c7574, gives the chapter's medians (strength
1.250 / 1.034 / 0.131 °C per sd, rho(f_ref, f_ins) +0.98, cosine -0.98).
ch03_depth_window: per entity, the deepest modelled depth with a value against the deepest borehole
observation, from scripts/icetemp_analysis/depth_window.py on alps_sensitivity_w0/1/2 (first design
point), run with the code and observation table of icetemp 43c7574 (2026-09-22). The current
observation table assigns boreholes differently (icetemp cf140d0) and gives other counts, so the
per-entity rows of that run are written out below (DEPTH_WINDOW).
"""
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUT = Path(__file__).parent
PREVIEW = Path(sys.argv[1]) if len(sys.argv) > 1 else None
RUNS = Path(__file__).resolve().parents[2] / 'data' / 'bayescal'
plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none'})

SURFACE, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
AFTER, BEFORE = '#2a78d6', '#eb6834'   # categorical slots 1-2, fixed order
AQUA = '#1baf7a'                        # categorical slot 3
NEUTRAL = '#bdbbb5'

PARAMS = ('refreeze_frac', 'insul_scale', 'advection_scale')
ENTITIES = [('Grenzgletscher@2505', 'a  tongue, 2505 m band (ice)'),
            ('Grenzgletscher@4005', 'b  4005 m band (firn)'),
            ('Grenzgletscher@4455', 'c  Colle Gnifetti, 4455 m band (firn)')]

# depth_window.py at icetemp 43c7574 on alps_sensitivity_w0/1/2, first design point:
# (entity, deepest observation [m], deepest modelled depth with a value [m], observations, of them below)
DEPTH_WINDOW = [
    ('Glacier de Tête Rousse@3135', 32.5, 14, 9, 7),
    ('Glacier de Tête Rousse@3145', 46.0, 19, 9, 6),
    ('Glacier de Tête Rousse@3155', 46.2, 19, 12, 9),
    ('Glacier de Tête Rousse@3165', 65.0, 24, 15, 12),
    ('Glacier de Tête Rousse@3175', 69.2, 24, 26, 21),
    ('Glacier de Tête Rousse@3185', 66.0, 19, 15, 13),
    ('Glacier de Tête Rousse@3195', 35.7, 14, 5, 5),
    ('Glacier de Tête Rousse@3205', 70.0, 14, 15, 15),
    ('Glacier des Bossons@3935', 18.0, 59, 1, 0),
    ('Glacier des Bossons@4005', 18.0, 59, 1, 0),
    ('Glacier des Bossons@4015', 18.0, 59, 1, 0),
    ('Glacier des Bossons@4125', 18.0, 59, 1, 0),
    ('Glacier des Bossons@4135', 18.0, 79, 1, 0),
    ('Glacier des Bossons@4195', 18.0, 54, 1, 0),
    ('Glacier des Bossons@4215', 18.0, 49, 1, 0),
    ('Glacier des Bossons@4245', 125.6, 44, 44, 29),
    ('Glacier des Bossons@4255', 103.1, 54, 25, 15),
    ('Glacier des Bossons@4265', 18.0, 44, 1, 0),
    ('Glacier des Bossons@4275', 20.0, 44, 4, 0),
    ('Glacier des Bossons@4285', 18.0, 39, 5, 0),
    ('Glacier des Bossons@4295', 18.0, 49, 1, 0),
    ('Glacier des Bossons@4305', 39.8, 49, 36, 0),
    ('Glacier des Bossons@4515', 12.0, 44, 1, 0),
    ('Gornergletscher@2525', 292.8, 179, 9, 3),
    ('Gornergletscher@2545', 318.4, 139, 9, 5),
    ('Gornergletscher@2555', 347.7, 139, 30, 11),
    ('Grenzgletscher@2015', 267.3, 79, 10, 7),
    ('Grenzgletscher@2505', 106.5, 159, 12, 0),
    ('Grenzgletscher@2525', 202.1, 179, 14, 4),
    ('Grenzgletscher@2605', 347.9, 79, 16, 11),
    ('Grenzgletscher@3895', 18.0, 54, 9, 0),
    ('Grenzgletscher@3955', 21.0, 59, 10, 0),
    ('Grenzgletscher@3985', 20.0, 49, 10, 0),
    ('Grenzgletscher@3995', 22.0, 49, 11, 0),
    ('Grenzgletscher@4005', 14.0, 49, 7, 0),
    ('Grenzgletscher@4095', 21.0, 54, 11, 0),
    ('Grenzgletscher@4105', 14.0, 54, 9, 0),
    ('Grenzgletscher@4155', 22.0, 44, 11, 0),
    ('Grenzgletscher@4205', 118.0, 44, 17, 4),
    ('Grenzgletscher@4255', 30.0, 44, 15, 0),
    ('Grenzgletscher@4285', 21.0, 34, 11, 0),
    ('Grenzgletscher@4295', 34.6, 44, 36, 0),
    ('Grenzgletscher@4305', 14.0, 44, 9, 0),
    ('Grenzgletscher@4335', 22.0, 39, 11, 0),
    ('Grenzgletscher@4405', 30.0, 34, 15, 0),
    ('Grenzgletscher@4445', 101.0, 39, 8, 6),
    ('Grenzgletscher@4455', 80.0, 44, 63, 22),
    ('Grenzgletscher@4465', 76.0, 34, 19, 14),
    ('Grenzgletscher@4475', 61.0, 34, 24, 6),
    ('Grenzgletscher@4485', 26.0, 24, 18, 1),
    ('Hintereisferner@2555', 11.1, 59, 12, 0),
    ('Hintereisferner@2615', 12.5, 59, 14, 0),
    ('Hintereisferner@2805', 10.0, 54, 5, 0),
    ('Hintereisferner@2955', 14.6, 39, 12, 0),
    ('Hintereisferner@3065', 14.8, 24, 8, 0),
    ('Hintereisferner@3155', 14.8, 34, 8, 0),
    ('Hohsaasgletscher@3275', 21.5, 44, 5, 0),
    ('Hohsaasgletscher@3285', 15.6, 29, 2, 0),
    ('Hohsaasgletscher@3295', 21.5, 34, 5, 0),
    ('Hohsaasgletscher@3335', 15.4, 39, 2, 0),
    ('Hohsaasgletscher@3345', 29.0, 39, 8, 0),
    ('Vadret dal Corvatsch@3295', 16.0, 9, 8, 4),
    ('Vadret dal Corvatsch@3315', 16.0, 9, 7, 4),
    ('Vadret dal Corvatsch@3325', 16.0, 9, 8, 4),
    ('Vadret dal Corvatsch@3345', 13.0, 14, 11, 0),
    ('Vadret dal Corvatsch@3365', 14.0, 14, 7, 0),
]


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


def depth_response(run, entity):
    """Notebook 02, sections B-C: per depth, d T / d (1 sd of each parameter), by least squares."""
    design = pd.read_csv(RUNS / run / 'design.csv')[list(PARAMS)]
    X, Y, depths = [], [], None
    for f in sorted((RUNS / run / 'training_runs').glob('design_*.output.pkl')):
        got = pickle.load(open(f, 'rb')).get(entity)
        if got is None:
            continue
        i = int(f.stem.split('_')[1].split('.')[0])
        depths = np.asarray(got[0], float) if depths is None else depths
        X.append(design.iloc[i].values)
        Y.append(np.asarray(got[1], float))
    X, Y = np.asarray(X), np.asarray(Y)
    keep = np.all(np.isfinite(Y), axis=0)
    Xs = (X - X.mean(axis=0)) / X.std(axis=0)
    coef, *_ = np.linalg.lstsq(np.column_stack([np.ones(len(Xs)), Xs]), Y[:, keep], rcond=None)
    return depths[keep], coef[1:].T, len(X)


def fig_depth_shapes():
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 8.6), facecolor=SURFACE,
                             gridspec_kw=dict(wspace=0.32, hspace=0.5))
    n_runs = set()
    curves = (('f_ref', 1, AFTER, '-', 'o'), ('f_ins', -1, BEFORE, (0, (4, 2)), 's'),
              ('f_adv', 1, AQUA, '-', '^'))
    for k, (entity, title) in enumerate(ENTITIES):
        z, J, n = depth_response('grenz_tuning3', entity)
        n_runs.add(n)
        cos = J[:, 0] @ J[:, 1] / np.linalg.norm(J[:, 0]) / np.linalg.norm(J[:, 1])
        for row, ax in enumerate(axes[:, k]):
            style(ax)
            ax.grid(color=GRID, lw=0.8)
            ax.set_axisbelow(True)
            ax.axvline(0, color=INK2, lw=0.8)
            xs = []
            for j, (_, sign, color, ls, mk) in enumerate(curves):
                x = sign * J[:, j]
                if row == 1:
                    x = x / np.linalg.norm(x)
                xs.append(x)
                ax.plot(x, z, color=color, lw=2, ls=ls, marker=mk, ms=4, mec=SURFACE, mew=1, zorder=4 - j)
            xs = np.concatenate(xs + [[0]])
            pad = 0.1 * (xs.max() - xs.min())
            ax.set_xlim(xs.min() - pad, xs.max() + pad)
            ax.set_ylim(z.max() * 1.04, 0)
            if row == 0:
                ax.set_title(title, fontsize=10, color=INK, loc='left')
                ax.set_xlabel('temperature change [°C per 1 sd]', color=INK)
            else:
                ax.set_xlabel('the same, each curve ÷ its length [–]', color=INK)
                ax.text(0.97, 0.03, f'cosine(f_ref, f_ins)\n= {cos:.2f}'.replace('-', '−'),
                        transform=ax.transAxes, ha='right', va='bottom', color=INK, fontsize=8.5,
                        bbox=dict(fc=SURFACE, ec='none', pad=2))
    axes[0, 0].set_ylabel('depth [m]', color=INK)
    axes[1, 0].set_ylabel('depth [m]', color=INK)
    fig.text(0.07, 0.925, 'Size and shape: how much 1 sd of each parameter changes the temperature at each depth',
             color=INK, fontsize=10, fontweight='bold')
    fig.text(0.07, 0.475, 'Shape only (what the cosine compares): f_ref and −f_ins nearly coincide on the firn bands',
             color=INK, fontsize=10, fontweight='bold')
    handles = [plt.Line2D([], [], color=AFTER, lw=2, marker='o', ms=5, mec=SURFACE, label='f_ref (refreezing)'),
               plt.Line2D([], [], color=BEFORE, lw=2, ls=(0, (4, 2)), marker='s', ms=5, mec=SURFACE,
                          label='− f_ins (insulation, sign flipped)'),
               plt.Line2D([], [], color=AQUA, lw=2, marker='^', ms=5, mec=SURFACE, label='f_adv (advection)')]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.995), ncol=3, fontsize=9,
               frameon=False, labelcolor=INK)
    fig.text(0.01, 0.01, f'Grenzgletscher sweep grenz_tuning3: {n_runs.pop()} runs, Latin hypercube over '
             '(f_ref, f_ins, f_adv), thermal spin-up on. Per depth, the runs\' temperatures\nregressed on the '
             'standardised parameters (notebooks/02_parameter_sensitivity.ipynb, section C). Each panel has its own '
             'x scale: on the tongue\nall three effects stay below 0.25 °C per sd. Output nodes every 1 m to 9 m, '
             'then every 5 m, so the top 10 m weigh heavily in the cosine.', color=INK2, fontsize=8)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.88, bottom=0.135)
    save(fig, 'ch03_depth_shapes')


def fig_depth_window():
    d = pd.DataFrame(DEPTH_WINDOW, columns=['entity', 'z_obs', 'z_model', 'n_obs', 'n_lost'])
    d['glacier'] = d.entity.str.split('@').str[0]
    n_obs, n_lost = d.n_obs.sum(), d.n_lost.sum()
    outrun = d.n_lost > 0
    assert (n_obs, n_lost, outrun.sum(), len(d)) == (776, 238, 25, 66)
    assert d.z_model.median() == 44 and d.z_obs.max() == 347.9

    fig, ax = plt.subplots(figsize=(8.6, 6.4), facecolor=SURFACE)
    style(ax)
    ax.set_xscale('log')
    ax.set_yscale('log')
    lim = (7, 420)
    ax.fill_between(lim, [lim[0], lim[0]], lim, color='#fbe3d8', lw=0, zorder=0)
    ax.plot(lim, lim, color=INK2, lw=1, zorder=1)
    ax.text(300, 255, '1 : 1', color=INK2, fontsize=8.5, ha='right', va='bottom', rotation=45)
    ax.text(330, 9.5, 'borehole deeper than the model column:\nthe readings below it are dropped',
            color=INK, fontsize=9, ha='right', va='bottom')
    size = lambda n: 12 + 5 * np.asarray(n, float)
    ax.scatter(d.z_obs[~outrun], d.z_model[~outrun], s=size(d.n_obs[~outrun]), color=NEUTRAL, alpha=0.8,
               edgecolor=SURFACE, lw=1.2, zorder=2)
    ax.scatter(d.z_obs[outrun], d.z_model[outrun], s=size(d.n_obs[outrun]), color=BEFORE, alpha=0.85,
               edgecolor=SURFACE, lw=1.2, zorder=3)
    tr = d[d.glacier == 'Glacier de Tête Rousse']
    ax.annotate(f'Glacier de Tête Rousse (8 bands): model {tr.z_model.min()}–{tr.z_model.max()} m,\n'
                f'boreholes to {tr.z_obs.max():.0f} m; {tr.n_lost.sum()} of {tr.n_obs.sum()} readings dropped',
                (66, 19), xytext=(85, 17.5), color=INK, fontsize=9, va='center',
                arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    g = d[d.entity == 'Grenzgletscher@2605'].iloc[0]
    ax.annotate(f'Grenzgletscher tongue (2605 m band):\ndeepest borehole {g.z_obs:.0f} m, model '
                f'{g.z_model} m;\n{g.n_lost} of {g.n_obs} readings dropped', (g.z_obs, g.z_model),
                xytext=(95, 30), color=INK, fontsize=9, va='center', ha='left',
                arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    ax.axhline(d.z_model.median(), color=INK2, lw=0.8, ls=(0, (2, 3)), zorder=1)
    ax.text(400, d.z_model.median() * 0.96, f'median {d.z_model.median():.0f} m', color=INK2,
            fontsize=8.5, va='top', ha='right')
    ax.text(8, 340, f'{n_obs - n_lost} of {n_obs} readings ({(n_obs - n_lost) / n_obs:.0%}) inside the '
            f'model column,\n{n_lost} below it, in {outrun.sum()} of {len(d)} entities', color=INK,
            fontsize=9.5, ha='left', va='top')
    ticks = [10, 20, 50, 100, 200, 400]
    for setter in (ax.set_xticks, ax.set_yticks):
        setter(ticks)
    ax.set_xticklabels([str(t) for t in ticks], color=INK2)
    ax.set_yticklabels([str(t) for t in ticks], color=INK2)
    ax.minorticks_off()
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_aspect('equal')
    ax.set_xlabel('deepest borehole reading of the entity [m]', color=INK)
    ax.set_ylabel('deepest modelled depth with a value [m]', color=INK)
    ax.set_title('Seven Alpine glaciers: where the model column stops short of the boreholes',
                 fontsize=10, color=INK, loc='left')
    handles = [plt.Line2D([], [], marker='o', ls='', ms=9, color=BEFORE, mec=SURFACE, label='loses readings'),
               plt.Line2D([], [], marker='o', ls='', ms=9, color=NEUTRAL, mec=SURFACE, label='keeps all'),
               plt.Line2D([], [], marker='o', ls='', ms=4, color=INK2, mec=SURFACE,
                          label='marker area ∝ number of readings')]
    ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.0, 0.86), fontsize=8.5, frameon=False,
              labelcolor=INK)
    fig.text(0.01, 0.012, 'One point per entity (glacier and 10 m band), 66 entities. scripts/icetemp_analysis/'
             'depth_window.py on alps_sensitivity_w0/1/2, first design point,\nwith the observation table '
             'of icetemp 43c7574 (2026-09-22). Modelled depths are the output nodes (…, 9, 14, 19, … m).',
             color=INK2, fontsize=8)
    fig.subplots_adjust(left=0.1, right=0.98, top=0.95, bottom=0.13)
    save(fig, 'ch03_depth_window')


if __name__ == '__main__':
    fig_depth_shapes()
    fig_depth_window()
    print('written:', sorted(p.name for p in OUT.glob('ch03_*.svg')))
