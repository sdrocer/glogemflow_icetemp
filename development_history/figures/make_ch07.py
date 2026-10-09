"""Figures for chapter 07 (the snow, the comparison and the calibration target).
Run from this folder: python make_ch07.py [preview_dir]

ch07_snow_resistance: snow logbooks (firnice_snowclim_03116.dat) of runs snowlog_w0 (before the snow
fix, GloGEM 3b33e4b) and snowfix_w0 (after it, c400a50), both at alps_r7 design 13.
ch07_hef_tongue: model minus boreholes at the four Hintereisferner tongue sites, same two runs, each
reading compared at its date (icetemp 2f603aa); numbers from Step 2 test (b) (GloGEM c400a50).
ch07_runs_rmse: pooled RMSE of each of the 32 design points in alps_r7 (re-matched with icetemp
2f603aa) and alps_s2 (snow fix), on the 769 common points; numbers from the Step 4 table of
scripts/icetemp_analysis/score_paired_ensembles.py (icetemp c9af84f), section "Every paired design point".
ch07_timing, ch07_snow_gap: schematics.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

OUT = Path(__file__).parent
PREVIEW = Path(sys.argv[1]) if len(sys.argv) > 1 else None
RUNS = Path(__file__).resolve().parents[2] / 'data' / 'bayescal'
LOGBOOK = 'training_runs/monthly/CentralEurope/PAST/firnice_temperature/firnice_snowclim_03116.dat'
plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none'})

SURFACE, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
AFTER, BEFORE = '#2a78d6', '#eb6834'   # categorical slots 1-2, fixed order; corrected state first, as in ch06
NEUTRAL = '#bdbbb5'                     # quantities shared by both runs
SNOW, ICE = '#f4f4f1', '#dde6ee'        # schematic fills

MONTHS = ['Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep']
ORDER = [10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9]   # glaciological year
SNO, R_USED, R_SERIES = 3, 10, 11                  # logbook columns (elev month firn_frac sno_mwe ...)

# model minus boreholes [C] at alps_r7 design 13: before (snowlog_w0) and after (snowfix_w0) the fix
HEF = [('2955 m', -4.77, -2.57), ('2805 m', -5.09, -2.74), ('2615 m', -5.09, -3.04),
       ('2555 m', -4.68, -2.87)]
HEF_ALL = ('all 41 points', -4.89, -2.82)

# pooled RMSE [C] per design point: (design point, alps_r7, alps_s2), 769 common points
RUNS_RMSE = [(0, 2.386, 2.220), (1, 2.491, 2.306), (2, 2.521, 2.349), (3, 2.347, 2.193),
             (4, 2.901, 2.751), (5, 2.561, 2.410), (6, 2.627, 2.460), (7, 2.596, 2.408),
             (8, 2.428, 2.286), (9, 2.402, 2.208), (10, 3.214, 3.065), (11, 4.357, 4.146),
             (12, 2.767, 2.616), (13, 2.327, 2.140), (14, 2.356, 2.165), (15, 2.422, 2.223),
             (16, 3.902, 3.773), (17, 2.768, 2.606), (18, 2.923, 2.757), (19, 2.731, 2.535),
             (20, 3.715, 3.583), (21, 2.944, 2.773), (22, 2.882, 2.721), (23, 2.583, 2.402),
             (24, 2.884, 2.752), (25, 2.525, 2.353), (26, 2.452, 2.318), (27, 2.684, 2.524),
             (28, 2.635, 2.482), (29, 3.067, 2.914), (30, 3.481, 3.322), (31, 2.615, 2.432)]


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


def key(ax, x, y, text, color, dx=0.35):
    """Direct label: a short line in the series colour, text in ink."""
    ax.plot([x, x + dx], [y, y], color=color, lw=2, solid_capstyle='round')
    ax.text(x + dx + 0.12, y, text, color=INK, va='center', fontsize=9)


def logbook(run, elev=2615):
    """One band's logbook rows, October to September."""
    d = np.loadtxt(RUNS / run / LOGBOOK, comments='#')
    d = d[d[:, 0] == elev]
    return np.array([d[d[:, 1] == m][0] for m in ORDER])


def fig_snow_resistance():
    before, after = logbook('snowlog_w0'), logbook('snowfix_w0')
    x = np.arange(12)
    fig, (ax, axs) = plt.subplots(2, 1, figsize=(9, 6.2), sharex=True, facecolor=SURFACE,
                                  gridspec_kw=dict(height_ratios=[3, 1.2], hspace=0.15))
    for a in (ax, axs):
        style(a)
        a.grid(axis='y', color=GRID, lw=0.8)
        a.set_axisbelow(True)

    ax.scatter(x, before[:, R_SERIES], s=150, facecolor='none', edgecolor=INK, lw=1.2, zorder=5)
    ax.plot(x, after[:, R_USED], color=AFTER, lw=2, marker='o', ms=6, mec=SURFACE, mew=1.5, zorder=4)
    ax.plot(x, before[:, R_USED], color=BEFORE, lw=2, marker='o', ms=6, mec=SURFACE, mew=1.5, zorder=3)
    jan = MONTHS.index('Jan')
    ax.annotate('', xy=(jan, before[jan, R_USED] + 0.3), xytext=(jan, before[jan, R_SERIES] - 0.45),
                arrowprops=dict(arrowstyle='<->', color=INK2, lw=1))
    ax.text(jan + 0.15, 4.4, 'Dec–Feb: 3–5× too little', color=INK, ha='left', va='center',
            fontsize=9)
    may = MONTHS.index('May')
    key(ax, may - 1.05, before[may, R_SERIES] + 0.95, 'after the fix (= layered)', AFTER)
    key(ax, may - 0.7, before[may, R_USED] + 0.85, 'before the fix', BEFORE)
    ax.set_ylim(0, 13.8)
    ax.set_ylabel('resistance [m² K W⁻¹]', color=INK)
    handles = [
        plt.Line2D([], [], marker='o', ls='', ms=11, mfc='none', mec=INK,
                   label='layered snow and ice (what it should be)'),
        plt.Line2D([], [], marker='o', color=AFTER, lw=2, ms=6, mec=SURFACE,
                   label='as the solver uses it, after the fix'),
        plt.Line2D([], [], marker='o', color=BEFORE, lw=2, ms=6, mec=SURFACE,
                   label='as the solver used it, before the fix'),
    ]
    ax.legend(handles=handles, loc='upper left', fontsize=8, frameon=False, labelcolor=INK)
    ax.set_title('Hintereisferner tongue, 2615 m: resistance from the surface to the first node '
                 'under the snow', fontsize=10, color=INK, loc='left')

    axs.bar(x, before[:, SNO], width=0.45, color=NEUTRAL, zorder=2)
    for m in ('Dec', 'Jan', 'Feb'):
        i = MONTHS.index(m)
        axs.text(i, before[i, SNO] + 0.04, f'{before[i, SNO]:.3f}', ha='center', va='bottom',
                 color=INK, fontsize=8)
    axs.set_ylim(0, 0.95)
    axs.set_yticks([0, 0.4, 0.8])
    axs.set_ylabel('snow\n[m w.e.]', color=INK)
    axs.text(-0.4, 0.86, 'snow at the start of the month (the same in both runs)', color=INK2,
             fontsize=8.5, va='top')
    axs.set_xticks(x)
    axs.set_xticklabels(MONTHS, color=INK2)
    axs.set_xlim(-0.6, 11.6)
    fig.text(0.01, 0.01, 'Monthly means, October 1971 – September 2025, alps_r7 design 13; '
             'snow logbooks of runs snowlog_w0 and snowfix_w0.', color=INK2, fontsize=8)
    fig.subplots_adjust(left=0.1, right=0.98, top=0.94, bottom=0.1)
    save(fig, 'ch07_snow_resistance')


def fig_hef_tongue():
    fig, ax = plt.subplots(figsize=(8.5, 3.6), facecolor=SURFACE)
    style(ax)
    ax.grid(axis='x', color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    rows = [(4.4 - i, *r) for i, r in enumerate(HEF)] + [(0, *HEF_ALL)]
    for y, label, b, a in rows:
        ax.plot([b, a], [y, y], color=GRID, lw=2, zorder=1)
        ax.scatter(b, y, s=70, color=BEFORE, edgecolor=SURFACE, lw=1.5, zorder=3)
        ax.scatter(a, y, s=70, color=AFTER, edgecolor=SURFACE, lw=1.5, zorder=4)
    _, _, b, a = rows[0]
    ax.text(b, 4.72, 'before', color=INK2, ha='center', va='bottom', fontsize=8.5)
    ax.text(a, 4.72, 'after', color=INK2, ha='center', va='bottom', fontsize=8.5)
    _, _, b, a = rows[-1]
    ax.text(b - 0.14, 0, f'{b:.2f}'.replace('-', '−'), color=INK, ha='right', va='center', fontsize=9)
    ax.text(a + 0.14, 0, f'{a:.2f}'.replace('-', '−'), color=INK, ha='left', va='center', fontsize=9)
    ax.axhline(0.7, color=GRID, lw=0.8)
    ax.axvline(0, color=INK2, lw=1)
    ax.text(-0.08, -0.55, 'model = boreholes', color=INK2, ha='right', va='center', fontsize=8.5)
    ax.set_yticks([r[0] for r in rows])
    ax.set_yticklabels([r[1] for r in rows], color=INK)
    ax.set_ylim(-0.8, 5.15)
    ax.set_xlim(-5.6, 0.4)
    ax.set_xlabel('model − boreholes [°C], each reading compared at its date', color=INK)
    handles = [plt.Line2D([], [], marker='o', ls='', ms=8, color=BEFORE, mec=SURFACE,
                          label='before the snow fix (snowlog_w0)'),
               plt.Line2D([], [], marker='o', ls='', ms=8, color=AFTER, mec=SURFACE,
                          label='after the snow fix (snowfix_w0)')]
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(0.56, 0.55), fontsize=8,
              frameon=False, labelcolor=INK)
    ax.set_title('Hintereisferner tongue: misfit at the four borehole sites (alps_r7 design 13)',
                 fontsize=10, color=INK, loc='left')
    fig.tight_layout()
    save(fig, 'ch07_hef_tongue')


def fig_runs_rmse():
    """Pooled RMSE of every design point before and after the snow fix, best run before it on top."""
    rows = sorted(RUNS_RMSE, key=lambda r: r[1])
    fig, ax = plt.subplots(figsize=(8.5, 7.6), facecolor=SURFACE)
    style(ax)
    ax.grid(axis='x', color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for y, (_, b, a) in enumerate(rows):
        ax.plot([a, b], [y, y], color=GRID, lw=2, zorder=1)
        ax.scatter(b, y, s=45, color=BEFORE, edgecolor=SURFACE, lw=1.2, zorder=3)
        ax.scatter(a, y, s=45, color=AFTER, edgecolor=SURFACE, lw=1.2, zorder=4)
    _, b, a = rows[0]
    ax.text(b, -0.75, 'before', color=INK2, ha='center', va='bottom', fontsize=8.5)
    ax.text(a, -0.75, 'after', color=INK2, ha='center', va='bottom', fontsize=8.5)
    ax.text(b + 0.06, 0, f'{b:.3f} → {a:.3f} °C, the best run in both', color=INK, ha='left',
            va='center', fontsize=8.5)
    ax.text(3.32, 9.2, 'lower in all 32 runs:\nby 0.164 °C in the median\n(0.129 to 0.210 °C)',
            color=INK, ha='left', va='top', fontsize=9)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([str(r[0]) for r in rows], color=INK2, fontsize=8)
    ax.set_ylim(len(rows) - 0.4, -1.5)
    ax.set_xlim(2.0, 4.5)
    ax.set_ylabel('design point, ordered by its RMSE before the fix', color=INK)
    ax.set_xlabel('RMSE over the 769 common points [°C]', color=INK)
    handles = [plt.Line2D([], [], marker='o', ls='', ms=8, color=BEFORE, mec=SURFACE,
                          label='before the snow fix (alps_r7)'),
               plt.Line2D([], [], marker='o', ls='', ms=8, color=AFTER, mec=SURFACE,
                          label='after the snow fix (alps_s2)')]
    ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.52, 0.86), fontsize=8.5,
              frameon=False, labelcolor=INK)
    ax.set_title('The snow fix in each of the 32 runs of the Alpine ensemble (same parameter sets)',
                 fontsize=10, color=INK, loc='left')
    fig.text(0.01, 0.01, 'Each reading compared at its date (icetemp 2f603aa); paired scores of '
             'alps_r7 and alps_s2 (icetemp c9af84f).', color=INK2, fontsize=8)
    fig.subplots_adjust(left=0.1, right=0.98, top=0.95, bottom=0.09)
    save(fig, 'ch07_runs_rmse')


def fig_timing():
    """Schematic: a July reading against the month-end profile and against the interpolated one."""
    fig, ax = plt.subplots(figsize=(9, 4.0), facecolor=SURFACE)
    style(ax)
    t = np.linspace(0, 92, 300)                  # days from 1 June
    curve = lambda d: -5 + 4 * np.sin(np.pi * d / 184)
    ax.plot(t, curve(t), color=INK2, lw=1.2, zorder=2)
    ax.text(12, curve(12) - 0.3, 'model temperature,\nnear the surface', color=INK2, ha='left',
            va='top', fontsize=8.5)
    ends = np.array([30, 61, 92])                # 30 June, 31 July, 31 August
    for e in ends:
        ax.axvline(e, color=GRID, lw=0.8, zorder=1)
    ax.scatter(ends, curve(ends), s=55, color=INK, edgecolor=SURFACE, lw=1.5, zorder=4)
    ax.plot(ends[:2], curve(ends[:2]), color=INK2, lw=1, ls=(0, (4, 3)), zorder=2)

    d = 40                                       # reading on 10 July
    f = (d - ends[0]) / (ends[1] - ends[0])
    t_new = curve(ends[0]) + f * (curve(ends[1]) - curve(ends[0]))
    ax.axvline(d, color=INK2, lw=1, ls=':', zorder=1)
    ax.scatter([d], [-5.55], marker='^', s=70, color=INK, zorder=4, clip_on=False)
    ax.text(d + 1, -5.55, 'reading on 10 July', color=INK, va='center', fontsize=9)
    ax.scatter([ends[1]], [curve(ends[1])], s=170, facecolor='none', edgecolor=BEFORE, lw=2.2,
               zorder=5)
    ax.annotate('before: every July reading was\ncompared with the profile of 31 July',
                xy=(ends[1] + 0.6, curve(ends[1]) - 0.2), xytext=(66, -3.7), color=INK,
                fontsize=9, va='center', arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    ax.scatter([d], [t_new], s=70, color=AFTER, edgecolor=SURFACE, lw=1.5, zorder=5)
    ax.annotate('now: the model on 10 July, on the straight\nline between the two month ends',
                xy=(d - 0.8, t_new + 0.1), xytext=(3, -1.25), color=INK, fontsize=9,
                va='center', arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    ax.plot([d + 1, ends[1] + 4.2], [t_new, t_new], color=INK2, lw=0.8, ls=':', zorder=1)
    ax.annotate('', xy=(ends[1] + 3.5, curve(ends[1])), xytext=(ends[1] + 3.5, t_new),
                arrowprops=dict(arrowstyle='<->', color=INK2, lw=1, shrinkA=0, shrinkB=0))
    ax.text(ends[1] + 4.6, (t_new + curve(ends[1])) / 2, 'error of the old\ncomparison',
            color=INK, va='center', fontsize=8.5)
    ax.text(ends[0] - 1, curve(ends[0]) - 0.3, 'profile written\n30 June', color=INK2,
            ha='right', va='top', fontsize=8.5)
    ax.text(ends[2] - 1, curve(ends[2]) - 0.3, 'profile written\n31 August', color=INK2,
            ha='right', va='top', fontsize=8.5)
    for mid, name in ((15, 'June'), (45.5, 'July'), (76.5, 'August')):
        ax.text(mid, 0.05, name, color=INK2, ha='center', va='bottom', fontsize=9)
    ax.set_xlim(0, 95)
    ax.set_ylim(-5.9, 0.45)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_ylabel('temperature near the surface', color=INK)
    ax.set_title('Comparing a reading with the model at its date (schematic)', fontsize=10,
                 color=INK, loc='left')
    fig.tight_layout()
    save(fig, 'ch07_timing')


def column(ax, x0, k_gap, title, note, heat_lw):
    """One ice-band column: snow filling the 0-1 m gap, ice below, nodes at 0, 1 and 2 m."""
    w = 1.0
    ax.add_patch(Rectangle((x0, 0), w, 1, fc=SNOW, ec=INK2, lw=0.8, hatch='.', zorder=1))
    ax.add_patch(Rectangle((x0, 1), w, 1.6, fc=ICE, ec=INK2, lw=0.8, zorder=1))
    ax.text(x0 + w * 0.64, 0.5, 'snow\n250 kg/m³', ha='center', va='center', fontsize=8.5,
            color=INK, bbox=dict(fc=SNOW, ec='none', pad=4))
    ax.text(x0 + w * 0.64, 1.45, 'ice', ha='center', va='center', fontsize=8.5, color=INK)
    ax.annotate('', xy=(x0 + 0.13, 0.12), xytext=(x0 + 0.13, 2.45),
                arrowprops=dict(arrowstyle='-|>', color=INK2, lw=heat_lw, mutation_scale=12 + 3 * heat_lw),
                zorder=3)
    for z, lab in ((0, 'node 0: surface,\nset by the air'), (1, 'node 1'), (2, 'node 2')):
        ax.scatter([x0 + w / 2], [z], s=40, color=INK, edgecolor=SURFACE, lw=1.5, zorder=4)
        ax.text(x0 - 0.1, z, lab, ha='right', va='center', fontsize=8.5, color=INK2)
    ax.text(x0 + w + 0.15, 0.5, k_gap, ha='left', va='center', fontsize=9, color=INK)
    ax.text(x0 + w + 0.15, 1.5, 'ice: 2.13 W/m/K', ha='left', va='center', fontsize=9, color=INK)
    ax.text(x0 - 1.0, -0.62, title, ha='left', va='bottom', fontsize=10, color=INK,
            fontweight='bold')
    ax.text(x0 - 1.0, 2.95, note, ha='left', va='top', fontsize=8.5, color=INK2)


def fig_snow_gap():
    fig, ax = plt.subplots(figsize=(10, 4.6), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    column(ax, 1.2, 'gap 0–1 m: (0.15 + 2.13) / 2\n= 1.14 W/m/K',
           'what the solver did',
           'the snow sat on node 0 and was averaged\nwith the ice of node 1: heat passed\n'
           'almost as through ice', heat_lw=4.5)
    column(ax, 6.0, 'gap 0–1 m: the snow alone\n= 0.15 W/m/K',
           'after the fix: snow and ice in series',
           'each gap adds the resistances\nof the snow and the ice it holds', heat_lw=1.2)
    ax.text(4.7, 2.3, 'arrows: heat flowing up to\nthe cold surface in winter', ha='center',
            va='center', fontsize=8.5, color=INK2)
    ax.set_xlim(0, 9.4)
    ax.set_ylim(3.05, -0.75)
    ax.axis('off')
    fig.tight_layout()
    save(fig, 'ch07_snow_gap')


if __name__ == '__main__':
    fig_snow_resistance()
    fig_hef_tongue()
    fig_runs_rmse()
    fig_timing()
    fig_snow_gap()
    print('written:', sorted(p.name for p in OUT.glob('ch07_*.svg')))
