"""Figures for chapter 06 (firn extent and the mass-balance calibration).
Run from this folder: python make_ch06.py [preview_dir]

Numbers from the firn line logs (firnice_firnline_<id>.dat) of runs firnlog_w0 (RGI6-numbered
calibration) and firnrgi7_w0 (RGI7 calibration), same thermal parameters, and the WGMS
Fluctuations of Glaciers overview (2017-06 release).
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt

OUT = Path(__file__).parent
PREVIEW = Path(sys.argv[1]) if len(sys.argv) > 1 else None
plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none'})

SURFACE, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
RGI7_COL, RGI6_COL = '#2a78d6', '#eb6834'   # categorical slots 1-2, fixed order

# ELA 2000-2019 [m]: each glacier's own RGI7 calibration, and the modelled mean of the two runs
GLACIERS = [  # (label, calibration, run with RGI6-numbered file, run with RGI7 file)
    ('Monte Rosa east face (01230)', 2827, 2385, 2826),
    ('Corvatsch (02103)', 3164, 3166, 3163),
    ('Tête Rousse (00777)', 3180, 3164, 3177),
    ('Hintereisferner (03116)', 3193, 3425, 3194),
    ('Hohsaas (01112)', 3213, 2884, 3212),
    ('Bossons (00773)', 3257, 3019, 3255),
    ('Gorner/Grenz (01225)', 3268, 2762, 3266),
]

YEARS = [1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986,
         1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001,
         2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2016]
WGMS = [2935, 3250, 2910, 2905, 2995, 2840, 2825, 2970, 2930, 2940, 3260, 3075, 2970, 3010, 3080,
        3070, 3130, 3080, 3115, 3260, 3155, 3050, 3145, 3080, 3100, 3050, 3160, 3105, 3050, 2955,
        3050, 3750, 3185, 3225, 3750, 3750, 3276, 3257, 3116, 3285, 3380, 3047, 2990, 3100]
HEF_RGI6 = [3115, 3315, 2995, 3115, 3065, 2875, 2975, 3045, 3075, 3085, 3295, 3235, 3005, 3235,
            3205, 3275, 3225, 3155, 3285, 3315, 3355, 3195, 3405, 3205, 3305, 3195, 3345, 3225,
            3245, 3065, 3455, 3685, 3305, 3385, 3525, 3445, 3465, 3425, 3425, 3325, 3375, 3285,
            3195, 3415]
HEF_RGI7 = [2875, 3135, 2765, 2855, 2885, 2635, 2745, 2845, 2835, 2905, 3105, 2985, 2765, 3025,
            2995, 3015, 2985, 2915, 3045, 3075, 3065, 3015, 3145, 2945, 3165, 3015, 3135, 3015,
            3035, 2835, 3235, 3455, 3095, 3205, 3275, 3245, 3245, 3205, 3205, 3145, 3185, 3045,
            2985, 3205]


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


def fig_calibration_mismatch():
    fig, ax = plt.subplots(figsize=(9, 4.2), facecolor=SURFACE)
    style(ax)
    ax.grid(axis='x', color=GRID, lw=0.8)
    for y, (label, cal, r6, r7) in enumerate(GLACIERS):
        ax.plot([r6, r7], [y, y], color=GRID, lw=2, zorder=1)
        ax.scatter(r6, y, s=64, color=RGI6_COL, edgecolor=SURFACE, lw=1.5, zorder=3)
        ax.scatter(r7, y, s=64, color=RGI7_COL, edgecolor=SURFACE, lw=1.5, zorder=4)
        ax.scatter(cal, y, s=150, facecolor='none', edgecolor=INK, lw=1.2, zorder=5)
    ax.set_yticks(range(len(GLACIERS)))
    ax.set_yticklabels([g[0] for g in GLACIERS], color=INK)
    ax.set_xlabel('mean ELA 2000–2019 [m a.s.l.]', color=INK)
    handles = [
        plt.Line2D([], [], marker='o', ls='', ms=11, mfc='none', mec=INK, label='own calibration (target)'),
        plt.Line2D([], [], marker='o', ls='', ms=8, color=RGI6_COL, mec=SURFACE,
                   label='run with the RGI6-numbered file'),
        plt.Line2D([], [], marker='o', ls='', ms=8, color=RGI7_COL, mec=SURFACE,
                   label='run with the RGI7 file'),
    ]
    ax.legend(handles=handles, loc='upper left', fontsize=8, frameon=False, labelcolor=INK)
    ax.set_title('Modelled ELA against each glacier’s own calibration', fontsize=10, color=INK,
                 loc='left')
    fig.tight_layout()
    save(fig, 'ch06_calibration_mismatch')


def fig_hef_ela():
    fig, ax = plt.subplots(figsize=(10, 4.2), facecolor=SURFACE)
    style(ax)
    ax.grid(axis='y', color=GRID, lw=0.8)
    ax.plot(YEARS, WGMS, color=INK, lw=2, marker='o', ms=4, label='WGMS (observed)')
    ax.plot(YEARS, HEF_RGI6, color=RGI6_COL, lw=2, label='model, RGI6-numbered file')
    ax.plot(YEARS, HEF_RGI7, color=RGI7_COL, lw=2, label='model, RGI7 file')
    for vals, name, col in ((HEF_RGI6, 'RGI6 file: mean +151 m', RGI6_COL),
                            (HEF_RGI7, 'RGI7 file: mean −68 m', RGI7_COL)):
        ax.text(YEARS[-1] + 0.8, vals[-1], name, color=INK, va='center', fontsize=9)
        ax.plot([YEARS[-1] + 0.2, YEARS[-1] + 0.7], [vals[-1]] * 2, color=col, lw=2)
    ax.set_xlim(1970.5, 2023)
    ax.set_ylabel('ELA [m a.s.l.]', color=INK)
    ax.legend(loc='upper left', fontsize=8, frameon=False, labelcolor=INK)
    ax.set_title('Hintereisferner: modelled and observed ELA (3750 m = no accumulation area)',
                 fontsize=10, color=INK, loc='left')
    fig.tight_layout()
    save(fig, 'ch06_hef_ela')


if __name__ == '__main__':
    fig_calibration_mismatch()
    fig_hef_ela()
    print('written:', sorted(p.name for p in OUT.glob('ch06_*.svg')))
