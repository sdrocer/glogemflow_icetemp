"""Figures for chapter 09 (snow and firn from data): the snowline calibration of 2026-10-16.
Run from this folder: python make_ch09.py [preview_dir]

ch09_snowline_definitions: schematic on Morteratsch's band hypsometry (RGI7 02216): her model snow
line (highest snow-free band, else the lowest band) against Sommer et al.'s statistic (10th
percentile of the elevation of the snow-covered area), on a spring day and a fully snow-covered day.
ch09_rmse_curves: RMSE of the daily snow line against the observed snow lines, by scanned c_prec, for
the 12 glaciers; scans 0.4-3.0 with all months (cal12_p10raw_wide) and May-September (cal12_seas_wide).
ch09_error_by_month: (a) the standard calibration of 2026-10-14 (monthly model, its own snow-line
definition, run mbout/regional), weighted bias by month with all days and without the days Sommer
marks as fully snow-covered; (b) after the snowline calibration (all months, p10; daily forward run
fwd_p10raw), weighted bias and share of the squared error by month.
ch09_toff_by_elevation: T_off of the 23 calibration glaciers before (standard monthly calibration,
refit_none) and after the snowline calibration (refit_p10raw: all months; refit_seas: May-September),
against the area-weighted mean elevation (toff_test.py tables).
ch09_cprec_bootstrap: c_prec that fits the May-September snow lines best (full-sample spline minimum over
0.4-3.0) and its 5-95 % range from 1000 bootstrap draws over the 20 hydrological years (cprec_bootstrap.py,
boot_result/cprec_bootstrap.csv; one snowline run per single c_prec, fixed target period); kept if the range
reaches neither end of the scan.
ch09_daily_vs_monthly: the daily and the monthly model on the same satellite snow lines (12 glaciers, 2000-2019,
May-September; scripts/icetemp_analysis/snowline/daily_vs_monthly.py): (a) weighted bias by month, monthly model at its
standard calibration (fwdB_none_fixed, interpolated between month ends) against the range of the daily runs over
c_prec 0.4-3.0 (boot/f*); (b) the summer's highest snow line (July-September, mean over the years) per glacier.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUT = Path(__file__).parent
PREVIEW = Path(sys.argv[1]) if len(sys.argv) > 1 else None
REPO = Path(__file__).resolve().parents[2]
SNOWLINES = REPO / 'data' / 'snowlines'
SNOWCAL = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week2/snowcal')
MBOUT = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week2/mbout/regional')
BANDS = Path('/itet-stor/jabeer/glogem/geometricdata/rgiv7/bands/centraleurope')
CAL = 'daily/CentralEurope/calibration/calibrate_m1_cID9_centraleurope_CentralEurope_snowcal_snowline.dat'
plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none'})

SURFACE, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
AFTER, BEFORE = '#2a78d6', '#eb6834'   # categorical slots 1-2: May-Sep / corrected first, as in ch06-07
NEUTRAL = '#bdbbb5'
SNOW, ICE = '#f4f4f1', '#dde6ee'

NAMES = {'00773': 'Bossons', '00774': 'Taconnaz', '01112': 'Hohsaas', '01127': 'Gruben', '01205': 'Lys',
         '01225': 'Gorner/Grenz', '01230': 'Monte Rosa E face', '02216': 'Morteratsch', '02596': 'Gr. Aletsch',
         '02611': 'Fiescher', '02821': 'Ortles (Alta)', '03116': 'Hintereis'}


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


def p10(snow, area, elev, step=10.):
    """GloGEM's snowline_p10: 10th percentile of the elevation of the snow-covered area."""
    js = np.where(snow & (area > 0))[0]
    if len(js) == 0:
        return elev[area > 0].max() + step / 2
    ca = np.cumsum(area[js]) / area[js].sum()
    k = int(np.argmax(ca >= 0.1))
    ca0 = ca[k - 1] if k > 0 else 0.
    return elev[js[k]] - step / 2 + step * (0.1 - ca0) / (ca[k] - ca0)


def fig_snowline_definitions():
    b = np.loadtxt(BANDS / '02216.dat', skiprows=5, usecols=(1, 3))
    elev, area = b[:, 0] + 5, b[:, 1]
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 5.2), sharey=True, facecolor=SURFACE)
    cases = [('a  spring day: snow above 2700 m', elev >= 2700),
             ('b  fully snow-covered day', np.ones_like(elev, bool))]
    for ax, (title, snow) in zip(axes, cases):
        style(ax)
        ax.barh(elev, area, height=10, color=list(np.where(snow, '#e6e5e0', '#bcd3ea')), edgecolor=list(np.where(snow, '#b3b1aa', '#7f9fbf')),
                lw=0.3, zorder=2)
        free = np.where(~snow & (area > 0))[0]
        hers = elev[free[-1]] if len(free) else elev[area > 0][0]
        som = p10(snow, area, elev)
        x1 = area.max() * 1.05
        for z, col, lab, va in ((hers, BEFORE, 'her model snow line', 'top'),
                                (som, AFTER, "Sommer's statistic (p10)", 'bottom')):
            ax.plot([0, x1], [z, z], color=col, lw=2, zorder=4)
            ax.text(x1, z, f' {lab}\n {z:.0f} m', color=col, va=va, fontsize=9)
        ax.set_title(title, loc='left', color=INK, fontsize=10)
        ax.set_xlabel('band area [km²]', color=INK)
        ax.set_xlim(0, area.max() * 1.9)
        ax.set_ylim(1880, 4000)
        ax.grid(axis='x', color=GRID, lw=0.6, zorder=0)
        if not len(free):
            ax.annotate('', xy=(x1 * 0.55, som), xytext=(x1 * 0.55, hers),
                        arrowprops=dict(arrowstyle='<->', color=INK2, lw=1))
            ax.text(x1 * 0.58, (som + hers) / 2, f'{som - hers:.0f} m: "error"\nthat no c_prec\ncan remove',
                    color=INK, va='center', fontsize=9)
    axes[0].set_ylabel('elevation [m a.s.l.]', color=INK)
    fig.text(0.01, 0.005, 'Morteratsch (RGI7 02216), 10 m bands; snow-covered bands grey, snow-free bands blue.\n'
             'Her rule: the highest snow-free band, else the lowest band. p10: 10 % of the snow-covered area lies below.',
             color=INK2, fontsize=8)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    save(fig, 'ch09_snowline_definitions')


def read_log(run):
    t = pd.read_csv(SNOWCAL / run / CAL, sep=r'\s+', dtype={'ID': str}).set_index('ID')
    cols = list(t.columns[6:])
    x = np.array([float(c) for c in cols])
    r = t[cols].astype(float).where(lambda v: v > -98)
    return t, x, r


def fig_rmse_curves():
    t_all, x, r_all = read_log('cal12_p10raw_wide')
    t_ms, _, r_ms = read_log('cal12_seas_wide')
    gain = 100 * (r_ms.iloc[:, 0] - r_ms.min(axis=1)) / r_ms.iloc[:, 0]
    order = gain.sort_values(ascending=False).index
    fig, axes = plt.subplots(3, 4, figsize=(11, 7.6), sharex=True, sharey=True, facecolor=SURFACE)
    for ax, gid in zip(axes.flat, order):
        style(ax)
        ax.grid(color=GRID, lw=0.6)
        for r, t, col, lab in ((r_all, t_all, BEFORE, 'all months'), (r_ms, t_ms, AFTER, 'May–September')):
            y = r.loc[gid] - r.loc[gid].min()
            ax.plot(x, y, color=col, lw=1.8, marker='o', ms=3.5, mec=SURFACE, mew=0.8, label=lab)
            ax.plot([t.loc[gid, 'CprecMin']], [0], marker='v', color=col, ms=7, mec=SURFACE, zorder=5)
        ax.set_title(f'{NAMES[gid]} ({gid})\nMay–Sep gain {gain[gid]:.0f} %', fontsize=9, color=INK, loc='left')
        ax.set_ylim(-5, 125)
    axes[0, 0].legend(fontsize=8, frameon=False, labelcolor=INK, loc='upper right')
    for ax in axes[-1]:
        ax.set_xlabel('c_prec', color=INK)
    for ax in axes[:, 0]:
        ax.set_ylabel('RMSE − its minimum [m]', color=INK)
    fig.text(0.01, 0.005, 'Each scan stops at the first rise of the RMSE; ▼: c_prec chosen (spline minimum). '
             '"Gain": how much better the best c_prec fits than 0.4 (May–Sep). Sorted by gain.',
             color=INK2, fontsize=8)
    fig.tight_layout(rect=(0, 0.025, 1, 1))
    save(fig, 'ch09_rmse_curves')


def wstats(res, sig):
    w = 1 / sig
    return np.sum(w * res) / w.sum(), np.sum(w * res ** 2)


def month_before():
    """2026-10-14 comparison: monthly model (standard calibration, own snow-line definition), May-Oct."""
    ny, y0 = 36, 1990
    mod = {}
    with open(MBOUT / 'monthly/CentralEurope/PAST/PAST_original/centraleurope_Snowline_mon_r1_CentralEurope_fair.dat') as f:
        next(f)
        for ln in f:
            p = ln.split()
            if p:
                mod[p[0]] = np.array(p[1:], float).reshape(ny, 12)
    rows = []
    for fn in sorted((SNOWLINES / 'tsla_csv').glob('tsla_RGI70-11.*.csv')):
        gid = fn.stem[-5:]
        d = pd.read_csv(fn, parse_dates=['date'])
        d = d[(d.date <= '2025-09-30') & d.month.between(5, 10) & (d.qa >= 0.666)]
        Y = d.date.dt.year + (d.date.dt.month >= 10)
        i = Y - y0
        ok = (i >= 0) & (i < ny)
        d, i = d[ok], i[ok]
        d['mod'] = mod[gid][i.values, ((d.date.dt.month - 10) % 12).values]
        d['full'] = d.flag.fillna('').str.contains('full_snow')
        rows.append(d)
    D = pd.concat(rows)
    D['res'] = D['mod'] - D.sla_m
    out = {}
    for m in range(5, 11):
        a = D[D.month == m]
        out[m] = (wstats(a.res.values, a.sigma_m.values)[0], wstats(a[~a.full].res.values, a[~a.full].sigma_m.values)[0])
    return out


def month_after():
    """Daily forward run with the all-months p10 calibration, against all raw observations 2000-2019."""
    mod = {}
    for ln in open(SNOWCAL / 'fwd_p10raw/daily/CentralEurope/PAST/PAST_original/centraleurope_Snowline_day_r1_CentralEurope_snowcal.dat'):
        p = ln.split()
        if p[0] == 'ID' or len(p) < 300:
            continue
        mod.setdefault(p[0], {})[int(p[1])] = np.array(p[3:], float)
    o = pd.read_csv(SNOWLINES / 'snowline_2000_2020/11_snowline_median_raw.dat', sep=r'\s+', skiprows=1, header=None,
                    names=['rid', 'sla', 'sig', 'year', 'doy'])
    o['gid'] = o.rid.str[-5:]
    k = 365 * (o.year - 2000 + 1) + o.doy - 274
    o = o[(k >= 0) & (k < 365 * 20)].copy()
    k = k[o.index]
    o['res'] = [mod[g][2000 + kk // 365][kk % 365] for g, kk in zip(o.gid, k)] - o.sla
    o['month'] = (pd.to_datetime('2001-01-01') + pd.to_timedelta(o.doy - 1, unit='D')).dt.month
    tot = np.sum(o.res ** 2 / o.sig)
    out = {}
    for m in (4, 5, 6, 7, 8, 9, 10, 11):
        a = o[o.month == m]
        bias, sse = wstats(a.res.values, a.sig.values)
        out[m] = (bias, sse / tot)
    return out


def fig_error_by_month():
    bef, aft = month_before(), month_after()
    names = {4: 'Apr', 5: 'May', 6: 'Jun', 7: 'Jul', 8: 'Aug', 9: 'Sep', 10: 'Oct', 11: 'Nov'}
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.4), facecolor=SURFACE, gridspec_kw={'width_ratios': [6, 8]})
    for a in (a1, a2):
        style(a)
        a.axhline(0, color=INK2, lw=0.8)
        a.grid(axis='y', color=GRID, lw=0.6)
    m1 = sorted(bef)
    xx = np.arange(len(m1))
    a1.bar(xx - 0.2, [bef[m][0] for m in m1], 0.4, color=NEUTRAL, label='all days', zorder=2)
    a1.bar(xx + 0.2, [bef[m][1] for m in m1], 0.4, color=BEFORE, label='without fully snow-covered days', zorder=2)
    a1.set_xticks(xx, [names[m] for m in m1])
    a1.set_ylabel('model − observed snow line [m]', color=INK)
    a1.set_title('a  standard calibration (14 Oct), monthly model:\n    the spring snow persists, with or without '
                 'fully covered days', loc='left', fontsize=9.5, color=INK)
    a1.legend(fontsize=8, frameon=False, labelcolor=INK, loc='lower right')
    m2 = sorted(aft)
    xx = np.arange(len(m2))
    a2.bar(xx, [aft[m][0] for m in m2], 0.55, color=AFTER, zorder=2, label='bias (left axis)')
    a2.set_xticks(xx, [names[m] for m in m2])
    a2.set_ylabel('model − observed snow line [m]', color=INK)
    a2.set_ylim(-300, 300)
    tw = a2.twinx()
    tw.plot(xx, [100 * aft[m][1] for m in m2], color=INK, marker='o', ms=5, lw=1.5, label='share of the squared error')
    tw.set_ylim(0, 40)
    tw.set_ylabel('share of the squared error [%]', color=INK)
    tw.tick_params(colors=INK2)
    for s in ('top',):
        tw.spines[s].set_visible(False)
    tw.spines['right'].set_color(GRID)
    a2.set_title('b  after the snowline calibration (all months, p10), daily model:\n    spring right on average; '
                 'October–November hold 45 % of the error', loc='left', fontsize=9.5, color=INK)
    h1, l1 = a2.get_legend_handles_labels()
    h2, l2 = tw.get_legend_handles_labels()
    a2.legend(h1 + h2, l1 + l2, fontsize=8, frameon=False, labelcolor=INK, loc='upper left')
    fig.text(0.01, 0.005, 'Weighted by 1/σ; 12 glaciers with Sommer snow lines (QA ≥ 4/6). (a) May–October, the '
             "model's own snow-line definition; (b) all observations, raw snow lines, 2000–2019.",
             color=INK2, fontsize=8)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, 'ch09_error_by_month')


def fig_toff_by_elevation():
    a = pd.read_csv(SNOWCAL / 'toff_test_p10raw.csv', dtype={'ID': str}).set_index('ID')
    s = pd.read_csv(SNOWCAL / 'toff_test_seas.csv', dtype={'ID': str}).set_index('ID')
    fig, ax = plt.subplots(figsize=(9, 5), facecolor=SURFACE)
    style(ax)
    ax.grid(color=GRID, lw=0.6)
    ax.axhline(0, color=INK2, lw=0.8)
    ax.axvline(3400, color=GRID, lw=1.2, ls='--')
    ax.text(3410, 8.6, 'above 3400 m: the test\'s "high glaciers"', color=INK2, fontsize=8)
    ch = a.changed_cprec | s.changed_cprec
    ax.scatter(a.zmean[~ch], a.T_off_before[~ch], s=36, color=NEUTRAL, zorder=3, label='no snow lines (unchanged)')
    ax.scatter(a.zmean[ch], a.T_off_before[ch], s=46, facecolor='none', edgecolor=INK, lw=1.2, zorder=4,
               label='standard calibration')
    ax.scatter(a.zmean[ch], a.T_off_after[ch], s=40, color=BEFORE, zorder=5, label='snowline, all months')
    ax.scatter(s.zmean[ch], s.T_off_after[ch], s=40, color=AFTER, zorder=6, label='snowline, May–September')
    for g in a.index[ch]:
        ax.plot([a.zmean[g]] * 3, [a.T_off_before[g], a.T_off_after[g], s.T_off_after[g]], color=GRID, lw=1, zorder=2)
    for g, dx, dy in (('00773', 30, 0), ('01205', 30, 0.2), ('01230', 30, 0)):
        ax.text(a.zmean[g] + dx, a.T_off_after[g] + dy, NAMES[g], color=INK, fontsize=8, va='center')
    for g in a.index[~ch]:
        if a.T_off_before[g] > 2:
            ax.text(a.zmean[g] + 30, a.T_off_before[g], g, color=INK2, fontsize=7.5, va='center')
    ax.set_xlabel('area-weighted mean elevation [m a.s.l.]', color=INK)
    ax.set_ylabel('T_off [°C]', color=INK)
    ax.legend(fontsize=8, frameon=False, labelcolor=INK, loc='upper left')
    fig.text(0.01, 0.005, '23 calibration glaciers, monthly re-fits (GLAMOS targets for 01309, 01039, 01150).\n'
             'Grey labels: the highest glaciers without snow lines (01565 Sphinxgrat, 01150 Alphubel, 01195 Colle Gnifetti '
             'area, 00786 Mont Blanc summit).', color=INK2, fontsize=8)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    save(fig, 'ch09_toff_by_elevation')


def fig_cprec_bootstrap():
    t = pd.read_csv(SNOWCAL / 'boot_result' / 'cprec_bootstrap.csv', dtype={'ID': str}).set_index('ID')
    t = t.sort_values('c_full')
    fig, ax = plt.subplots(figsize=(8.5, 4.6), facecolor=SURFACE)
    style(ax)
    ax.grid(axis='x', color=GRID, lw=0.6)
    ax.axvspan(0.4, 3.0, color='#f3f2ee', zorder=0)
    ax.axvline(1.3, color=INK2, lw=1, ls='--', zorder=1)
    ax.text(1.32, len(t) - 0.4, 'standard calibration: 1.3\n(its lower bound)', color=INK2, fontsize=8, va='top')
    for i, (gid, r) in enumerate(t.iterrows()):
        col = AFTER if r.kept else NEUTRAL
        ax.plot([r.c_p05, r.c_p95], [i, i], color=col, lw=5, solid_capstyle='round', zorder=2)
        ax.plot([r.c_full], [i], marker='o', ms=7, color=INK if r.kept else INK2, mec=SURFACE, zorder=3)
    ax.set_yticks(range(len(t)), [f'{NAMES[g]} ({g})' for g in t.index])
    ax.set_xlim(0.3, 3.1)
    ax.set_xlabel('c_prec that fits the snow lines best', color=INK)
    handles = [plt.Line2D([], [], color=AFTER, lw=5, label='90 % range, kept'),
               plt.Line2D([], [], color=NEUTRAL, lw=5, label='90 % range reaches the end of the scan: not kept'),
               plt.Line2D([], [], marker='o', color=INK, lw=0, label='best fit (all years)')]
    ax.legend(handles=handles, fontsize=8, frameon=False, labelcolor=INK, loc='lower right')
    fig.text(0.01, 0.005, 'May–September snow lines (Sommer et al., raw, p10), one daily run per c_prec 0.4–3.0 (shaded), '
             '1000 draws over the 20 years.', color=INK2, fontsize=8)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, 'ch09_cprec_bootstrap')


def fig_daily_vs_monthly():
    import runpy
    script = REPO / 'scripts' / 'icetemp_analysis' / 'snowline' / 'daily_vs_monthly.py'
    argv, sys.argv = sys.argv, [str(script)]   # the script reads its own arguments
    try:
        g = runpy.run_path(str(script))
    finally:
        sys.argv = argv
    o, D, t, wb = g['o'], g['D'], g['t'], g['wb']
    months, lab = [5, 6, 7, 8, 9], ['May', 'Jun', 'Jul', 'Aug', 'Sep']
    bm, lo, hi = [], [], []
    for mo in months:
        s = o[o.month == mo]
        bm.append(wb(s.mon.values, s.sla_raw_m.values, s.w.values)[0])
        bd = [wb(s[c].values, s.sla_raw_m.values, s.w.values)[0] for c in D]
        lo.append(min(bd)); hi.append(max(bd))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4.8), facecolor=SURFACE, gridspec_kw=dict(width_ratios=[1, 1.15]))
    style(a1)
    a1.grid(axis='y', color=GRID, lw=0.6)
    a1.axhline(0, color=INK2, lw=0.8)
    a1.fill_between(range(5), lo, hi, color=AFTER, alpha=0.22, lw=0, label='daily, range over c_prec 0.4–3.0')
    a1.plot(range(5), bm, color=BEFORE, lw=2, marker='o', ms=8, mec=SURFACE, mew=2, label='monthly, standard calibration')
    a1.annotate(f'+{bm[-1]:.0f} m', (4, bm[-1]), xytext=(-8, 0), textcoords='offset points', ha='right', va='center',
                color=INK, fontsize=9)
    a1.set_xticks(range(5), lab)
    a1.set_ylabel('snow line, model minus image [m]', color=INK)
    a1.set_title('(a) bias by month', loc='left', color=INK, fontsize=10)
    a1.legend(fontsize=8, frameon=False, labelcolor=INK, loc='upper left')
    style(a2)
    a2.grid(axis='x', color=GRID, lw=0.6)
    tt = t.sort_values('mon_max')
    for i, (gid, r) in enumerate(tt.iterrows()):
        a2.plot([r.day_lo, r.day_hi], [i, i], color=AFTER, alpha=0.45, lw=6, solid_capstyle='round', zorder=2)
        a2.plot([r.mon_max], [i], marker='o', ms=8, color=BEFORE, mec=SURFACE, mew=2, lw=0, zorder=3)
    a2.set_yticks(range(len(tt)), [f'{NAMES[gid]}' for gid in tt.index])
    a2.set_xlabel('highest snow line of the summer [m a.s.l.]', color=INK)
    a2.set_title("(b) the summer's highest snow line", loc='left', color=INK, fontsize=10)
    handles = [plt.Line2D([], [], color=AFTER, alpha=0.45, lw=6, label='daily, range over c_prec'),
               plt.Line2D([], [], marker='o', color=BEFORE, lw=0, ms=7, label='monthly')]
    a2.legend(handles=handles, fontsize=8, frameon=False, labelcolor=INK, loc='upper left')
    fig.text(0.01, 0.005, '12 glaciers with Sommer snow lines, hydrological years 2000–2019, May–September, '
             f'{len(o)} readings weighted by 1/σ.\nMonthly model: interpolated between month ends to each image date. '
             'Daily runs: c_prec fixed at 0.4, 0.6, ... 3.0, melt fitted to the same geodetic target.',
             color=INK2, fontsize=8)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    save(fig, 'ch09_daily_vs_monthly')


if __name__ == '__main__':
    fig_snowline_definitions()
    fig_rmse_curves()
    fig_error_by_month()
    fig_toff_by_elevation()
    fig_cprec_bootstrap()
    fig_daily_vs_monthly()
