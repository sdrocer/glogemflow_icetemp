"""How well do the snow lines fix c_prec? A bootstrap over years (plans of 2026-10-16).

Daily (A1, afternoon plan): for every scanned c_prec there is one snowline run with that single value (DDF
and T_off fitted as in the scan) and a daily forward run writing the p10 snow line (boot/c<c>, boot/f<c>).
Monthly (the monthly snow-line test, evening plan): for every c_prec a monthly calibration with the c_prec of
the 12 snow-line glaciers fixed at that value and a monthly forward run 1999-2025 (mtest/c<c>, mtest/f<c>);
the monthly snow line is interpolated between month ends to each reading's date, and only the readings up to
each glacier-year's highest one are used (the rise of the snow line; the September drop after it is out of
reach of a monthly model, chapter 9).
Per glacier:
  RMSE(c) = sqrt(sum w r^2 / sum w), w = 1/sigma, r = model - observed snow line
  full-sample c_prec: minimum of a natural cubic spline through RMSE(c) on 1000 points (the whole range,
                      not the scan's first rise)
  bootstrap: the hydrological years resampled with replacement (1000 draws, seed fixed), the spline minimum
             of each draw -> 5-95 % range
  kept (rule fixed in the plan): the range reaches neither end of the scan (0.4, 3.0)
Writes cprec_bootstrap[_monthly].csv and, for the kept glaciers, a calibrate_cprec_file. Daily: its Cprec is
carried over to the monthly model by precipitation amount (P_daily(1)/P_monthly(1), as cprec_transfer.py);
monthly: the value itself. Only ID and Cprec are read by GloGEM, the other columns are placeholders.

Usage: python cprec_bootstrap.py BOOT_DIR OBS_FILE PRECIP_DAILY PRECIP_MONTHLY OUT_DIR [--n 1000]
       python cprec_bootstrap.py monthly MTEST_DIR TSLA_DIR OUT_DIR [--n 1000] [--years 2000 2025] [--all-readings]
"""
import argparse
import calendar
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

SUB, CC = 'centraleurope', 'CentralEurope_snowcal'
SNOW_MULTIPLIER = 1.2
TSLA = ['00773', '00774', '01112', '01127', '01205', '01225', '01230', '02216', '02596', '02611', '02821', '03116']
Y0_MON, NY_MON = 1999, 27   # monthly forward runs: hydrological years 1999-2025


def read_snowline_day(fwd):
    fn = f'{fwd}/daily/CentralEurope/PAST/PAST_original/{SUB}_Snowline_day_r1_{CC}.dat'
    mod = {}
    for ln in open(fn):
        p = ln.split()
        if p[0] == 'ID' or len(p) < 300:
            continue
        mod.setdefault(p[0], {})[int(p[1])] = np.array(p[3:], float)
    return mod


def read_snowline_mon(fwd):
    fn = f'{fwd}/monthly/CentralEurope/PAST/PAST_original/{SUB}_Snowline_mon_r1_CentralEurope_fair.dat'
    mod = {}
    for ln in open(fn).readlines()[1:]:
        p = ln.split()
        if p:
            mod[p[0]] = np.array(p[1:], float).reshape(NY_MON, 12)
    return mod


def month_end(arr, year, month):
    i = (year + 1 if month >= 10 else year) - Y0_MON
    return arr[i, (month - 10) % 12] if 0 <= i < NY_MON else np.nan


def interp(arr, t):
    prev = t.replace(day=1) - pd.Timedelta(days=1)
    a, b = month_end(arr, prev.year, prev.month), month_end(arr, t.year, t.month)
    return a + (b - a) * t.day / calendar.monthrange(t.year, t.month)[1]


def spline_min(x, y):
    xf = np.linspace(x[0], x[-1], 1000)
    return xf[int(np.argmin(CubicSpline(x, y, bc_type='natural')(xf)))]


def bootstrap(R, w, yi, ny, cs, rng, n):
    """R (nc, nobs) residuals, w weights, yi year index per reading -> RMSE(c), full-sample c, draws."""
    sw2 = R ** 2 * w                             # per observation, per c
    S = np.zeros((len(cs), ny))                  # per-year sums, so a draw is a weighted sum over years
    W = np.zeros(ny)
    for j in range(ny):
        sel = yi == j
        S[:, j] = sw2[:, sel].sum(axis=1)
        W[j] = w[sel].sum()
    rmse = np.sqrt(S.sum(axis=1) / W.sum())
    draws = np.empty(n)
    for b in range(n):
        cnt = np.bincount(rng.integers(0, ny, ny), minlength=ny)
        draws[b] = spline_min(cs, np.sqrt((S * cnt).sum(axis=1) / (W * cnt).sum()))
    return rmse, spline_min(cs, rmse), draws


def summary(gid, n_obs, cs, rmse, c_full, draws):
    lo, hi = np.percentile(draws, [5, 95])
    at_end = np.mean((draws <= cs[0] + 1e-9) | (draws >= cs[-1] - 1e-9))
    kept = bool(lo > cs[0] + 1e-9 and hi < cs[-1] - 1e-9)
    return dict(ID=gid, n_obs=n_obs, c_full=c_full, c_p05=lo, c_p95=hi, share_draws_at_end=at_end,
                kept=kept, rmse_min=rmse.min(), rmse_at_04=rmse[0], **{f'rmse_{c:.1f}': v for c, v in zip(cs, rmse)})


def write_fixed(t, col, fn):
    with open(fn, 'w') as f:
        f.write('ID        Ba         Bw     Area     ELA   AAR    dBdz   Bt    DDFsnow  DDFice   Cprec   T_off  Flag\n')
        for gid, r in t[t.kept].iterrows():
            f.write(f'{gid} {0:8.3f} {0:8.3f} {0:10.3f} {0:5d} {0:5.1f} {0:8.3f} {0:8.3f} {0:6.3f} {0:6.3f} '
                    f'{r[col]:8.4f} {0:6.2f} {1:2d}\n')


def annual(run, tres, var, y0=2001, y1=2019):
    d = pd.read_csv(f'{run}/{tres}/CentralEurope/PAST/PAST_original/{SUB}_{var}_r1_{CC}.dat', sep=r'\s+',
                    dtype={'ID': str}).set_index('ID')
    d.columns = d.columns.astype(int)
    return d.loc[:, y0:y1].mean(axis=1)


def rising_readings(tsla_dir, y0, y1, all_readings=False):
    """Sommer raw snow lines (qa >= 0.666, May-September), per glacier-year up to the highest reading."""
    out = []
    for gid in TSLA:
        d = pd.read_csv(f'{tsla_dir}/tsla_RGI70-11.{gid}.csv', parse_dates=['date'])
        d = d[(d.qa >= 0.666) & d.month.between(5, 9) & d.hydro_year.between(y0, y1)].copy()
        if not all_readings:
            top = d.loc[d.groupby('hydro_year').sla_raw_m.idxmax()].set_index('hydro_year').date
            d = d[d.date <= d.hydro_year.map(top)]
        d['gid'] = gid
        out.append(d)
    return pd.concat(out, ignore_index=True)


def main_daily(argv):
    ap = argparse.ArgumentParser(description='daily bootstrap (A1)')
    ap.add_argument('boot')
    ap.add_argument('obs')
    ap.add_argument('precip_daily')
    ap.add_argument('precip_monthly')
    ap.add_argument('out_dir')
    ap.add_argument('--n', type=int, default=1000)
    a = ap.parse_args(argv)

    fwd = sorted(glob.glob(f'{a.boot}/f*'), key=lambda p: float(re.sub(r'^f', '', os.path.basename(p))))
    cs = np.array([float(os.path.basename(p)[1:]) for p in fwd])
    mods = [read_snowline_day(p) for p in fwd]
    t0 = 1999   # first run year (one spin-up year before the target period)
    o = pd.read_csv(a.obs, sep=r'\s+', skiprows=1, header=None, names=['rid', 'sla', 'sig', 'year', 'doy'])
    o['gid'] = o.rid.str[-5:]
    k = 365 * (o.year - t0 + 1) + o.doy - 274
    o['hy'], o['col'] = t0 + k // 365, k % 365
    o = o[o.hy.between(2000, 2019)].reset_index(drop=True)
    o['w'] = 1 / o.sig

    rng = np.random.default_rng(20261016)
    years = np.arange(2000, 2020)
    rows = []
    for gid, og in o.groupby('gid'):
        R = np.array([[m[gid][h][c] for h, c in zip(og.hy, og.col)] for m in mods]) - og.sla.values  # (nc, nobs)
        rmse, c_full, draws = bootstrap(R, og.w.values, np.searchsorted(years, og.hy.values), len(years), cs, rng, a.n)
        rows.append(summary(gid, len(og), cs, rmse, c_full, draws))
    t = pd.DataFrame(rows).set_index('ID')

    p_d = annual(a.precip_daily, 'daily', 'Accumulation_sfc') / SNOW_MULTIPLIER + annual(a.precip_daily, 'daily', 'Rain_sfc')
    p_m = annual(a.precip_monthly, 'monthly', 'Accumulation_sfc') / SNOW_MULTIPLIER + annual(a.precip_monthly, 'monthly', 'Rain_sfc')
    t['ratio'] = (p_d / p_m).reindex(t.index)
    t['c_monthly'] = t.c_full * t.ratio
    os.makedirs(a.out_dir, exist_ok=True)
    t.to_csv(f'{a.out_dir}/cprec_bootstrap.csv', float_format='%.4f')
    write_fixed(t, 'c_monthly', f'{a.out_dir}/cprec_fixed_bootstrap_{CC}.dat')
    pd.set_option('display.width', 200)
    print(t[['n_obs', 'c_full', 'c_p05', 'c_p95', 'share_draws_at_end', 'kept', 'rmse_min', 'rmse_at_04', 'ratio',
             'c_monthly']].round(3).to_string())
    print(f'kept {int(t.kept.sum())} of {len(t)}')


def main_monthly(argv):
    ap = argparse.ArgumentParser(description='monthly bootstrap (the monthly snow-line test)')
    ap.add_argument('mtest')
    ap.add_argument('tsla_dir')
    ap.add_argument('out_dir')
    ap.add_argument('--n', type=int, default=1000)
    ap.add_argument('--years', type=int, nargs=2, default=(2000, 2025))
    ap.add_argument('--all-readings', action='store_true', help='all May-September readings, not only the rise')
    a = ap.parse_args(argv)

    fwd = sorted((p for p in glob.glob(f'{a.mtest}/f*') if os.path.isdir(p)), key=lambda p: float(os.path.basename(p)[1:]))
    cs = np.array([float(os.path.basename(p)[1:]) for p in fwd])
    mods = [read_snowline_mon(p) for p in fwd]
    o = rising_readings(a.tsla_dir, *a.years, all_readings=a.all_readings)
    o['w'] = 1 / o.sigma_m

    rng = np.random.default_rng(20261019)
    years = np.arange(a.years[0], a.years[1] + 1)
    rows = []
    for gid, og in o.groupby('gid'):
        R = np.array([[interp(m[gid], d) for d in og.date] for m in mods]) - og.sla_raw_m.values
        rmse, c_full, draws = bootstrap(R, og.w.values, np.searchsorted(years, og.hydro_year.values), len(years), cs,
                                        rng, a.n)
        rows.append(summary(gid, len(og), cs, rmse, c_full, draws))
    t = pd.DataFrame(rows).set_index('ID')
    tag = '_all' if a.all_readings else ''
    os.makedirs(a.out_dir, exist_ok=True)
    t.to_csv(f'{a.out_dir}/cprec_bootstrap_monthly{tag}.csv', float_format='%.4f')
    write_fixed(t, 'c_full', f'{a.out_dir}/cprec_fixed_monthly{tag}_CentralEurope_fair.dat')
    pd.set_option('display.width', 200)
    print(t[['n_obs', 'c_full', 'c_p05', 'c_p95', 'share_draws_at_end', 'kept', 'rmse_min', 'rmse_at_04']]
          .round(3).to_string())
    print(f'kept {int(t.kept.sum())} of {len(t)}; readings {len(o)}, years {a.years[0]}-{a.years[1]}'
          f'{", all May-September" if a.all_readings else ", up to each year`s highest"}')


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'monthly':
        main_monthly(sys.argv[2:])
    else:
        main_daily(sys.argv[1:])
