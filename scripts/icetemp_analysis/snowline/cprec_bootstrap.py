"""A1 of the plan of 2026-10-16 (afternoon): how well do the snow lines fix c_prec? A bootstrap over years.

For every scanned c_prec there is one snowline run with that single value (DDF and T_off fitted as in the
scan) and a daily forward run writing the p10 snow line (boot/c<c>, boot/f<c>). Per glacier:
  RMSE(c) = sqrt(sum w r^2 / sum w), w = 1/sigma, r = model - observed snow line (May-September, raw)
  full-sample c_prec: minimum of a natural cubic spline through RMSE(c) on 1000 points (the whole range,
                      not the scan's first rise)
  bootstrap: the hydrological years 2000-2019 resampled with replacement (1000 draws, seed fixed), the
             spline minimum of each draw -> 5-95 % range
  kept (rule fixed in the plan): the range reaches neither end of the scan (0.4, 3.0)
Writes cprec_bootstrap.csv and, for the kept glaciers, a calibrate_cprec_file whose Cprec is carried over
by precipitation amount (P_daily(1)/P_monthly(1), as cprec_transfer.py); only ID and Cprec are read by GloGEM,
the other columns are placeholders.

Usage: python cprec_bootstrap.py BOOT_DIR OBS_FILE PRECIP_DAILY PRECIP_MONTHLY OUT_DIR [--n 1000]
"""
import argparse
import glob
import os
import re

import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

SUB, CC = 'centraleurope', 'CentralEurope_snowcal'
SNOW_MULTIPLIER = 1.2


def read_snowline_day(fwd):
    fn = f'{fwd}/daily/CentralEurope/PAST/PAST_original/{SUB}_Snowline_day_r1_{CC}.dat'
    mod = {}
    for ln in open(fn):
        p = ln.split()
        if p[0] == 'ID' or len(p) < 300:
            continue
        mod.setdefault(p[0], {})[int(p[1])] = np.array(p[3:], float)
    return mod


def spline_min(x, y):
    xf = np.linspace(x[0], x[-1], 1000)
    return xf[int(np.argmin(CubicSpline(x, y, bc_type='natural')(xf)))]


def annual(run, tres, var, y0=2001, y1=2019):
    d = pd.read_csv(f'{run}/{tres}/CentralEurope/PAST/PAST_original/{SUB}_{var}_r1_{CC}.dat', sep=r'\s+',
                    dtype={'ID': str}).set_index('ID')
    d.columns = d.columns.astype(int)
    return d.loc[:, y0:y1].mean(axis=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('boot')
    ap.add_argument('obs')
    ap.add_argument('precip_daily')
    ap.add_argument('precip_monthly')
    ap.add_argument('out_dir')
    ap.add_argument('--n', type=int, default=1000)
    a = ap.parse_args()

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
        w = og.w.values
        sw2 = (R ** 2 * w)                            # per observation, per c
        # per-year sums, so a bootstrap draw is a weighted sum over years
        yi = np.searchsorted(years, og.hy.values)
        S = np.zeros((len(cs), len(years)))
        W = np.zeros(len(years))
        for j in range(len(years)):
            sel = yi == j
            S[:, j] = sw2[:, sel].sum(axis=1)
            W[j] = w[sel].sum()
        rmse = np.sqrt(S.sum(axis=1) / W.sum())
        c_full = spline_min(cs, rmse)
        draws = np.empty(a.n)
        for b in range(a.n):
            cnt = np.bincount(rng.integers(0, len(years), len(years)), minlength=len(years))
            draws[b] = spline_min(cs, np.sqrt((S * cnt).sum(axis=1) / (W * cnt).sum()))
        lo, hi = np.percentile(draws, [5, 95])
        at_end = np.mean((draws <= cs[0] + 1e-9) | (draws >= cs[-1] - 1e-9))
        kept = bool(lo > cs[0] + 1e-9 and hi < cs[-1] - 1e-9)
        rows.append(dict(ID=gid, n_obs=len(og), c_full=c_full, c_p05=lo, c_p95=hi, share_draws_at_end=at_end,
                         kept=kept, rmse_min=rmse.min(), rmse_at_04=rmse[0],
                         **{f'rmse_{c:.1f}': v for c, v in zip(cs, rmse)}))
    t = pd.DataFrame(rows).set_index('ID')

    p_d = annual(a.precip_daily, 'daily', 'Accumulation_sfc') / SNOW_MULTIPLIER + annual(a.precip_daily, 'daily', 'Rain_sfc')
    p_m = annual(a.precip_monthly, 'monthly', 'Accumulation_sfc') / SNOW_MULTIPLIER + annual(a.precip_monthly, 'monthly', 'Rain_sfc')
    t['ratio'] = (p_d / p_m).reindex(t.index)
    t['c_monthly'] = t.c_full * t.ratio
    os.makedirs(a.out_dir, exist_ok=True)
    t.to_csv(f'{a.out_dir}/cprec_bootstrap.csv', float_format='%.4f')
    with open(f'{a.out_dir}/cprec_fixed_bootstrap_{CC}.dat', 'w') as f:
        f.write('ID        Ba         Bw     Area     ELA   AAR    dBdz   Bt    DDFsnow  DDFice   Cprec   T_off  Flag\n')
        for gid, r in t[t.kept].iterrows():
            f.write(f'{gid} {0:8.3f} {0:8.3f} {0:10.3f} {0:5d} {0:5.1f} {0:8.3f} {0:8.3f} {0:6.3f} {0:6.3f} '
                    f'{r.c_monthly:8.4f} {0:6.2f} {1:2d}\n')
    pd.set_option('display.width', 200)
    print(t[['n_obs', 'c_full', 'c_p05', 'c_p95', 'share_draws_at_end', 'kept', 'rmse_min', 'rmse_at_04', 'ratio',
             'c_monthly']].round(3).to_string())
    print(f'kept {int(t.kept.sum())} of {len(t)}')


if __name__ == '__main__':
    main()
