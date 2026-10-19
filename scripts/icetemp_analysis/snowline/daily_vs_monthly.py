"""Daily versus monthly GloGEM snow lines on the same satellite readings (book chapter 9; plan of 16 Oct, evening).

Monthly: fwdB_none_fixed (standard monthly calibration, hydrological years 1999-2025), Snowline_mon interpolated
between month ends to each reading's date. Daily: boot/f<c> (c_prec fixed at each scan value 0.4-3.0, melt fitted
to the same geodetic target, 1999-2019), the p10 snow line of the reading's day.
Readings: Sommer TSLA sla_raw_m, qa >= 0.666, May-September, hydrological years 2000-2019, weights 1/sigma.
Prints (1) the weighted bias and RMSE by month, (2) the summer's highest snow line (July-September).

Usage: python daily_vs_monthly.py [SNOWCAL_DIR] [TSLA_DIR]
"""
import calendar
import glob
import os
import sys

import numpy as np
import pandas as pd

S = sys.argv[1] if len(sys.argv) > 1 else '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week2/snowcal'
T = sys.argv[2] if len(sys.argv) > 2 else \
    '/home/jabeer/projects/glogemflow_development/glogemflow_icetemp/data/snowlines/tsla_csv'
TSLA = ['00773', '00774', '01112', '01127', '01205', '01225', '01230', '02216', '02596', '02611', '02821', '03116']
Y0_MON, NY_MON = 1999, 27   # monthly run: hydrological years 1999-2025
Y0_DAY = 1999               # daily runs: hydrological years 1999-2019, 1999 the spin-up year
YEARS = (2000, 2019)


def monthly(fn):
    out = {}
    for ln in open(fn).readlines()[1:]:
        p = ln.split()
        if p:
            out[p[0]] = np.array(p[1:], float).reshape(NY_MON, 12)
    return out


def month_end(arr, year, month):
    i = (year + 1 if month >= 10 else year) - Y0_MON
    return arr[i, (month - 10) % 12] if 0 <= i < NY_MON else np.nan


def interp(arr, t):
    prev = t.replace(day=1) - pd.Timedelta(days=1)
    a, b = month_end(arr, prev.year, prev.month), month_end(arr, t.year, t.month)
    return a + (b - a) * t.day / calendar.monthrange(t.year, t.month)[1]


def daily(fwd):
    m = {}
    for ln in open(f'{fwd}/daily/CentralEurope/PAST/PAST_original/'
                   'centraleurope_Snowline_day_r1_CentralEurope_snowcal.dat'):
        p = ln.split()
        if p[0] == 'ID' or len(p) < 300:
            continue
        m.setdefault(p[0], {})[int(p[1])] = np.array(p[3:], float)
    return m


def wb(mod, obs, w):
    ok = np.isfinite(mod)
    r, w = mod[ok] - obs[ok], w[ok]
    return np.sum(w * r) / w.sum(), np.sqrt(np.sum(w * r * r) / w.sum())


mon = monthly(f'{S}/fwdB_none_fixed/monthly/CentralEurope/PAST/PAST_original/'
              'centraleurope_Snowline_mon_r1_CentralEurope_fair.dat')
runs = sorted(glob.glob(f'{S}/boot/f*'), key=lambda p: float(os.path.basename(p)[1:]))
cs = [float(os.path.basename(p)[1:]) for p in runs]
day = [daily(p) for p in runs]

obs = []
for gid in TSLA:
    d = pd.read_csv(f'{T}/tsla_RGI70-11.{gid}.csv', parse_dates=['date'])
    d = d[(d.qa >= 0.666) & d.month.between(5, 9) & d.hydro_year.between(*YEARS)].copy()
    d['gid'] = gid
    obs.append(d)
o = pd.concat(obs, ignore_index=True)
o['w'] = 1 / o.sigma_m
o['mon'] = [interp(mon[g], t) for g, t in zip(o.gid, o.date)]
doy = o.date.dt.dayofyear.values.copy()
doy[o.date.dt.is_leap_year.values & (doy > 59)] -= 1          # the runs' 365-day calendar
k = 365 * (o.date.dt.year.values - Y0_DAY + 1) + doy - 274     # day index from 1 Oct of the first run year
o['col'] = k % 365
assert ((Y0_DAY + k // 365) == o.hydro_year).all()
D = []
for c, m in zip(cs, day):
    o[f'd{c:.1f}'] = [m[g][h][j] for g, h, j in zip(o.gid, o.hydro_year, o.col)]
    D.append(f'd{c:.1f}')

print(f'{o.gid.nunique()} glaciers, {len(o)} readings, hydrological years {YEARS[0]}-{YEARS[1]}')
print('month     n  bias_mon  bias_day(min..max)  wrmse_mon  wrmse_day(min..max)')
for mo, s in o.groupby('month'):
    bm, rm = wb(s.mon.values, s.sla_raw_m.values, s.w.values)
    bd, rd = zip(*[wb(s[c].values, s.sla_raw_m.values, s.w.values) for c in D])
    print(f'{mo:5d} {len(s):5d} {bm:8.0f}  {min(bd):7.0f}..{max(bd):5.0f}  {rm:9.0f}  {min(rd):7.0f}..{max(rd):5.0f}')

# the summer's highest snow line, July-September
rows = []
for g in TSLA:
    yrs = range(YEARS[0], YEARS[1] + 1)
    m_max = np.mean([mon[g][y - Y0_MON, 9:12].max() for y in yrs])
    d_max = [np.mean([m[g][y][273:365].max() for y in yrs]) for m in day]
    s = o[(o.gid == g) & (o.month >= 7)]
    n = s.groupby('hydro_year').size()
    grp = s[s.hydro_year.isin(n[n >= 3].index)].groupby('hydro_year')
    rows.append(dict(gid=g, mon_max=m_max, day_lo=min(d_max), day_hi=max(d_max), n_years=grp.ngroups,
                     obs_max=grp.sla_raw_m.max().mean(), mon_on_obs=grp.mon.max().mean(),
                     day_on_obs_lo=min(grp[c].max().mean() for c in D),
                     day_on_obs_hi=max(grp[c].max().mean() for c in D)))
t = pd.DataFrame(rows).set_index('gid')
print(t.round(0).to_string())
inside = int(((t.mon_max >= t.day_lo) & (t.mon_max <= t.day_hi)).sum())
print(f'monthly minus the middle of the daily range: {(t.mon_max - (t.day_lo + t.day_hi) / 2).mean():.0f} m; '
      f'monthly inside the daily range on {inside} of {len(t)} glaciers')
