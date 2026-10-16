"""T_off test of the snowline recalibration, fixed in the plan of 2026-10-16 before any result:
does the recalibration raise T_off, which made high firn too warm?

Compared, before and after (two monthly calibrations of the same glaciers and targets):
  median T_off, share with T_off > 2 C, share at both bounds (c_prec at the lower bound of its
  phase-1 range or fixed at or below it, and DDF at 4.5); all glaciers and those with an area-weighted
  mean elevation above 3400 m; T_off against mean elevation (slope per 100 m).
Flag: the median T_off of the high glaciers rises, or more of them reach both bounds.

Usage: python toff_test.py BEFORE_FINAL_CAL AFTER_FINAL_CAL OUT_CSV [--c1-low 1.3] [--high 3400]
"""
import argparse

import numpy as np
import pandas as pd

BANDS = '/itet-stor/jabeer/glogem/geometricdata/rgiv7/bands/centraleurope'
DDF_MAX = 4.5


def mean_elev(gid):
    b = np.loadtxt(f'{BANDS}/{gid}.dat', skiprows=5, usecols=(1, 2, 3))
    return float(np.sum((b[:, 0] + b[:, 1]) / 2 * b[:, 2]) / b[:, 2].sum())


def read(fn):
    return pd.read_csv(fn, sep=r'\s+', dtype={'ID': str}).set_index('ID')


def summary(t, tag, c1_low, sel):
    x = t[sel]
    both = (x[f'Cprec_{tag}'] <= c1_low + 0.005) & (x[f'DDFsnow_{tag}'] >= DDF_MAX - 0.005)
    slope = np.polyfit(x.zmean, x[f'T_off_{tag}'], 1)[0] * 100 if len(x) > 2 else np.nan
    return dict(n=len(x), median_toff=x[f'T_off_{tag}'].median(), share_toff_gt2=(x[f'T_off_{tag}'] > 2).mean(),
                share_both_bounds=both.mean(), n_both_bounds=int(both.sum()), toff_per_100m=slope)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('before')
    ap.add_argument('after')
    ap.add_argument('out_csv')
    ap.add_argument('--c1-low', type=float, default=1.3)
    ap.add_argument('--high', type=float, default=3400.)
    a = ap.parse_args()
    b, f = read(a.before), read(a.after)
    cols = ['Cprec', 'DDFsnow', 'T_off', 'Flag']
    t = b[cols].add_suffix('_before').join(f[cols].add_suffix('_after'), how='inner')
    t['zmean'] = [mean_elev(g) for g in t.index]
    t['changed_cprec'] = (t.Cprec_after - t.Cprec_before).abs() > 0.0005
    t.to_csv(a.out_csv, float_format='%.3f')
    rows = []
    for name, sel in (('all', t.zmean > -1), (f'mean elevation > {a.high:.0f} m', t.zmean > a.high),
                      ('c_prec changed', t.changed_cprec)):
        for tag in ('before', 'after'):
            rows.append(dict(sample=name, run=tag, **summary(t, tag, a.c1_low, sel)))
    s = pd.DataFrame(rows)
    pd.set_option('display.width', 200)
    print(t.round(3).sort_values('zmean').to_string())
    print()
    print(s.round(3).to_string(index=False))
    hi = s[s['sample'].str.startswith('mean elevation')].set_index('run')
    flag = (hi.loc['after', 'median_toff'] > hi.loc['before', 'median_toff'] + 1e-9
            or hi.loc['after', 'n_both_bounds'] > hi.loc['before', 'n_both_bounds'])
    print('\nFLAG (high glaciers: median T_off rises or more at both bounds):', bool(flag))


if __name__ == '__main__':
    main()
