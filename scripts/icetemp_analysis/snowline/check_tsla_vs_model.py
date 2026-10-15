"""Check that the prepared TSLA files reproduce yesterday's model-vs-observation scores
(snowline/out/A_snowline_stats_qa0.67.csv, regional run): weights 1/sigma, months May-Oct,
model = end-of-month snowline of the observation month, observations 2000-01-01..2025-09-30.

Yesterday's filter QA_FLAG >= 0.67 kept only QA 5/6 and 1 (QA values are k/6), so the
reproduction uses qa >= 0.67 on the prepared files; extra variants are reported alongside.
"""
import calendar
import glob
import os

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
DATA = f'{REPO}/data/snowlines'
OUT = DATA

REF = f'{DATA}/reference/A_snowline_stats_qa0.67.csv'
MOD = ('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week2/mbout/regional/monthly/CentralEurope/'
       'PAST/PAST_original/centraleurope_Snowline_mon_r1_CentralEurope_fair.dat')
Y0, NY = 1990, 36   # first hydrological year; columns Oct..Sep, end-of-month state


def read_model():
    out = {}
    with open(MOD) as f:
        next(f)
        for ln in f:
            p = ln.split()
            if p:
                out[p[0]] = np.array(p[1:], float).reshape(NY, 12)
    return out


def model_month(arr, year, month):
    Y = year + 1 if month >= 10 else year
    i = Y - Y0
    return arr[i, (month - 10) % 12] if 0 <= i < NY else np.nan


def model_interp(arr, t):
    nd = calendar.monthrange(t.year, t.month)[1]
    prev = t.replace(day=1) - pd.Timedelta(days=1)
    a, b = model_month(arr, prev.year, prev.month), model_month(arr, t.year, t.month)
    return a + (b - a) * t.day / nd


def wstats(m, o, s):
    ok = np.isfinite(m) & np.isfinite(o)
    w, d = 1 / s[ok], m[ok] - o[ok]
    return np.sqrt(np.sum(w * d * d) / w.sum()), np.sum(w * d) / w.sum(), int(ok.sum())


def main():
    mod = read_model()
    ref = pd.read_csv(REF, dtype={'gid': str})
    ref = ref[ref.run == 'regional'].set_index('gid')
    rows = []
    for fn in sorted(glob.glob(f'{OUT}/tsla_csv/tsla_RGI70-11.*.csv')):
        gid = fn[-9:-4]
        d = pd.read_csv(fn, parse_dates=['date'])
        arr = mod[gid]
        d['mod'] = [model_month(arr, t.year, t.month) for t in d.date]
        d['mod_int'] = [model_interp(arr, t) for t in d.date]
        base = (d.date <= '2025-09-30') & d.month.between(5, 10)
        rec = dict(gid=gid)
        for tag, sel, col in [('repr', base & (d.qa >= 0.67), 'mod'),
                              ('repr_int', base & (d.qa >= 0.67), 'mod_int'),
                              ('qa4of6', base & (d.qa >= 0.666), 'mod'),
                              ('cal_AprNov_qa4of6', d.hydro_year.between(2000, 2019) & (d.qa >= 0.666), 'mod')]:
            r, b, n = wstats(d[col].values[sel], d.sla_m.values[sel], d.sigma_m.values[sel])
            rec.update({f'{tag}_wrmse': r, f'{tag}_bias': b, f'{tag}_n': n})
        rr = ref.loc[gid]
        rec.update(ref_wrmse=rr.rmse_MJJASO, ref_bias=rr.bias_MJJASO, ref_n=int(rr.n_MJJASO),
                   ref_int_wrmse=rr.int_rmse_MJJASO, ref_int_bias=rr.int_bias_MJJASO)
        rec['d_wrmse'] = rec['repr_wrmse'] - rec['ref_wrmse']
        rec['d_bias'] = rec['repr_bias'] - rec['ref_bias']
        rec['d_int_wrmse'] = rec['repr_int_wrmse'] - rec['ref_int_wrmse']
        rec['match_1m'] = bool(abs(rec['d_wrmse']) < 1 and abs(rec['d_bias']) < 1 and rec['repr_n'] == rec['ref_n']
                               and abs(rec['d_int_wrmse']) < 1)
        rows.append(rec)
    t = pd.DataFrame(rows)
    t.to_csv(f'{OUT}/check_vs_yesterday.csv', index=False, float_format='%.2f')
    pd.set_option('display.width', 250, 'display.max_columns', 40)
    print(t.round(1).to_string(index=False))
    print('all match to 1 m:', t.match_1m.all())


if __name__ == '__main__':
    main()
