"""Carry the c_prec of a daily snowline calibration into the monthly model, by precipitation amount.

c_prec multiplies the reanalysis precipitation. The daily and monthly inputs can differ (daily
amounts below P_threshold are dropped, and the two paths may pick different grid points), so the
same factor can put different amounts on the glacier. The monthly factor is the daily one times
P_daily(1) / P_monthly(1): the glacier-wide precipitation at c_prec = 1 over the target years
(snowfall without the snow multiplier, plus rain), from one forward run in each mode.

Writes a GloGEM calibration file for calibrate_cprec_file (only ID and Cprec are read; the other
columns are copied from the daily calibration) and a table of the transfer.

Usage: python cprec_transfer.py DAILY_CAL_DIR PRECIP_DAILY_DIR PRECIP_MONTHLY_DIR OUT_DIR
         [--catchment CentralEurope_snowcal] [--years 2001 2019] [--snow-multiplier 1.2]
"""
import argparse
import os

import pandas as pd

REGION, SUB = 'CentralEurope', 'centraleurope'


def annual(run, tres, var, cc, y0, y1):
    fn = f'{run}/{tres}/{REGION}/PAST/PAST_original/{SUB}_{var}_r1_{cc}.dat'
    d = pd.read_csv(fn, sep=r'\s+', dtype={'ID': str}).set_index('ID')
    d.columns = d.columns.astype(int)
    return d.loc[:, y0:y1].mean(axis=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('daily_cal')
    ap.add_argument('precip_daily')
    ap.add_argument('precip_monthly')
    ap.add_argument('out_dir')
    ap.add_argument('--catchment', default='CentralEurope_snowcal')
    ap.add_argument('--years', type=int, nargs=2, default=(2001, 2019))
    ap.add_argument('--snow-multiplier', type=float, default=1.2)
    a = ap.parse_args()
    cc, (y0, y1) = a.catchment, a.years

    p = {}
    for tres, run in (('daily', a.precip_daily), ('monthly', a.precip_monthly)):
        p[tres] = (annual(run, tres, 'Accumulation_sfc', cc, y0, y1) / a.snow_multiplier
                   + annual(run, tres, 'Rain_sfc', cc, y0, y1))
    cal = f'{a.daily_cal}/daily/{REGION}/calibration'
    head = open(f'{cal}/calibrate_m1_cID9_{SUB}_final_era5_{cc}.dat').readline()
    fin = pd.read_csv(f'{cal}/calibrate_m1_cID9_{SUB}_final_era5_{cc}.dat', sep=r'\s+', dtype={'ID': str})
    log = pd.read_csv(f'{cal}/calibrate_m1_cID9_{SUB}_{cc}_snowline.dat', sep=r'\s+', dtype={'ID': str})
    t = fin.merge(log[['ID', 'Status', 'Nscan', 'CprecMin', 'Bound', 'RMSEfinal']], on='ID')
    t = t[t.Status == 2].copy()
    t['P_daily_1'] = t.ID.map(p['daily'])
    t['P_monthly_1'] = t.ID.map(p['monthly'])
    t['ratio'] = t.P_daily_1 / t.P_monthly_1
    t['Cprec_daily'] = t.Cprec
    t['Cprec_monthly'] = t.Cprec * t.ratio
    if t.ratio.isna().any():
        raise SystemExit(f'no precipitation for {list(t.ID[t.ratio.isna()])}')

    os.makedirs(a.out_dir, exist_ok=True)
    with open(f'{a.out_dir}/cprec_fixed_{cc}.dat', 'w') as f:
        f.write(head)
        for _, r in t.iterrows():
            f.write(f'{r.ID} {r.Ba:8.3f} {r.Bw:8.3f} {r.Area:10.3f} {int(r.ELA):5d} {r.AAR:5.1f} {r.dBdz:8.3f} '
                    f'{r.Bt:8.3f} {r.DDFsnow:6.3f} {r.DDFice:6.3f} {r.Cprec_monthly:8.4f} {r.T_off:6.2f} {int(r.Flag):2d}\n')
    cols = ['ID', 'Nscan', 'Bound', 'RMSEfinal', 'Cprec_daily', 'P_daily_1', 'P_monthly_1', 'ratio', 'Cprec_monthly',
            'DDFsnow', 'T_off', 'Flag']
    t[cols].to_csv(f'{a.out_dir}/cprec_transfer_{cc}.csv', index=False, float_format='%.4f')
    print(t[cols].round(3).to_string(index=False))


if __name__ == '__main__':
    main()
