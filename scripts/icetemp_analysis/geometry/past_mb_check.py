"""Test 1c, part 1 (plan of 20 Oct): why is GloGEM's mass balance before 2000 too positive: too much winter snow or too
little summer melt? The 8 GLAMOS glaciers observed since before 1972 (mass-balance run with the stage-2 thermal
settings and calibration, week3/geom/mb_glamoslong) against their GLAMOS fixed-date series: annual balance Ba, winter
balance Bw (only years with a measured winter balance), summer balance Bs = Ba - Bw, by period.

Usage: python scripts/icetemp_analysis/geometry/past_mb_check.py [OUT_CSV]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

FIX = Path('/itet-stor/jabeer/glazioarch/GlacioBaseData/MassBalance/glacier')
RUN = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3/geom/mb_glamoslong/monthly/CentralEurope/PAST/'
           'PAST_original')
GLACIERS = {'gietro': '00879', 'schwarzberg': '01084', 'hohlaub': '01144', 'allalin': '01182', 'clariden': '01838',
            'silvretta': '01949', 'gries': '02438', 'aletsch': '02596'}
FLAG = {'clariden': 'RGI7 outline holds 43 % of the GLAMOS glacier'}
MEAS_W = (1, 2, 4, 5)       # GLAMOS methods with a measured winter balance (as compare_snow.py)
PERIODS = ((1972, 1985), (1986, 1999), (2000, 2019), (2020, 2025))


def glamos(name):
    rows = []
    for line in open(FIX / f'{name}_fix.dat', encoding='latin-1'):
        p = line.split()
        if line.startswith('#') or len(p) < 9:
            continue
        val = lambda s: np.nan if s.lower() == 'nan' else float(s.rstrip('.')) / 1000
        try:
            rows.append(dict(year=int(p[4][:4]), meth=int(p[0]), Bw=val(p[5]), Ba=val(p[6])))
        except ValueError:
            continue
    d = pd.DataFrame(rows).drop_duplicates('year').set_index('year')
    d.loc[~d.meth.isin(MEAS_W), 'Bw'] = np.nan
    return d


def model(var):
    d = pd.read_csv(RUN / f'centraleurope_{var}_r1_CentralEurope_glamoslong.dat', sep=r'\s+', dtype={'ID': str})
    d = d.set_index('ID')
    d.columns = d.columns.astype(int)
    return d


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else str(RUN.parents[3] / 'past_mb_check.csv')
    Ba, Bw = model('Annual_Balance_sfc'), model('Winter_balance_sfc')
    rows = []
    for name, gid in GLACIERS.items():
        o = glamos(name)
        for y0, y1 in PERIODS:
            yrs = [y for y in range(y0, y1 + 1) if y in o.index and y in Ba.columns]
            wy = [y for y in yrs if np.isfinite(o.Bw.get(y, np.nan))]
            r = dict(glacier=name, rgi7=gid, period=f'{y0}-{y1}', n=len(yrs), n_winter=len(wy), flag=FLAG.get(name, ''),
                     Ba_obs=o.Ba.loc[yrs].mean(), Ba_mod=Ba.loc[gid, yrs].mean())
            if wy:
                r.update(Bw_obs=o.Bw.loc[wy].mean(), Bw_mod=Bw.loc[gid, wy].mean(),
                         Bs_obs=(o.Ba - o.Bw).loc[wy].mean(), Bs_mod=(Ba.loc[gid, wy] - Bw.loc[gid, wy]).mean())
            rows.append(r)
    T = pd.DataFrame(rows)
    for c in ('Ba', 'Bw', 'Bs'):
        T[f'd{c}'] = T[f'{c}_mod'] - T[f'{c}_obs']
    T.to_csv(out, index=False, float_format='%.3f')
    pd.set_option('display.width', 250)
    print(T.drop(columns='flag').round(2).to_string(index=False))
    M = T[T.glacier != 'clariden'].groupby('period')[['n', 'Ba_obs', 'Ba_mod', 'dBa', 'n_winter', 'dBw', 'dBs']].mean()
    print('\nmean of the 7 glaciers (Clariden out: outline), model minus observed [m w.e. a-1]:')
    print(M.round(2).to_string())


if __name__ == '__main__':
    main()
