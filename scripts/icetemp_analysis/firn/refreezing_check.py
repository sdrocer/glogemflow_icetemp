"""B of the missing-heat plan: how much meltwater refreezes on ice and on firn bands, with the mass balance's default
(parametrised) refreezing and with GloGEM's full refreezing model (repaired, GloGEMflow 01dbb65).

Two monthly forward runs 1999-2025 of the 23 calibration glaciers, each on its own calibration: C0 = fwdB_none_fixed
(parametrised; refit_none_fixed) and B = fwdRF2 (full; refit_rffull2). Per glacier and hydrological year 2000-2025, the
area-weighted refreezing [m w.e. a-1] on ice bands and on firn bands, a band being firn where the mean balance of the
5 years up to then is > 0 (the mass balance's own mask rule, from its belev file).

Usage: python refreezing_check.py [OUT_CSV] [C0_RUN] [B_RUN]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

S = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week2/snowcal')
SUB = 'monthly/CentralEurope/PAST/PAST_original'
BANDS = Path('/itet-stor/jabeer/glogem/geometricdata/rgiv7/bands/centraleurope')
Y0, Y1 = 2000, 2025


def table(run, kind, gid):
    f = run / SUB / ('refr_elevation' if kind == 'refrelev' else 'mb_elevation') / f'{kind}_{gid}.dat'
    d = pd.read_csv(f, sep=r'\s+').set_index('Elev')
    d.columns = d.columns.astype(int)
    return d.where(d > -98)


def area(gid):
    b = np.atleast_2d(np.loadtxt(BANDS / f'{gid}.dat', skiprows=5))
    return pd.Series(b[:, 3], index=((b[:, 1] + b[:, 2]) / 2).round().astype(int))


def summary(run, gid):
    r, bl, a = table(run, 'refrelev', gid), table(run, 'belev', gid), area(gid)
    out = []
    for y in range(Y0, Y1 + 1):
        if y not in r.columns or y not in bl.columns:
            continue
        m5 = bl.loc[:, [c for c in bl.columns if y - 4 <= c <= y]].mean(axis=1)
        w = a.reindex(r.index).fillna(0)
        ok = r[y].notna() & (w > 0)
        for kind, sel in (('ice', ok & (m5 <= 0)), ('firn', ok & (m5 > 0))):
            out.append(dict(gid=gid, year=y, kind=kind, area=w[sel].sum(),
                            refr=(r.loc[sel, y] * w[sel]).sum() / w[sel].sum() if w[sel].sum() > 0 else np.nan))
    return out


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else str(S / 'refreezing_check.csv')
    runs = dict(C0=S / (sys.argv[2] if len(sys.argv) > 2 else 'fwdB_none_fixed'),
                B=S / (sys.argv[3] if len(sys.argv) > 3 else 'fwdRF2'))
    gids = sorted(p.stem.split('_')[1] for p in (runs['C0'] / SUB / 'refr_elevation').glob('refrelev_*.dat'))
    rows = []
    for lab, run in runs.items():
        for g in gids:
            rows += [dict(r, run=lab) for r in summary(run, g)]
    T = pd.DataFrame(rows)
    T.to_csv(out, index=False, float_format='%.4f')
    P = T.groupby(['gid', 'run', 'kind']).refr.mean().unstack(['run', 'kind'])
    pd.set_option('display.width', 200)
    print(f'mean refreezing {Y0}-{Y1} [m w.e. a-1, as written], ice and firn bands')
    print(P.round(3).to_string())
    print('all glaciers (plain mean):', P.mean().round(3).to_dict())


if __name__ == '__main__':
    main()
