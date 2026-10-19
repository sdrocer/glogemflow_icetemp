"""Observed band balances for the perfect-firn test (plan of 2026-10-19, step 2): firnobs_<RGI7 id>.dat for GloGEM's
firnice_firn_obs_dir (procedures/read/read_firn_obs.pro).

From the band tables of band_compare.py (GLAMOS distributed balance averaged over the RGI7 10 m bands): every band and
hydrological year with at least 3 GLAMOS cells gets its observed balance. Year 0 holds the ledger's start balance, the
mean of the band's first ten observed years (as check 5, ledger_vs_firnmodel.py), only for bands observed by the run's
first year (RUN0); later-observed glaciers (Fiescher, Gorner, Morteratsch) keep GloGEM's start. Band-years without an
observation keep GloGEM's own balance in the run.

Usage: python write_firn_obs.py [OUT_DIR]   (FIRNMODEL_DIR: band_tables_v2/, as band_compare.py)
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from band_compare import GL, MINCELLS  # noqa: E402

FM = os.environ.get('FIRNMODEL_DIR', '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week2/firnmodel')
RUN0 = 1972   # first hydrological year of the thermal runs (window Oct 1971 - Sep 2025)


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3/firnobs'
    os.makedirs(out, exist_ok=True)
    for g, gid in GL.items():
        t = pd.read_csv(f'{FM}/band_tables_v2/annual_bands_{g}_{gid}.csv')
        t = t[(t.glamos_cells >= MINCELLS) & np.isfinite(t.mb_mwe)].sort_values(['band_mid', 'year'])
        first = t.groupby('band_mid').year.min()
        start = t.groupby('band_mid').mb_mwe.apply(lambda s: s.iloc[:10].mean())[first <= RUN0]
        with open(f'{out}/firnobs_{gid}.dat', 'w') as f:
            f.write(f'# {g} ({gid}): year (0 = ledger start balance), band elevation [m], GLAMOS band balance [m w.e.]\n')
            for z, b in start.items():
                f.write(f'{0:5d} {z:8.1f} {b:9.4f}\n')
            for r in t.itertuples():
                f.write(f'{int(r.year):5d} {r.band_mid:8.1f} {r.mb_mwe:9.4f}\n')
        print(f'{gid} {g}: {t.band_mid.nunique()} bands ({len(start)} with a start balance), {len(t)} band-years, '
              f'{t.year.min()}-{t.year.max()}')


if __name__ == '__main__':
    main()
