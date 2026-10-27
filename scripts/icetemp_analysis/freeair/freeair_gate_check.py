"""Gate checks of the free-air switch (GloGEM freeair-thermal; week3/launch_freeair_gate.sh, tag fag).

1. air: the logged band air of the free-air run (firnice_tairth_<id>.dat) equals the free-air table at the band
   heights and the glacier's grid cell (nodes.json of the orientation investigation), every month.
2. surface: in the run with Gilbert's cooling, the last month's surface (firnice_temp_validation_<id>.dat, column
   T_sfc, written to 3 decimals) equals min(0, 0.95 * logged air + DT) for every band.
The byte comparisons (switch off against o4; mass-balance files of the free-air runs against the switch-off run) are
compare_run_files.py.

Usage: python freeair_gate_check.py FA_RUN FAC_RUN [--dt -1.4] [--slope 0.95]   (f_sur F: --slope 1 --dt -F)
"""
import argparse
import json
from pathlib import Path

import numpy as np
import xarray as xr

BAYESCAL = Path('/home/jabeer/projects/glogemflow_development/glogemflow_icetemp/data/bayescal')
SUB = Path('training_runs/monthly/CentralEurope/PAST/firnice_temperature')
TABLE = Path('/scratch-fourth/jabeer/tmp_icetemp_trees/week3/freeair/freeair_CentralEurope.nc')
NODES = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3/orientation_investigation/evidence/nodes.json')
BAL_MONTH = 10


def read_band_table(path):
    lines = path.read_text().splitlines()
    bands = np.array([t for t in lines[1].split() if t not in ('year', 'month')], float)
    return bands, np.loadtxt(lines[2:], ndmin=2)


def check_air(run, tab, nodes):
    worst = 0.0
    for f in sorted((BAYESCAL / run / SUB).glob('firnice_tairth_*.dat')):
        gid = f.stem.split('_')[-1]
        nd = nodes[gid]
        bands, d = read_band_table(f)
        T = tab.tfa.sel(latitude=nd['lat'], longitude=nd['lon']).values.astype(float)     # (time, height)
        h = tab.height.values
        cal_year = d[:, 0].astype(int) - (d[:, 1] >= BAL_MONTH)
        tix = (cal_year - int(tab.year.values[0])) * 12 + d[:, 1].astype(int) - 1
        ref = np.array([np.interp(bands, h, T[t]) for t in tix])
        got = d[:, 2:]
        ok = got > -500                                                                     # noval rows (no band)
        diff = np.abs(np.where(ok, got - ref, 0.0)).max()
        worst = max(worst, diff)
        print(f'air  {gid}: cell {nd["lon"]:.2f} E {nd["lat"]:.2f} N, {len(bands)} bands, {ok.sum()} values, '
              f'max |logged - table| {diff:.2e} C')
    return worst


def check_surface(run, dt, slope=0.95):
    worst = 0.0
    for f in sorted((BAYESCAL / run / SUB).glob('firnice_temp_validation_*.dat')):
        gid = f.stem.split('_')[-1]
        lines = f.read_text().splitlines()
        cols = next(ln for ln in lines if ln.startswith('# elev')).lstrip('#').split()
        v = np.array([ln.split() for ln in lines if ln.strip() and not ln.startswith('#')
                      and len(ln.split()) == len(cols)], float)
        js, je = cols.index('T_sfc'), cols.index('elev')
        bands, d = read_band_table(f.parent / f'firnice_tairth_{gid}.dat')
        # rows are hydrological year x calendar month: the run's last month is September of the last label
        last = d[(d[:, 0] == d[:, 0].max()) & (d[:, 1] == BAL_MONTH - 1)][0, 2:]
        surf = v[:, js]
        elev = v[:, je]
        k = np.array([int(np.argmin(np.abs(bands - e))) for e in elev])
        expect = np.minimum(0.0, slope * last[k] + dt)
        diff = np.abs(surf - expect).max()
        worst = max(worst, diff)
        print(f'surf {gid}: {len(surf)} bands, max |surface - min(0, {slope} air {dt:+.1f})| {diff:.2e} C')
    return worst


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('fa_run')
    ap.add_argument('fac_run')
    ap.add_argument('--dt', type=float, default=-1.4)
    ap.add_argument('--slope', type=float, default=0.95, help='1.0 for f_sur (surface = air - f_sur)')
    a = ap.parse_args()
    tab = xr.open_dataset(TABLE)
    nodes = json.load(open(NODES))['nodes']
    wa = check_air(a.fa_run, tab, nodes)
    ws = check_surface(a.fac_run, a.dt, a.slope)
    print(f'air: worst {wa:.2e} C ({"PASS" if wa < 0.01 else "FAIL"} at 0.01); '
          f'surface: worst {ws:.2e} C ({"PASS" if ws < 0.01 else "FAIL"} at 0.01)')


if __name__ == '__main__':
    main()
