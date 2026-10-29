"""Check of GloGEM's air_temperature = 'freeair' switch (free-air mass-balance plan, Step 1).

A past run with write_hypsometry_files = 'y' writes, per glacier, temp_<id>.dat: each band's initial surface
height and the decade mean of the mass balance's band air (tgs, hydrological years ending in the decade year).
With the switch on, that air must equal the free-air table at the band height plus the glacier's T_off.
The table is interpolated as freeair_at.pro does (linear in height, clamped at the table's ends). Each glacier's
climate cell and ERA5 ground (hclim) come from Step 0's band_change.csv (GloGEM's cell rule, computed separately);
the cells must be those of the run log ('[free air] cell ...'). The first snapshot (one year, divided by 120, with
the months before ERA5) is skipped.

With a blend depth (Step 2), bands below the ground get (1 - w) * free air + w * today's air, w = (hclim - z) / depth
in [0, 1]; today's air is read from the same snapshots of a 2 m run on the same calibration (its air is linear in
the months, so the decade means blend the same way).

Run: python -I mb_freeair_check.py <run log> <hypsometry folder> <calibration .dat> <free-air table .nc>
                                   <band_change.csv> [<hypsometry folder of the 2 m run> <blend depth, m>]
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

TOL = 2e-5   # written with f12.5; the table is float32


def band_air_expected(tab, lon, lat, heights, months):
    """Mean over `months` (indices into the table time) of the free air at `heights`."""
    node = tab.tfa.sel(longitude=lon, latitude=lat, method='nearest')
    assert abs(float(node.longitude) - lon) < 1e-3 and abs(float(node.latitude) - lat) < 1e-3
    fa = node.values.astype(np.float64)                      # (time, height)
    h0, dh = float(tab.height[0]), float(tab.height[1] - tab.height[0])
    nh = tab.height.size
    x = np.clip((heights - h0) / dh, 0.0, nh - 1.0)
    k = np.minimum(x.astype(int), nh - 2)
    w = x - k
    sub = fa[months]                                         # (n, height)
    return (sub[:, k] * (1 - w) + sub[:, k + 1] * w).mean(axis=0)


def main(log, hypso_dir, cal_file, table, cells_csv, ref_dir=None, depth=None):
    logged = {(float(a), float(b)) for a, b in
              re.findall(r'\[free air\] cell\s+([-\d.]+) E\s+([-\d.]+) N', Path(log).read_text())}
    assert logged, 'no free-air cell in the log'
    cells = pd.read_csv(cells_csv, dtype={'id': str}).groupby('id')[['node_lon', 'node_lat', 'hclim']].first()
    cal = {l.split()[0]: float(l.split()[11]) for l in Path(cal_file).read_text().splitlines()[1:] if l.strip()}
    tab = xr.open_dataset(table)
    tyear, tmon = tab.year.values.astype(int), tab.month.values.astype(int)
    worst, used = 0.0, set()
    for f in sorted(Path(hypso_dir).glob('temp_*.dat')):
        gid = f.stem.split('_')[1]
        lon, lat, hclim = cells.loc[gid]
        used.add((round(lon, 2), round(lat, 2)))
        rows = np.loadtxt(f, skiprows=1, ndmin=2)
        heights, snaps = rows[:, 0], rows[:, 1:]
        if depth is not None:
            ref = np.loadtxt(Path(ref_dir) / f.name, skiprows=1, ndmin=2)
            assert np.array_equal(ref[:, 0], heights)
            w = np.clip((hclim - heights) / float(depth), 0.0, 1.0)
        for j in range(1, snaps.shape[1]):
            y_end = 1940 + 10 * j                            # snapshots at hydrological years 1940, 1950, ...
            if y_end > 2025:
                break
            hyd = tyear + (tmon >= 10)
            months = np.where((hyd > y_end - 10) & (hyd <= y_end))[0]
            assert months.size == 120
            exp = band_air_expected(tab, lon, lat, heights, months) + cal[gid]
            ok = snaps[:, j] > -90                           # bands without ice carry the no-value
            if depth is not None:
                exp = (1 - w) * exp + w * ref[:, 1 + j]
                ok &= ref[:, 1 + j] > -90
            d = np.abs(snaps[ok, j] - exp[ok]).max() if ok.any() else 0.0
            worst = max(worst, d)
            below = int((ok & (heights < hclim)).sum())
            print(f'{gid} {y_end}: {ok.sum():3d} bands ({below} below the ground), max |model - expected| = {d:.2e} K')
    assert used == logged, f'cells: Step 0 {sorted(used)}, run log {sorted(logged)}'
    print(f'worst {worst:.2e} K: ' + ('PASS' if worst < TOL else 'FAIL'))
    return worst < TOL


if __name__ == '__main__':
    sys.exit(0 if main(*sys.argv[1:8]) else 1)
