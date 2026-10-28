"""Check of GloGEM's air_temperature = 'freeair' switch (free-air mass-balance plan, Step 1).

A past run with write_hypsometry_files = 'y' writes, per glacier, temp_<id>.dat: each band's initial surface
height and the decade mean of the mass balance's band air (tgs, hydrological years ending in the decade year).
With the switch on, that air must equal the free-air table at the band height plus the glacier's T_off.
The table is interpolated as freeair_at.pro does (linear in height, clamped at the table's ends); the climate
cell is read from the run log ('[free air] cell ...'). The first snapshot (one year, divided by 120, with the
months before ERA5) is skipped.

Run: python -I mb_freeair_check.py <run log> <hypsometry folder> <calibration .dat> <free-air table .nc>
"""
import re
import sys
from pathlib import Path

import numpy as np
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


def main(log, hypso_dir, cal_file, table):
    cells = re.findall(r'\[free air\] cell\s+([-\d.]+) E\s+([-\d.]+) N', Path(log).read_text())
    assert cells, 'no free-air cell in the log'
    assert len(set(cells)) == 1, f'more than one cell in the log: {set(cells)} (the check assumes one)'
    lon, lat = map(float, cells[0])
    cal = {l.split()[0]: float(l.split()[11]) for l in Path(cal_file).read_text().splitlines()[1:] if l.strip()}
    tab = xr.open_dataset(table)
    tyear, tmon = tab.year.values.astype(int), tab.month.values.astype(int)
    worst = 0.0
    for f in sorted(Path(hypso_dir).glob('temp_*.dat')):
        gid = f.stem.split('_')[1]
        rows = np.loadtxt(f, skiprows=1, ndmin=2)
        heights, snaps = rows[:, 0], rows[:, 1:]
        for j in range(1, snaps.shape[1]):
            y_end = 1940 + 10 * j                            # snapshots at hydrological years 1940, 1950, ...
            if y_end > 2025:
                break
            hyd = tyear + (tmon >= 10)
            months = np.where((hyd > y_end - 10) & (hyd <= y_end))[0]
            assert months.size == 120
            exp = band_air_expected(tab, lon, lat, heights, months) + cal[gid]
            ok = snaps[:, j] > -90                           # bands without ice carry the no-value
            d = np.abs(snaps[ok, j] - exp[ok]).max() if ok.any() else 0.0
            worst = max(worst, d)
            print(f'{gid} {y_end}: {ok.sum():3d} bands, max |model - (table + T_off)| = {d:.2e} K')
    print(f'worst {worst:.2e} K: ' + ('PASS' if worst < TOL else 'FAIL'))
    return worst < TOL


if __name__ == '__main__':
    sys.exit(0 if main(*sys.argv[1:5]) else 1)
