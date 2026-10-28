"""Where and when does the mass balance's air change with air_temperature = 'freeair'? (free-air mass-balance plan,
Step 0; no model runs.)

For every CentralEurope glacier band (RGI7, GloGEM's inventory and band files): today's band air without T_off,
T_2m(cell) + (z - hclim) * dtdz[month], and the free air at the band height, per calendar month, averaged over
the calibration window (hydrological years 2000-2019, Oct 1999 - Sep 2019). T_off cancels in the difference.

The climate cell follows GloGEM: 0.1 deg cells with edges at multiples of 0.1 deg minus 0.05 deg, in single
precision (setup_grid_loops.pro, glogem.pro), the cell's climate point
from the mean position of every glacier of the region in it (read_climate_series.pro), the nearest .mdi grid point
in latitude and longitude separately (downscale_gcmdata_monthly.pro). Glaciers with volume_ini > 0 and area inside
size_range, as GloGEM selects them.

Writes <out>/band_change.csv (one row per band) and prints summaries. Run:
python -I freeair_band_change.py [out_dir]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

GEOM = Path('/itet-stor/jabeer/glogem/geometricdata/rgiv7')
CLIM = Path('/itet-stor/jabeer/glogem/climatedata/reanalysis/monthly/era5/CentralEurope')
TABLE = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3/freeair/freeair_CentralEurope.nc')
OUT = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/mbfreeair/step0')
GRID_STEP = 0.1
SIZE_RANGE = (0.002, 100000.0)
HYD_YEARS = (2000, 2019)
BAL_MONTH = 10


def tokens(path):
    """All numbers after the header line, as IDL's free-format readf sees them."""
    with open(path) as f:
        f.readline()
        return np.array(f.read().split(), dtype=np.float64)


def read_clim():
    t = tokens(CLIM / 'clim_CentralEurope.mdi')
    nt, nlon, nlat = (int(v) for v in t[:3])
    i = 5                                                    # nvar: two numbers
    rtime = t[i:i + nt]; i += nt
    rlon = t[i:i + nlon]; i += nlon
    rlat = t[i:i + nlat]; i += nlat
    relev = t[i:i + nlon * nlat].reshape(nlon, nlat); i += nlon * nlat
    rtemp = t[i:i + nt * nlon * nlat].reshape(nt, nlon, nlat)
    ryear = np.trunc(rtime)
    rmon = np.floor(((rtime + 0.04) - ryear) * 12 + 0.6 - (ryear - ryear[0]) / 2800. * 12 + 0.5)
    g = tokens(CLIM / 'tgrad_CentralEurope.mdi')
    nm, gl, ga = (int(v) for v in g[:3])
    j = 4 + nm
    assert gl == nlon and ga == nlat and np.allclose(g[j:j + gl], rlon) and np.allclose(g[j + gl:j + gl + ga], rlat)
    j += gl + ga
    rtg = g[j:j + nm * gl * ga].reshape(nm, gl, ga)
    rlon = np.where(rlon >= 180, rlon - 360, rlon)           # read_climatepast_monthly.pro
    return dict(rlon=rlon, rlat=rlat, relev=relev, rtemp=rtemp, ryear=ryear.astype(int), rmon=rmon.astype(int),
                rtg=rtg)


def read_inventory():
    rows = []
    with open(GEOM / 'files' / 'thick_centraleurope.dat') as f:
        f.readline()
        for line in f:
            a = line.replace(';', ' ').split()
            rows.append((a[0], float(a[1]), float(a[2]), float(a[3]), float(a[4])))
    inv = pd.DataFrame(rows, columns=['id', 'lon', 'lat', 'area', 'volume'])
    inv['run'] = (inv.volume > 0) & (inv.area > SIZE_RANGE[0]) & (inv.area < SIZE_RANGE[1])
    return inv


def climate_points(inv, clim):
    """GloGEM's climate grid point of every glacier (indices into rlon, rlat)."""
    gs = np.float32(GRID_STEP)                               # grid_step = 0.1 is single precision in IDL
    edge = lambda x0, i: float(x0 + np.float32(i) * gs)      # lon0[0] + gx*grid_step, single precision
    # setup_grid_loops.pro: cell edges at multiples of grid_step minus half a step
    lo_edge = lambda v: np.float32(np.floor(v / float(gs))) * gs - gs / np.float32(2)
    hi_edge = lambda v: np.float32(np.trunc(v / float(gs))) * gs + gs / np.float32(2) + np.float32(2) * gs
    lon0 = (lo_edge(inv.lon.min()), hi_edge(inv.lon.max()))
    lat0 = (lo_edge(inv.lat.min()), hi_edge(inv.lat.max()))
    ngx = max(int((lon0[1] - lon0[0]) / gs), 1)
    ngy = max(int((lat0[1] - lat0[0]) / gs), 1)
    ilon = np.full(len(inv), -1); ilat = np.full(len(inv), -1)
    run = inv.run.values
    for gx in range(ngx):
        lo = (edge(lon0[0], gx), edge(lon0[0], gx + 1))
        inx = (inv.lon.values >= lo[0]) & (inv.lon.values < lo[1])
        if not (inx & run).any():
            continue
        for gy in range(ngy):
            la = (edge(lat0[0], gy), edge(lat0[0], gy + 1))
            cell = inx & (inv.lat.values >= la[0]) & (inv.lat.values < la[1]) & run
            if not cell.any():
                continue
            gmid = (inv.lat.values[cell].mean(), inv.lon.values[cell].mean())
            ilat[cell] = np.argmin(np.abs(clim['rlat'] - gmid[0]))
            ilon[cell] = np.argmin(np.abs(clim['rlon'] - gmid[1]))
    assert (ilon[run] >= 0).all(), 'a glacier fell outside every cell'
    return ilon, ilat


def read_bands(gid):
    b = np.loadtxt(GEOM / 'bands' / 'centraleurope' / f'{gid}.dat', skiprows=5, ndmin=2)
    return b[:, 1] + 5.0, b[:, 3]                            # process_hypsometry_data.pro: elev = hstart + 5


def main(out=OUT):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    clim = read_clim()
    inv = read_inventory()
    ilon, ilat = climate_points(inv, clim)
    hyd = clim['ryear'] + (clim['rmon'] >= BAL_MONTH)
    win = (hyd >= HYD_YEARS[0]) & (hyd <= HYD_YEARS[1])
    tab = xr.open_dataset(TABLE)
    thyd = tab.year.values.astype(int) + (tab.month.values.astype(int) >= BAL_MONTH)
    twin = (thyd >= HYD_YEARS[0]) & (thyd <= HYD_YEARS[1])
    tmon = tab.month.values.astype(int)
    h0, dh, nh = float(tab.height[0]), float(tab.height[1] - tab.height[0]), tab.height.size
    rows, nodes = [], {}
    for k in np.where(inv.run.values)[0]:
        key = (ilon[k], ilat[k])
        if key not in nodes:
            lon, lat = clim['rlon'][key[0]], clim['rlat'][key[1]]
            t2 = clim['rtemp'][:, key[0], key[1]]
            node = tab.tfa.sel(longitude=lon, latitude=lat, method='nearest')
            assert abs(float(node.longitude) - lon) < 1e-3 and abs(float(node.latitude) - lat) < 1e-3, (lon, lat)
            fa = node.values.astype(np.float64)              # (time, height)
            nodes[key] = dict(
                lon=lon, lat=lat, hclim=clim['relev'][key], dtdz=clim['rtg'][:, key[0], key[1]],
                t2=np.array([t2[win & (clim['rmon'] == m)].mean() for m in range(1, 13)]),
                fa=np.array([fa[twin & (tmon == m)].mean(axis=0) for m in range(1, 13)]))   # (12, height)
        nd = nodes[key]
        z, area = read_bands(inv.id[k])
        x = np.clip((z - h0) / dh, 0.0, nh - 1.0)
        kk = np.minimum(x.astype(int), nh - 2); w = x - kk
        tfa = nd['fa'][:, kk] * (1 - w) + nd['fa'][:, kk + 1] * w                     # (12, nband)
        t2m = nd['t2'][:, None] + (z[None, :] - nd['hclim']) * nd['dtdz'][:, None]
        for b in range(z.size):
            rows.append([inv.id[k], z[b], area[b], nd['lon'], nd['lat'], nd['hclim'],
                         *t2m[:, b], *tfa[:, b]])
    cols = (['id', 'z', 'area', 'node_lon', 'node_lat', 'hclim'] + [f't2m_{m:02d}' for m in range(1, 13)]
            + [f'tfa_{m:02d}' for m in range(1, 13)])
    df = pd.DataFrame(rows, columns=cols)
    df.to_csv(out / 'band_change.csv', index=False, float_format='%.4f')
    d = df[[f'tfa_{m:02d}' for m in range(1, 13)]].values - df[[f't2m_{m:02d}' for m in range(1, 13)]].values
    a = df.area.values
    print(f'{inv.run.sum()} glaciers, {len(df)} bands, {a.sum():.1f} km2, {len(nodes)} climate points')
    below = df.z.values < df.hclim.values
    print(f'below ERA5 ground: {a[below].sum():.1f} km2 ({100 * a[below].sum() / a.sum():.1f} % of the area), '
          f'{df.id[below].nunique()} glaciers')
    print('area-weighted free air minus 2 m air [K], hydrological years 2000-2019')
    print('height        area   ' + ' '.join(f'{m:>5d}' for m in range(1, 13)) + '   year')
    for lo in range(1000, 4800, 250):
        s = (df.z.values >= lo) & (df.z.values < lo + 250)
        if a[s].sum() < 0.5:
            continue
        mm = (d[s] * a[s, None]).sum(axis=0) / a[s].sum()
        print(f'{lo:4d}-{lo + 250:4d} {a[s].sum():7.1f}   ' + ' '.join(f'{v:5.1f}' for v in mm) + f'  {mm.mean():5.2f}')
    mm = (d * a[:, None]).sum(axis=0) / a.sum()
    print('all          ' + f'{a.sum():7.1f}   ' + ' '.join(f'{v:5.1f}' for v in mm) + f'  {mm.mean():5.2f}')


if __name__ == '__main__':
    main(*sys.argv[1:2])
