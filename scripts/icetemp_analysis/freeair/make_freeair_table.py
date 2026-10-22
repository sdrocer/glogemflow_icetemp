"""ERA5 free-air temperature at fixed heights for every grid point of the CentralEurope box, 1940-2025.

At each 0.25 deg grid point and month, T is interpolated linearly in geopotential height (z / g) between the two
pressure levels around each target height (500-5000 m, every 100 m). Heights below the 1000 hPa level or above
400 hPa would be NaN; none occur for this grid. Levels below ERA5's ground are ERA5's own extrapolations.

Read by GloGEM (firnice_thermal_air = 'freeair', setup_firnice_freeair.pro), so the file keeps a node's whole
series contiguous: dims (latitude, longitude, time, height); IDL sees them reversed.

Run: python -I make_freeair_table.py [out.nc]
"""
import sys
from pathlib import Path

import numpy as np
import xarray as xr

SRC = Path('/itet-stor/jabeer/glogem/climatedata/reanalysis/monthly/era5_pressure_levels/CentralEurope')
FILES = ['era5_pl_monthly_alps_1940_1983.nc', 'era5_pl_monthly_alps_1984_2025.nc']
OUT = Path('/scratch-fourth/jabeer/tmp_icetemp_trees/week3/freeair/freeair_CentralEurope.nc')
HEIGHTS = np.arange(500.0, 5000.0 + 1.0, 100.0)
G = 9.80665


def load():
    ds = xr.concat([xr.open_dataset(SRC / f) for f in FILES], dim='valid_time')
    ds = ds.sortby('pressure_level', ascending=False)          # 1000 -> 400 hPa: height increases with index
    return ds


def interp_heights(zlev, tlev, heights):
    """zlev, tlev: (nlev, ...) with zlev increasing along axis 0; returns (..., nheight)."""
    zl = np.moveaxis(zlev, 0, -1)[..., None, :]                # (..., 1, nlev)
    tl = np.moveaxis(tlev, 0, -1)[..., None, :]
    k = (zl < heights[:, None]).sum(axis=-1)                   # (..., nheight): first level at or above h
    ok = (k > 0) & (k < zlev.shape[0])
    k = np.clip(k, 1, zlev.shape[0] - 1)
    z0 = np.take_along_axis(zl[..., 0, :], k - 1, axis=-1)
    z1 = np.take_along_axis(zl[..., 0, :], k, axis=-1)
    t0 = np.take_along_axis(tl[..., 0, :], k - 1, axis=-1)
    t1 = np.take_along_axis(tl[..., 0, :], k, axis=-1)
    out = t0 + (t1 - t0) * (heights - z0) / (z1 - z0)
    return np.where(ok, out, np.nan)


def main(out=OUT):
    ds = load()
    time = ds.valid_time.values
    nt, nlat, nlon = time.size, ds.latitude.size, ds.longitude.size
    tfa = np.empty((nlat, nlon, nt, HEIGHTS.size), dtype=np.float32)
    zmin = np.empty((nlat, nlon, nt), dtype=np.float32)
    for i0 in range(0, nt, 60):
        sl = slice(i0, min(i0 + 60, nt))
        z = ds.z.isel(valid_time=sl).values / G                 # (t, lev, lat, lon)
        t = ds.t.isel(valid_time=sl).values - 273.15
        assert np.all(np.diff(z, axis=1) > 0), 'levels not monotonic in height'
        v = interp_heights(np.moveaxis(z, 1, 0), np.moveaxis(t, 1, 0), HEIGHTS)   # (t, lat, lon, h)
        tfa[:, :, sl, :] = np.moveaxis(v, 0, 2)
        zmin[:, :, sl] = np.moveaxis(z[:, 0], 0, 2)
    if np.isnan(tfa).any():
        raise SystemExit(f'{int(np.isnan(tfa).sum())} heights outside the levels')
    yr = time.astype('datetime64[Y]').astype(int) + 1970
    mo = time.astype('datetime64[M]').astype(int) % 12 + 1
    o = xr.Dataset(
        {'tfa': (('latitude', 'longitude', 'time', 'height'), tfa,
                 {'units': 'degC', 'long_name': 'ERA5 free-air temperature at height (monthly mean)'}),
         'z1000': (('latitude', 'longitude', 'time'), zmin,
                   {'units': 'm', 'long_name': 'geopotential height of the 1000 hPa level'}),
         'year': (('time',), yr.astype(np.int32)), 'month': (('time',), mo.astype(np.int32))},
        coords={'latitude': ds.latitude.values, 'longitude': ds.longitude.values,
                'time': np.arange(nt, dtype=np.int32), 'height': HEIGHTS},
        attrs={'source': 'ERA5 monthly averaged reanalysis on pressure levels (CDS), ' + ', '.join(FILES),
               'method': 'linear in geopotential height between pressure levels, per grid point and month',
               'levels_hPa': ' '.join(str(int(p)) for p in ds.pressure_level.values)})
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    o.to_netcdf(out, format='NETCDF4_CLASSIC')
    print(out, tfa.shape, 'years', yr[0], yr[-1])


if __name__ == '__main__':
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else OUT)
