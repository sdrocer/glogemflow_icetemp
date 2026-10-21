"""Orientation factors for every band and part of a glacier: step 1c of the Step 5 plan (approved 21 Oct 2026).

Each glacierised 10 m band of a glacier is split into three parts:
  G  gentle, slope below 10 deg;
  N  slope 10 deg or more, aspect 270-90 deg;
  S  slope 10 deg or more, aspect 90-270 deg.
F of a part = June-September energy-balance melt of its pixels (area mean) over that of flat open ground at the band,
the same energy balance as the borehole tables (orientation_eb.py, write_orient_tables.py): psi0 per band such that
the flat-ground melt equals the band's own degree-day thermal melt (June-September 1981-2020, at least 30 mm), kept
when the hourly air temperature is shifted by dT = -4..+6 K; central strength (c 10 W m-2 K-1, Jungfraujoch weather).

Pixels: EU-DEM 25 m inside GloGEM's glacier mask (a pixel belongs to the glacier when its centre falls in a mask
cell), slope and aspect from the EU-DEM (aspect to true north), terrain horizons from EU-DEM (72 azimuths, as week
1), each pixel joined to the nearest glacierised model band by its EU-DEM elevation. Clear-sky radiation hourly for
June-September, vectorised over pixels, sun positions at the glacier centre and the air mass of the band (the
physics of orientation_eb.clear_sky_hourly). An empty part takes the band's mean, a band without pixels the nearest
band's factors (both flagged). The bands, their monthly air temperature and degree-day melt come from a run with
firnice_write_melth = 'y' (firnice_tairth_<id>.dat, firnice_melth_<id>.dat; glacierised bands from
firnice_snowclim_<id>.dat).

Writes OUT/C_G, OUT/C_N, OUT/C_S/orient_<id>.dat (one file per glacier of the run), OUT/parts.csv (per glacier,
band and part: pixels, area share, mean slope and aspect, F at dT 0, flags) and OUT/parts.json (inputs, counts).

Usage: PYTHONPATH=src python scripts/icetemp_analysis/orientation/part_factors.py MELTH_RUN_DIR [OUT] [--ids ...]
           [--workers N]
"""
import argparse
import hashlib
import json
import sys
from multiprocessing import get_context
from pathlib import Path

import numpy as np
import pandas as pd
from pyproj import Transformer
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import exposure_check as ec  # noqa: E402
import orientation_eb as oe  # noqa: E402
import write_orient_tables as wot  # noqa: E402

W3 = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3')
MASKS = Path('/itet-stor/jabeer/glogem/geometricdata/rgiv7/glmask/centraleurope')
ZONES = W3 / 'orientation_investigation/data/region11_grid_utm_zone.csv'
FIRNICE = Path('training_runs/monthly/CentralEurope/PAST/firnice_temperature')
GENTLE = np.radians(10.0)
PARTS = ('G', 'N', 'S')
JJAS_DAYS = slice(151, 273)          # 1 June - 30 September (0-based day of a 365-day year)
DTS = wot.DTS
C, ALB, TTIP, DMIN = 10.0, wot.ALB, wot.TTIP, wot.DMIN
CHUNK = 40


def read_grid(f):
    with open(f) as fh:
        txt = fh.read().split('\n', 6)
    hd = {line.split()[0].lower(): float(line.split()[1]) for line in txt[:6]}
    return hd, np.array(txt[6].split(), float).reshape(int(hd['nrows']), int(hd['ncols']))


def read_band_table(path):
    """firnice_melth/tairth file -> (band elevations, DataFrame year, month, then one column per band)."""
    lines = path.read_text().splitlines()
    head = [t for t in lines[1].split() if t not in ('year', 'month')]   # older files repeat the two names
    bands = np.array(head, float)
    d = np.loadtxt(lines[2:])
    return bands, d


def glacierised_bands(run, gid):
    d = np.loadtxt(run / FIRNICE / f'firnice_snowclim_{gid}.dat', comments='#', ndmin=2)
    return np.unique(d[:, 0])


def band_inputs(run, gid, bands_g):
    """Per glacierised band: monthly air temperature [NY, 12] (calendar 1981-2020) and June-September degree-day
    melt (mm, 1981-2020 mean) from the melth run."""
    eb, T = read_band_table(run / FIRNICE / f'firnice_tairth_{gid}.dat')
    em, M = read_band_table(run / FIRNICE / f'firnice_melth_{gid}.dat')
    if not np.array_equal(eb, em):
        raise SystemExit(f'{gid}: tairth and melth bands differ')
    lab, mon = T[:, 0].astype(int), T[:, 1].astype(int)
    out = {}
    for b in bands_g:
        k = int(np.argmin(np.abs(eb - b)))
        if abs(eb[k] - b) > 1:
            raise SystemExit(f'{gid}: band {b} not in the melth file')
        Tm = np.empty((oe.NY, 12))
        dd = 0.0
        for iy, Y in enumerate(range(oe.Y0, oe.Y1 + 1)):
            for m in range(1, 13):
                row = np.where((lab == (Y + 1 if m >= 10 else Y)) & (mon == m))[0]
                if len(row) != 1:
                    raise SystemExit(f'{gid}: no row for {Y}-{m} in the tairth file')
                Tm[iy, m - 1] = T[row[0], 2 + k]
                if 6 <= m <= 9:
                    dd += M[row[0], 2 + k]
        if not np.isfinite(Tm).all() or (Tm < -500).any():
            raise SystemExit(f'{gid}: band {b} has no thermal air temperature in 1981-2020')
        out[b] = dict(Tm=Tm, dd_jjas=1000.0 * dd / oe.NY)
    return out


def glacier_pixels(gid, eu, zone):
    """EU-DEM 25 m pixels whose centre falls in a mask cell: lon, lat, z, slope, aspect (rad, from true north)."""
    hd, m = read_grid(MASKS / f'gl_{gid}.grid')
    cs_m, nr, nc = hd['cellsize'], int(hd['nrows']), int(hd['ncols'])
    utm = f'EPSG:326{int(zone)}'
    to_ll = Transformer.from_crs(utm, 'EPSG:4326', always_xy=True)
    ll_to_utm = Transformer.from_crs('EPSG:4326', utm, always_xy=True)
    ll_to_3035 = Transformer.from_crs('EPSG:4326', 'EPSG:3035', always_xy=True)
    from_3035 = Transformer.from_crs('EPSG:3035', 'EPSG:4326', always_xy=True)
    x0, y0 = hd['xllcorner'], hd['yllcorner']
    cx = np.array([x0, x0 + nc * cs_m, x0, x0 + nc * cs_m])
    cy = np.array([y0, y0, y0 + nr * cs_m, y0 + nr * cs_m])
    X3, Y3 = ll_to_3035.transform(*to_ll.transform(cx, cy))
    a, ex0, ey0, cs = eu.block(min(X3) - 100, max(X3) + 100, min(Y3) - 100, max(Y3) + 100)
    dzr, dzc = np.gradient(a, cs)
    slope = np.arctan(np.hypot(dzc, -dzr))
    asp_g = np.arctan2(-dzc, dzr)
    rr, cc = np.mgrid[0:a.shape[0], 0:a.shape[1]]
    X, Y = ex0 + cc * cs, ey0 - rr * cs
    lon, lat = from_3035.transform(X.ravel(), Y.ravel())
    ux, uy = ll_to_utm.transform(lon, lat)
    ci = np.floor((ux - x0) / cs_m).astype(int)
    ri = np.floor((y0 + nr * cs_m - uy) / cs_m).astype(int)
    inside = (ci >= 0) & (ci < nc) & (ri >= 0) & (ri < nr)
    sel = np.zeros(inside.shape, bool)
    sel[inside] = m[ri[inside], ci[inside]] == 1
    sel &= np.isfinite(a.ravel()) & np.isfinite(slope.ravel())
    lonc, latc = np.mean(np.asarray(lon)[sel]), np.mean(np.asarray(lat)[sel])
    gam = ec.north_bearing('EPSG:3035', lonc, latc)
    return dict(lon=np.asarray(lon)[sel], lat=np.asarray(lat)[sel], z=a.ravel()[sel], slope=slope.ravel()[sel],
                aspect=np.mod(asp_g.ravel()[sel] - gam, 2 * np.pi), latc=latc, lonc=lonc, mask_cell=cs_m)


class Sun:
    """Sun positions at one latitude for June-September, 10-min steps, and per band the clear-sky beam."""

    def __init__(self, lat):
        dec, e0 = ec.spencer(ec.DAYS)
        h, az = ec.sun_geometry(lat, dec)
        self.h, self.az = h[JJAS_DAYS], az[JJAS_DAYS]
        self.e0 = e0[JJAS_DAYS]
        self.up = self.h > 0
        azd = np.degrees(self.az)
        self.b0 = np.floor(azd / 5.0).astype(int) % 72
        self.b1 = (self.b0 + 1) % 72
        self.wb = (azd - 5.0 * np.floor(azd / 5.0)) / 5.0
        self.sh, self.ch = np.sin(self.h), np.cos(self.h)

    def beam(self, z, tau=ec.TAU):
        m = np.where(self.up, ec.air_mass(np.where(self.up, self.h, 0.1), z), 0)
        tm = np.where(self.up, tau ** m, 0)
        sp = ec.S0 * self.e0[:, None]
        return sp * tm, 0.3 * (1 - tm) * sp * self.sh * self.up, sp * tm * self.sh * self.up

    def direct(self, beam, slope, aspect, hor):
        """Hour means [P, days, 24] of clear-sky direct on the planes, with their terrain horizons."""
        cosi = (np.sin(self.h)[None] * np.cos(slope)[:, None, None]
                + np.cos(self.h)[None] * np.sin(slope)[:, None, None] * np.cos(self.az[None] - aspect[:, None, None]))
        hh = (1 - self.wb)[None] * hor[:, self.b0] + self.wb[None] * hor[:, self.b1]
        lit = self.up[None] & (cosi > 0) & (self.h[None] > hh)
        s = beam[None] * np.maximum(cosi, 0) * lit
        return s.reshape(s.shape[0], s.shape[1], 24, 6).mean(3)


def hour_means(a):
    return a.reshape(a.shape[0], 24, 6).mean(2)


def eb_jjas(Th, S, p0, c=C):
    """June-September energy-balance melt, mm w.e. per year; Th, S [..., NY, days, 24]."""
    psi = np.where(Th >= TTIP, p0 + c * Th, p0)
    return np.maximum((1 - ALB) * S + psi, 0).sum(axis=(-3, -2, -1)) * 3600 / oe.LF / oe.NY


def one_glacier(args):
    gid, run, zone, w, sig, dif = args
    eu = ec.EuDem()
    bands_g = glacierised_bands(run, gid)
    inp = band_inputs(run, gid, bands_g)
    px = glacier_pixels(gid, eu, zone)
    band_of = bands_g[np.argmin(np.abs(px['z'][:, None] - bands_g[None, :]), axis=1)]
    part = np.where(px['slope'] < GENTLE, 'G', np.where((px['aspect'] >= np.pi / 2) & (px['aspect'] < 1.5 * np.pi),
                                                        'S', 'N'))
    hor = np.empty((len(px['z']), 72))
    for i in range(0, len(hor), 800):
        hor[i:i + 800], _ = eu.horizons(px['lon'][i:i + 800], px['lat'][i:i + 800])
    svf = ec.sky_view(px['slope'], px['aspect'], hor)
    sun = Sun(px['latc'])
    sig_j = sig.reshape(oe.NY, 365)[:, JJAS_DAYS, None]
    dif_j = dif.reshape(oe.NY, 365)[:, JJAS_DAYS, None]
    rows, tables = [], {p: [] for p in PARTS}
    for b in bands_g:
        Th = oe.Forcing(inp[b]['Tm'], w, sig, dif).Th[:, JJAS_DAYS, :]
        beam, fdif, fdir = sun.beam(b)
        S_flat = sig_j * hour_means(fdir)[None] + dif_j * hour_means(fdif)[None]
        target = max(inp[b]['dd_jjas'], DMIN)
        p0 = brentq(lambda x: eb_jjas(Th, S_flat, x) - target, -900, 400, xtol=0.02)
        flat = np.array([eb_jjas(Th + dT, S_flat, p0) for dT in DTS])
        sel = np.where(band_of == b)[0]
        melt = np.zeros((len(sel), len(DTS)))
        for i in range(0, len(sel), CHUNK):
            j = sel[i:i + CHUNK]
            sdir = sun.direct(beam, px['slope'][j], px['aspect'][j], hor[j])
            S = sig_j[None] * sdir[:, None] + dif_j[None] * svf[j, None, None, None] * hour_means(fdif)[None, None]
            for k, dT in enumerate(DTS):
                melt[i:i + CHUNK, k] = eb_jjas(Th[None] + dT, S, p0)
        F_band = melt.mean(axis=0) / flat if len(sel) else None
        for p in PARTS:
            q = part[sel] == p
            if q.any():
                F, flag = melt[q].mean(axis=0) / flat, ''
            elif len(sel):
                F, flag = F_band, 'empty part: band mean'
            else:
                F, flag = None, 'no pixels: nearest band'
            rows.append(dict(rgi=gid, band=b, part=p, n_pix=int(q.sum()), share=float(q.mean()) if len(sel) else 0.0,
                             slope=float(np.degrees(px['slope'][sel][q]).mean()) if q.any() else np.nan,
                             aspect=float(np.degrees(px['aspect'][sel][q]).mean()) if q.any() else np.nan,
                             psi0=p0, dd_jjas=inp[b]['dd_jjas'], flag=flag,
                             **{f'F{dT:+.0f}': (np.nan if F is None else float(F[k])) for k, dT in enumerate(DTS)}))
    R = pd.DataFrame(rows)
    fcols = [f'F{dT:+.0f}' for dT in DTS]
    have = R.dropna(subset=[fcols[0]])
    for i in R.index[R[fcols[0]].isna()]:
        near = have[(have.part == R.at[i, 'part'])]
        k = (near.band - R.at[i, 'band']).abs().idxmin()
        R.loc[i, fcols] = near.loc[k, fcols].to_numpy()
    print(f'{gid}: {len(px["z"])} pixels, {len(bands_g)} bands, parts G/N/S '
          f'{(part == "G").mean():.2f}/{(part == "N").mean():.2f}/{(part == "S").mean():.2f}', flush=True)
    return R


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('run_dir')
    ap.add_argument('out', nargs='?', default=str(W3 / 'orient_parts'))
    ap.add_argument('--ids', default=None)
    ap.add_argument('--workers', type=int, default=12)
    a = ap.parse_args()
    run, out = Path(a.run_dir), Path(a.out)
    ids = a.ids.split(',') if a.ids else sorted(set((W3 / 'orient_glaciers.txt').read_text().split()))
    zones = pd.read_csv(ZONES, dtype={'glacier_id': str}).set_index('glacier_id').utm_zone.to_dict()
    w = oe.read_station(wot.INPUTS['jungfraujoch'])
    sig, dif, _ = oe.sky_weights(w)
    jobs = [(g, run, zones.get(g, 32), w, sig, dif) for g in ids]
    with get_context('fork').Pool(a.workers) as pool:
        R = pd.concat(pool.map(one_glacier, jobs), ignore_index=True)
    fcols = [f'F{dT:+.0f}' for dT in DTS]
    for p in PARTS:
        d = out / f'C_{p}'
        d.mkdir(parents=True, exist_ok=True)
        for g in ids:
            B = [(int(round(r['band'])), [r[c] for c in fcols]) for _, r in R[(R.rgi == g) & (R.part == p)].iterrows()]
            wot.write(d / f'orient_{g}.dat', wot.table_text(g, B, f'C, part {p}'))
    R.to_csv(out / 'parts.csv', index=False, float_format='%.5f')
    meta = dict(run=str(run), glaciers=ids, gentle_deg=10.0, strength='C', dT=DTS, n_rows=len(R),
                flags=R.flag.value_counts().to_dict(),
                inputs={k: dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest())
                        for k, p in dict(jungfraujoch=wot.INPUTS['jungfraujoch'], zones=ZONES).items()})
    wot.write(out / 'parts.json', json.dumps(meta, indent=1) + '\n')
    print(json.dumps({k: meta[k] for k in ('n_rows', 'flags')}, indent=1))


if __name__ == '__main__':
    main()
