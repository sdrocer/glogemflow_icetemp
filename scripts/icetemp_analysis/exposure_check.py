"""Exposure check: does sunshine predict wet vs dry 10-20 m firn at high sites beyond band air temperature?

PRE-REGISTRATION (written 2026-10-13 before any radiation value or test statistic was computed;
identical text in <OUT_DIR>/preregistration.md)

Sites: all rows of forcing_check/forcing_check.csv with band_model >= 3300 (76 entities, 16 RGI7
glaciers: 00245, 00773, 00774, 00786, 01112, 01127, 01150, 01205, 01225, 01230, 01565, 02103,
02216, 02596, 02611, 02821). Wet = T_meas > -1.0 C. Air temperature = T_air_mb_notoff.
Massifs (leave-one-out units) and glaciers (permutation blocks = RGI7 glacier_id) as in MASSIFS.
Holes of an entity: borehole.csv rows with the same base glacier name and 10 m band; site sunshine
of the entity = mean over its holes.

Radiation: Spencer (1971) sun, 365 days, 10-min steps of local solar time, S0 = 1361 W m-2.
Direct = S0 E0 tau^m, tau = 0.75, m = Kasten & Young (1989) air mass x p/p0 (standard atmosphere),
times max(cos incidence, 0), only with the sun above the terrain horizon. Diffuse = isotropic,
horizontal 0.3 (1 - tau^m) S0 E0 sin(h) (Campbell & Norman 1998 eq. 11.13) times the sky-view factor
(Dozier & Frew 1990; effective horizon = max(terrain, the plane's own horizon)). No terrain
reflection. Horizon: EU-DEM 25 m, 72 azimuths (5 deg, true north), 100 m - 10 km (25 m steps to
2 km, 50 m beyond), bilinear, earth curvature with refraction 0.13, origin = EU-DEM height at the point.
Site slope/aspect: least-squares plane over DEM cells within 75 m of the hole, first complete DEM of
RGE ALTI 5 m, swissALTI3D 10 m, EU-DEM 25 m; aspect rotated to true north. Reference = flat unshaded
surface at the same latitude and elevation. Ratio = surface / reference, summed over JJA (test
variable) and over the year (reported). Band sunshine: RGI7 glacier pixels with glacier-DEM height
in [band-5, band+5) (widened to +-25 m if < 5 pixels), slope/aspect by np.gradient on the glacier DEM,
horizon from EU-DEM per pixel, all pixels (no subsampling); band ratio = mean of pixel ratios.

Test: A = logit(wet) ~ T_air; B = logit(wet) ~ T_air + JJA ratio; standardised with the training
fold's mean/SD; L2 logistic regression C = 1. Leave-one-massif-out predictions, pooled AUC.
Gain = AUC_B - AUC_A. Permutation p: 2000 block permutations of sunshine between glaciers
(seed 20261013; donor values drawn without replacement if the donor has enough, else all of them plus
draws with replacement), p = (1 + #{gain_perm >= gain}) / 2001. For site and for band sunshine.
PASS (fixed): gain >= 0.05 AND p < 0.05 with BAND sunshine.

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/exposure_check.py
"""
import glob
import json
import logging
import pickle
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window
from pyproj import CRS, Transformer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

REPO = Path(__file__).resolve().parents[2]
FORCING = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/forcing_check/forcing_check.csv')
CACHE = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/cache/datahandler.pkl')
BOREHOLES = REPO / 'glenglat/data/borehole.csv'
OUT_DIR = Path(sys.argv[1] if len(sys.argv) > 1
               else '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week1/exposure')
HM = Path('/itet-stor/jabeer/glazioarch/GeoBaseData/HeightModels')
SWISS = HM / 'swissALTI3D/r2024/swissALTI3D_2024_10m_LV95_LN02.tif'
RGE_DIR = HM / ('France/RGEALTI_2-0_5M_ASC_LAMB93-IGN69_D074_2020-10-15/RGEALTI/'
                '1_DONNEES_LIVRAISON_2021-10-00009/RGEALTI_MNT_5M_ASC_LAMB93_IGN69_D074')
EUDEM = HM / 'eudem/eu-dem_etrs89.vrt'
RGI = Path('/itet-stor/jabeer/glogem/geometricdata/rgiv7')
GL_CRS = 'EPSG:32632'

# Fixed settings
MIN_BAND = 3300
WET_T = -1.0
S0 = 1361.0
TAU = 0.75
DT_MIN = 10
PLANE_R = 75.0
AZ = np.arange(0, 360, 5.0)
RAY_D = np.concatenate([np.arange(100, 2000, 25.0), np.arange(2000, 10000.01, 50.0)])
EARTH_R = 6371000.0
REFRACTION = 0.13
BAND_HALF, BAND_HALF_WIDE, MIN_PIX = 5, 25, 5
C_REG = 1.0
N_PERM = 2000
SEED = 20261013
PASS_GAIN, PASS_P = 0.05, 0.05
JJA = [5, 6, 7]  # month indices 0-based
MASSIFS = {
    'Mont Blanc': ['00773', '00774', '00777', '00786'],
    'Monte Rosa': ['01205', '01225', '01230'],
    'Mischabel-Saas': ['01112', '01127', '01150'],
    'Bernese Alps': ['02596', '02611', '01565', '01265'],
    'Valais west and Diablerets': ['01039', '01309'],
    'Titlis': ['02669'],
    'Ecrins': ['00245'],
    'Bernina': ['02103', '02216'],
    'Ortles': ['02821'],
    'Oetztal': ['03116'],
}
MASSIF_OF = {g: m for m, gs in MASSIFS.items() for g in gs}

log = logging.getLogger('exposure_check')
DAYS = np.arange(1, 366)
MONTH_OF_DAY = pd.to_datetime('2001-01-01') + pd.to_timedelta(DAYS - 1, 'D')
MONTH_OF_DAY = MONTH_OF_DAY.month.to_numpy() - 1


# ---------------------------------------------------------------- sun and radiation

def spencer(days):
    g = 2 * np.pi * (days - 1) / 365.0
    dec = (0.006918 - 0.399912 * np.cos(g) + 0.070257 * np.sin(g) - 0.006758 * np.cos(2 * g)
           + 0.000907 * np.sin(2 * g) - 0.002697 * np.cos(3 * g) + 0.00148 * np.sin(3 * g))
    e0 = (1.000110 + 0.034221 * np.cos(g) + 0.001280 * np.sin(g) + 0.000719 * np.cos(2 * g)
          + 0.000077 * np.sin(2 * g))
    return dec, e0


def sun_geometry(lat_deg, dec, dt_min=DT_MIN):
    """Sun elevation and azimuth (rad, azimuth clockwise from north), [days, steps]."""
    t = (np.arange(0, 1440, dt_min) + dt_min / 2) / 60.0
    w = np.radians(15 * (t - 12))[None, :]
    phi = np.radians(lat_deg)
    d = dec[:, None]
    sinh = np.sin(phi) * np.sin(d) + np.cos(phi) * np.cos(d) * np.cos(w)
    h = np.arcsin(np.clip(sinh, -1, 1))
    az = np.arctan2(-np.cos(d) * np.sin(w), np.sin(d) * np.cos(phi) - np.cos(d) * np.sin(phi) * np.cos(w))
    return h, np.mod(az, 2 * np.pi)


def air_mass(h, z):
    hd = np.degrees(h)
    m = 1.0 / (np.sin(h) + 0.50572 * np.maximum(hd + 6.07995, 1e-6) ** -1.6364)
    return m * (1 - 2.25577e-5 * z) ** 5.25588


def sky_view(slope, aspect, hor):
    """Dozier & Frew (1990). slope, aspect [P] rad; hor [P, 72] terrain horizon elevation, rad."""
    phi = np.radians(AZ)[None, :]
    s, a = slope[:, None], aspect[:, None]
    self_h = np.arctan(-np.tan(s) * np.cos(phi - a))
    e = np.maximum(hor, self_h)
    H = np.pi / 2 - e
    f = np.cos(s) * np.sin(H) ** 2 + np.sin(s) * np.cos(phi - a) * (H - np.sin(H) * np.cos(H))
    return f.mean(axis=1)


def radiation(lat, z, slope, aspect, hor, tau=TAU, s0=S0, sun=None):
    """Monthly totals [P, 12] (J m-2): surface direct, flat direct, flat diffuse; and sky-view
    factors (as pre-registered, and with the terrain horizon floored at 0 deg).
    Surface total = direct + svf * flat diffuse. lat, z, slope, aspect: [P]; hor: [P, 72] rad
    or None (no terrain)."""
    P = len(lat)
    dec, e0 = spencer(DAYS) if sun is None else sun
    sdir, fdir, fdif = np.zeros((P, 12)), np.zeros((P, 12)), np.zeros((P, 12))
    hh_all = hor if hor is not None else np.full((P, len(AZ)), -np.pi / 2)
    svf = sky_view(slope, aspect, hh_all)
    svf0 = sky_view(slope, aspect, np.maximum(hh_all, 0))
    for p in range(P):
        h, az = sun_geometry(lat[p], dec)
        up = h > 0
        m = np.where(up, air_mass(np.where(up, h, 0.1), z[p]), 0)
        tm = np.where(up, tau ** m, 0)
        sp = s0 * e0[:, None]
        beam_n = sp * tm
        diff_h = 0.3 * (1 - tm) * sp * np.sin(h) * up
        cosi = (np.sin(h) * np.cos(slope[p]) + np.cos(h) * np.sin(slope[p]) * np.cos(az - aspect[p]))
        lit = up & (cosi > 0)
        if hor is not None:
            hh = np.interp(np.degrees(az), np.append(AZ, 360.0), np.append(hor[p], hor[p][0]))
            lit &= h > hh
        k = DT_MIN * 60
        np.add.at(sdir[p], MONTH_OF_DAY, (beam_n * np.maximum(cosi, 0) * lit).sum(1) * k)
        np.add.at(fdir[p], MONTH_OF_DAY, (beam_n * np.sin(h) * up).sum(1) * k)
        np.add.at(fdif[p], MONTH_OF_DAY, diff_h.sum(1) * k)
    return sdir, fdir, fdif, svf, svf0


def ratios(sdir, fdir, fdif, svf):
    """JJA ratio, annual ratio, monthly ratios [P, 12]."""
    surf = sdir + svf[:, None] * fdif
    flat = fdir + fdif
    return surf[:, JJA].sum(1) / flat[:, JJA].sum(1), surf.sum(1) / flat.sum(1), surf / flat


# ---------------------------------------------------------------- DEMs

_tf = {}


def tf(dst):
    if dst not in _tf:
        _tf[dst] = Transformer.from_crs('EPSG:4326', dst, always_xy=True)
    return _tf[dst]


def north_bearing(crs, lon, lat):
    """Grid azimuth of true north (rad) at lon/lat."""
    t = tf(crs)
    x1, y1 = t.transform(lon, lat)
    x2, y2 = t.transform(lon, lat + 1e-4)
    return np.arctan2(x2 - x1, y2 - y1)


class EuDem:
    def __init__(self):
        self.src = rasterio.open(EUDEM)
        self.crs = 'EPSG:3035'

    def block(self, xmin, xmax, ymin, ymax):
        T = self.src.transform
        c0, c1 = int(np.floor((xmin - T.c) / T.a)), int(np.ceil((xmax - T.c) / T.a))
        r0, r1 = int(np.floor((ymax - T.f) / T.e)), int(np.ceil((ymin - T.f) / T.e))
        a = self.src.read(1, window=Window(c0, r0, c1 - c0, r1 - r0)).astype(float)
        a[(a < -1000) | (a > 9000)] = np.nan
        x0 = T.c + (c0 + 0.5) * T.a
        y0 = T.f + (r0 + 0.5) * T.e
        return a, x0, y0, T.a

    def horizons(self, lon, lat):
        """Terrain horizon elevation [P, 72] (rad) and EU-DEM height at the points."""
        x, y = tf(self.crs).transform(np.asarray(lon), np.asarray(lat))
        x, y = np.atleast_1d(x), np.atleast_1d(y)
        pad = RAY_D[-1] + 200
        a, x0, y0, cs = self.block(x.min() - pad, x.max() + pad, y.min() - pad, y.max() + pad)
        gam = np.array([north_bearing(self.crs, lo, la) for lo, la in zip(np.atleast_1d(lon), np.atleast_1d(lat))])

        def bil(px, py):
            fc = (px - x0) / cs
            fr = (py - y0) / -cs
            c, r = np.floor(fc).astype(int), np.floor(fr).astype(int)
            dc, dr = fc - c, fr - r
            v = (a[r, c] * (1 - dc) * (1 - dr) + a[r, c + 1] * dc * (1 - dr)
                 + a[r + 1, c] * (1 - dc) * dr + a[r + 1, c + 1] * dc * dr)
            return v

        z0 = bil(x, y)
        th = np.radians(AZ)[None, :, None] + gam[:, None, None]
        d = RAY_D[None, None, :]
        px = x[:, None, None] + d * np.sin(th)
        py = y[:, None, None] + d * np.cos(th)
        zz = bil(px, py) - d ** 2 / (2 * EARTH_R) * (1 - REFRACTION)
        ang = np.arctan((zz - z0[:, None, None]) / d)
        hor = np.nanmax(ang, axis=2)
        return hor, z0

    def disc(self, x, y, r):
        a, x0, y0, cs = self.block(x - r - 50, x + r + 50, y - r - 50, y + r + 50)
        return _disc_from_block(a, x0, y0, cs, x, y, r)


def _disc_from_block(a, x0, y0, cs, x, y, r):
    rr, cc = np.mgrid[0:a.shape[0], 0:a.shape[1]]
    X, Y = x0 + cc * cs, y0 - rr * cs
    m = (X - x) ** 2 + (Y - y) ** 2 <= r ** 2
    if m.sum() < 6:
        return None
    z = a[m]
    if not np.all(np.isfinite(z)):
        return None
    return X[m], Y[m], z, cs


class SwissAlti:
    crs = 'EPSG:2056'

    def __init__(self):
        self.src = rasterio.open(SWISS)

    def disc(self, x, y, r):
        T, src = self.src.transform, self.src
        pad = r + 30
        c0, c1 = int(np.floor((x - pad - T.c) / T.a)), int(np.ceil((x + pad - T.c) / T.a))
        r0, r1 = int(np.floor((y + pad - T.f) / T.e)), int(np.ceil((y - pad - T.f) / T.e))
        if c0 < 0 or r0 < 0 or c1 > src.width or r1 > src.height:
            return None
        a = src.read(1, window=Window(c0, r0, c1 - c0, r1 - r0)).astype(float)
        a[(a < -1000) | (a > 9000)] = np.nan
        return _disc_from_block(a, T.c + (c0 + 0.5) * T.a, T.f + (r0 + 0.5) * T.e, T.a, x, y, r)


class RgeAlti:
    crs = 'EPSG:2154'
    CS = 5.0

    def __init__(self):
        self.tiles = {}
        for f in glob.glob(str(RGE_DIR / 'RGEALTI_FXX_*_MNT_LAMB93_IGN69.asc')):
            mt = re.search(r'FXX_(\d{4})_(\d{4})_MNT', f)
            X, Y = int(mt.group(1)) * 1000, int(mt.group(2)) * 1000
            self.tiles[(X, Y)] = f
        self._cache = {}

    def _tile(self, key):
        if key not in self._cache:
            with open(self.tiles[key]) as fh:
                hd = {k.lower(): float(v) for k, v in (fh.readline().split() for _ in range(6))}
            a = pd.read_csv(self.tiles[key], skiprows=6, sep=r'\s+', header=None,
                            dtype=np.float32, engine='c').to_numpy(dtype=float)
            a[a == hd['nodata_value']] = np.nan
            xc0 = hd['xllcorner'] + hd['cellsize'] / 2
            yc0 = hd['yllcorner'] + hd['cellsize'] * (hd['nrows'] - 0.5)
            self._cache[key] = (a, xc0, yc0)
        return self._cache[key]

    def disc(self, x, y, r):
        cs = self.CS
        pad = r + 15
        # 5 m lattice of cell centres (multiples of 5 m) covering the disc
        xs = np.arange(np.floor((x - pad) / cs) * cs, x + pad, cs)
        ys = np.arange(np.ceil((y + pad) / cs) * cs, y - pad, -cs)
        a = np.full((len(ys), len(xs)), np.nan)
        for (X, Y) in self.tiles:
            if X - 5 > xs[-1] or X + 5000 < xs[0] or Y + 5 < ys[-1] or Y - 5000 > ys[0]:
                continue
            t, xc0, yc0 = self._tile((X, Y))
            ci = np.round((xs - xc0) / cs).astype(int)
            ri = np.round((yc0 - ys) / cs).astype(int)
            okc, okr = (ci >= 0) & (ci < t.shape[1]), (ri >= 0) & (ri < t.shape[0])
            if okc.any() and okr.any():
                sub = t[np.ix_(ri[okr], ci[okc])]
                blk = a[np.ix_(okr, okc)]
                a[np.ix_(okr, okc)] = np.where(np.isfinite(sub), sub, blk)
        return _disc_from_block(a, xs[0], ys[0], cs, x, y, r)


def plane_fit(X, Y, z):
    A = np.column_stack([np.ones_like(X), X - X.mean(), Y - Y.mean()])
    coef, *_ = np.linalg.lstsq(A, z, rcond=None)
    b, c = coef[1], coef[2]
    resid = z - A @ coef
    return np.arctan(np.hypot(b, c)), np.arctan2(-b, -c), float(np.sqrt(np.mean(resid ** 2)))


def read_ascii_grid(f):
    with open(f) as fh:
        hd = {k.lower(): float(v) for k, v in (fh.readline().split() for _ in range(6))}
    a = np.loadtxt(f, skiprows=6)
    return hd, a


# ---------------------------------------------------------------- validation

def validate():
    res = {}
    lat = 46.0
    # (1) integrator vs FAO-56 eq. 21 extraterrestrial radiation (tau = 1, no atmosphere)
    for doy in (80, 172, 355):
        J = np.array([doy])
        dr = 1 + 0.033 * np.cos(2 * np.pi * J / 365)
        dec = 0.409 * np.sin(2 * np.pi * J / 365 - 1.39)
        phi = np.radians(lat)
        ws = np.arccos(-np.tan(phi) * np.tan(dec))
        ra = 24 * 60 / np.pi * 0.0820 * dr * (ws * np.sin(phi) * np.sin(dec) + np.cos(phi) * np.cos(dec) * np.sin(ws))
        h, _ = sun_geometry(lat, dec)
        num = (1367.0 * dr[:, None] * np.maximum(np.sin(h), 0)).sum() * DT_MIN * 60 / 1e6
        res[f'Ra_doy{doy}'] = {'fao56_MJ': float(ra[0]), 'numerical_MJ': float(num),
                               'rel_err': float(num / ra[0] - 1)}
    # (2) flat clear-sky daily totals, sea level and 3500 m, against FAO-56 Rso = (0.75 + 2e-5 z) Ra
    dec, e0 = spencer(DAYS)
    for z in (0.0, 3500.0):
        sd_, fd_, fdf_, _, _ = radiation(np.array([lat]), np.array([z]), np.zeros(1), np.zeros(1), None)
        f = fd_ + fdf_
        days_m = np.bincount(MONTH_OF_DAY)
        jun = f[0, 5] / days_m[5] / 1e6
        J = 172
        dr = 1 + 0.033 * np.cos(2 * np.pi * J / 365)
        dd = 0.409 * np.sin(2 * np.pi * J / 365 - 1.39)
        phi = np.radians(lat)
        ws = np.arccos(-np.tan(phi) * np.tan(dd))
        ra = 24 * 60 / np.pi * 0.0820 * dr * (ws * np.sin(phi) * np.sin(dd) + np.cos(phi) * np.cos(dd) * np.sin(ws))
        res[f'flat_june_mean_z{int(z)}'] = {'model_MJ_per_day': float(jun),
                                            'fao56_Rso_solstice_MJ': float((0.75 + 2e-5 * z) * ra),
                                            'annual_GJ': float(f[0].sum() / 1e9)}
    # zenith sun instantaneous, sea level
    m = 1.0
    res['zenith_sea_level_Wm2'] = {'beam': S0 * TAU ** m, 'diffuse': 0.3 * (1 - TAU ** m) * S0,
                                   'global': S0 * TAU ** m + 0.3 * (1 - TAU ** m) * S0}
    # solstice noon elevation
    h, _ = sun_geometry(lat, spencer(np.array([172]))[0])
    res['noon_elev_doy172_deg'] = float(np.degrees(h.max()))
    # (3) sky view of a tilted plane in flat surroundings (terrain horizon 0) = (1 + cos beta)/2
    sl = np.radians([0, 20, 40, 60])
    sv = sky_view(sl, np.radians([180, 90, 0, 270]), np.zeros((4, len(AZ))))
    res['svf_tilted'] = {f'{int(np.degrees(b))}deg': [float(v), float((1 + np.cos(b)) / 2)]
                         for b, v in zip(sl, sv)}
    # flat with uniform 20 deg horizon: cos^2(20)
    sv = sky_view(np.zeros(1), np.zeros(1), np.full((1, len(AZ)), np.radians(20)))
    res['svf_flat_20deg_horizon'] = [float(sv[0]), float(np.cos(np.radians(20)) ** 2)]
    # south 30-degree slope at 46 N: annual ratio, no terrain
    for name, a in (('south30', 180.0), ('north30', 0.0)):
        sd_, fd_, fdf_, _, sv0 = radiation(np.array([lat]), np.array([3500.0]), np.radians([30.0]),
                                           np.radians([a]), None)
        rj, ra, _ = ratios(sd_, fd_, fdf_, sv0)
        res[f'{name}_flat_surroundings_ratio_annual_JJA'] = [float(ra[0]), float(rj[0])]
    return res


# ---------------------------------------------------------------- sites and bands

def band_of(e):
    return int(np.floor(e / 10) * 10 + 5)


def load_entities():
    fc = pd.read_csv(FORCING, dtype={'glacier_id': str})
    fc = fc[fc.band_model >= MIN_BAND].copy()
    dh = pickle.load(open(CACHE, 'rb'))
    gm = {g.glacier_name: g for g in dh['glaciers']}
    fc['base_name'] = [gm[n].base_glacier_name for n in fc.glacier_name]
    fc['entity_band'] = [gm[n].elevation_band for n in fc.glacier_name]
    fc['wet'] = (fc.T_meas > WET_T).astype(int)
    fc['massif'] = fc.glacier_id.map(MASSIF_OF)
    assert fc.massif.notna().all(), fc[fc.massif.isna()]
    return fc.reset_index(drop=True)


def site_sunshine(ent, eu, dems, sun):
    bh = pd.read_csv(BOREHOLES)
    bh = bh[bh.elevation.notna()].copy()
    bh['band'] = bh.elevation.apply(band_of)
    rows = []
    for _, e in ent.iterrows():
        holes = bh[(bh.glacier_name == e.base_name) & (bh.band == e.entity_band)]
        assert len(holes), e.glacier_name
        for _, hrow in holes.iterrows():
            lon, lat = float(hrow.longitude), float(hrow.latitude)
            got = None
            for name, dem in dems:
                x, y = tf(dem.crs).transform(lon, lat)
                d = dem.disc(x, y, PLANE_R)
                if d is not None:
                    got = (name, dem, d)
                    break
            name, dem, (X, Y, z, cs) = got
            slope, asp_g, rms = plane_fit(X, Y, z)
            aspect = np.mod(asp_g - north_bearing(dem.crs, lon, lat), 2 * np.pi)
            rows.append(dict(glacier_name=e.glacier_name, hole_id=hrow.id, source_id=hrow.source_id,
                             lat=lat, lon=lon, elevation=float(hrow.elevation), dem=name,
                             cellsize=cs, n_cells=len(z), dem_z_mean=float(z.mean()),
                             plane_rms=rms, slope_deg=float(np.degrees(slope)),
                             aspect_deg=float(np.degrees(aspect))))
    holes = pd.DataFrame(rows)
    hor, z0 = np.zeros((len(holes), len(AZ))), np.zeros(len(holes))
    gid = holes.glacier_name.map(dict(zip(ent.glacier_name, ent.glacier_id)))
    for g in gid.unique():
        i = np.where(gid == g)[0]
        hor[i], z0[i] = eu.horizons(holes.lon.to_numpy()[i], holes.lat.to_numpy()[i])
    holes['eudem_z'] = z0
    holes['max_horizon_deg'] = np.degrees(hor.max(1))
    sd_, fd_, fdf_, svf, svf0 = radiation(holes.lat.to_numpy(), holes.elevation.to_numpy(),
                                          np.radians(holes.slope_deg.to_numpy()),
                                          np.radians(holes.aspect_deg.to_numpy()), hor, sun=sun)
    holes['svf'], holes['svf_floor0'] = svf, svf0
    rj, ra, rm = ratios(sd_, fd_, fdf_, svf)
    holes['ratio_jja'], holes['ratio_annual'] = rj, ra
    holes['ratio_jja_floor0'] = ratios(sd_, fd_, fdf_, svf0)[0]
    holes['direct_share_jja'] = sd_[:, JJA].sum(1) / (sd_[:, JJA].sum(1) + svf * fdf_[:, JJA].sum(1))
    for k in range(12):
        holes[f'ratio_m{k + 1:02d}'] = rm[:, k]
    return holes, hor


def band_sunshine(ent, sun):
    rows, pix_all = [], []
    eu = EuDem()
    for gid, g in ent.groupby('glacier_id'):
        hd, z = read_ascii_grid(RGI / f'dem/centraleurope/dem_{gid}.grid')
        _, m = read_ascii_grid(RGI / f'glmask/centraleurope/gl_{gid}.grid')
        cs, nr, nc = hd['cellsize'], int(hd['nrows']), int(hd['ncols'])
        dzr, dzc = np.gradient(z, cs)
        dzdx, dzdy = dzc, -dzr
        slope = np.arctan(np.hypot(dzdx, dzdy))
        asp_g = np.arctan2(-dzdx, -dzdy)
        rr, cc = np.mgrid[0:nr, 0:nc]
        X = hd['xllcorner'] + (cc + 0.5) * cs
        Y = hd['yllcorner'] + (nr - rr - 0.5) * cs
        inv = Transformer.from_crs(GL_CRS, 'EPSG:4326', always_xy=True)
        lonc, latc = inv.transform(X[m == 1].mean(), Y[m == 1].mean())
        gam = north_bearing(GL_CRS, lonc, latc)
        for band in sorted(g.band_model.unique()):
            half = BAND_HALF
            sel = (m == 1) & (z >= band - half) & (z < band + half)
            if sel.sum() < MIN_PIX:
                half = BAND_HALF_WIDE
                sel = (m == 1) & (z >= band - half) & (z < band + half)
            lon, lat = inv.transform(X[sel], Y[sel])
            hor, _ = eu.horizons(lon, lat)
            sl, asp = slope[sel], np.mod(asp_g[sel] - gam, 2 * np.pi)
            sd_, fd_, fdf_, svf, svf0 = radiation(np.asarray(lat), z[sel], sl, asp, hor, sun=sun)
            rj, ra, _ = ratios(sd_, fd_, fdf_, svf)
            rj0 = ratios(sd_, fd_, fdf_, svf0)[0]
            rows.append(dict(glacier_id=gid, band_model=band, band_half_width=half, n_pix=int(sel.sum()),
                             cellsize=cs, band_ratio_jja=float(rj.mean()), band_ratio_annual=float(ra.mean()),
                             band_ratio_jja_floor0=float(rj0.mean()),
                             band_ratio_jja_sd=float(rj.std()), band_slope_mean=float(np.degrees(sl).mean()),
                             band_svf_mean=float(svf.mean()),
                             band_south_frac=float(np.mean(np.cos(asp) < 0))))
            pix_all.append(pd.DataFrame(dict(glacier_id=gid, band_model=band, x=X[sel], y=Y[sel], z=z[sel],
                                             slope_deg=np.degrees(sl), aspect_deg=np.degrees(asp),
                                             svf=svf, ratio_jja=rj, ratio_annual=ra)))
            log.info('band %s@%d: %d px (+-%d m), JJA %.3f', gid, band, sel.sum(), half, rj.mean())
    return pd.DataFrame(rows), pd.concat(pix_all, ignore_index=True)


# ---------------------------------------------------------------- test

def lomo_probs(X, y, massif, C=C_REG):
    p = np.zeros(len(y))
    for mname in np.unique(massif):
        te = massif == mname
        tr = ~te
        mu, sd = X[tr].mean(0), X[tr].std(0)
        sd[sd == 0] = 1
        clf = LogisticRegression(C=C, max_iter=10000)
        clf.fit((X[tr] - mu) / sd, y[tr])
        p[te] = clf.predict_proba((X[te] - mu) / sd)[:, 1]
    return p


def block_permute(vals, glacier, rng):
    gl = np.unique(glacier)
    donors = rng.permutation(gl)
    out = np.empty_like(vals)
    for g, d in zip(gl, donors):
        idx = np.where(glacier == g)[0]
        dv = vals[glacier == d]
        if len(dv) >= len(idx):
            out[idx] = rng.choice(dv, len(idx), replace=False)
        else:
            out[idx] = np.concatenate([rng.permutation(dv), rng.choice(dv, len(idx) - len(dv), replace=True)])
    return out


def run_test(ent, col, C=C_REG, n_perm=N_PERM):
    y = ent.wet.to_numpy()
    T = ent.T_air_mb_notoff.to_numpy()
    s = ent[col].to_numpy()
    massif = ent.massif.to_numpy()
    glacier = ent.glacier_id.to_numpy()
    pA = lomo_probs(T[:, None], y, massif, C)
    pB = lomo_probs(np.column_stack([T, s]), y, massif, C)
    aucA, aucB = roc_auc_score(y, pA), roc_auc_score(y, pB)
    gain = aucB - aucA
    rng = np.random.default_rng(SEED)
    null = np.empty(n_perm)
    for k in range(n_perm):
        sp = block_permute(s, glacier, rng)
        null[k] = roc_auc_score(y, lomo_probs(np.column_stack([T, sp]), y, massif, C)) - aucA
    p = (1 + np.sum(null >= gain)) / (1 + n_perm)
    # in-sample fit for interpretation
    Xs = np.column_stack([T, s])
    Xs = (Xs - Xs.mean(0)) / Xs.std(0)
    clf = LogisticRegression(C=C, max_iter=10000).fit(Xs, y)
    return dict(variable=col, C=C, auc_without=float(aucA), auc_with=float(aucB), gain=float(gain),
                p_perm=float(p), null_mean=float(null.mean()), null_q95=float(np.quantile(null, 0.95)),
                coef_std_Tair=float(clf.coef_[0, 0]), coef_std_sun=float(clf.coef_[0, 1]),
                spearman_sun_Tair=float(pd.Series(s).corr(pd.Series(T), method='spearman'))), pA, pB, null


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s',
                        handlers=[logging.FileHandler(OUT_DIR / 'exposure_check.log', mode='w'),
                                  logging.StreamHandler()])
    val = validate()
    (OUT_DIR / 'validation.json').write_text(json.dumps(val, indent=2))
    log.info('validation: %s', json.dumps(val))

    ent = load_entities()
    log.info('%d entities, %d glaciers, %d massifs, %d wet', len(ent), ent.glacier_id.nunique(),
             ent.massif.nunique(), ent.wet.sum())
    sun = spencer(DAYS)
    eu = EuDem()
    dems = [('RGEALTI5', RgeAlti()), ('swissALTI3D10', SwissAlti()), ('EUDEM25', eu)]
    holes, hor = site_sunshine(ent, eu, dems, sun)
    holes.to_csv(OUT_DIR / 'site_holes.csv', index=False)
    np.save(OUT_DIR / 'site_horizons_rad.npy', hor)
    agg = holes.groupby('glacier_name').agg(
        n_holes=('hole_id', 'size'), site_ratio_jja=('ratio_jja', 'mean'),
        site_ratio_annual=('ratio_annual', 'mean'), site_ratio_jja_floor0=('ratio_jja_floor0', 'mean'),
        site_ratio_jja_sd=('ratio_jja', 'std'),
        site_slope=('slope_deg', 'mean'), site_aspect=('aspect_deg', lambda a: float(np.degrees(
            np.arctan2(np.sin(np.radians(a)).mean(), np.cos(np.radians(a)).mean())) % 360)),
        site_svf=('svf', 'mean'), site_dem=('dem', lambda d: '+'.join(sorted(set(d))))).reset_index()
    ent = ent.merge(agg, on='glacier_name', how='left')
    bands, pix = band_sunshine(ent, sun)
    pix.to_csv(OUT_DIR / 'band_pixels.csv', index=False)
    ent = ent.merge(bands, on=['glacier_id', 'band_model'], how='left')
    assert ent.site_ratio_jja.notna().all() and ent.band_ratio_jja.notna().all()

    results = {}
    deviations = [
        'Validation case (3) re-specified: (1 + cos beta)/2 is the sky view of a tilted plane in flat '
        'surroundings (terrain horizon 0); with only the plane\'s own horizon the pre-registered formula '
        'correctly gives 1. Main test unchanged (pre-registered formula, horizon not floored).',
        'Added non-decisive sensitivity: terrain horizon floored at 0 deg for the sky view (sky below the '
        'horizontal excluded), columns *_floor0.',
    ]
    probs = {}
    for col in ('band_ratio_jja', 'site_ratio_jja', 'band_ratio_annual', 'site_ratio_annual'):
        r, pA, pB, null = run_test(ent, col)
        results[col] = r
        probs[col] = (pA, pB)
        np.save(OUT_DIR / f'null_{col}.npy', null)
        log.info('%s', r)
    for col in ('band_ratio_jja_floor0', 'site_ratio_jja_floor0'):
        r, *_ = run_test(ent, col)
        results[col] = r
        log.info('%s', r)
    for col in ('band_ratio_jja', 'site_ratio_jja'):
        r, *_ = run_test(ent, col, C=1e6)
        results[f'{col}_C1e6'] = r
        log.info('%s', r)
    ent['p_without'] = probs['band_ratio_jja'][0]
    ent['p_with_band'] = probs['band_ratio_jja'][1]
    ent['p_with_site'] = probs['site_ratio_jja'][1]
    ent.to_csv(OUT_DIR / 'entities.csv', index=False)
    b = results['band_ratio_jja']
    verdict = 'PASS' if (b['gain'] >= PASS_GAIN and b['p_perm'] < PASS_P) else 'FAIL'
    out = dict(n_entities=len(ent), n_glaciers=int(ent.glacier_id.nunique()),
               n_massifs=int(ent.massif.nunique()), n_wet=int(ent.wet.sum()), n_holes=len(holes),
               dem_use=holes.dem.value_counts().to_dict(),
               widened_bands=bands[bands.band_half_width > BAND_HALF][['glacier_id', 'band_model', 'n_pix']]
               .to_dict('records'),
               results=results, verdict=verdict, deviations=deviations)
    (OUT_DIR / 'results.json').write_text(json.dumps(out, indent=2, default=str))
    log.info('VERDICT %s', verdict)


if __name__ == '__main__':
    main()
