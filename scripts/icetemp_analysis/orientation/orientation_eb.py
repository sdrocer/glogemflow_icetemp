"""Hourly surface energy balance behind the orientation factor F (orientation investigation of 20-21 Oct 2026,
designer/verify_eb.py, revise_F.py, revise_sens.py; plan of 21 Oct, step 1b).

Q = (1 - albedo) S + psi, psi = psi0 + c T for T >= T_tip, else psi0 (Giesen & Oerlemans 2012, TC 6, eq. 3);
melt = sum of positive Q / L_f. S is clear-sky radiation on the plane (Spencer sun, Kasten & Young air mass,
transmissivity 0.75, terrain horizon; exposure_check.py of week 1), weighted by Jungfraujoch's daily relative
sunshine (direct beam) and clearness (diffuse, times the sky view). The hourly air temperature is the band's monthly
temperature plus Jungfraujoch's daily anomaly and a sine diurnal cycle of the daily range, 1981-2020.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import exposure_check as ec  # noqa: E402

LF = 334000.0
Y0, Y1 = 1981, 2020
NY = Y1 - Y0 + 1
MON_LEN = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
DOY_M = np.repeat(np.arange(12), MON_LEN)
JJAS = [5, 6, 7, 8]
JFJ = dict(lat=46.547556, z=3571.0)


def clear_sky_hourly(lat, z, slope, aspect, hor, tau=ec.TAU):
    """Hour means [365, 24]: clear-sky direct on the plane (terrain horizon), flat direct, flat diffuse, sun up."""
    dec, e0 = ec.spencer(ec.DAYS)
    h, az = ec.sun_geometry(lat, dec)
    up = h > 0
    m = np.where(up, ec.air_mass(np.where(up, h, 0.1), z), 0)
    tm = np.where(up, tau ** m, 0)
    sp = ec.S0 * e0[:, None]
    beam = sp * tm
    fdif = 0.3 * (1 - tm) * sp * np.sin(h) * up
    fdir = beam * np.sin(h) * up
    cosi = np.sin(h) * np.cos(slope) + np.cos(h) * np.sin(slope) * np.cos(az - aspect)
    lit = up & (cosi > 0)
    if hor is not None:
        hh = np.interp(np.degrees(az), np.append(ec.AZ, 360.0), np.append(hor, hor[0]))
        lit &= h > hh
    sdir = beam * np.maximum(cosi, 0) * lit
    r = lambda a: a.reshape(365, 24, 6).mean(2)  # noqa: E731
    return r(sdir), r(fdir), r(fdif), up.reshape(365, 24, 6).mean(2)


def read_station(path):
    """MeteoSwiss daily series (ogd-smn_*_d_historical.csv), 1981-2020 without 29 February."""
    d = pd.read_csv(path, sep=';', encoding='latin1')
    d['t'] = pd.to_datetime(d.reference_timestamp, format='%d.%m.%Y %H:%M')
    d = d[(d.t.dt.year >= Y0) & (d.t.dt.year <= Y1)].set_index('t')
    idx = pd.date_range(f'{Y0}-01-01', f'{Y1}-12-31', freq='D')
    idx = idx[~((idx.month == 2) & (idx.day == 29))]
    d = d.reindex(idx)
    out = pd.DataFrame(index=idx)
    out['T'] = d.tre200d0.interpolate(limit=3)
    out['R'] = (d.tre200dx - d.tre200dn).interpolate(limit=3)
    out['G'] = d.gre000d0
    out['SD'] = d.sre000d0
    out['doy'] = np.tile(np.arange(365), NY)
    out['year'], out['month'] = out.index.year, out.index.month
    for c in ('G', 'SD', 'R'):
        out[c] = out[c].fillna(out.groupby('doy')[c].transform('mean'))
    out['T_anom'] = out['T'] - out.groupby(['year', 'month'])['T'].transform('mean')
    return out


def sky_weights(w, st=JFJ):
    """Daily direct-beam fraction sig and diffuse scaling dif of the station, and its clear-sky flat direct."""
    sd_j, fd_j, fi_j, up_j = clear_sky_hourly(st['lat'], st['z'], 0.0, 0.0, None)
    d = w.doy.to_numpy()
    gclr = (fd_j + fi_j).mean(1)
    sdmax = up_j.sum(1) * 60.0
    k = np.clip(w.G.to_numpy() / gclr[d], 0.02, 1.15)
    sig = np.clip(w.SD.to_numpy() / sdmax[d], 0.0, 1.0)
    G = k * gclr[d]
    over = sig * fd_j.mean(1)[d] > 0.95 * G
    sig = np.where(over, 0.95 * G / fd_j.mean(1)[d], sig)
    dif = (G - sig * fd_j.mean(1)[d]) / fi_j.mean(1)[d]
    return sig, dif, fd_j


def concentrate(sig, fd_st):
    """[NY, 365, 24] beam weights for sunny spells: full clear-sky beam in the hours nearest noon until the day's
    direct energy on flat ground at the station equals sig x its clear-sky daily direct; the last hour fractional."""
    doy = np.tile(np.arange(365), NY)
    out = np.zeros((len(sig), 24))
    order = np.argsort(-fd_st, axis=1)
    for i in range(len(sig)):
        dd = doy[i]
        tot = sig[i] * fd_st[dd].sum()
        acc = 0.0
        for hh in order[dd]:
            b = fd_st[dd, hh]
            if b <= 0 or acc >= tot:
                break
            take = min(1.0, (tot - acc) / b)
            out[i, hh] = take
            acc += take * b
    return out.reshape(NY, 365, 24)


class Forcing:
    """Hourly air temperature Th [NY, 365, 24] of one band and the all-sky weights of the station."""

    def __init__(self, Tm, w, sig, dif, conc=None):
        hrs = np.arange(24) + 0.5
        Td = Tm[:, DOY_M] + w.T_anom.to_numpy().reshape(NY, 365)
        R = w.R.to_numpy().reshape(NY, 365)
        self.Th = Td[:, :, None] + 0.5 * R[:, :, None] * np.sin(2 * np.pi * (hrs - 9.0) / 24.0)[None, None, :]
        self.sig = sig.reshape(NY, 365)[:, :, None] if conc is None else conc
        self.dif = dif.reshape(NY, 365)[:, :, None]

    def shortwave(self, sdir, fdif, svf):
        return self.sig * sdir[None] + self.dif * svf * fdif[None]


def monthly_melt(Th, S, psi0, c, alb=0.80, ttip=1.0):
    """Climatological monthly energy-balance melt [12], mm w.e. per month."""
    psi = np.where(Th >= ttip, psi0 + c * Th, psi0)
    q = np.maximum((1 - alb) * S + psi, 0).sum(2) * 3600 / LF
    out = np.zeros(12)
    np.add.at(out, DOY_M, q.sum(0))
    return out / NY


def jjas(m):
    return float(m[JJAS].sum())


def band_T(F, ent):
    """Monthly band air temperature [NY, 12] (GloGEM ERA5, no t_offset, ERA5 lapse rate), 1981-2020."""
    f = F[(F.entity == ent) & (F.year >= Y0) & (F.year <= Y1)].sort_values(['year', 'month'])
    return f['T'].to_numpy().reshape(NY, 12)
