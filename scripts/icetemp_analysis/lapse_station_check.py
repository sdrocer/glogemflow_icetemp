"""Lapse-rate check at high Swiss stations: ERA5 air temperature extrapolated the way GloGEM does
it, with ERA5's monthly lapse rate dtdz ('reanalysis', mass balance) or with the thermal floor
min(dtdz, -0.0065 K/m) ('floor', thermal model), against measured monthly mean air temperature.

Station data: MeteoSwiss open data (data.geo.admin.ch), homogenised NBCN series (ths200m0) where
one exists, else the SwissMetNet monthly mean (tre200m0). Files are downloaded beforehand into
<out_dir>/raw and read here as data only.

Decision rule (fixed before any result): the version with the smaller |mean annual bias| at the
stations above 3000 m (Jungfraujoch, Corvatsch) wins; if the two differ by less than 0.3 C there,
keep the floor.

Usage: python scripts/icetemp_analysis/lapse_station_check.py [out_dir]
"""
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from forcing_check import read_climate, band_air

OUT_DIR = Path(sys.argv[1] if len(sys.argv) > 1
               else '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week1/lapse_stations')
RAW = OUT_DIR / 'raw'
PERIOD = (1981, 2020)
DECISION_STATIONS = ['JUN', 'COV']
TIE_MARGIN = 0.3
SEASONS = {'DJF': (12, 1, 2), 'MAM': (3, 4, 5), 'JJA': (6, 7, 8), 'SON': (9, 10, 11)}
VERSIONS = {'reanalysis': 'mb_notoff', 'floor': 'th_notoff'}
STATIONS = {  # abbr: (file, parameter)
    'JUN': ('ogd-nbcn_jun_m.csv', 'ths200m0'),
    'COV': ('ogd-smn_cov_m.csv', 'tre200m0'),
    'WFJ': ('ogd-smn_wfj_m.csv', 'tre200m0'),
    'SAE': ('ogd-nbcn_sae_m.csv', 'ths200m0'),
}
CROSSCHECK = {'JUN': ('ogd-smn_jun_m.csv', 'tre200m0'), 'SAE': ('ogd-smn_sae_m.csv', 'tre200m0')}
BASE_URL = 'https://data.geo.admin.ch/ch.meteoschweiz.{coll}/{abbr}/{file}'

log = logging.getLogger('lapse_station_check')


def read_meta():
    m = pd.read_csv(RAW / 'ogd-smn_meta_stations.csv', sep=';', encoding='latin1', dtype=str)
    m = m.set_index('station_abbr')
    out = {}
    for s in STATIONS:
        r = m.loc[s]
        out[s] = dict(name=r['station_name'], elev=float(r['station_height_masl']),
                      lat=float(r['station_coordinates_wgs84_lat']),
                      lon=float(r['station_coordinates_wgs84_lon']),
                      exposition=r['station_exposition_en'], data_since=r['station_data_since'])
    return out


def read_station(fname, param):
    d = pd.read_csv(RAW / fname, sep=';', dtype=str)
    t = pd.to_datetime(d['reference_timestamp'], format='%d.%m.%Y %H:%M')
    v = pd.to_numeric(d[param], errors='coerce')
    s = pd.Series(v.values, index=pd.MultiIndex.from_arrays([t.dt.year, t.dt.month],
                                                            names=['year', 'month']))
    return s[(s.index.get_level_values(0) >= PERIOD[0]) & (s.index.get_level_values(0) <= PERIOD[1])]


def era5_at(clim, ilat, ilon, elev, months):
    idx = [clim['index'][ym] for ym in months]
    T = clim['rtemp'][idx, ilon, ilat]
    mon = np.array([m for _, m in months])
    hclim = float(clim['relev'][ilon, ilat])
    air = band_air(T, mon, elev, hclim, clim['rtg'][:, ilon, ilat], 0.0)
    return {v: pd.Series(air[k], index=pd.MultiIndex.from_tuples(months, names=['year', 'month']))
            for v, k in VERSIONS.items()}, hclim, T


def stats(model, obs):
    ok = obs.notna()
    diff = (model - obs)[ok]
    mon = diff.index.get_level_values(1)
    out = {'annual': float(diff.groupby(mon).mean().mean())}   # mean over calendar months
    for s, ms in SEASONS.items():
        out[s] = float(diff[np.isin(mon, ms)].mean())
    am = model[ok] - model[ok].groupby(mon).transform('mean')
    ao = obs[ok] - obs[ok].groupby(mon).transform('mean')
    out['rms_anom'] = float(np.sqrt(((am - ao) ** 2).mean()))
    out['n_months'] = int(ok.sum())
    return out


def main():
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    clim = read_climate()
    meta = read_meta()
    months = [(y, m) for y in range(PERIOD[0], PERIOD[1] + 1) for m in range(1, 13)]
    rlat, rlon = clim['rlat'], clim['rlon']
    rows, summary = [], {'period': PERIOD, 'stations': {}}
    for s, (fname, param) in STATIONS.items():
        md = meta[s]
        obs = read_station(fname, param).reindex(pd.MultiIndex.from_tuples(months))
        coll = 'ogd-nbcn' if fname.startswith('ogd-nbcn') else 'ogd-smn'
        il0 = int(np.argmin(np.abs(rlat - md['lat'])))
        io0 = int(np.argmin(np.abs(rlon - md['lon'])))
        st = dict(meta=md, source=BASE_URL.format(coll=coll, abbr=s.lower(), file=fname),
                  parameter=param, n_missing=int(obs.isna().sum()), nodes=[])
        for dla in (-1, 0, 1):
            for dlo in (-1, 0, 1):
                il, io = il0 + dla, io0 + dlo
                ser, hclim, T = era5_at(clim, il, io, md['elev'], months)
                node = dict(nearest=(dla == 0 and dlo == 0), lat=float(rlat[il]), lon=float(rlon[io]),
                            hclim=hclim, dtdz_mean=float(clim['rtg'][:, io, il].mean()))
                for v, mod in ser.items():
                    r = stats(mod, obs)
                    node[v] = r
                    rows.append(dict(station=s, version=v, **{k: node[k] for k in
                                     ('nearest', 'lat', 'lon', 'hclim', 'dtdz_mean')}, **r))
                # lapse rate the station implies between the node and the station, per season
                ok = obs.notna().values
                mon = np.array([m for _, m in months])
                dT = obs.values[ok] - T[ok]
                node['implied_lapse'] = {k: float(dT[np.isin(mon[ok], ms)].mean() / (md['elev'] - hclim))
                                         for k, ms in [('annual', tuple(range(1, 13)))] + list(SEASONS.items())}
                node['dtdz_season'] = {k: float(clim['rtg'][[m - 1 for m in ms], io, il].mean())
                                       for k, ms in SEASONS.items()}
                st['nodes'].append(node)
        if s in CROSSCHECK:
            f2, p2 = CROSSCHECK[s]
            o2 = read_station(f2, p2).reindex(pd.MultiIndex.from_tuples(months))
            st['crosscheck_raw_minus_homog'] = float((o2 - obs).mean())
        summary['stations'][s] = st
    pd.DataFrame(rows).to_csv(OUT_DIR / 'lapse_station_stats.csv', index=False, float_format='%.4f')

    # decision
    nearest = {s: next(n for n in st['nodes'] if n['nearest']) for s, st in summary['stations'].items()}
    absb = {v: float(np.mean([abs(nearest[s][v]['annual']) for s in DECISION_STATIONS])) for v in VERSIONS}
    per = {s: {v: nearest[s][v]['annual'] for v in VERSIONS} for s in DECISION_STATIONS}
    diff = abs(absb['reanalysis'] - absb['floor'])
    winner = 'floor' if diff < TIE_MARGIN else min(absb, key=absb.get)
    # sensitivity: the same rule on every node of the 3x3 block (paired by position)
    sens = []
    for k in range(9):
        a = {v: float(np.mean([abs(summary['stations'][s]['nodes'][k][v]['annual'])
                               for s in DECISION_STATIONS])) for v in VERSIONS}
        d = abs(a['reanalysis'] - a['floor'])
        sens.append('floor' if d < TIE_MARGIN else min(a, key=a.get))
    summary['decision'] = dict(mean_abs_annual_bias=absb, per_station_annual_bias=per,
                               difference=diff, tie_margin=TIE_MARGIN, winner=winner,
                               winner_by_node_position=sens)
    (OUT_DIR / 'lapse_station_summary.json').write_text(json.dumps(summary, indent=1))

    for s, st in summary['stations'].items():
        n = nearest[s]
        log.info(f"\n{s} {st['meta']['name']} {st['meta']['elev']:.0f} m ({st['parameter']}, "
                 f"missing {st['n_missing']}/{len(months)}) node {n['lat']:.2f}N {n['lon']:.2f}E "
                 f"hclim {n['hclim']:.0f} m, dE {st['meta']['elev'] - n['hclim']:.0f} m, "
                 f"mean dtdz {n['dtdz_mean'] * 1000:.2f} K/km")
        for v in VERSIONS:
            r = n[v]
            log.info(f"  {v:10s} annual {r['annual']:+.2f}  " + '  '.join(
                f"{k} {r[k]:+.2f}" for k in SEASONS) + f"  rms_anom {r['rms_anom']:.2f}")
        log.info('  implied lapse K/km: ' + '  '.join(f"{k} {x * 1000:.2f}" for k, x in n['implied_lapse'].items())
                 + ' | ERA5 dtdz: ' + '  '.join(f"{k} {x * 1000:.2f}" for k, x in n['dtdz_season'].items()))
        log.info('  3x3 annual bias (hclim): ' + '; '.join(
            f"{x['hclim']:.0f} m R {x['reanalysis']['annual']:+.2f} F {x['floor']['annual']:+.2f}"
            for x in st['nodes']))
        if 'crosscheck_raw_minus_homog' in st:
            log.info(f"  raw minus homogenised station mean: {st['crosscheck_raw_minus_homog']:+.2f}")
    log.info(f"\nDecision: mean |annual bias| JUN+COV {absb}; difference {diff:.2f} -> {winner}; "
             f"by node position: {sens}")


if __name__ == '__main__':
    main()
