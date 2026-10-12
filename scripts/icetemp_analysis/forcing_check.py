"""Step 0 forcing check: measured 10-20 m firn temperature against GloGEM's band air temperature.

In cold, dry firn the temperature at 10-20 m is close to the mean annual air temperature. This
rebuilds the band air temperature of the ce25_pred_w0 training run from the ERA5 files the way
GloGEM builds it (climate node of the 0.1-degree loop cell, monthly lapse rate dtdz, calibrated
T_off), averages it over the 120 months before each point's first reading and compares it with the
measured 10-20 m mean, for {mass-balance lapse rate dtdz, thermal lapse rate} x {with T_off,
without}. The site selections are fixed below, before any result is computed. Validation: the
rebuilt thermal surface temperature against the T_surf the run logged (firnice_snowclim_<id>.dat).

Usage: python scripts/icetemp_analysis/forcing_check.py [out_dir]
"""
import json
import logging
import pickle
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
RUN = REPO / 'data/bayescal/ce25_pred_w0/training_runs'
FIRNICE = RUN / 'monthly/CentralEurope/PAST/firnice_temperature'
CALI = RUN / 'monthly/CentralEurope/calibration/calibrate_m1_cID9_centraleurope_final_era5.dat'
LOOKUP = RUN / 'glenglat_borehole_elevations.dat'
GLOGEM_DATA = Path('/itet-stor/jabeer/glogem')
CLIM = GLOGEM_DATA / 'climatedata/reanalysis/monthly/era5/CentralEurope/clim_CentralEurope.mdi'
TGRAD = CLIM.with_name('tgrad_CentralEurope.mdi')
INV = GLOGEM_DATA / 'geometricdata/rgiv7/files/thick_centraleurope.dat'
INV_HF = GLOGEM_DATA / 'geometricdata/rgiv7/files_HF/thick_centraleurope.dat'
CATCHMENT = GLOGEM_DATA / 'data/catchments/RGI11_CentralEurope_bayescal.dat'
CACHE = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/cache/datahandler.pkl')
OUT_DIR = Path(sys.argv[1] if len(sys.argv) > 1
               else '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/forcing_check')

# Selections, fixed before any result was computed
DEPTH_RANGE = (10.0, 20.0)   # m, points compared with the air temperature
WINDOW_MONTHS = 120          # months before a point's first reading month
SEL_A_MIN_BAND = 4400        # (a) model band at or above [m a.s.l.]
SEL_B_MAX_T = -5.0           # (b) measured 10-20 m mean below [C]
TOFF_MIN = 2.0               # T_off subset: T_off at or above [C]
VALID_TOL = 0.01             # C, rebuilt vs logged thermal surface temperature

# GloGEM settings of this run
F32 = np.float32
GRID_STEP = F32(0.1)            # settings.pro grid_step: a single-precision literal in IDL
SIZE_RANGE = (0.0, 100000.0)    # glogem.pro:183, catchment runs
THERMAL_LAPSE = -0.0065         # initialise_firnicetemp_spinup.pro: dtdz_icetemp = dtdz < -0.0065
RUN_MONTHS = [(y, m) for y in range(1971, 2026) for m in range(1, 13)
              if (1971, 10) <= (y, m) <= (2025, 9)]          # tran = [1972, 2025], bal_month = 10

VARIANTS = {'mb_toff': 'dtdz, with T_off', 'mb_notoff': 'dtdz, without T_off',
            'th_toff': 'thermal lapse, with T_off', 'th_notoff': 'thermal lapse, without T_off'}
# the thermal surface as the model applies it (capped at 0 C each month), reported alongside
APPLIED = {'surf_toff': 'applied surface, with T_off', 'surf_notoff': 'applied surface, without T_off'}

log = logging.getLogger('forcing_check')


class IdlReader:
    """IDL free-format readf: a call fills n values from as many lines as it needs and drops the
    rest of its last line."""

    def __init__(self, path):
        self.lines, self.i = Path(path).read_text().splitlines(), 0

    def skip(self):
        self.i += 1

    def read(self, n):
        vals = []
        while len(vals) < n:
            vals += self.lines[self.i].split()
            self.i += 1
        return np.array(vals[:n], dtype=float)

    def scalar(self):
        return int(self.read(1)[0])


def read_climate():
    """clim_ and tgrad_ files as read_climatepast_monthly.pro reads them, with its month index."""
    r = IdlReader(CLIM)
    r.skip()
    ntime, nlons, nlats = r.scalar(), r.scalar(), r.scalar()
    r.read(2)                                                   # nvar
    rtime, rlon, rlat = r.read(ntime), r.read(nlons), r.read(nlats)
    relev = np.array([r.read(nlats) for _ in range(nlons)])
    rtemp = np.array([r.read(nlats) for _ in range(ntime * nlons)]).reshape(ntime, nlons, nlats)
    t = IdlReader(TGRAD)
    t.skip()
    nmon, tlons, tlats = t.scalar(), t.scalar(), t.scalar()
    t.read(1)                                                   # nvar
    t.read(nmon)
    tlon, tlat = t.read(tlons), t.read(tlats)
    rtg = np.array([t.read(tlats) for _ in range(nmon * tlons)]).reshape(nmon, tlons, tlats)
    rlon = np.where(rlon >= 180, rlon - 360, rlon)
    tlon = np.where(tlon >= 180, tlon - 360, tlon)

    # era5 month hack (read_climatepast_monthly.pro:64-70); the float literals and the INT/FLOAT
    # term are single precision in IDL, round() is half away from zero (x > 0 here)
    ryear = np.trunc(rtime).astype(int)
    x = (((rtime + float(F32(0.04))) - ryear) * 12 + float(F32(0.6))
         - (F32(ryear - ryear[0]) / F32(2800.0) * F32(12)).astype(float))
    rmon = np.floor(x + 0.5).astype(int)
    index = {}
    for i, ym in enumerate(zip(ryear.tolist(), rmon.tolist())):
        index.setdefault(ym, i)                                 # where(...)[0]: first match
    calendar = [(1940 + i // 12, i % 12 + 1) for i in range(ntime)]
    log.info(f'climate: {ntime} months, rtime {rtime[0]:.4f}..{rtime[-1]:.4f} (fractional year at '
             f'the start of each month); (ryear, rmon) = calendar month for every entry: '
             f'{list(zip(ryear.tolist(), rmon.tolist())) == calendar}, unique: {len(index) == ntime}, '
             f'closest rounding margin {np.min(np.abs(x - np.floor(x) - 0.5)):.3f}')
    same_grid = (tlon.shape == rlon.shape and tlat.shape == rlat.shape
                 and np.allclose(tlon, rlon) and np.allclose(tlat, rlat))
    log.info(f'tgrad grid identical to clim grid: {same_grid} (GloGEM indexes rtg with the clim '
             f'indices either way); dtdz range {rtg.min():.5f}..{rtg.max():.5f} K/m')
    if not same_grid:
        log.warning('tgrad grid differs from the clim grid; dtdz follows the code (clim indices)')
    return dict(rlat=rlat, rlon=rlon, relev=relev, rtemp=rtemp, rtg=rtg, index=index)


def read_inventory(path):
    """thick_<region>.dat as read_glacier_inventory_batchfile.pro: id, then 19 numbers per row."""
    ids, vals = [], []
    for line in path.read_text().splitlines()[1:]:
        a = line.split()
        ids.append(a[0].split(';')[0])
        vals.append([float(s.rstrip(';')) for s in a[1:20]])
    return ids, np.array(vals)


def climate_nodes(run_ids, clim):
    """Climate node of each run glacier: nearest to the mean position of every inventory glacier in
    its loop cell (setup_grid_loops.pro:19-22, glogem.pro:314-343, downscale_gcmdata_monthly.pro:21-42)."""
    ids, inv = read_inventory(INV)
    ids_hf, inv_hf = read_inventory(INV_HF)
    if ids_hf != ids:
        raise SystemExit('files and files_HF inventories are not in the same order')
    lon, lat, vol = inv[:, 0], inv[:, 1], inv[:, 3]
    area = inv_hf[:, 2]                        # RGI7: areas from files_HF
    run = np.array([i in run_ids for i in ids])
    ok = (vol > 0) & (area > SIZE_RANGE[0]) & (area < SIZE_RANGE[1])
    log.info(f'inventory: {len(ids)} glaciers, {int(ok.sum())} pass volume > 0 and the size range; '
             f'{int(run.sum())} run glaciers')

    def bounds(v):                             # loop extent from the run glaciers, single precision
        lo = F32(F32(np.floor(v.min() / float(GRID_STEP))) * GRID_STEP) - GRID_STEP / F32(2)
        hi = (F32(F32(np.trunc(v.max() / float(GRID_STEP))) * GRID_STEP)
              + GRID_STEP / F32(2) + F32(2) * GRID_STEP)
        return lo, hi

    lon0, lat0 = bounds(lon[run]), bounds(lat[run])
    ngx = max(int(np.trunc((lon0[1] - lon0[0]) / GRID_STEP)), 1)
    ngy = max(int(np.trunc((lat0[1] - lat0[0]) / GRID_STEP)), 1)
    log.info(f'loop cells: lon0 {lon0[0]:.5f}..{lon0[1]:.5f} ({ngx}), lat0 {lat0[0]:.5f}..'
             f'{lat0[1]:.5f} ({ngy}), half-open, edges in single precision')
    rlat, rlon, relev = clim['rlat'], clim['rlon'], clim['relev']
    nodes = {}
    for gx in range(ngx):
        lo = [float(lon0[0] + F32(gx + k) * GRID_STEP) for k in (0, 1)]
        for gy in range(ngy):
            la = [float(lat0[0] + F32(gy + k) * GRID_STEP) for k in (0, 1)]
            cell = (lat >= la[0]) & (lat < la[1]) & (lon >= lo[0]) & (lon < lo[1]) & ok
            batch = np.flatnonzero(cell & run)
            if len(batch) == 0:
                continue
            gc = np.flatnonzero(cell)
            gmid = (lat[gc].mean(), lon[gc].mean())
            il, io = int(np.argmin(np.abs(rlat - gmid[0]))), int(np.argmin(np.abs(rlon - gmid[1])))
            for k in batch:
                nodes[ids[k]] = dict(ilat=il, ilon=io, lat=float(rlat[il]), lon=float(rlon[io]),
                                     hclim=float(relev[io, il]), cell=[gx, gy], cell_lon=lo,
                                     cell_lat=la, n_cell=len(gc), mean_lat=float(gmid[0]),
                                     mean_lon=float(gmid[1]), glacier_lat=float(lat[k]),
                                     glacier_lon=float(lon[k]),
                                     dtdz=clim['rtg'][:, io, il].tolist())
    missing = sorted(set(run_ids) - set(nodes))
    if missing:
        raise SystemExit(f'run glaciers without a loop cell: {missing}')
    return nodes


def read_toff(run_ids):
    """T_off per glacier: column 12 of the calibration file (read_calibration_params.pro:30-38), at
    the numerically nearest ID (apply_calibration_params.pro:15). The flag fill touches columns
    9-11 only."""
    r = IdlReader(CALI)
    r.skip()
    n = len(r.lines) - 1
    da = r.read(13 * n).reshape(n, 13)
    out = {}
    for gid in run_ids:
        k = int(np.argmin(np.abs(float(gid) - da[:, 0])))
        if da[k, 0] != float(gid):
            log.warning(f'{gid}: not in the calibration file, nearest ID {da[k, 0]:05.0f} used')
        out[gid] = float(da[k, 11])
    return out


def read_lookup():
    """glacier_id -> borehole band elevations in profile order (setup_firnice_profile_from_glenglat.pro)."""
    out = {}
    for line in LOOKUP.read_text().splitlines():
        parts = line.split()
        if parts and not parts[0].startswith('#'):
            out.setdefault(parts[0], [float(p.split(':')[0]) for p in parts[1:]])
    return out


def read_snowclim(gid):
    """(elev, month, T_surf) rows, fixed format (i6, i6, f10.2, 3f7.3, 4f8.2, 2f8.3)."""
    rows = [(int(ln[0:6]), int(ln[6:12]), float(ln[43:51]))
            for ln in (FIRNICE / f'firnice_snowclim_{gid}.dat').read_text().splitlines()
            if ln.strip() and not ln.startswith('#')]
    return np.array(rows, dtype=float)


def model_band(g, lookup, bands):
    """Band the run compared the borehole with: the 'Point elevation' of its profile file, else
    the nearest glacierized band (prepare_output_firnicetemp.pro:106-114)."""
    near = int(bands[np.argmin(np.abs(bands - g.elevation_band))])
    elevs = lookup.get(g.glacier_id, [])
    if g.elevation_band in elevs:
        k = elevs.index(g.elevation_band) + 1
        head = (FIRNICE / f'temp_ID{k}_{g.glacier_id}.dat').open().readline()
        band = int(re.search(r'Point elevation\s+(-?\d+)', head).group(1))
        if band != near:
            log.warning(f'{g.glacier_name}: profile header {band} m, nearest glacierized band {near} m')
        return band, f'temp_ID{k}'
    return near, 'nearest band'


def band_air(T, mon, E, hclim, dtdz, toff):
    """Monthly band air temperature variants: tg (glogem.pro:530) with the mass-balance lapse rate,
    and the thermal surface (firnice_temperature_model.pro:237) before and after its 0 C cap."""
    g = dtdz[mon - 1]
    mb = T + (E - hclim) * g
    th = mb + max(E - hclim, 0.0) * (np.minimum(g, THERMAL_LAPSE) - g)
    return {'mb_toff': mb + toff, 'mb_notoff': mb, 'th_toff': th + toff, 'th_notoff': th,
            'surf_toff': np.minimum(0.0, th + toff), 'surf_notoff': np.minimum(0.0, th)}


def node_series(clim, node, months):
    idx = [clim['index'].get(ym) for ym in months]
    if None in idx:
        raise SystemExit(f'months missing from the climate record: {months[idx.index(None)]}')
    return clim['rtemp'][idx, node['ilon'], node['ilat']], np.array([m for _, m in months])


def validate(clim, nodes, toff):
    """Rebuilt thermal surface temperature, mean per calendar month over the run, against T_surf."""
    per_glacier, n_rows = {}, 0
    for gid in sorted(nodes):
        nd = nodes[gid]
        T, mon = node_series(clim, nd, RUN_MONTHS)
        dtdz = np.array(nd['dtdz'])
        worst = 0.0
        for E, m, logged in read_snowclim(gid):
            surf = band_air(T, mon, E, nd['hclim'], dtdz, toff[gid])['surf_toff']
            worst = max(worst, abs(surf[mon == int(m)].mean() - logged))
            n_rows += 1
        per_glacier[gid] = worst
        log.info(f'validation {gid}: max |rebuilt - logged T_surf| {worst:.4f} C')
    vmax = max(per_glacier.values())
    return dict(max_abs_diff=vmax, n_rows=n_rows, n_glaciers=len(per_glacier),
                tolerance=VALID_TOL, passed=bool(vmax <= VALID_TOL), per_glacier=per_glacier)


def month_add(ym, k):
    i = (ym // 100) * 12 + ym % 100 - 1 + k
    return (i // 12, i % 12 + 1)


def entity_rows(glaciers, clim, nodes, toff):
    """Per entity: the measured 10-20 m mean and the band air temperature over the 120 months
    before each point's first reading month, both averaged over the points with the obs weights."""
    lookup = read_lookup()
    bands = {gid: np.unique(read_snowclim(gid)[:, 0]) for gid in nodes}
    rows, points, skipped = [], {}, []
    for g in glaciers:
        sel = (g.depths >= DEPTH_RANGE[0]) & (g.depths <= DEPTH_RANGE[1])
        if g.glacier_id not in nodes or not sel.any():
            skipped.append((g.glacier_name, g.glacier_id,
                            'no 10-20 m point' if g.glacier_id in nodes else 'not a run glacier'))
            continue
        nd, gid = nodes[g.glacier_id], g.glacier_id
        E, source = model_band(g, lookup, bands[gid])
        w = g.weights[sel]
        pts = []
        oms = [om for om, s in zip(g.obs_months, sel) if s]
        for d, t, om in zip(g.depths[sel], g.T_obs[sel], oms):
            first = int(min(om[0]))
            months = [month_add(first, k) for k in range(-WINDOW_MONTHS, 0)]
            T, mon = node_series(clim, nd, months)
            air = {k: float(v.mean()) for k, v in
                   band_air(T, mon, E, nd['hclim'], np.array(nd['dtdz']), toff[gid]).items()}
            pts.append(dict(depth=float(d), T_obs=float(t), first_month=first,
                            window=[100 * months[0][0] + months[0][1],
                                    100 * months[-1][0] + months[-1][1]], **air))
        T_meas = float(np.average(g.T_obs[sel], weights=w))
        row = dict(glacier_name=g.glacier_name, glacier_id=gid, band_model=E,
                   band_entity=int(g.elevation_band), band_source=source, T_off=toff[gid],
                   node_lat=nd['lat'], node_lon=nd['lon'], hclim=nd['hclim'],
                   n_points=int(sel.sum()), depth_min=float(g.depths[sel].min()),
                   depth_max=float(g.depths[sel].max()),
                   first_month_min=min(p['first_month'] for p in pts),
                   first_month_max=max(p['first_month'] for p in pts),
                   labels='/'.join(sorted(set(map(str, g.firn_label[sel])))), T_meas=T_meas)
        for k in [*VARIANTS, *APPLIED]:
            row[f'T_air_{k}'] = float(np.average([p[k] for p in pts], weights=w))
        for k in [*VARIANTS, *APPLIED]:
            row[f'd_{k}'] = T_meas - row[f'T_air_{k}']
        row.update(sel_a=E >= SEL_A_MIN_BAND, sel_b=T_meas < SEL_B_MAX_T,
                   toff_ge_2=toff[gid] >= TOFF_MIN)
        rows.append(row)
        points[g.glacier_name] = pts
    for name, gid, why in skipped:
        log.info(f'not compared: {name} ({gid}): {why}')
    return pd.DataFrame(rows), points, skipped


def summarise(df):
    groups = {f'(a) band >= {SEL_A_MIN_BAND} m': df.sel_a,
              f'(a) and T_off >= {TOFF_MIN:g} C': df.sel_a & df.toff_ge_2,
              f'(b) 10-20 m mean < {SEL_B_MAX_T:g} C': df.sel_b,
              f'(b) and T_off >= {TOFF_MIN:g} C': df.sel_b & df.toff_ge_2,
              f'T_off >= {TOFF_MIN:g} C, all entities': df.toff_ge_2}
    out = []
    for name, m in groups.items():
        sub = df[m]
        rec = dict(group=name, n_entities=len(sub), n_glaciers=int(sub.glacier_id.nunique()),
                   glacier_ids=sorted(sub.glacier_id.unique()), entities=sub.glacier_name.tolist())
        for k in [*VARIANTS, *APPLIED]:
            d = sub[f'd_{k}'].to_numpy()
            rec[k] = dict(mean=float(d.mean()) if len(d) else None,
                          rms=float(np.sqrt(np.mean(d ** 2))) if len(d) else None)
        out.append(rec)
    return out


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    log.setLevel(logging.INFO)
    fh = logging.FileHandler(OUT_DIR / 'forcing_check.log', mode='w')
    fh.setFormatter(logging.Formatter('%(message)s'))
    sh = logging.StreamHandler()
    sh.setLevel(logging.WARNING)
    log.addHandler(fh)
    log.addHandler(sh)
    log.info(f'forcing_check {datetime.now():%Y-%m-%d %H:%M}; run {RUN}')
    log.info(f'selections (fixed beforehand): depth {DEPTH_RANGE} m, window {WINDOW_MONTHS} months, '
             f'(a) band >= {SEL_A_MIN_BAND} m, (b) mean < {SEL_B_MAX_T} C, T_off >= {TOFF_MIN} C')

    run_ids = sorted(p.stem.split('_')[-1] for p in FIRNICE.glob('firnice_snowclim_*.dat'))
    if CATCHMENT.exists():
        cat = {ln.strip()[9:14] for ln in CATCHMENT.read_text().splitlines()[1:] if ln.strip()}
        if cat != set(run_ids):
            log.warning(f'catchment file differs from the run output ({sorted(cat ^ set(run_ids))}); '
                        f'the run output defines the glacier set')
    clim = read_climate()
    nodes = climate_nodes(run_ids, clim)
    toff = read_toff(run_ids)
    for gid in run_ids:
        nd = nodes[gid]
        log.info(f'{gid}: cell {nd["cell"]} ({nd["n_cell"]} glaciers, mean {nd["mean_lat"]:.4f}N '
                 f'{nd["mean_lon"]:.4f}E) -> node {nd["lat"]:.2f}N {nd["lon"]:.2f}E, hclim '
                 f'{nd["hclim"]:.1f} m, T_off {toff[gid]:.2f} C, dtdz {min(nd["dtdz"]):.5f}..'
                 f'{max(nd["dtdz"]):.5f} K/m')

    valid = validate(clim, nodes, toff)
    log.info(f'VALIDATION: max |rebuilt - logged thermal surface temperature| = '
             f'{valid["max_abs_diff"]:.4f} C over {valid["n_rows"]} band-months of '
             f'{valid["n_glaciers"]} glaciers (tolerance {VALID_TOL} C, file precision 0.01): '
             f'{"PASS" if valid["passed"] else "FAIL"}')

    with open(CACHE, 'rb') as f:
        glaciers = pickle.load(f)['glaciers']
    df, points, skipped = entity_rows(glaciers, clim, nodes, toff)
    summary = summarise(df)
    df.to_csv(OUT_DIR / 'forcing_check.csv', index=False, float_format='%.4f')
    meta = dict(created=f'{datetime.now():%Y-%m-%dT%H:%M:%S}', run=str(RUN), climate=str(CLIM),
                tgrad=str(TGRAD), inventory=str(INV), inventory_hf=str(INV_HF),
                calibration=str(CALI), cache=str(CACHE), run_months='1971-10..2025-09',
                selections=dict(depth_range=DEPTH_RANGE, window_months=WINDOW_MONTHS,
                                sel_a_min_band=SEL_A_MIN_BAND, sel_b_max_T=SEL_B_MAX_T,
                                toff_min=TOFF_MIN),
                variants={**VARIANTS, **APPLIED},
                difference='measured 10-20 m firn temperature minus band air temperature, C')
    entities = df.to_dict(orient='records')
    for e in entities:
        e['points'] = points[e['glacier_name']]
    with open(OUT_DIR / 'forcing_check.json', 'w') as f:
        json.dump(dict(meta=meta, validation=valid, nodes=nodes, T_off=toff, summary=summary,
                       entities=entities, not_compared=skipped), f, indent=1,
                  default=lambda o: o.item() if hasattr(o, 'item') else str(o))

    lines = [f'Validation: max |rebuilt - logged T_surf| {valid["max_abs_diff"]:.4f} C over '
             f'{valid["n_rows"]} band-months, {valid["n_glaciers"]} glaciers '
             f'({"PASS" if valid["passed"] else "FAIL"}, tolerance {VALID_TOL} C)',
             f'{len(df)} entities compared ({len(skipped)} not compared); '
             f'measured 10-20 m minus air temperature, mean / RMS over entities [C]:',
             f'{"group":34s} {"n":>3s} {"gl":>3s}' + ''.join(f'  {k:>13s}' for k in VARIANTS)]
    for rec in summary:
        cells = ''.join(f'  {rec[k]["mean"]:+6.2f} /{rec[k]["rms"]:5.2f}' if rec['n_entities']
                        else f'  {"-":>13s}' for k in VARIANTS)
        lines.append(f'{rec["group"]:34s} {rec["n_entities"]:3d} {rec["n_glaciers"]:3d}{cells}')
    lines.append(f'outputs: {OUT_DIR}/forcing_check.{{json,csv,log}}')
    for ln in lines:
        log.info(ln)
        print(ln)
    if not valid['passed']:
        sys.exit(1)


if __name__ == '__main__':
    main()
