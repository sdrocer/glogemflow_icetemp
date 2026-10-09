"""Step 5a: observation and model rows for the three calibration targets on the 7 Alpine test
glaciers, and the helpers the 5a analyses share. Analysis only: no package code is changed.

Targets (one number compared per ...):
  T1   entity (glacier, 10 m band): data-weighted mean over its points, as calibrator.py pools
  T2   pooled point: an entity's 1 m depth bin (DataHandler.calibration_glaciers)
  T3   pooled point and calendar month (yyyymm): the mean of that month's readings; a profile
       dated only to a range of months is spread evenly over them (data._reading_months), so one
       such reading fills up to 25 rows with the same value
  T3r  T3 with each range-dated reading kept whole: one row per point and month for readings
       dated to a month (as in T3), one row per point and date range for the others, the model
       averaged over the range's months at mid-month (as T2 averages it)

Observations: DataHandler(region='CentralEurope', include_estimated=True, include_unassessed=True).
The raw readings are captured during its own load() (a hook on _resolve_unknown_firn receives each
entity's readings just before they are pooled), so filtering, pooling bins and the months of
range-dated profiles are DataHandler's own; re-pooling them reproduces every entity exactly.

Model: T2 from the re-matched ensembles (rescored/<ens>_step5/design_XXXX.pkl); T3 and T3r from
the archived raw monthly profiles through runner.match_model_to_obs itself, each row passed as a
point of its own (the profile file of runner.parse_training_output; model_at_depth at the
point's depth; linear in time between the end of the previous month and the end of the row's
month at the readings' mean day fraction; the run's first month keeps its month-end value); T1
pooled from T2. Arrays are (design points x rows), NaN where the model does not resolve the depth.
Residual = model - observed.

Build the cache and run the checks (from the repo root; about 3 min for the borehole data, then
about 1 min):
    PYTHONPATH=src python scripts/icetemp_analysis/target_analysis.py [--rebuild] [--workers 6]
Use (with scripts/icetemp_analysis on sys.path):
    import target_analysis as ta
    obs = ta.load_obs()                # tables T1, T2, T3, T3r, readings; glaciers; meta
    mod = ta.load_model('alps_s2')     # design, Y2, Y3, Y3r, T1 arrays, common-row masks
    t2 = ta.target('T2', obs, mod)     # rows, y, Y over the common rows; residual Y - y
    C = ta.target_error_cov(t2, obs)   # D + S (August model error)
    t3 = ta.target('T3r', obs, mod)    # T3 with month-to-month error (Stage 2), date spread in D:
    C = (ta.target_error_cov(t3, obs, seasonal=ta.seasonal_cov)
         + np.diag(ta.date_spread_var(t3.rows, obs, mod['Y3'][13])))
    ta.misfit_ratio(t3.Y[13] - t3.y, C, t3.rows)   # ratio, level and month-to-month parts
Sensitivities of the error model (Step 5a review): a fitted error model,
    em = ta.fit_error_model(t2.Y[13] - t2.y, t2.rows)        # variance, length, offset, nugget (ML)
    C = ta.target_error_cov(t2, obs, **{k: em[k] for k in ('variance', 'length', 'offset', 'nugget')})
and part of the measurement error shared between a point's T3 rows: error_cov(rows, shared_d=0.5).
"""
import argparse
import contextlib
import io
import json
import multiprocessing as mp
import pickle
import re
import sys
import tarfile
import tempfile
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import score_paired_ensembles as spe  # noqa: E402
from icetemp.calibration import data as D  # noqa: E402
from icetemp.calibration.calibrator import (S_MODEL_LENGTH_M, S_MODEL_VARIANCE,  # noqa: E402
                                            _weighted_mean_and_var)
from icetemp.calibration.priors import Priors  # noqa: E402
from icetemp.calibration.runner import (PROFILE_ARCHIVE_DIR, match_model_to_obs,  # noqa: E402
                                        model_at_depth)
from icetemp.io import read_profile_file  # noqa: E402

REPO = HERE.parents[1]
BAYESCAL = REPO / 'data' / 'bayescal'
SCRATCH = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees')
CACHE = SCRATCH / 'step5a' / 'cache'
TMP = SCRATCH / 'step5a' / 'stage1' / 'tmp'          # archives are unpacked here, one design at a time
RESCORED = SCRATCH / 'rescored'                      # <ens>_step5/design_XXXX.pkl
ENSEMBLES = ('alps_s2', 'alps_r7')                   # Step 4 (snow fix) and its reference
DATA_PKL = CACHE / 'datahandler.pkl'
OBS_PKL = CACHE / 'obs_rows.pkl'
MODEL_PKL = {e: CACHE / f'model_rows_{e}.pkl' for e in ENSEMBLES}
CHECKS_JSON = CACHE / 'checks.json'

PARAMS, SHORT = spe.PARAMS, spe.SHORT                # refreeze_frac/f_ref, insul_scale/f_ins, ...
TEST_GLACIERS = spe.TEST_GLACIERS
F_REF_MIN = spe.F_REF_MIN                            # response surfaces leave out f_ref < 0.1
DEPTH_BINS = spe.DEPTH_BINS                          # (lo, hi] m: 0-5, 5-15, 15-30, > 30
BEST_RUN = 13                                        # best run of both ensembles (Step 4)
SEASON = {12: 'DJF', 1: 'DJF', 2: 'DJF', 3: 'MAM', 4: 'MAM', 5: 'MAM', 6: 'JJA', 7: 'JJA',
          8: 'JJA', 9: 'SON', 10: 'SON', 11: 'SON'}
SEASONS = ('DJF', 'MAM', 'JJA', 'SON')
PROFILE_RE = re.compile(r'temp_ID\d+_(\w+)\.dat')
# the scorer's numbers after the tie-break fix (icetemp 436f3ed): design 13 on the 771 points
# finite in all 32 runs of both ensembles
EXPECTED_COMMON = 771
EXPECTED_BEST_RMSE = {'alps_s2': 2.129, 'alps_r7': 2.319}
TOL = 1e-9


# ── small helpers ─────────────────────────────────────────────────────────────────────────────
def _write(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.part')
    with open(tmp, 'wb') as f:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    tmp.replace(path)


def _read(path):
    with open(path, 'rb') as f:
        return pickle.load(f)


def rmse(R, axis=-1):
    """Root mean square along `axis` (NaN-free input)."""
    return np.sqrt(np.mean(np.asarray(R) ** 2, axis=axis))


def prev_month(ym):
    """yyyymm of the month before."""
    return ym - 1 if ym % 100 > 1 else ym - 89


def season_of(ym):
    return SEASON[int(ym) % 100]


def depth_band(depth):
    """Scorer depth band of each depth: '0-5 m', '5-15 m', '15-30 m', '> 30 m' (lo < z <= hi)."""
    z = np.asarray(depth, dtype=float)
    out = np.full(z.shape, '', dtype=object)
    for lo, hi, name in DEPTH_BINS:
        out[(z > lo) & (z <= hi)] = name
    return out


# ── parameters and design ─────────────────────────────────────────────────────────────────────
def prior_sd(campaign=False):
    """Prior sd of (f_ref, f_ins, f_adv): 0.2887, 0.3259, 0.2444. campaign=True: f_ref on the
    campaign's U[0.1, 1] (sd 0.2598)."""
    pri = Priors(bounds_override={'refreeze_frac': (0.1, 1.0)} if campaign else None)
    return np.array([pri.refreeze_frac.std(), pri.insul_scale.std(), pri.advection_scale.std()])


def load_design(ens='alps_s2'):
    """Design of an ensemble: DataFrame indexed by design point, columns f_ref, f_ins, f_adv and
    keep (f_ref >= 0.1, the points response surfaces use). Every worker's design.csv must agree."""
    files = sorted(BAYESCAL.glob(f'{ens}_w*/design.csv'))
    if not files:
        raise FileNotFoundError(f'no design.csv under {BAYESCAL}/{ens}_w*')
    d = pd.read_csv(files[0])[PARAMS]
    for f in files[1:]:
        if not np.array_equal(pd.read_csv(f)[PARAMS].to_numpy(), d.to_numpy()):
            raise RuntimeError(f'{f} differs from {files[0]}')
    d.columns = SHORT
    d.index.name = 'design'
    d['keep'] = d.f_ref >= F_REF_MIN
    return d


def to_sd_units(theta, ref=None, campaign=False):
    """Parameters in prior sd units, relative to `ref` (default 0)."""
    theta = np.atleast_2d(np.asarray(theta, dtype=float))
    ref = np.zeros(theta.shape[1]) if ref is None else np.asarray(ref, dtype=float)
    return (theta - ref) / prior_sd(campaign)[:theta.shape[1]]


# ── borehole data ─────────────────────────────────────────────────────────────────────────────
def _point_index(depths, g):
    """Pooled point (index into g.depths) of each reading: DataHandler's 1 m bins, ordered by mean
    depth as data._pool_by_depth_bin orders them."""
    bins = np.round(depths / D.DEPTH_BIN_WIDTH).astype(int)
    uniq, inv = np.unique(bins, return_inverse=True)
    mean_depth = np.array([depths[inv == j].mean() for j in range(len(uniq))])
    order = np.argsort(mean_depth)
    if not np.array_equal(mean_depth[order], g.depths):
        raise RuntimeError(f'{g.glacier_name}: bins of the captured readings do not give its points')
    k_of_bin = np.empty(len(uniq), dtype=int)
    k_of_bin[order] = np.arange(len(uniq))
    return k_of_bin[inv]


def _reading_months(r, ym_lo, ym_hi):
    """Per reading: (months, weights, day fractions), as DataHandler.load builds them."""
    return [D._reading_months(a, b, f, ym_lo, ym_hi)
            for a, b, f in zip(r['ym0'].to_numpy(), r['ym1'].to_numpy(), r['frac'].to_numpy())]


def _repool_matches(g, r, ym_lo, ym_hi):
    """True if pooling the captured readings again gives entity g exactly."""
    out = D._pool_by_depth_bin(r['depth'].to_numpy(float), r['temperature'].to_numpy(float),
                               r['firn_label'].to_numpy(object), months=_reading_months(r, ym_lo, ym_hi))
    depths, T, is_firn, n, lab, months = out
    same = (np.array_equal(depths, g.depths) and np.array_equal(T, g.T_obs)
            and np.array_equal(is_firn, g.is_firn) and np.array_equal(n, g.n_readings)
            and np.array_equal(lab.astype(str), np.asarray(g.firn_label).astype(str)))
    return same and all(all(np.array_equal(np.asarray(a), np.asarray(b)) for a, b in zip(m1, m2))
                        for m1, m2 in zip(months, g.obs_months))


def _load_with_capture():
    """DataHandler.load() with each entity's readings captured as it pools them. Returns
    (dh, captured [(glacier name, readings DataFrame, resolved firn labels)], seconds, summary)."""
    dh = D.DataHandler(region='CentralEurope', include_estimated=True, include_unassessed=True)
    captured = []
    resolve = dh._resolve_unknown_firn

    def capture(name, labels, gdf):
        out = resolve(name, labels, gdf)
        captured.append((name, gdf.copy(), np.asarray(out, dtype=object).copy()))
        return out

    dh._resolve_unknown_firn = capture          # called once per entity, right before pooling
    t0 = time.time()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', FutureWarning)
        dh.load()
    load_s = time.time() - t0
    del dh._resolve_unknown_firn
    with contextlib.redirect_stdout(io.StringIO()):
        summary = dh.summary()
    if len(captured) != len(dh.calibration_glaciers):
        raise RuntimeError(f'captured {len(captured)} reading sets for {len(dh.calibration_glaciers)} entities')
    return dh, captured, load_s, summary


def _readings_table(glaciers, captured, ym_lo, ym_hi):
    """The captured readings as one table (DataHandler order within each entity), with each
    reading's pooled point k. Returns (table, number of entities that re-pool exactly)."""
    frames, n_repooled = [], 0
    for g, (name, gdf, lab) in zip(glaciers, captured):
        if name != g.base_glacier_name:
            raise RuntimeError(f'capture order: {name} against {g.glacier_name}')
        t0 = pd.to_datetime(gdf['date_meas'], errors='coerce')
        t1 = pd.to_datetime(gdf['date_end'], errors='coerce').fillna(t0)   # as DataHandler.load
        r = pd.DataFrame({
            'entity': g.glacier_name,
            'borehole_id': gdf['borehole_id'].astype(str).to_numpy(object),
            'profile_id': gdf['profile_id'].astype(str).to_numpy(object),
            'bh_label': gdf['label'].to_numpy(object),
            'source_id': gdf['source_id'].to_numpy(object),
            'depth': gdf['depth'].to_numpy(float),
            'temperature': gdf['temperature'].to_numpy(float),
            'equilibrium': gdf['equilibrium'].astype(str).to_numpy(object),
            'date_meas': t0.dt.strftime('%Y-%m-%d').to_numpy(object),
            'date_end': t1.dt.strftime('%Y-%m-%d').to_numpy(object),
            'ym0': gdf['_ym0'].to_numpy(int),
            'ym1': gdf['_ym1'].to_numpy(int),
            'frac': gdf['_frac'].to_numpy(float),
            'mass_balance_area': gdf['mass_balance_area'].to_numpy(object),
            'firn_label': lab,
        })
        r['range_dated'] = r.ym0 != r.ym1
        r['k'] = _point_index(r['depth'].to_numpy(float), g)
        n_repooled += _repool_matches(g, r, ym_lo, ym_hi)
        frames.append(r)
    return pd.concat(frames, ignore_index=True), n_repooled


def load_datahandler(rebuild=False):
    """DataHandler(CentralEurope, estimated and unassessed included), loaded once (about 3 min) and
    cached in CACHE/datahandler.pkl. Every entity's captured readings are re-pooled and compared
    with DataHandler's points (meta n_entities_repooled_exactly). Returns a dict:
      glaciers   all calibration entities (GlacierCalibrationData), DataHandler order
      readings   raw readings of the 7 test glaciers' entities (DataFrame, DataHandler order
                 within each entity): entity, borehole_id, profile_id, bh_label (glenglat
                 borehole label), source_id, depth, temperature, equilibrium (profile),
                 date_meas, date_end (the profile's dates; date_end = date_meas if blank), ym0,
                 ym1 (first and last month of the reading's date range), frac (fraction of the
                 month elapsed at the reading), range_dated (ym0 != ym1), mass_balance_area,
                 firn_label (resolved), k (pooled point)
      systems    entity -> glacier system (data.glacier_systems over all entities: the folds)
      meta       window, counts, summary text, load time"""
    if DATA_PKL.exists() and not rebuild:
        return _read(DATA_PKL)
    dh, captured, load_s, summary = _load_with_capture()
    gl = dh.calibration_glaciers
    readings, n_repooled = _readings_table(gl, captured, dh.ym_lo, dh.ym_hi)
    n_all = len(readings)
    test = {g.glacier_name for g in gl if g.base_glacier_name in TEST_GLACIERS}
    readings = readings[readings.entity.isin(test)].reset_index(drop=True)   # keep the 7 test glaciers

    glenglat = Path(dh.glenglat_dir)
    meta = dict(
        built=datetime.now().isoformat(timespec='seconds'), load_seconds=round(load_s, 1),
        ym_lo=int(dh.ym_lo), ym_hi=int(dh.ym_hi), n_entities=len(gl), n_readings_all=n_all,
        n_readings=len(readings), n_outside_window=int(dh.n_outside_window),
        n_firn_from_wgms=int(dh.n_firn_from_wgms), n_entities_repooled_exactly=int(n_repooled), summary=summary,
        glenglat_mtimes={p.name: p.stat().st_mtime for p in
                         (glenglat / 'borehole.csv', glenglat / 'profile.csv', glenglat / 'measurement.csv')},
    )
    out = dict(glaciers=gl, readings=readings, systems=D.glacier_systems(gl), meta=meta)
    _write(DATA_PKL, out)
    return out


def system_short(label):
    """Short name of a glacier system label (data.glacier_systems)."""
    if 'Grenzgletscher' in label or 'Gornergletscher' in label:
        return 'Monte Rosa'
    return {'Glacier des Bossons': 'Bossons', 'Glacier de Tête Rousse': 'Tête Rousse',
            'Hohsaasgletscher': 'Hohsaas', 'Vadret dal Corvatsch': 'Corvatsch'}.get(label, label)


def test_entities(glaciers, ensembles=ENSEMBLES):
    """Entities of the 7 test glaciers present in the re-matched outputs, in the scorer's order
    (glacier name, then elevation)."""
    known = {g.glacier_name: g for g in glaciers}
    names = set()
    for ens in ensembles:
        for f in sorted((RESCORED / f'{ens}_step5').glob('design_*.pkl')):
            names |= set(_read(f))
    return sorted((known[n] for n in names if n in known and known[n].base_glacier_name in TEST_GLACIERS),
                  key=lambda g: (g.base_glacier_name, g.elevation))


# ── observation rows ──────────────────────────────────────────────────────────────────────────
def _majority_label(lab):
    """Firn label of a set of readings, voted as data._pool_by_depth_bin votes."""
    firn = np.mean(lab != 'ice') >= 0.5
    side = lab[(lab != 'ice') == firn]
    vals, counts = np.unique(side.astype(str), return_counts=True)
    return vals[np.argmax(counts)]


def _join(values):
    return '/'.join(sorted({str(v) for v in values if isinstance(v, str) and v and v != 'nan'}))


def build_obs(data, entities):
    """Observation tables of the test entities. Returns a dict:
      T2   one row per pooled point: t2, e, entity, glacier, rgi, elev, band, system,
           system_short, k, depth, obs, sigma, weight (depth weight), firn_label, is_firn,
           equilibrium, n_readings, n_months, n_seasons, n_years, seasons, ym_first, ym_last,
           range_dated, in_window (all its readings' date ranges inside the run window),
           bh_label, row_label
      T3   one row per (pooled point, month): t3, t2, e, entity, ..., k, depth (the point's,
           where the model is read), depth_readings (mean depth of the month's readings), ym,
           year, month, season, tmon (months since year 0), t (decimal year at the readings'
           mean date), obs (weighted mean of the month's readings), count (the month's weight
           in the point: its readings, a range-dated one split evenly over its months),
           n_readings, n_range (range-dated readings among them), range_only, frac (readings'
           mean fraction of the month elapsed), sigma, weight, firn_label (the point's),
           firn_label_month (the month's readings), equilibrium, bh_label, date_first,
           date_last, row_label
      T3r  as T3, but a range-dated reading is not split: one row per (point, month) for readings
           dated to a month and per (point, date range) for range-dated ones (kind 'month' or
           'range'); months and fracs hold the months the model is averaged over (range: every
           month at 0.5, equal weights); count = readings in the row; obs their mean; frac, tmon
           and t at the row's mean date (for a range its centre); ym_first, ym_last, n_months
      T1   one row per entity: e, entity, glacier, rgi, elev, band, system, system_short,
           equilibrium, n_points, obs_all (weighted mean over all its points; over the common
           rows use t1_pool)
      readings  the raw readings of the test entities (load_datahandler columns plus e, t2)
      glaciers  the test entities (GlacierCalibrationData), the order of e
      checks    largest differences of the T3 months, weights and day fractions from
                DataHandler's obs_months"""
    systems = data['systems']
    ym_lo, ym_hi = data['meta']['ym_lo'], data['meta']['ym_hi']
    R_all = data['readings']
    t1_rows, t2_rows, t3_rows, t3r_rows, rd = [], [], [], [], []
    dmax = dict(month_lists=0, count=0.0, frac=0.0)
    for e, g in enumerate(entities):
        r = R_all[R_all.entity == g.glacier_name].reset_index(drop=True)
        months = _reading_months(r, ym_lo, ym_hi)
        z, T, kk = r.depth.to_numpy(float), r.temperature.to_numpy(float), r.k.to_numpy(int)
        lab, bh = r.firn_label.to_numpy(object), r.bh_label.to_numpy(object)
        d0, d1, rng = r.date_meas.to_numpy(object), r.date_end.to_numpy(object), r.range_dated.to_numpy()
        frac_r = r.frac.to_numpy(float)
        inwin = ((r.ym0 >= ym_lo) & (r.ym1 <= ym_hi)).to_numpy()
        base = dict(e=e, entity=g.glacier_name, glacier=g.base_glacier_name, rgi=g.glacier_id,
                    elev=float(g.elevation), band=int(g.elevation_band), system=systems[g.glacier_name],
                    system_short=system_short(systems[g.glacier_name]), equilibrium=g.equilibrium)
        t1_rows.append(dict(base, n_points=int(g.n_obs),
                            obs_all=float(np.sum(g.weights * g.T_obs) / np.sum(g.weights))))
        t2_start = len(t2_rows)
        rd.append(r.assign(e=e, t2=t2_start + kk))
        for k in range(g.n_obs):
            idx = np.flatnonzero(kk == k)
            acc, groups = {}, {}
            for i in idx:                     # readings in DataHandler order, months in order
                for ym, w, f in zip(*months[i]):
                    a = acc.setdefault(int(ym), [0.0, 0.0, 0.0, 0.0, []])
                    a[0] += w
                    a[1] += w * f
                    a[2] += w * T[i]
                    a[3] += w * z[i]
                    a[4].append(i)
                key = (('range',) + tuple(int(y) for y in months[i][0]) if rng[i]
                       else ('month', int(months[i][0][0])))
                groups.setdefault(key, []).append(i)
            yms = sorted(acc)
            g_yms, g_w, g_f = (np.asarray(x) for x in g.obs_months[k])
            if not np.array_equal(np.asarray(yms), g_yms):
                dmax['month_lists'] += 1
            else:
                dmax['count'] = max(dmax['count'], float(np.max(np.abs(np.array([acc[y][0] for y in yms]) - g_w))))
                dmax['frac'] = max(dmax['frac'], float(np.max(np.abs(np.array([acc[y][1] / acc[y][0] for y in yms])
                                                                     - g_f))))
            t2 = len(t2_rows)
            seasons = sorted({season_of(y) for y in yms}, key=SEASONS.index)
            point = dict(base, t2=t2, k=k, depth=float(g.depths[k]), sigma=float(g.sigma[k]),
                         weight=float(g.weights[k]), firn_label=str(g.firn_label[k]))
            t2_rows.append(dict(point, obs=float(g.T_obs[k]), is_firn=bool(g.is_firn[k]),
                                n_readings=int(g.n_readings[k]), n_months=len(yms), n_seasons=len(seasons),
                                n_years=len({y // 100 for y in yms}), seasons='/'.join(seasons),
                                ym_first=int(yms[0]), ym_last=int(yms[-1]), range_dated=bool(rng[idx].any()),
                                in_window=bool(inwin[idx].all()), bh_label=_join(bh[idx]),
                                row_label=f'{g.glacier_name} {g.depths[k]:.1f} m'))
            for j, ym in enumerate(yms):
                w, wf, wT, wz, ii = acc[ym]
                ii = np.asarray(ii)
                # weights and fractions taken from obs_months itself (equal; checked above)
                cnt = float(g_w[j]) if len(g_w) == len(yms) else w
                frac = float(g_f[j]) if len(g_f) == len(yms) else wf / w
                t3_rows.append(dict(point, depth_readings=wz / w, ym=int(ym), year=int(ym // 100),
                                    month=int(ym % 100), season=season_of(ym),
                                    tmon=int(ym // 100) * 12 + int(ym % 100) - 1,
                                    t=int(ym // 100) + (int(ym % 100) - 1 + frac) / 12.0,
                                    obs=wT / w, count=cnt, n_readings=len(ii), n_range=int(rng[ii].sum()),
                                    range_only=bool(rng[ii].all()), frac=frac,
                                    firn_label_month=str(_majority_label(lab[ii])), bh_label=_join(bh[ii]),
                                    date_first=min(map(str, d0[ii])), date_last=max(map(str, d1[ii])),
                                    row_label=f'{g.glacier_name} {g.depths[k]:.1f} m {ym // 100}-{ym % 100:02d}'))
            for key in sorted(groups, key=lambda q: (q[1], q[0], len(q))):
                ii = np.asarray(groups[key])
                m = np.asarray(key[1:], dtype=int)
                fr = np.array([np.mean(frac_r[ii])]) if key[0] == 'month' else np.full(len(m), 0.5)
                tm = (m // 100) * 12 + m % 100 - 1 + fr        # months since year 0 at each date
                t3r_rows.append(dict(point, kind=key[0], months=m, fracs=fr, ym_first=int(m[0]), ym_last=int(m[-1]),
                                     n_months=len(m), season=season_of(m[0]) if len(m) == 1 else '',
                                     frac=float(fr.mean()), tmon=float(tm.mean()), t=float(tm.mean() / 12.0),
                                     obs=float(T[ii].mean()), count=len(ii), depth_readings=float(z[ii].mean()),
                                     firn_label_month=str(_majority_label(lab[ii])), bh_label=_join(bh[ii]),
                                     date_first=min(map(str, d0[ii])), date_last=max(map(str, d1[ii])),
                                     row_label=f'{g.glacier_name} {g.depths[k]:.1f} m {m[0] // 100}-{m[0] % 100:02d}'
                                               + ('' if len(m) == 1 else f'..{m[-1] // 100}-{m[-1] % 100:02d}')))
    T3 = pd.DataFrame(t3_rows)
    T3.insert(0, 't3', np.arange(len(T3)))
    T3r = pd.DataFrame(t3r_rows)
    T3r.insert(0, 't3r', np.arange(len(T3r)))
    return dict(T1=pd.DataFrame(t1_rows), T2=pd.DataFrame(t2_rows), T3=T3, T3r=T3r,
                readings=pd.concat(rd, ignore_index=True), glaciers=list(entities), checks=dmax,
                meta=dict(data['meta'], built_obs=datetime.now().isoformat(timespec='seconds')))


# ── model rows ────────────────────────────────────────────────────────────────────────────────
def archive_path(ens, i):
    """Archived raw profiles of design point i: net_scratch/icetemp_profiles/<ens>_w<k>/design_XXXX.tar.gz."""
    hits = [p for p in PROFILE_ARCHIVE_DIR.glob(f'{ens}_w*/design_{i:04d}.tar.gz')
            if re.fullmatch(re.escape(ens) + r'_w\d+', p.parent.name)]
    if len(hits) != 1:
        raise FileNotFoundError(f'{ens} design {i}: {len(hits)} archives found')
    return hits[0]


def select_profile(folder, g, parsed=None):
    """(path, (elev, depths, df)) of the profile file runner.parse_training_output compares entity
    g with: the band nearest its elevation; on a band edge its own band; (None, None) if none."""
    parsed = {} if parsed is None else parsed
    own_band = getattr(g, 'elevation_band', g.elevation)
    best_path, best_diff, best = None, (np.inf, np.inf), None
    for cand in sorted(Path(folder).glob(f'temp_ID*_{g.glacier_id}.dat')):
        if cand not in parsed:
            parsed[cand] = read_profile_file(cand)
        elev_m = parsed[cand][0]
        diff = ((abs(elev_m - g.elevation), abs(elev_m - own_band)) if elev_m == elev_m
                else (np.inf, np.inf))
        if diff < best_diff:
            best_diff, best_path, best = diff, cand, parsed[cand]
    return best_path, best


def profile_arrays(df):
    """(node depths, temperatures [months x nodes], calendar yyyymm -> row) of a profile file,
    as runner.match_model_to_obs sets them up (rows are month ends)."""
    depth_cols = [c for c in df.columns if c not in ('year', 'month')]
    zgrid = np.asarray(depth_cols, dtype=float)
    T = df[depth_cols].to_numpy(dtype=float)
    hyd_year = df['year'].to_numpy(dtype=int)
    month = df['month'].to_numpy(dtype=int)
    cal_year = hyd_year - (month >= month[0])
    row_of = {int(y) * 100 + int(m): r for r, (y, m) in enumerate(zip(cal_year, month))}
    return zgrid, T, row_of


def at_date(series, row_of, ym, f):
    """Model series at month ym, fraction f of it elapsed: linear between the end of the previous
    month and the end of ym; month-end value in the run's first month; NaN outside the run.
    (A re-implementation of match_model_to_obs for one month, used as a check.)"""
    if ym not in row_of:
        return np.nan
    v = series[row_of[ym]]
    if prev_month(ym) in row_of:
        v = (1.0 - f) * series[row_of[prev_month(ym)]] + f * v
    return v


def rows_at_dates(df, depths, months, fracs):
    """Model at rows given by depth, months and day fractions (equal weights over a row's months),
    through runner.match_model_to_obs with each row as a point of its own."""
    p = SimpleNamespace(depths=np.asarray(depths, dtype=float),
                        obs_months=[(np.atleast_1d(np.asarray(m, dtype=int)), np.ones(np.size(m)),
                                     np.atleast_1d(np.asarray(f, dtype=float))) for m, f in zip(months, fracs)])
    return match_model_to_obs(df, p)[1][:len(p.depths)]


def _match_design(task):
    """Model rows of one design point from its archived profiles: T2 again (match_model_to_obs),
    the basal row, T3 and T3r (rows_at_dates), and T3 once more through at_date (check)."""
    ens, i, ents, starts, r3, r3r, n2 = task
    tgz = archive_path(ens, i)
    Y2 = np.full(n2, np.nan)
    Y3, Y3_check = np.full(len(r3['e']), np.nan), np.full(len(r3['e']), np.nan)
    Y3r = np.full(len(r3r['e']), np.nan)
    basal, basal_z = np.full(len(ents), np.nan), np.full(len(ents), np.nan)
    used, spans = [], set()
    gids = {g.glacier_id for g in ents}
    TMP.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=TMP) as tmp:
        with tarfile.open(tgz) as tar:
            keep = [m for m in tar.getmembers()
                    if (mm := PROFILE_RE.fullmatch(m.name)) and mm.group(1) in gids]
            tar.extractall(tmp, members=keep, filter='data')
        parsed = {}
        for e, g in enumerate(ents):
            path, prof = select_profile(tmp, g, parsed)
            if prof is None or prof[2].empty:
                used.append((None, np.nan))
                continue
            elev, _, df = prof
            used.append((path.name, float(elev)))
            n = g.n_obs
            d_out, v_out = match_model_to_obs(df, g)
            Y2[starts[e]:starts[e] + n] = v_out[:n]
            if len(v_out) > n:
                basal[e], basal_z[e] = v_out[n], d_out[n]
            zgrid, T, row_of = profile_arrays(df)
            spans.add((min(row_of), max(row_of), len(row_of)))
            j3 = np.flatnonzero(r3['e'] == e)
            if j3.size:
                Y3[j3] = rows_at_dates(df, g.depths[r3['k'][j3]], r3['ym'][j3], r3['frac'][j3])
                for k in np.unique(r3['k'][j3]):
                    series = model_at_depth(zgrid, T, g.depths[k])
                    for j in j3[r3['k'][j3] == k]:
                        Y3_check[j] = at_date(series, row_of, int(r3['ym'][j]), float(r3['frac'][j]))
            jr = np.flatnonzero(r3r['e'] == e)
            if jr.size:
                Y3r[jr] = rows_at_dates(df, g.depths[r3r['k'][jr]], [r3r['months'][j] for j in jr],
                                        [r3r['fracs'][j] for j in jr])
    return dict(design=i, Y2=Y2, Y3=Y3, Y3_check=Y3_check, Y3r=Y3r, basal=basal, basal_depth=basal_z,
                used=used, spans=sorted(spans), archive=str(tgz))


def _rescored_T2(ens, n_des, entities, starts, n2):
    """T2 model values from rescored/<ens>_step5 (design points x points); the depths must be the
    observation points'."""
    Y2 = np.full((n_des, n2), np.nan)
    for i in range(n_des):
        out = _read(RESCORED / f'{ens}_step5' / f'design_{i:04d}.pkl')
        for e, g in enumerate(entities):
            o = out.get(g.glacier_name)
            if o is None:
                continue
            d, T = np.asarray(o[0], float), np.asarray(o[1], float)
            if not (len(T) == len(d) >= g.n_obs and np.allclose(d[:g.n_obs], g.depths, rtol=0, atol=1e-6)):
                raise RuntimeError(f'{ens} design {i}: {g.glacier_name} matched at other depths')
            Y2[i, starts[e]:starts[e] + g.n_obs] = T[:g.n_obs]
    return Y2


def t1_pool(T2, mask=None):
    """T1 pooling over the rows in `mask` (default all): (A, ents) with A (len(ents) x len(T2)),
    each entity's row its depth weights over its masked points normalised to 1 (calibrator.py's
    pooling matrix), ents the entity indices e that keep at least one point. T1 = Y2 @ A.T."""
    mask = np.ones(len(T2), bool) if mask is None else np.asarray(mask, bool)
    e, w = T2['e'].to_numpy(int), T2['weight'].to_numpy(float)
    ents = np.unique(e[mask])
    A = np.zeros((len(ents), len(T2)))
    for a, ent in enumerate(ents):
        m = mask & (e == ent)
        A[a, m] = w[m] / w[m].sum()
    return A, ents


def build_model(ens, obs, workers=6):
    """Model arrays of one ensemble (dict):
      design    load_design(ens); theta (n x 3)
      Y2        model at the T2 rows (design points x rows) from rescored/<ens>_step5, NaN where
                unresolved
      Y3, Y3r   model at the T3 / T3r rows from the archived profiles
      Y2_fresh, Y3_check  T2 re-matched here and T3 through at_date (checks only)
      basal, basal_depth  deepest resolved node per entity (match_model_to_obs's last value)
      common2   T2 rows finite in every design point; common3 / common3r: T3 / T3r rows of those
                points; finite3 / finite3r: rows finite in every design point
      A1, ents1, Y1, obs1  T1 over common2 (t1_pool): Y1 = Y2 @ A1.T, obs1 = A1 @ obs
      band_used per entity: profile file and band elevation (design 0); same in every design
      spans     (first yyyymm, last yyyymm, months) of the profile files read"""
    design = load_design(ens)
    n_des = len(design)
    T2, T3, T3r, ents = obs['T2'], obs['T3'], obs['T3r'], obs['glaciers']
    starts = T2.groupby('e').t2.min().reindex(range(len(ents))).to_numpy(int)
    Y2 = _rescored_T2(ens, n_des, ents, starts, len(T2))
    r3 = dict(e=T3.e.to_numpy(int), k=T3.k.to_numpy(int), ym=T3.ym.to_numpy(int), frac=T3.frac.to_numpy(float))
    r3r = dict(e=T3r.e.to_numpy(int), k=T3r.k.to_numpy(int), months=list(T3r.months), fracs=list(T3r.fracs))
    tasks = [(ens, i, ents, starts, r3, r3r, len(T2)) for i in range(n_des)]
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context('fork')) as ex:
        res = sorted(ex.map(_match_design, tasks), key=lambda r: r['design'])
    stack = lambda key: np.vstack([r[key] for r in res])
    Y3, Y3r = stack('Y3'), stack('Y3r')
    common2 = np.isfinite(Y2).all(axis=0)
    A1, ents1 = t1_pool(T2, common2)
    used0 = res[0]['used']
    band_used = pd.DataFrame({
        'entity': [g.glacier_name for g in ents], 'band': [g.elevation_band for g in ents],
        'file': [u[0] for u in used0], 'band_used': [u[1] for u in used0],
        'same_in_all': [all(r['used'][e] == used0[e] for r in res) for e in range(len(ents))]})
    return dict(
        ens=ens, design=design, theta=design[SHORT].to_numpy(float), Y2=Y2, Y3=Y3, Y3r=Y3r,
        Y2_fresh=stack('Y2'), Y3_check=stack('Y3_check'), basal=stack('basal'), basal_depth=stack('basal_depth'),
        common2=common2, common3=np.isin(T3.t2.to_numpy(), np.flatnonzero(common2)),
        common3r=np.isin(T3r.t2.to_numpy(), np.flatnonzero(common2)),
        finite3=np.isfinite(Y3).all(axis=0), finite3r=np.isfinite(Y3r).all(axis=0),
        A1=A1, ents1=ents1, Y1=Y2[:, common2] @ A1[:, common2].T, obs1=A1 @ T2.obs.to_numpy(float),
        band_used=band_used, spans=sorted({s for r in res for s in r['spans']}),
        archives=[r['archive'] for r in res], source=str(RESCORED / f'{ens}_step5'),
        built=datetime.now().isoformat(timespec='seconds'))


# ── loading the cache ─────────────────────────────────────────────────────────────────────────
def load_obs():
    """Cached observation tables (build_obs); build_cache() first."""
    if not OBS_PKL.exists():
        raise FileNotFoundError(f'{OBS_PKL} missing: run target_analysis.py to build the cache')
    return _read(OBS_PKL)


def load_model(ens='alps_s2'):
    """Cached model arrays of one ensemble (build_model); build_cache() first."""
    if not MODEL_PKL[ens].exists():
        raise FileNotFoundError(f'{MODEL_PKL[ens]} missing: run target_analysis.py to build the cache')
    return _read(MODEL_PKL[ens])


def common_rows(models, obs=None):
    """(common2, common3, common3r): T2 rows finite in every design point of every model in
    `models`, and the T3 and T3r rows of those points."""
    obs = load_obs() if obs is None else obs
    c2 = np.logical_and.reduce([m['common2'] for m in models])
    pts = np.flatnonzero(c2)
    return c2, np.isin(obs['T3'].t2.to_numpy(), pts), np.isin(obs['T3r'].t2.to_numpy(), pts)


def target(name, obs, mod, common=None):
    """One target over the common rows: SimpleNamespace(name, rows (table, index reset), y
    (observations), Y (model, design points x rows), mask (rows of the full table), A (T1 only:
    pooling matrix over the T2 table)). common: (common2, common3, common3r) from common_rows,
    default this ensemble's own. Residual = Y - y."""
    c2, c3, c3r = common if common is not None else (mod['common2'], mod['common3'], mod['common3r'])
    T2 = obs['T2']
    if name == 'T1':
        A, ents = t1_pool(T2, c2)
        return SimpleNamespace(name=name, rows=obs['T1'].iloc[ents].reset_index(drop=True),
                               y=A[:, c2] @ T2.obs.to_numpy(float)[c2], Y=mod['Y2'][:, c2] @ A[:, c2].T,
                               mask=np.isin(np.arange(len(obs['T1'])), ents), A=A)
    table, Y, mask = {'T2': (T2, mod['Y2'], c2), 'T3': (obs['T3'], mod['Y3'], c3),
                      'T3r': (obs['T3r'], mod['Y3r'], c3r)}[name]
    return SimpleNamespace(name=name, rows=table[mask].reset_index(drop=True), y=table.obs.to_numpy(float)[mask],
                           Y=Y[:, mask], mask=mask, A=None)


# ── groups, folds and the error model ─────────────────────────────────────────────────────────
def group_masks(rows):
    """Borehole groups of score_paired_ensembles.borehole_groups: [(name, mask)], 'all common
    points' first. Rows need rgi, band, glacier and (for the Colle Gnifetti depth split) depth."""
    rows = rows if 'depth' in rows else rows.assign(depth=np.nan)
    return spe.borehole_groups(rows.reset_index(drop=True))


def group_of(rows):
    """One borehole group name per row (the scorer's main groups; 'other' if none)."""
    out = np.full(len(rows), 'other', dtype=object)
    for name, m in group_masks(rows)[1:]:
        if not name.startswith('  ') and name != 'other':
            out[m] = name
    return out


def fold_labels(rows):
    """(fold index per row, fold names): one fold per glacier system (6 on the 7 test glaciers)."""
    names = sorted(rows['system'].unique())
    return rows['system'].map({n: i for i, n in enumerate(names)}).to_numpy(int), names


def obs_var(rows):
    """Observation error variance per row: glenglat sigma^2 (x3 sigma for 'estimated' profiles).
    DataHandler gives every point of an entity the same sigma (one value per borehole entry), and
    each row gets the full sigma^2, independent between rows; a point's T3 rows therefore imply a
    smaller measurement error for the point's mean than its single T2 row does."""
    return rows['sigma'].to_numpy(float) ** 2


def obs_cov(rows, shared=0.0):
    """Observation error D: diag(obs_var) (default). shared (0-1): the fraction of sigma^2 that is
    common to all rows of one pooled point (a systematic error per borehole entry), the rest
    independent per row; shared=1 gives a point's mean the measurement error of its T2 row.
    A sensitivity only."""
    s = rows['sigma'].to_numpy(float)
    if not shared:
        return np.diag(s ** 2)
    t2 = rows['t2'].to_numpy(int)
    return (1.0 - shared) * np.diag(s ** 2) + shared * np.where(t2[:, None] == t2[None, :], np.outer(s, s), 0.0)


def s_model_cov(rows, variance=S_MODEL_VARIANCE, length=S_MODEL_LENGTH_M, seasonal=None, offset=0.0, nugget=0.0):
    """Model-error covariance S between T2, T3 or T3r rows (calibrator.py's August term):
    variance * exp(-|dz| / length) between rows of one entity, 0 across entities. Rows of one
    point in different months share its depth, so this term correlates them fully.
    seasonal: optional callable(rows) -> (n x n) array, added within each entity (the seasonal
    term between months, measured in Stage 2).
    offset: variance of a misfit shared by all rows of an entity (an entity offset, as a random
    effect); nugget: variance shared by the rows of one pooled point only. Both 0 in the August
    model; they enter the fitted error model (fit_error_model), a sensitivity."""
    z, e = rows['depth'].to_numpy(float), rows['e'].to_numpy(int)
    same = e[:, None] == e[None, :]
    S = np.where(same, variance * np.exp(-np.abs(z[:, None] - z[None, :]) / length) + offset, 0.0)
    if nugget:
        t2 = rows['t2'].to_numpy(int)
        S = S + np.where(t2[:, None] == t2[None, :], nugget, 0.0)
    if seasonal is not None:
        S = S + np.where(same, seasonal(rows), 0.0)
    return S


def error_cov(rows, variance=S_MODEL_VARIANCE, length=S_MODEL_LENGTH_M, seasonal=None, offset=0.0, nugget=0.0,
              shared_d=0.0):
    """D + S: observation error plus model error (no emulator term, no per-glacier offset unless
    `offset` is given)."""
    return obs_cov(rows, shared_d) + s_model_cov(rows, variance, length, seasonal, offset, nugget)


def target_error_cov(tgt, obs, variance=S_MODEL_VARIANCE, length=S_MODEL_LENGTH_M, seasonal=None, offset=0.0,
                     nugget=0.0, shared_d=0.0):
    """D + S of a target (target()): error_cov of its rows; for T1 pooled from the T2 rows,
    A (D + S) A^T, as calibrator.py pools it."""
    if tgt.name == 'T1':
        return tgt.A @ error_cov(obs['T2'], variance, length, offset=offset, nugget=nugget) @ tgt.A.T
    return error_cov(tgt.rows, variance, length, seasonal, offset, nugget, shared_d)


# ── a fitted error model (sensitivity) ────────────────────────────────────────────────────────
ERROR_MODEL_AUGUST = dict(variance=S_MODEL_VARIANCE, length=S_MODEL_LENGTH_M, offset=0.0, nugget=0.0)


def error_model_loglik(r, rows, variance, length, offset=0.0, nugget=0.0, mask=None):
    """Gaussian log-likelihood of T2 residuals r under D + S(variance, length, offset, nugget),
    block-diagonal by entity; mask: rows to use (default all)."""
    mask = np.ones(len(r), bool) if mask is None else np.asarray(mask, bool)
    e = rows['e'].to_numpy(int)
    z, s2 = rows['depth'].to_numpy(float), rows['sigma'].to_numpy(float) ** 2
    ll = 0.0
    for ent in np.unique(e[mask]):
        j = np.flatnonzero(mask & (e == ent))
        C = (np.diag(s2[j] + nugget) + variance * np.exp(-np.abs(z[j, None] - z[None, j]) / length) + offset)
        L = np.linalg.cholesky(C)
        w = np.linalg.solve(L, r[j])
        ll += -0.5 * (w @ w) - np.log(np.diag(L)).sum() - 0.5 * len(j) * np.log(2 * np.pi)
    return float(ll)


def fit_error_model(r, rows, mask=None, starts=((3.5, 15.0, 1.0, 0.1), (1.5, 30.0, 2.0, 0.3), (2.5, 60.0, 3.0, 0.05))):
    """Maximum-likelihood (variance, length, offset, nugget) of D + S for T2 residuals r (one run;
    optionally only the rows in mask), D fixed at glenglat's sigma^2. Variance 0.01-50 C^2, length
    1-500 m, offset and nugget 1e-4-50 C^2; L-BFGS-B on their logarithms from several starts.
    Returns dict(variance, length, offset, nugget, loglik, loglik_august)."""
    from scipy.optimize import minimize
    lo, hi = np.log([1e-2, 1.0, 1e-4, 1e-4]), np.log([50.0, 500.0, 50.0, 50.0])
    f = lambda x: -error_model_loglik(r, rows, *np.exp(np.clip(x, lo, hi)), mask=mask)
    best = None
    for s in starts:
        res = minimize(f, np.log(s), method='L-BFGS-B', bounds=list(zip(lo, hi)))
        if best is None or res.fun < best.fun:
            best = res
    v, l, o, n = np.exp(np.clip(best.x, lo, hi))
    return dict(variance=float(v), length=float(l), offset=float(o), nugget=float(n), loglik=float(-best.fun),
                loglik_august=error_model_loglik(r, rows, S_MODEL_VARIANCE, S_MODEL_LENGTH_M, mask=mask))


# ── month-to-month model error for T3 rows (Stage 2) ──────────────────────────────────────────
# Fitted to the within-point semivariogram of the T3 month rows of alps_s2 (pair likelihood; median
# over the 29 runs with f_ref >= 0.1; the six Colle Gnifetti melt-spike rows left out).
SEASONAL = dict(
    a=3.0, zs=3.5,                        # within one borehole: variance a exp(-z / zs) [C^2]
    ice=0.25, high=5.4, high_band=4000.0,  # variance factors: ice points; firn of entities at band >= 4000 m
    w=0.27, tau=7.0,                      # lag correlation w cos(2 pi dt / 12) + (1 - w) exp(-dt / tau) [months]
    u0=0.45,                              # rows of different boreholes of one point [C^2]
    beta=0.035, ld=50.0,                  # slow drift [C / yr] and its depth-correlation length [m]
    lz=(2.5, 9.0, 20.0), lens=(1.2, 4.3, 200.0),  # depth-correlation length [m] at these depths
)


def _row_times(rows):
    """Per row: (times [months since year 0, at the readings' date], weights summing to 1).
    T3: its month; T3r: its months (equal weights); other rows: one time 0."""
    if 'kind' in rows:
        out = []
        for m, f in zip(rows['months'], rows['fracs']):
            m = np.asarray(m, dtype=int)
            out.append(((m // 100) * 12 + m % 100 - 1 + np.asarray(f, dtype=float), np.full(len(m), 1.0 / len(m))))
        return out
    if 'tmon' in rows and 'frac' in rows:
        t = rows['tmon'].to_numpy(float) + rows['frac'].to_numpy(float)
        return [(np.array([x]), np.ones(1)) for x in t]
    return [(np.zeros(1), np.ones(1)) for _ in range(len(rows))]


def _depth_corr(z, lz, lens):
    """Matern-3/2 correlation in depth with a length that varies with depth (log-linear between
    the knots lz, constant beyond); Paciorek & Schervish (2006), positive definite for any profile."""
    lv = np.exp(np.interp(z, lz, np.log(lens)))
    L2 = 0.5 * (lv[:, None] ** 2 + lv[None, :] ** 2)
    s = np.sqrt(3.0) * np.abs(z[:, None] - z[None, :]) / np.sqrt(L2)
    return np.sqrt(np.outer(lv, lv) / L2) * (1.0 + s) * np.exp(-s)


def seasonal_cov(rows, params=None):
    """Month-to-month model error between T3 or T3r rows: the `seasonal` term of s_model_cov.
    Zero across entities and for points with a single row. Three parts, smooth in depth
    (_depth_corr), each centred on the point's count-weighted mean, so that a point's pooled value
    (T2) and the entity means (T1) are unchanged and S keeps carrying the level:
      within one borehole  v(z) [w cos(2 pi dt / 12) + (1 - w) exp(-dt / tau)], v = a exp(-z / zs),
                           x ice for ice points, x high for firn of entities at band >= high_band;
                           a row spanning several months (T3r range rows) averages over them
      between boreholes    u0 for rows of different boreholes (labels split at '/')
      slow drift           (beta dt)^2 / 2 between rows dt years apart: a misfit trend per entity
    params: overrides of SEASONAL."""
    p = dict(SEASONAL, **(params or {}))
    e, t2 = rows['e'].to_numpy(int), rows['t2'].to_numpy(int)
    z = rows['depth'].to_numpy(float)
    cnt = rows['count'].to_numpy(float) if 'count' in rows else np.ones(len(rows))
    ice = rows['firn_label'].to_numpy(object) == 'ice'
    high = ~ice & (rows['band'].to_numpy(float) >= p['high_band'])
    v = p['a'] * np.exp(-z / p['zs']) * np.where(ice, p['ice'], np.where(high, p['high'], 1.0))
    times = _row_times(rows)
    labels = [str(b).split('/') for b in rows['bh_label']]
    K = np.zeros((len(rows), len(rows)))
    for ent in np.unique(e):
        idx = np.flatnonzero(e == ent)
        pts = t2[idx]
        if len(np.unique(pts)) == len(idx):
            continue
        T = np.concatenate([times[i][0] for i in idx])
        Wm = np.zeros((len(idx), len(T)))
        k = 0
        for a, i in enumerate(idx):
            Wm[a, k:k + len(times[i][0])] = times[i][1]
            k += len(times[i][0])
        dt = np.abs(T[:, None] - T[None, :])
        rho = Wm @ (p['w'] * np.cos(2 * np.pi * dt / 12.0) + (1 - p['w']) * np.exp(-dt / p['tau'])) @ Wm.T
        zi = z[idx]
        corr = _depth_corr(zi, p['lz'], p['lens'])
        names = sorted({x for i in idx for x in labels[i]})
        B = np.zeros((len(idx), len(names)))
        for a, i in enumerate(idx):
            for x in labels[i]:
                B[a, names.index(x)] += 1.0 / len(labels[i])
        tm = Wm @ T / 12.0
        tm = tm - tm.mean()
        Kb = (np.sqrt(np.outer(v[idx], v[idx])) * rho + p['u0'] * (B @ B.T)) * corr \
            + p['beta'] ** 2 * _depth_corr(zi, (0.0,), (p['ld'],)) * np.outer(tm, tm)
        C = np.eye(len(idx))
        for pt in np.unique(pts):
            j = np.flatnonzero(pts == pt)
            C[np.ix_(j, j)] -= cnt[idx][j] / cnt[idx][j].sum()
        K[np.ix_(idx, idx)] = C @ Kb @ C.T
    return K


def date_spread_var(rows, obs, y3):
    """Per T3r row: variance of the model over the months a range-dated reading may have been taken
    in [C^2] (y3: one design's model at every row of obs['T3']); 0 for rows dated to a month.
    An addition to D for T3r."""
    T3 = obs['T3']
    lut = dict(zip(zip(T3.t2.to_numpy(int), T3.ym.to_numpy(int)), np.asarray(y3, dtype=float)))
    return np.array([np.var([lut[(t2, int(m))] for m in months]) if kind == 'range' else 0.0
                     for kind, t2, months in zip(rows['kind'], rows['t2'], rows['months'])])


def misfit_ratio(r, Sigma, rows):
    """r^T Sigma^-1 r / n (Sigma block-diagonal by entity, as the error model is) and its split into
    a level part (GLS point means: one dimension per pooled point) and a month-to-month part (the
    remaining n - points dimensions), each divided by its number of dimensions; ~1 if honest."""
    e, t2 = rows['e'].to_numpy(int), rows['t2'].to_numpy(int)
    tot = lvl = 0.0
    n_pts = 0
    for ent in np.unique(e):
        idx = np.flatnonzero(e == ent)
        L = np.linalg.cholesky(Sigma[np.ix_(idx, idx)])
        pts, inv = np.unique(t2[idx], return_inverse=True)
        Z = np.zeros((len(idx), len(pts)))
        Z[np.arange(len(idx)), inv] = 1.0
        wr = np.linalg.solve(L, r[idx])
        wZ = np.linalg.solve(L, Z)
        u = wZ.T @ wr
        tot += wr @ wr
        lvl += u @ np.linalg.solve(wZ.T @ wZ, u)
        n_pts += len(pts)
    n = len(r)
    return dict(n=n, points=n_pts, ratio=tot / n, level=lvl / n_pts,
                month_to_month=(tot - lvl) / (n - n_pts) if n > n_pts else np.nan)


# ── checks and counts ─────────────────────────────────────────────────────────────────────────
def point_means(rows, values):
    """Count-weighted mean over each point's rows (T3 or T3r; last axis = rows), NaN if any row is.
    Returns (n_sets x T2 points)."""
    t2, w = rows.t2.to_numpy(int), rows['count'].to_numpy(float)
    n2 = int(t2.max()) + 1
    den = np.bincount(t2, weights=w, minlength=n2)
    return np.vstack([np.bincount(t2, weights=w * v, minlength=n2) for v in np.atleast_2d(values)]) / den


def _t3_obs_groupby(obs):
    """T3 observations recomputed from the readings with a pandas group-by: (t2, ym) -> weighted
    mean, against the T3 table. Returns (rows missing either side, largest difference)."""
    lo, hi = obs['meta']['ym_lo'], obs['meta']['ym_hi']
    R = obs['readings']
    parts = []
    for t2, a, b, f, T in zip(R.t2.to_numpy(), R.ym0.to_numpy(), R.ym1.to_numpy(), R.frac.to_numpy(),
                              R.temperature.to_numpy()):
        yms, ws, _ = D._reading_months(a, b, f, lo, hi)
        parts += [(t2, ym, w, w * T) for ym, w in zip(yms, ws)]
    E = pd.DataFrame(parts, columns=['t2', 'ym', 'w', 'wT']).groupby(['t2', 'ym'])[['w', 'wT']].sum()
    J = obs['T3'].set_index(['t2', 'ym']).join(E, how='outer')
    missing = int(J.obs.isna().sum() + J.w.isna().sum())
    return missing, float(np.nanmax(np.abs(J.wT / J.w - J.obs)))


def run_checks(obs, models):
    """Consistency checks of the cache; returns a dict (also written to CACHE/checks.json)."""
    T1, T2, T3, T3r = obs['T1'], obs['T2'], obs['T3'], obs['T3r']
    missing, d_gb = _t3_obs_groupby(obs)
    out = dict(datahandler=dict(
        entities=obs['meta']['n_entities'], repooled_exactly=obs['meta']['n_entities_repooled_exactly'],
        test_entities=len(obs['glaciers']), test_readings=len(obs['readings']),
        t3_month_lists_differing=obs['checks']['month_lists'], t3_count_max_diff=obs['checks']['count'],
        t3_frac_max_diff=obs['checks']['frac'], t3_obs_groupby_missing_rows=missing, t3_obs_groupby_max_diff=d_gb))
    inwin = T2.in_window.to_numpy()
    for name, rows in (('T3', T3), ('T3r', T3r)):
        d = np.abs(point_means(rows, rows.obs.to_numpy())[0] - T2.obs.to_numpy())
        out[f'T2_obs_vs_{name}'] = dict(points=len(T2), in_window=int(inwin.sum()),
                                        max_diff_in_window=float(d[inwin].max()), max_diff_all=float(d.max()),
                                        passes=bool(d[inwin].max() <= TOL))
    for ens, m in models.items():
        Y2 = m['Y2']
        fin2 = np.isfinite(Y2)
        for name, rows, Y in (('T3', T3, m['Y3']), ('T3r', T3r, m['Y3r'])):
            M = point_means(rows, Y)
            finM = np.isfinite(M)
            both = fin2 & finM & inwin[None, :]
            dm = np.abs(Y2 - M)[both]
            out[f'{ens}_T2_model_vs_{name}'] = dict(
                values_compared=int(both.sum()), max_diff=float(dm.max()), passes=bool(dm.max() <= TOL),
                t2_finite_some_rows_nan=int((fin2 & ~finM).sum()), t2_nan_all_rows_finite=int((~fin2 & finM).sum()))
        same_nan = np.array_equal(np.isnan(m['Y3']), np.isnan(m['Y3_check']))
        d_chk = np.abs(m['Y3'] - m['Y3_check'])[np.isfinite(m['Y3'])]
        out[f'{ens}_T3_matcher_vs_at_date'] = dict(nan_pattern_equal=bool(same_nan), max_diff=float(d_chk.max()),
                                                   passes=bool(same_nan and d_chk.max() <= TOL))
        same_fresh = np.array_equal(np.isnan(Y2), np.isnan(m['Y2_fresh']))
        d_fresh = np.abs(Y2 - m['Y2_fresh'])[np.isfinite(Y2) & np.isfinite(m['Y2_fresh'])]
        out[f'{ens}_T2_rescored_vs_archives'] = dict(nan_pattern_equal=bool(same_fresh),
                                                     max_diff=float(d_fresh.max()),
                                                     passes=bool(same_fresh and d_fresh.max() <= TOL))
        # T1 two ways: the pooling matrix, and calibrator._weighted_mean_and_var entity by entity
        c2 = m['common2']
        e_of, w_of, y_of = T2.e.to_numpy(int), T2.weight.to_numpy(float), T2.obs.to_numpy(float)
        loop_obs, loop_mod = [], []
        for ent in m['ents1']:
            sel = c2 & (e_of == ent)
            zero = np.zeros(sel.sum())
            loop_obs.append(_weighted_mean_and_var(y_of[sel], w_of[sel], zero)[0])
            loop_mod.append([_weighted_mean_and_var(Y2[i, sel], w_of[sel], zero)[0] for i in range(len(Y2))])
        d1o = float(np.abs(np.asarray(loop_obs) - m['obs1']).max())
        d1m = float(np.abs(np.asarray(loop_mod).T - m['Y1']).max())
        t1 = target('T1', obs, m)
        d1t = float(max(np.abs(t1.y - m['obs1']).max(), np.abs(t1.Y - m['Y1']).max()))
        out[f'{ens}_T1_vs_T2'] = dict(entities=len(m['ents1']), entities_all=len(T1), max_diff_obs=d1o,
                                      max_diff_model=d1m, max_diff_target_accessor=d1t,
                                      passes=bool(max(d1o, d1m, d1t) <= 1e-12))
        out[f'{ens}_files'] = dict(entities_same_file_in_all_designs=int(m['band_used'].same_in_all.sum()),
                                   off_band=[f'{r.entity} -> {r.band_used:.0f}' for r in m['band_used'].itertuples()
                                             if r.band_used != r.band],
                                   spans=[list(s) for s in m['spans']])
        out[f'{ens}_common'] = dict(common2=int(c2.sum()), common3=int(m['common3'].sum()),
                                    common3r=int(m['common3r'].sum()), finite3=int(m['finite3'].sum()),
                                    finite3r=int(m['finite3r'].sum()),
                                    finite3_outside_common_points=int((m['finite3'] & ~m['common3']).sum()),
                                    common3_all_finite=bool(np.isfinite(m['Y3'][:, m['common3']]).all()),
                                    common3r_all_finite=bool(np.isfinite(m['Y3r'][:, m['common3r']]).all()))
    c2, c3, c3r = common_rows(list(models.values()), obs)
    y = T2.obs.to_numpy()
    rm = {}
    for ens, m in models.items():
        r_joint = rmse(m['Y2'][BEST_RUN, c2] - y[c2])
        r_own = rmse(m['Y2'][BEST_RUN, m['common2']] - y[m['common2']])
        rm[ens] = dict(joint_rows=int(c2.sum()), rmse_joint=round(float(r_joint), 4),
                       own_rows=int(m['common2'].sum()), rmse_own=round(float(r_own), 4),
                       expected=EXPECTED_BEST_RMSE[ens],
                       best_run_joint=int(np.argmin(rmse(m['Y2'][:, c2] - y[c2][None, :], axis=1))),
                       passes=bool(c2.sum() == EXPECTED_COMMON and abs(r_joint - EXPECTED_BEST_RMSE[ens]) < 5e-4))
    out['T2_rmse_design13'] = rm
    # for reference: RMSE (unweighted over each target's rows) at design 13 and the best run
    out['rmse_by_target'] = {}
    for ens, m in models.items():
        for name in ('T1', 'T2', 'T3', 'T3r'):
            t = target(name, obs, m, (c2, c3, c3r))
            r = rmse(t.Y - t.y[None, :], axis=1)
            out['rmse_by_target'][f'{ens} {name}'] = dict(rows=len(t.y), design13=round(float(r[BEST_RUN]), 4),
                                                          best_run=int(np.argmin(r)), best=round(float(r.min()), 4))
    out['counts'] = counts(obs, c2, c3, c3r)
    dh = out['datahandler']
    out['all_pass'] = bool(
        dh['repooled_exactly'] == dh['entities'] and dh['t3_month_lists_differing'] == 0
        and max(dh['t3_count_max_diff'], dh['t3_frac_max_diff'], dh['t3_obs_groupby_max_diff']) <= TOL
        and dh['t3_obs_groupby_missing_rows'] == 0
        and out['T2_obs_vs_T3']['passes'] and out['T2_obs_vs_T3r']['passes']
        and all(out[f'{e}_{k}']['passes'] for e in models
                for k in ('T2_model_vs_T3', 'T2_model_vs_T3r', 'T3_matcher_vs_at_date', 'T2_rescored_vs_archives',
                          'T1_vs_T2'))
        and all(v['passes'] for v in rm.values()))
    return out


def counts(obs, common2, common3, common3r):
    """Rows per target, depth band, season and firn label, over the common rows."""
    T1, T2, T3, T3r = obs['T1'], obs['T2'], obs['T3'], obs['T3r']
    A1, ents1 = t1_pool(T2, common2)
    out = dict(T1=dict(all=len(T1), common=len(ents1)), T2=dict(all=len(T2), common=int(common2.sum())),
               T3=dict(all=len(T3), common=int(common3.sum())), T3r=dict(all=len(T3r), common=int(common3r.sum())))
    t2c, t3c, t3rc = T2[common2], T3[common3], T3r[common3r]
    bands = [b for _, _, b in DEPTH_BINS]
    out['T2_by_depth'] = {b: int((depth_band(t2c.depth) == b).sum()) for b in bands}
    out['T2_by_firn_label'] = {k: int(v) for k, v in t2c.firn_label.value_counts().items()}
    out['T2_mixing'] = dict(several_months=int((t2c.n_months > 1).sum()),
                            several_seasons=int((t2c.n_seasons > 1).sum()),
                            several_years=int((t2c.n_years > 1).sum()), range_dated=int(t2c.range_dated.sum()),
                            not_in_window=int((~t2c.in_window).sum()))
    out['T2_single_season'] = {s: int(((t2c.n_seasons == 1) & (t2c.seasons == s)).sum()) for s in SEASONS}
    zb, se = depth_band(t3c.depth), t3c.season.to_numpy()
    out['T3_by_depth_season'] = {b: {s: int(((zb == b) & (se == s)).sum()) for s in SEASONS} for b in bands}
    out['T3_by_firn_label'] = {k: int(v) for k, v in t3c.firn_label.value_counts().items()}
    per_point = t3c.groupby('t2').size()
    out['T3_rows_per_point'] = {lab: int(((per_point >= lo) & (per_point <= hi)).sum())
                                for lab, lo, hi in (('1', 1, 1), ('2', 2, 2), ('3-5', 3, 5), ('6-12', 6, 12),
                                                    ('13-60', 13, 60), ('> 60', 61, 10 ** 9))}
    out['T3_range_only_rows'] = int(t3c.range_only.sum())
    out['T3_rows_mixing_range_and_dated'] = int(((t3c.n_range > 0) & ~t3c.range_only).sum())
    mo = t3rc[t3rc.kind == 'month']
    zr, sr = depth_band(mo.depth), mo.season.to_numpy()
    out['T3r_by_kind'] = {k: int(v) for k, v in t3rc.kind.value_counts().items()}
    out['T3r_month_rows_by_depth_season'] = {b: {s: int(((zr == b) & (sr == s)).sum()) for s in SEASONS}
                                             for b in bands}
    zrr = depth_band(t3rc[t3rc.kind == 'range'].depth)
    out['T3r_range_rows_by_depth'] = {b: int((zrr == b).sum()) for b in bands}
    out['by_group'] = {name: dict(T1=int(np.isin(ents1, T2.e[m & common2]).sum()), T2=int((m & common2).sum()),
                                  T3=int(np.isin(t3c.t2, np.flatnonzero(m & common2)).sum()),
                                  T3r=int(np.isin(t3rc.t2, np.flatnonzero(m & common2)).sum()))
                       for name, m in group_masks(T2)[1:] if not name.startswith('  ')}
    out['by_system'] = {s: dict(T1=int(np.isin(ents1, T2.e[(T2.system_short == s).to_numpy() & common2]).sum()),
                                T2=int((t2c.system_short == s).sum()), T3=int((t3c.system_short == s).sum()),
                                T3r=int((t3rc.system_short == s).sum()))
                        for s in sorted(T2.system_short.unique())}
    return out


# ── build ─────────────────────────────────────────────────────────────────────────────────────
def build_cache(rebuild=False, workers=6):
    """Build (or reuse) every cache file and run the checks. Returns (obs, models, checks)."""
    t0 = time.time()
    data = load_datahandler(rebuild=rebuild)
    print(f'DataHandler: {data["meta"]["n_entities"]} entities, {data["meta"]["n_readings"]} readings of the '
          f'test glaciers ({time.time() - t0:.0f} s)', flush=True)
    if OBS_PKL.exists() and not rebuild:
        obs = load_obs()
    else:
        obs = build_obs(data, test_entities(data['glaciers']))
        _write(OBS_PKL, obs)
    print(f'observation rows: T1 {len(obs["T1"])}, T2 {len(obs["T2"])}, T3 {len(obs["T3"])}, '
          f'T3r {len(obs["T3r"])} ({time.time() - t0:.0f} s)', flush=True)
    models = {}
    for ens in ENSEMBLES:
        if MODEL_PKL[ens].exists() and not rebuild:
            models[ens] = load_model(ens)
        else:
            models[ens] = build_model(ens, obs, workers=workers)
            _write(MODEL_PKL[ens], models[ens])
        print(f'{ens}: model rows T2 {models[ens]["Y2"].shape}, T3 {models[ens]["Y3"].shape}, '
              f'T3r {models[ens]["Y3r"].shape} ({time.time() - t0:.0f} s)', flush=True)
    checks = run_checks(obs, models)
    CHECKS_JSON.write_text(json.dumps(checks, indent=1, default=str), encoding='utf-8')
    return obs, models, checks


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--rebuild', action='store_true', help='rebuild every cache file')
    ap.add_argument('--workers', type=int, default=6, help='processes reading the archives')
    args = ap.parse_args()
    _, _, checks = build_cache(rebuild=args.rebuild, workers=args.workers)
    print(json.dumps(checks, indent=1, default=str))
    print(f'cache: {CACHE}; all checks pass: {checks["all_pass"]}')


if __name__ == '__main__':
    main()
