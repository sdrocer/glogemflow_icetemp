"""(d) model values with the depth measured from the ice surface: a reading z m below the ice
surface is read at z + h in the model grid, h = the snow depth the thermal step used (start-of-
month snow of that month, z_mass), on ice bands only (firn bands: snow and firn share one column).
The value at a date is linear in time between the end of the previous month (read at z + h_prev)
and the end of the reading's month (z + h_month), as runner.match_model_to_obs. Variants:
'shallow' shifts only points <= 3 m, 'all' every depth. Every T3 and T3r row of all 32 alps_s2
runs from the archives; h = 0 must reproduce the cached Y3 / Y3r exactly. Saves snow_shift.pkl."""
import os
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import sys, time, pickle, tarfile, tempfile, multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/checks')
import numpy as np, pandas as pd
import c3lib as C
ta = C.ta
from icetemp.calibration.runner import model_at_depth
obs = ta.load_obs(); mod = ta.load_model('alps_s2')
SR = pickle.load(open(C.OUT / 'snow_rows.pkl', 'rb'))['S']          # aligned with obs['T3'] (all 3426 rows)
T3, T3r = obs['T3'], obs['T3r']
ents = obs['glaciers']
TMPD = C.OUT / 'tmp'
# per-month snow depth h per (entity, ym): from the T3 bookkeeping (h_m for the row's month,
# h_p for the previous month); collect a lookup (e, ym) -> h for every month any row touches
H = {}
for r in SR.itertuples():
    H[(r.e, int(r.ym))] = r.h_m
    H[(r.e, ta.prev_month(int(r.ym)))] = r.h_p
# T3r range rows: their months; snow for those months from snowbh again via the T3 rows of the same
# point (every month of a range row is a T3 row of that point), so H covers them


def values(series_at, ym_list, frac_list, row_of):
    out = []
    for ym, f in zip(ym_list, frac_list):
        if ym not in row_of:
            out.append(np.nan); continue
        v = series_at(ym)[row_of[ym]]
        pm = ta.prev_month(ym)
        if pm in row_of:
            v = (1.0 - f) * series_at(pm)[row_of[pm]] + f * v
        out.append(v)
    return np.asarray(out, float)


def work(task):
    i, variant = task
    tgz = ta.archive_path('alps_s2', i)
    gids = {g.glacier_id for g in ents}
    Y3 = np.full(len(T3), np.nan); Y3r = np.full(len(T3r), np.nan)
    TMPD.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=TMPD) as tmp:
        with tarfile.open(tgz) as tar:
            keep = [m for m in tar.getmembers() if (mm := ta.PROFILE_RE.fullmatch(m.name)) and mm.group(1) in gids]
            tar.extractall(tmp, members=keep, filter='data')
        parsed = {}
        for e, g in enumerate(ents):
            path, prof = ta.select_profile(tmp, g, parsed)
            if prof is None or prof[2].empty:
                continue
            zgrid, T, row_of = ta.profile_arrays(prof[2])
            cache = {}

            def at(z):
                key = round(float(z), 9)
                if key not in cache:
                    cache[key] = model_at_depth(zgrid, T, z)
                return cache[key]

            def shift(k, ym):
                z = float(g.depths[k])
                if variant == 'none' or (variant == 'shallow' and z > 3.0):
                    return z
                return z + H.get((e, ym), 0.0)

            for j in np.flatnonzero(T3.e.to_numpy() == e):
                k, ym, f = int(T3.k.iat[j]), int(T3.ym.iat[j]), float(T3.frac.iat[j])
                Y3[j] = values(lambda yy: at(shift(k, yy)), [ym], [f], row_of)[0]
            for j in np.flatnonzero(T3r.e.to_numpy() == e):
                k = int(T3r.k.iat[j])
                ms, fs = [int(x) for x in T3r.months.iat[j]], [float(x) for x in T3r.fracs.iat[j]]
                v = values(lambda yy: at(shift(k, yy)), ms, fs, row_of)
                Y3r[j] = float(v.mean()) if np.isfinite(v).all() else np.nan
    return i, variant, Y3, Y3r


if __name__ == '__main__':
    t0 = time.time()
    tasks = [(i, v) for v in ('none', 'shallow', 'all') for i in range(32)]
    with ProcessPoolExecutor(max_workers=8, mp_context=mp.get_context('fork')) as ex:
        res = list(ex.map(work, tasks))
    out = {}
    for v in ('none', 'shallow', 'all'):
        out[v] = dict(Y3=np.vstack([r[2] for r in sorted(res, key=lambda r: r[0]) if r[1] == v]),
                      Y3r=np.vstack([r[3] for r in sorted(res, key=lambda r: r[0]) if r[1] == v]))
    d3 = np.nanmax(np.abs(out['none']['Y3'] - mod['Y3'])); d3r = np.nanmax(np.abs(out['none']['Y3r'] - mod['Y3r']))
    same_nan = (np.array_equal(np.isnan(out['none']['Y3']), np.isnan(mod['Y3'])),
                np.array_equal(np.isnan(out['none']['Y3r']), np.isnan(mod['Y3r'])))
    print(f'h = 0 reproduces the cache: Y3 max diff {d3:.2e}, Y3r {d3r:.2e}, NaN patterns equal {same_nan} [{time.time()-t0:.0f} s]')
    # T2 from T3 by count-weighted means
    for v in ('none', 'shallow', 'all'):
        out[v]['Y2'] = ta.point_means(T3, out[v]['Y3'])
    d2 = np.nanmax(np.abs(out['none']['Y2'][:, mod['common2']] - mod['Y2'][:, mod['common2']]))
    print(f'T2 pooled from the h = 0 T3 rows vs cached Y2 on the common points: max diff {d2:.2e}')
    for v in ('shallow', 'all'):
        dd = out[v]['Y3'] - mod['Y3']
        ch = np.abs(dd).max(axis=0) > 1e-9
        print(f'{v}: T3 rows changed {int(ch.sum())}; mean change at design 13 over changed rows {np.nanmean(dd[13, ch]):+.3f} C, '
              f'RMS {np.sqrt(np.nanmean(dd[13, ch] ** 2)):.3f}; NaN rows now {int(np.isnan(out[v]["Y3"]).any(axis=0).sum())} '
              f'(was {int(np.isnan(mod["Y3"]).any(axis=0).sum())})')
    pickle.dump(out, open(C.OUT / 'snow_shift.pkl', 'wb'))
