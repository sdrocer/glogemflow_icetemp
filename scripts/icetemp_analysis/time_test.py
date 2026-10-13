"""Time test: measured vs modelled temperature change at repeat borehole sites.

Phase 'pairs' builds the repeat-site pairs from glenglat alone (no model values) and writes
pairs.csv. Phase 'model' reads pairs.csv, computes the measured change at 15/20/50 m and the
modelled change at one band (the earlier hole's) over the same months, for four runs, and writes
results.csv and summary.md.

Usage: PYTHONPATH=src python scripts/icetemp_analysis/time_test.py {pairs|model} [out_dir]
"""
import itertools
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
GLENGLAT = REPO / 'glenglat/data'
EXTRA_SOURCES = ['carturan2023-ortles']     # glenglat subfolder outside the specified data
GLIMS_LOOKUP = REPO / 'data/rgi7_centraleurope_glims_lookup.csv'
RUNS = ['ce25_pred_w0', 'ce25_F_w0', 'ce25_H_w0', 'ce25_FH_w0']
FIRNICE = 'training_runs/monthly/CentralEurope/PAST/firnice_temperature'
OUT = Path(sys.argv[2] if len(sys.argv) > 2
           else '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week1/time_test')

# Pre-registered rule
DEPTHS = [15.0, 20.0, 50.0]
MIN_YEARS = 3                    # calendar years between visit midpoints
MAX_DIST = 50.0                  # m, rule (b)
MAX_DZ = 20.0                    # m, rule (b)
EQ_OK = {'true', 'estimated', ''}
BOX = dict(lat=(43.0, 48.5), lon=(5.0, 17.0))
WINDOW = (197110, 202509)        # model run: hydrological years 1972-2025
NOISE_DIST, NOISE_DZ = 300.0, 40.0   # same-campaign hole pairs for the site noise
MAX_BAND_GAP = 50                # m, hole band to nearest modelled band
CG_01230_BANDS = (4455, 4465)    # sensitivity: all Colle Gnifetti holes on outline 01230
N_BOOT, SEED = 2000, 1
PUBLISHED = [  # (label, site, depth, year1, year2, published change)
    ('Hoelzle et al. 2011, Colle Gnifetti 20 m 2000-2008', 'Colle Gnifetti', 20.0, 2000, 2008, 1.3),
    ('Hoelzle et al. 2011, Grenzgletscher 1991-2008', 'Grenzgletscher 4250 m', 20.0, 1991, 2008, 6.8),
    ('Vincent et al. 2020, Col du Dome 50 m 1994-2017', 'Col du Dome', 50.0, 1994, 2017, 1.5),
]


def band_of(elev):
    return int(np.floor(elev / 10.0) * 10 + 5)


def ym_index(ym):
    return (ym // 100) * 12 + ym % 100 - 1


def site_of(row):
    """Site label for grouping and reporting."""
    name, label, z = row.glacier_name, str(row.label), row.elevation
    if name == 'Grenzgletscher':
        if label.startswith('CG'):
            return 'Colle Gnifetti'
        if label.startswith('SJ'):
            return 'Seserjoch'
        if label.startswith('GG') and abs(z - 4250) <= 15:
            return 'Grenzgletscher 4250 m'
        return f'Grenzgletscher {int(round(z, -2))} m'
    if name == 'Glacier des Bossons' and z >= 4200:
        return 'Col du Dome'
    return name


# ---------------------------------------------------------------------------------------------
# glenglat
# ---------------------------------------------------------------------------------------------
def load_glenglat(extra=True):
    bh = pd.read_csv(GLENGLAT / 'borehole.csv')
    bh = bh[bh.latitude.between(*BOX['lat']) & bh.longitude.between(*BOX['lon'])].copy()
    sources = [('main', GLENGLAT), ('beer2026', GLENGLAT / 'beer2026')]
    if extra:
        sources += [(s, GLENGLAT / s) for s in EXTRA_SOURCES]
    prof, meas = [], []
    for tag, d in sources:
        p = pd.read_csv(d / 'profile.csv', dtype={'equilibrium': str}, keep_default_na=False,
                        na_values={'date_min': [''], 'date_max': ['']})
        p['data_source'] = tag
        prof.append(p)
        meas.append(pd.read_csv(d / 'measurement.csv'))
    prof = pd.concat(prof, ignore_index=True)
    meas = pd.concat(meas, ignore_index=True)
    if prof.duplicated(['borehole_id', 'id']).any():
        raise SystemExit('duplicate (borehole_id, profile id) across glenglat folders')
    prof = prof[prof.borehole_id.isin(bh.id)]
    meas = meas[meas.borehole_id.isin(bh.id)]
    look = pd.read_csv(GLIMS_LOOKUP)
    look['outline'] = look.rgi7_glogem.str.split('.').str[-1]
    bh = bh.merge(look[['glims_id', 'outline']], on='glims_id', how='left')
    bh['band'] = bh.elevation.map(band_of)
    bh['site'] = bh.apply(site_of, axis=1)
    return bh.set_index('id', drop=False), prof, meas


def profile_months(p):
    """(months yyyymm, weights, fractions) of each profile, as the calibration dates readings
    (data.py: up to 31 days at the midpoint, longer as the range of months), inside the run."""
    t0 = pd.to_datetime(p.date_min, errors='coerce')
    t1 = pd.to_datetime(p.date_max, errors='coerce').fillna(t0)
    tm = t0 + (t1 - t0) / 2
    long = (t1 - t0).dt.days > 31
    out = []
    for a, b, m, lg in zip(t0, t1, tm, long):
        if pd.isna(a):
            out.append(None)
            continue
        if not lg:
            ym = m.year * 100 + m.month
            ent = [(ym, 1.0, m.day / m.days_in_month)] if WINDOW[0] <= ym <= WINDOW[1] else []
        else:
            lo = max(ym_index(a.year * 100 + a.month), ym_index(WINDOW[0]))
            hi = min(ym_index(b.year * 100 + b.month), ym_index(WINDOW[1]))
            yms = [(i // 12) * 100 + i % 12 + 1 for i in range(lo, hi + 1)]
            ent = [(ym, 1.0 / len(yms), 0.5) for ym in yms]
        out.append(ent or None)
    return out, tm


def eligible_profiles(prof, meas):
    p = prof[prof.equilibrium.isin(EQ_OK)].copy()
    p['months'], p['t_mid'] = profile_months(p)
    p = p[p.months.notna()].copy()
    p['year'] = p.t_mid.dt.year
    zr = meas.groupby(['borehole_id', 'profile_id']).depth.agg(['min', 'max'])
    p = p.join(zr, on=['borehole_id', 'id'], how='inner')
    return p


def visits(p):
    """(hole, depth) -> {year: [profile ids spanning the depth]}"""
    out = {}
    for d in DEPTHS:
        sp = p[(p['min'] <= d) & (p['max'] >= d)]
        for (h, y), g in sp.groupby(['borehole_id', 'year']):
            out.setdefault((h, d), {})[int(y)] = sorted(g.id.astype(int).tolist())
    return out


def dist_m(a, b):
    la1, lo1, la2, lo2 = map(np.radians, (a.latitude, a.longitude, b.latitude, b.longitude))
    h = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return float(2 * 6371000 * np.arcsin(np.sqrt(h)))


def build_pairs(bh, p):
    vis = visits(p)
    holes = sorted({h for h, _ in vis})
    src = p.groupby('borehole_id').data_source.agg(lambda s: '/'.join(sorted(set(s))))
    rows = []

    def add(rule, h1, y1, h2, y2, d, dd, dz):
        a, b = bh.loc[h1], bh.loc[h2]
        rows.append(dict(
            rule=rule, site=a.site, glacier=a.glacier_name, depth=d,
            hole1=h1, label1=a.label, elev1=a.elevation, band1=a.band, outline1=a.outline,
            glims1=a.glims_id, year1=y1, profiles1=' '.join(map(str, vis[(h1, d)][y1])),
            hole2=h2, label2=b.label, elev2=b.elevation, band2=b.band, outline2=b.outline,
            glims2=b.glims_id, year2=y2, profiles2=' '.join(map(str, vis[(h2, d)][y2])),
            dt_years=y2 - y1, dist_m=round(dd, 1), dz_m=round(dz, 1),
            extra_source=any(s in EXTRA_SOURCES for s in (src[h1] + '/' + src[h2]).split('/'))))

    for h in holes:                                   # (a) same hole, earliest vs latest visit
        for d in DEPTHS:
            ys = sorted(vis.get((h, d), {}))
            if len(ys) >= 2 and ys[-1] - ys[0] >= MIN_YEARS:
                add('a', h, ys[0], h, ys[-1], d, 0.0, 0.0)
    for h1, h2 in itertools.combinations(holes, 2):  # (b) neighbouring holes
        a, b = bh.loc[h1], bh.loc[h2]
        if a.glacier_name != b.glacier_name:
            continue
        dd, dz = dist_m(a, b), b.elevation - a.elevation
        if dd > MAX_DIST or abs(dz) > MAX_DZ:
            continue
        for d in DEPTHS:
            va, vb = vis.get((h1, d), {}), vis.get((h2, d), {})
            best = None
            for ya, yb in itertools.product(sorted(va), sorted(vb)):
                sep = abs(yb - ya)
                if sep >= MIN_YEARS and (best is None or sep > best[0]):
                    best = (sep, ya, yb)
            if best is None:
                continue
            _, ya, yb = best
            if ya <= yb:
                add('b', h1, ya, h2, yb, d, dd, dz)
            else:
                add('b', h2, yb, h1, ya, d, dd, -dz)
    df = pd.DataFrame(rows).sort_values(['site', 'depth', 'year1', 'hole1', 'hole2'])
    df.insert(0, 'pair_id', [f'P{i:03d}' for i in range(1, len(df) + 1)])
    return df


def phase_pairs():
    OUT.mkdir(parents=True, exist_ok=True)
    bh, prof, meas = load_glenglat()
    p = eligible_profiles(prof, meas)
    pairs = build_pairs(bh, p)
    pairs.to_csv(OUT / 'pairs.csv', index=False)
    print(f'{len(pairs)} pair-depths written to {OUT / "pairs.csv"}')
    print(pairs.groupby(['site', 'rule', 'depth']).size().to_string())


# ---------------------------------------------------------------------------------------------
# measured values
# ---------------------------------------------------------------------------------------------
_PROFILES = {}


def profile_T(meas, h, pid, d):
    """Measured temperature at depth d, linear within the profile, NaN outside it."""
    if not _PROFILES:
        g = meas.groupby(['borehole_id', 'profile_id', 'depth']).temperature.mean()
        for (b, k), s in g.groupby(level=[0, 1]):
            _PROFILES[(int(b), int(k))] = (s.index.get_level_values(2).to_numpy(float),
                                           s.to_numpy(float))
    z, t = _PROFILES.get((int(h), int(pid)), (np.array([]), np.array([])))
    if len(z) == 0 or d < z[0] or d > z[-1]:
        return np.nan
    return float(np.interp(d, z, t))


def site_noise(bh, p, meas):
    """RMS difference at each depth between holes of one campaign (same glacier, same calendar
    year, within NOISE_DIST m and NOISE_DZ m), each hole the mean of its profiles that year."""
    vis = visits(p)
    rows = []
    for d in DEPTHS:
        by_year = {}
        for (h, dd), v in vis.items():
            if dd != d:
                continue
            for y, pids in v.items():
                by_year.setdefault(y, []).append((h, pids))
        for y, lst in by_year.items():
            for (h1, p1), (h2, p2) in itertools.combinations(sorted(lst), 2):
                a, b = bh.loc[h1], bh.loc[h2]
                if a.glacier_name != b.glacier_name:
                    continue
                if dist_m(a, b) > NOISE_DIST or abs(a.elevation - b.elevation) > NOISE_DZ:
                    continue
                t1 = np.nanmean([profile_T(meas, h1, k, d) for k in p1])
                t2 = np.nanmean([profile_T(meas, h2, k, d) for k in p2])
                rows.append(dict(depth=d, year=y, site=a.site, hole1=h1, label1=a.label,
                                 hole2=h2, label2=b.label, diff=t2 - t1,
                                 dist_m=round(dist_m(a, b), 1),
                                 dz_m=b.elevation - a.elevation))
    df = pd.DataFrame(rows)
    noise = df.groupby('depth')['diff'].agg(lambda x: float(np.sqrt(np.mean(x ** 2))))
    return noise, df


# ---------------------------------------------------------------------------------------------
# model values
# ---------------------------------------------------------------------------------------------
class ModelFiles:
    def __init__(self, run):
        self.dir = REPO / 'data/bayescal' / run / FIRNICE
        self.cache = {}
        self.bands = {}

    def bands_of(self, gid):
        if gid not in self.bands:
            out = {}
            for f in sorted(self.dir.glob(f'temp_ID*_{gid}.dat')):
                head = f.open().readline()
                out[int(re.search(r'Point elevation\s+(-?\d+)', head).group(1))] = f
            self.bands[gid] = out
        return self.bands[gid]

    def file_for(self, gid, band):
        b = self.bands_of(gid)
        if not b:
            return None, None
        near = min(b, key=lambda e: (abs(e - band), e))
        return near, b[near]

    def read(self, path):
        """Monthly profiles keyed by calendar yyyymm (end of month); checks the calendar."""
        if path in self.cache:
            return self.cache[path]
        lines = path.read_text().splitlines()
        z = np.array([float(c) for c in lines[1].split()[2:]])
        raw = np.array([[float(v) for v in ln.split()] for ln in lines[2:] if ln.strip()])
        hyd, mon = raw[:, 0].astype(int), raw[:, 1].astype(int)
        T = raw[:, 2:2 + len(z)]
        T = np.where(T <= -98.0, np.nan, T)
        cal = hyd - (mon >= 10)
        ym = cal * 100 + mon
        idx = np.array([ym_index(v) for v in ym])
        if not (mon[0] == 10 and ym[0] == WINDOW[0] and ym[-1] == WINDOW[1]
                and np.all(np.diff(idx) == 1)):
            raise SystemExit(f'{path}: calendar check failed ({ym[:3]}..{ym[-3:]})')
        out = (z, T, {int(v): r for r, v in enumerate(ym)})
        self.cache[path] = out
        return out


def model_series(z, T, d):
    """Linear between the two nodes around d, NaN if either is below the modelled column."""
    k = int(np.searchsorted(z, d))
    if k < len(z) and z[k] == d:
        return T[:, k]
    if k == 0 or k >= len(z):
        return np.full(T.shape[0], np.nan)
    f = (d - z[k - 1]) / (z[k] - z[k - 1])
    return (1 - f) * T[:, k - 1] + f * T[:, k]


def model_at(mf, path, d, entries_list):
    """Mean over profiles of the model on each profile's dates (runner.match_model_to_obs:
    linear between the end of the previous month and the end of the month, at the fraction)."""
    z, T, row = mf.read(path)
    s = model_series(z, T, d)
    vals = []
    for entries in entries_list:
        v, w = [], []
        for ym, wt, fr in entries:
            if ym not in row:
                continue
            x = s[row[ym]]
            prev = ym - 1 if ym % 100 > 1 else ym - 89
            if prev in row:
                x = (1 - fr) * s[row[prev]] + fr * x
            v.append(x)
            w.append(wt)
        v = np.asarray(v)
        vals.append(float(np.sum(v * w) / np.sum(w)) if len(v) and np.isfinite(v).all() else np.nan)
    return float(np.mean(vals)) if vals and np.isfinite(vals).all() else np.nan


# ---------------------------------------------------------------------------------------------
# statistics
# ---------------------------------------------------------------------------------------------
def slope0(x, y):
    return float(np.sum(x * y) / np.sum(x * x)) if len(x) else np.nan


def boot_ci(df, ycol, rng):
    sites = df.site.unique()
    out = []
    for _ in range(N_BOOT):
        pick = rng.choice(sites, len(sites), replace=True)
        s = pd.concat([df[df.site == k] for k in pick])
        out.append(slope0(s.dT_meas.to_numpy(), s[ycol].to_numpy()))
    return np.nanpercentile(out, [2.5, 97.5])


def sign3(v, thr):
    return np.where(v > thr, 1, np.where(v < -thr, -1, 0))


def run_stats(df, ycol, noise, rng):
    df = df[np.isfinite(df[ycol]) & np.isfinite(df.dT_meas)]
    x, y = df.dT_meas.to_numpy(), df[ycol].to_numpy()
    thr = df.depth.map(noise).to_numpy()
    rec = dict(n=len(df), n_sites=df.site.nunique(), ols0=slope0(x, y))
    if len(df) >= 3:
        ts = stats.theilslopes(y, x)
        rec.update(theil=ts.slope, theil_lo=ts.low_slope, theil_hi=ts.high_slope,
                   theil_icpt=ts.intercept)
        rec['boot_lo'], rec['boot_hi'] = boot_ci(df, ycol, rng)
        w = df.groupby('site').site.transform(lambda s: 1.0 / len(s)).to_numpy()
        rec['ols0_site_w'] = float(np.sum(w * x * y) / np.sum(w * x * x))
    sm, sy = sign3(x, thr), sign3(y, thr)
    rec['sign3'] = float(np.mean(sm == sy))
    big = sm != 0
    rec['n_big'] = int(big.sum())
    rec['sign_big'] = float(np.mean(np.sign(y[big]) == np.sign(x[big]))) if big.any() else np.nan
    rec['mean_meas'], rec['mean_mod'] = float(x.mean()), float(y.mean())
    rec['rms_diff'] = float(np.sqrt(np.mean((y - x) ** 2)))
    return rec


def phase_model():
    pairs = pd.read_csv(OUT / 'pairs.csv', dtype={'outline1': str, 'outline2': str})
    for k in (1, 2):                             # unlabelled holes: glenglat id
        pairs[f'label{k}'] = pairs[f'label{k}'].fillna('id' + pairs[f'hole{k}'].astype(str))
    bh, prof, meas = load_glenglat()
    p = eligible_profiles(prof, meas)
    months = {(int(r.borehole_id), int(r.id)): r.months for r in p.itertuples()}
    noise, noise_df = site_noise(bh, p, meas)
    noise_df.to_csv(OUT / 'site_noise_pairs.csv', index=False, float_format='%.3f')
    mfs = {r: ModelFiles(r) for r in RUNS}

    rows = []
    for r in pairs.itertuples():
        p1 = [int(k) for k in str(r.profiles1).split()]
        p2 = [int(k) for k in str(r.profiles2).split()]
        T1 = float(np.mean([profile_T(meas, r.hole1, k, r.depth) for k in p1]))
        T2 = float(np.mean([profile_T(meas, r.hole2, k, r.depth) for k in p2]))
        e1 = [months[(r.hole1, k)] for k in p1]
        e2 = [months[(r.hole2, k)] for k in p2]
        row = dict(pair_id=r.pair_id, rule=r.rule, site=r.site, depth=r.depth,
                   label1=r.label1, hole1=r.hole1, year1=r.year1,
                   label2=r.label2, hole2=r.hole2, year2=r.year2, dt_years=r.dt_years,
                   extra_source=r.extra_source, T1_meas=T1, T2_meas=T2, dT_meas=T2 - T1,
                   band1=r.band1, band2=r.band2, outline1=r.outline1, outline2=r.outline2)
        # model: outline and band of the earlier hole, else of the later one
        use = None
        for gid, band in ((r.outline1, r.band1), (r.outline2, r.band2)):
            if isinstance(gid, str):
                near, path = mfs[RUNS[0]].file_for(gid, band)
                if path is not None and abs(near - band) <= MAX_BAND_GAP:
                    use = (gid, band, near)
                    break
        if use is None:
            row['model_note'] = 'no model output for either hole'
        else:
            gid, band, near = use
            row.update(model_outline=gid, model_band=near,
                       model_note='' if (gid, band) == (r.outline1, r.band1)
                       else 'later hole used')
            for run in RUNS:
                _, path = mfs[run].file_for(gid, band)
                m1, m2 = model_at(mfs[run], path, r.depth, e1), model_at(mfs[run], path, r.depth, e2)
                row[f'T1_{run}'], row[f'T2_{run}'], row[f'dT_{run}'] = m1, m2, m2 - m1
            if r.site == 'Colle Gnifetti':      # sensitivity: outline 01230 for every CG hole
                b30 = min(CG_01230_BANDS, key=lambda e: (abs(e - r.band1), e))
                row['cg01230_band'] = b30
                for run in RUNS:
                    path = mfs[run].bands_of('01230')[b30]
                    row[f'dTcg_{run}'] = (model_at(mfs[run], path, r.depth, e2)
                                          - model_at(mfs[run], path, r.depth, e1))
        rows.append(row)
    res = pd.DataFrame(rows)
    for run in RUNS:                             # CG sensitivity column filled for other sites
        res[f'dTcg_{run}'] = res[f'dTcg_{run}'].where(res.site == 'Colle Gnifetti',
                                                      res[f'dT_{run}'])
    res.to_csv(OUT / 'results.csv', index=False, float_format='%.3f')
    write_summary(res, noise, noise_df)


def fmt(v, f='+.2f'):
    return '-' if v is None or (isinstance(v, float) and not np.isfinite(v)) else format(v, f)


def write_summary(res, noise, noise_df):
    rng = np.random.default_rng(SEED)
    main = res[~res.extra_source]
    L = ['# Time test: measured vs modelled temperature change at repeat sites', '',
         f'Pairs (pair x depth) in pairs.csv: {len(res)}; specified data: {len(main)}; '
         f'with model values (pred_w0): {int(np.isfinite(main["dT_ce25_pred_w0"]).sum())}.',
         f'Depths {DEPTHS} m; >= {MIN_YEARS} calendar years; rule (b) <= {MAX_DIST:g} m, '
         f'+-{MAX_DZ:g} m, same glenglat glacier_name.', '',
         '## Site noise (RMS difference between holes of one campaign)', '',
         '| depth [m] | RMS [C] | n hole pairs | sites |', '|---|---|---|---|']
    for d in DEPTHS:
        s = noise_df[noise_df.depth == d]
        L.append(f'| {d:g} | {noise.get(d, np.nan):.2f} | {len(s)} | '
                 + ', '.join(f'{k} ({n})' for k, n in s.site.value_counts().items()) + ' |')
    L += ['', 'Sign classes use this threshold for measured and modelled change alike '
          '(|change| <= noise = no change).', '']

    def table(df, ycols, title):
        out = [f'## {title}', '',
               '| run | n | sites | OLS slope (0) | site-bootstrap 95% | site-weighted OLS | '
               'Theil-Sen [95%] (intercept) | sign agree 3-class | sign agree, |meas|>noise (n) | '
               'mean meas / mod | RMS(mod-meas) |',
               '|---|---|---|---|---|---|---|---|---|---|---|']
        for run, yc in ycols:
            s = run_stats(df, yc, noise, rng)
            out.append(
                f'| {run} | {s["n"]} | {s["n_sites"]} | {fmt(s["ols0"])} | '
                f'{fmt(s.get("boot_lo"))}..{fmt(s.get("boot_hi"))} | {fmt(s.get("ols0_site_w"))} | '
                f'{fmt(s.get("theil"))} [{fmt(s.get("theil_lo"))}, {fmt(s.get("theil_hi"))}] '
                f'({fmt(s.get("theil_icpt"))}) | {s["sign3"]:.0%} | '
                f'{fmt(s["sign_big"], ".0%")} ({s["n_big"]}) | '
                f'{s["mean_meas"]:+.2f} / {s["mean_mod"]:+.2f} | {s["rms_diff"]:.2f} |')
        return out + ['']

    L += table(main, [(r, f'dT_{r}') for r in RUNS], 'Pooled, all depths (specified data)')
    for d in DEPTHS:
        L += table(main[main.depth == d], [(r, f'dT_{r}') for r in RUNS], f'Depth {d:g} m')
    L += table(main[main.rule == 'a'], [(r, f'dT_{r}') for r in RUNS], 'Rule (a) only, same hole')
    L += table(main[main.site != 'Glacier de Tête Rousse'],
               [(r, f'dT_{r}') for r in RUNS], 'Without Tete Rousse (temperate)')
    L += table(main, [(r, f'dTcg_{r}') for r in RUNS],
               'Sensitivity: all Colle Gnifetti holes on outline 01230')
    if res.extra_source.any():
        L += table(res, [(r, f'dT_{r}') for r in RUNS],
                   'Sensitivity: with the Ortles subfolder (carturan2023-ortles)')

    L += ['## Post-hoc sensitivities (chosen after seeing the results above)', '']
    L += table(main[main.site != 'Grenzgletscher 4250 m'], [(r, f'dT_{r}') for r in RUNS],
               'Post-hoc: without Grenzgletscher 4250 m')
    L += table(main[main.year1 >= 1985], [(r, f'dT_{r}') for r in RUNS],
               'Post-hoc: earlier reading 1985 or later (spin-up)')
    cg = noise_df[noise_df.label1.astype(str).str.match(r'CG(08|13)-')
                  & noise_df.label2.astype(str).str.match(r'CG(08|13)-')]
    L += ['Post-hoc: site noise from CG08-1/2/3 and CG13-1..4 alone: '
          + ', '.join(f'{d:g} m {np.sqrt(np.mean(g["diff"] ** 2)):.2f} C (n={len(g)})'
                      for d, g in cg.groupby('depth')), '']

    L += ['## Per site (pred_w0)', '',
          '| site | n | years | mean meas | mean mod (pred / F / H / FH) | level bias pred '
          '(mod - meas, earlier reading) | outlines used |',
          '|---|---|---|---|---|---|---|']
    for site, g in res.groupby('site'):
        g = g[np.isfinite(g.dT_ce25_pred_w0)]
        if g.empty:
            continue
        mods = ' / '.join(f'{g[f"dT_{r}"].mean():+.2f}' for r in RUNS)
        L.append(f'| {site} | {len(g)} | {g.year1.min()}-{g.year2.max()} | '
                 f'{g.dT_meas.mean():+.2f} | {mods} | '
                 f'{(g.T1_ce25_pred_w0 - g.T1_meas).mean():+.2f} | '
                 f'{", ".join(sorted(set(g.model_outline.dropna())))} |')
    L += ['', '## Largest disagreements (pred_w0, |mod - meas|)', '',
          '| pair | site | depth | holes | years | meas | mod | band |', '|---|---|---|---|---|---|---|---|']
    g = res[np.isfinite(res.dT_ce25_pred_w0)].copy()
    g['err'] = (g.dT_ce25_pred_w0 - g.dT_meas).abs()
    for r in g.sort_values('err', ascending=False).head(12).itertuples():
        L.append(f'| {r.pair_id} | {r.site} | {r.depth:g} | {r.label1} -> {r.label2} | '
                 f'{r.year1}-{r.year2} | {r.dT_meas:+.2f} | {r.dT_ce25_pred_w0:+.2f} | '
                 f'{r.model_outline}/{r.model_band:.0f} |')

    L += ['', '## Published numbers', '',
          '| reference | pairs here (holes, years) | meas here | mod pred / F / H / FH | published |',
          '|---|---|---|---|---|']
    for lab, site, d, y1, y2, pub in PUBLISHED:
        g = res[(res.site == site) & (res.depth == d) & (res.year1 == y1) & (res.year2 == y2)
                & np.isfinite(res.dT_ce25_pred_w0)]
        if g.empty:
            L.append(f'| {lab} | none | - | - | {pub:+.1f} |')
            continue
        holes = '; '.join(f'{a}->{b}' for a, b in zip(g.label1, g.label2))
        mods = ' / '.join(f'{g[f"dT_{r}"].mean():+.2f}' for r in RUNS)
        L.append(f'| {lab} | {holes} | {g.dT_meas.mean():+.2f} | {mods} | {pub:+.1f} |')

    nomod = res[res.model_note.fillna('') != '']
    L += ['', '## Model mapping notes', '']
    for r in nomod.itertuples():
        L.append(f'- {r.pair_id} {r.site} {r.label1}->{r.label2} {r.depth:g} m: {r.model_note}')
    nan = res[(res.model_note.fillna('') == '') & ~np.isfinite(res.dT_ce25_pred_w0)]
    for r in nan.itertuples():
        L.append(f'- {r.pair_id} {r.site} {r.label1}->{r.label2} {r.depth:g} m: '
                 f'depth below the modelled column at band {r.model_band:.0f}')
    L += ['', '## Caveats', '',
          '- All repeat sites lie on glaciers used to develop the model: this tests the trend, '
          'not the level.',
          '- The thermal spin-up settles to the 1942-71 climate; early readings (1976-80) may '
          'carry a spin-up bias in the modelled change.',
          '- Rule (b) pairs mix spatial and temporal differences; the site noise bounds the '
          'spatial part.',
          '- Profile files are labelled with the hydrological year (Oct-Dec rows belong to the '
          'previous calendar year); converted and checked continuous Oct 1971-Sep 2025.']
    (OUT / 'summary.md').write_text('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    {'pairs': phase_pairs, 'model': phase_model}[sys.argv[1]]()
