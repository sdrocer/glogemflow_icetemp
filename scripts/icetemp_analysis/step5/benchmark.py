"""Step 5, rules 8-10 of the rules plan (27 Oct): the benchmark, each massif predicted from the other six.

  folds  for each massif with readings below 15 m: rules 1 and 4 (step5_common.fit_rule1_4) on all other entities
         (also those without readings below 15 m, which are not scored), from the C6b start of the main fit and the
         config's fixed parameters -> OUT/folds.json (with the training set and commit), OUT/fold_thetas.csv (the
         real GloGEM runs of step 4.2, tag s5fK, K the fold's number) and each fold's joint grid posterior.
  score  every method on the left-out massifs' readings below 15 m (the points of the config's fit exclusions left
         out here too, forcing_fair_test.scored_entities); every run's parts checked against their design row and the
         config (assemble_training.check_part), every prediction finite and at the reading's depth:
           step5      the real run at each fold's theta (five parts merged, tag s5fK)
           step5_emu  the emulator at each fold's theta (reported; the decision uses step5 once its runs exist)
           grid       f_ref chosen from 0.03-1.00 at f_ins 1.017 by M1 on the other massifs (forcing_fair_test.lomo,
                      runs s5g00-07)
           start      the starting values f_ref 0.55, f_ins 1.0 (run s5s00)
           elev_line  T = a + b * band elevation, least squares on the other massifs' readings below 15 m, each
                      outline weighted equally
           air_line   T = a + b * the band's mean free-air temperature October 1971 - September 2025 (nearest grid
                      point of the free-air table, linear in height), fitted the same way
           nearest5   the mean of the 5 bands of the other massifs nearest in elevation (each band its mean below 15 m)
           elev_class classes only: two elevation thresholds (temperate below the lower, polythermal between, cold
                      above) that classify most bands of the other massifs right (ties: the lowest thresholds)
           suter      classes only: cold above the lower limit of cold firn of Suter et al. (2001), 3400 m on north and
                      4150 m on south faces, 3775 m on east, west and gentle parts (the hole's part); temperate below
         Measures: M1 (forcing_fair_test.m1) overall, in the zone 3000-3800 m and above 3800 m; regime classes per
         band (counts per measured class and the confusion table, the hole's tolerance) and per site (right when
         more than half its bands are); frozen beds for the run methods, those close to melting (measured at -2 C or
         warmer) apart from the very cold ones; for step5 the share of readings inside their 90 % predictive range
         (the fold's grid posterior, its error model's variance for an unseen site and the t of each reading).
         Uncertainty: whole rule sites resampled (N_BOOT times; the n-th copy of a site joins the n-th copies of the
         other sites of its outline, so a draw of every site once gives M1 itself); the grid choice and the fold
         results stay fixed. Adoption (rule 9): the 90th percentile of the resampled M1(step5) - M1(grid) at most
         0.1 C, overall and in the zone, and 80-97 % of the readings inside their 90 % ranges; the gate and the
         emulator agreement of fit_step5.py are reported, not part of it. Rule 10 with rule 3: the map's uncertainty
         per parameter is the largest of the result's own (step5_common.result_sd), the sd of the fold results, and,
         for each site decisive in check 3 (site_checks.csv), its shift and its fit's sd. Also reported, not used:
         (n - 1) times the fold sd, about how far one massif's own value lies from the result (for a mean exactly).

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/step5/benchmark.py folds|score
           [CONFIG] [--fit DIR] [--out DIR] [--stand-in] [--test]
"""
import argparse
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'orientation'))
import assemble_training as at  # noqa: E402
import step5_common as sc  # noqa: E402
from step5_common import ft  # noqa: E402
from icetemp.calibration import thermal_structure as ts  # noqa: E402
from icetemp.calibration.emulator import BASAL  # noqa: E402

W3 = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3')
HOLE_PARTS = W3 / 'orient_parts5' / 'hole_parts.csv'
FREEAIR = W3 / 'freeair' / 'freeair_CentralEurope.nc'
GRID = (0.03, 0.10, 0.25, 0.40, 0.55, 0.70, 0.85, 1.00)
GRID_FINS, START = 1.017, (0.55, 1.0)
SUTER = dict(N=3400.0, S=4150.0, other=3775.0)
N_NEAREST, N_BOOT, SEED = 5, 4000, 20261028
ADOPT_MARGIN, COVER = 0.1, (0.80, 0.97)
CLASSES = ft.CLASSES
# run names before _oCP_w0; STAND_IN (--stand-in, a test of the code only): the free-air base runs fb1 for all three
RUNS = dict(grid='s5g{k:02d}_FRHS_r{fr:.2f}_i1.017_a1.00_fa_fs1.4', start='s5s00_FRHS_r0.55_i1.000_a1.00_fa_fs1.4',
            fold='s5f{k}_*_fa_fs1.4', fold_tag='s5f')
# the redo on the free-air mass balance (config air_temperature = 'freeair'; plan of 29 Oct)
RUNS_MBFA = dict(grid='s6g{k:02d}_FRHS_r{fr:.2f}_i1.017_a1.00_mbfa_fa_fs1.4',
                 start='s6s00_FRHS_r0.55_i1.000_a1.00_mbfa_fa_fs1.4', fold='s6f{k}_*_mbfa_fa_fs1.4', fold_tag='s6f')
STAND_IN = dict(grid='fb1_FRHS_r{fr:.2f}_i1.017_a1.00_fa', start='fb1_FRHS_r0.55_i1.017_a1.00_fa',
                fold='fb1_FRHS_r0.55_i1.017_a1.00_fa')
PARTS = 'GNESW'
CHECK_RUNS = True   # off with --stand-in, whose runs are at other points


# ── runs ───────────────────────────────────────────────────────────────────────────────────
def run_names(config):
    return RUNS_MBFA if config.air_temperature == 'freeair' and CHECK_RUNS else RUNS


def one_run(pattern, theta, config, missing_ok=False):
    """The merged run matching pattern; its five part runs ran theta with the config's thermal settings."""
    hits = sorted(ft.BAYESCAL.glob(pattern + '_oCP_w0'))
    if missing_ok and not hits:
        return None
    if len(hits) != 1:
        raise SystemExit(f'{pattern}: {len(hits)} run folders, expected 1')
    name = hits[0].name
    if CHECK_RUNS:
        for p in PARTS:
            at.check_part(ft.BAYESCAL / f'{name[:-len("_oCP_w0")]}_oC{p}_w0', np.asarray(theta, float), config)
    return name


def grid_runs(config):
    return {fr: one_run(run_names(config)['grid'].format(k=k, fr=fr), (fr, GRID_FINS, 1.0), config)
            for k, fr in enumerate(GRID)}


def start_run(config):
    return one_run(run_names(config)['start'], (*START, 1.0), config)


def fold_run(f, config):
    return one_run(run_names(config)['fold'].format(k=f['fold']), f['theta'], config, missing_ok=True)


def emu_rows(emu):
    """(entity, point index) -> the emulator's output row, depth rows only."""
    return {(n, int(k)): r for r, ((n, d), k) in enumerate(zip(emu.x_index_used(), emu.x_points_used()))
            if d != BASAL and k is not None}


def output(run):
    return pickle.load(open(ft.output_path(run, 'hole_year'), 'rb'))


# ── points ─────────────────────────────────────────────────────────────────────────────────
def frame(E, glaciers, sites):
    """One row per scored point: entity, its index k in the entity's points, hole, year, outline, massif, band, site."""
    by = {g.glacier_name: g for g in glaciers}
    rows = []
    for name, e in E.items():
        g = by[name]
        for j, k in enumerate(e['idx']):
            rows.append(dict(entity=name, k=int(k), hole=str(g.point_hole[k]), year=int(g.point_year[k]),
                             rgi=e['rgi'], massif=e['massif'], band=e['band'], site=sites[name],
                             depth=float(e['depth'][j]), obs=float(e['obs'][j]), sig2=float(g.sigma[k]) ** 2))
    return pd.DataFrame(rows)


def add_run(P, col, run):
    out = output(run)
    vals = []
    for e, k, d in zip(P.entity, P.k, P.depth):
        if out.get(e) is None:
            vals.append(np.nan)
            continue
        dep, val = np.asarray(out[e][0], float), np.asarray(out[e][1], float)
        if not np.isclose(dep[k], d):
            raise SystemExit(f'{run}: {e} point {k} is at {dep[k]} m in the output, {d} m in the data')
        vals.append(float(val[k]))
    P[col] = vals


def massifs_with_points(P):
    return sorted(P.massif.unique())


# ── folds ──────────────────────────────────────────────────────────────────────────────────
def grid_posterior(cal, n=61):
    """Grid points (full theta) and normalised weights, as BayesianCalibrator.grid_moments."""
    eff = cal._effective_bounds()
    free = cal._free_param_names()
    axes = [np.linspace(*eff[name], n) for name in free]
    pts = np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1).reshape(-1, len(free))
    lp = np.array([cal.log_posterior(cal._expand_theta(p)) for p in pts])
    ok = np.isfinite(lp)
    w = np.exp(lp[ok] - lp[ok].max())
    return np.array([cal._expand_theta(p) for p in pts[ok]]), w / w.sum()


def cmd_folds(pl, F, out, head):
    glaciers = pl._calibration_glaciers()
    E = ft.scored_entities(glaciers, ts.load_borehole_table(), pl.config.fit_exclusions)
    massifs = sorted({e['massif'] for e in E.values()})
    print(f'{len(massifs)} massifs with readings below 15 m: {massifs}', flush=True)
    folds, rows = [], []
    for k, m in enumerate(massifs):
        train = [g for g in glaciers if sc.massif_of(g) != m]
        print(f'fold {k} ({m}): {len(glaciers) - len(train)} entities left out', flush=True)
        r = sc.fit_rule1_4(pl, sizes=F['rounds'][0]['sizes'], glaciers=train, fixed_params=F['config_fixed_params'])
        pts, w = grid_posterior(r['cal'])
        folds.append(dict(fold=k, massif=m, n_left_out=len(glaciers) - len(train), theta=r['theta'].tolist(),
                          names=list(r['moments']['names']), sd=r['moments']['sd'].tolist(), rule4=r['rule4'],
                          fixed=r['fixed'], sizes=sc.only_sizes(r['sizes']), rounds=r['rounds'],
                          second_estimate_gain=r['second_estimate_gain']))
        with open(out / f'fold{k}_posterior.pkl', 'wb') as f:
            pickle.dump(dict(points=pts, weights=w, sizes=sc.only_sizes(r['sizes'])), f)
        rows.append(dict(fold=k, massif=m, tag=f"{run_names(pl.config)['fold_tag']}{k}", **dict(zip(('refreeze_frac', 'insul_scale', 'advection_scale'),
                                                                     r['theta'].tolist()))))
    rec = dict(commit=head, training={k: pl.training[k] for k in ('n_runs', 'sha256')}, folds=folds)
    (out / 'folds.json').write_text(json.dumps(rec, indent=1, default=float))
    pd.DataFrame(rows).to_csv(out / 'fold_thetas.csv', index=False)
    print(pd.DataFrame(rows).to_string(index=False))


# ── simple rules ───────────────────────────────────────────────────────────────────────────
def wls_line(x, y, unit):
    w = 1.0 / pd.Series(unit).map(pd.Series(unit).value_counts()).to_numpy(float)
    X = np.column_stack([np.ones(len(x)), x])
    beta = np.linalg.solve(X.T @ (X * w[:, None]), X.T @ (w * y))
    return beta


def free_air_means(glaciers, path=FREEAIR, y0=1971, y1=2025):
    """entity -> mean free-air temperature at its band, October y0 - September y1, nearest grid point."""
    import xarray as xr
    ds = xr.open_dataset(path)
    ym = ds['year'].to_numpy() * 100 + ds['month'].to_numpy()
    sel = (ym >= y0 * 100 + 10) & (ym <= y1 * 100 + 9)
    out = {}
    for g in glaciers:
        col = ds['tfa'].sel(latitude=g.latitude, longitude=g.longitude, method='nearest').isel(time=np.flatnonzero(sel))
        out[g.glacier_name] = float(col.interp(height=float(g.elevation_band)).mean())
    return out


def thresholds(bands, cls):
    """(z1, z2) classifying most of these bands right: temperate below z1, polythermal in [z1, z2), cold above."""
    cand = np.r_[-np.inf, np.unique(bands) + 5.0, np.inf]
    best, arg = -1, None
    for i, z1 in enumerate(cand):
        for z2 in cand[i:]:
            pred = np.where(bands < z1, 'temperate', np.where(bands < z2, 'polythermal', 'cold'))
            n = int((pred == cls).sum())
            if n > best:
                best, arg = n, (float(z1), float(z2))
    return arg


def entity_part(P, parts):
    return {e: pd.Series([parts.get(h, 'G') for h in d.hole]).mode().iloc[0] for e, d in P.groupby('entity')}


# ── measures ───────────────────────────────────────────────────────────────────────────────
def classes_of(P, col, E):
    return ft.entity_classes(P, col, E)


def class_measures(obs_cls, mod_cls, names, site_of):
    counts = ft.class_counts(obs_cls, mod_cls, names)
    conf = pd.crosstab(pd.Series([obs_cls[n]['hole'] for n in names], name='measured'),
                       pd.Series([mod_cls[n]['hole'] for n in names], name='predicted'))
    right = {n: obs_cls[n]['hole'] == mod_cls[n]['hole'] for n in names}
    sites = pd.Series({n: site_of[n] for n in names})
    per_site = {s: bool(np.mean([right[n] for n in sites.index[sites == s]]) > 0.5) for s in sorted(set(sites))}
    return dict(counts={c: f'{a} of {b}' for c, (a, b) in counts.items()},
                bands_right=int(sum(right.values())), bands=len(names),
                sites_right=int(sum(per_site.values())), sites=len(per_site), confusion=conf.to_dict())


def run_measures(P, col):
    zone = P.band.between(*ft.ZONE)
    return dict(M1=ft.m1(P, col), zone_M1=ft.m1(P[zone], col), above_M1=ft.m1(P[P.band >= 3800], col))


def frozen_beds(E, chosen, base_of, part_of_hole):
    """ps.frozen_beds_parts with run bases per chosen key (each massif its run)."""
    import hole_parts as hp_mod
    import parts_score as ps
    hp_mod.SCHEME['name'] = 'GNESW'
    ps.CFG.update(letters='GNESW')
    ps.GRID = sorted(set(chosen.values()), key=str)
    ps.PART_RUN = lambda p, key: f'{base_of[key]}_oC{p}_w0'
    FB, _ = ps.frozen_beds_parts(E, chosen, part_of_hole)
    F = FB[np.isfinite(FB['T_CP'].astype(float))]
    near = F.frozen_meas & ~F.deep_cold
    found = F.frozen_meas & F['frozen_CP']
    return dict(found=int(found.sum()), of=int(F.frozen_meas.sum()),
                near_melting_found=int((found & near).sum()), near_melting_of=int(near.sum()),
                very_cold_found=int((found & F.deep_cold).sum()), very_cold_of=int((F.frozen_meas & F.deep_cold).sum()),
                temperate_called_frozen=int((~F.frozen_meas & F['frozen_CP']).sum())), FB


def coverage(P, folds, emu, out_dir, w_max_depth, level=0.90):
    """Share of the left-out readings inside their 90 % predictive range, per fold's grid posterior."""
    row_of = emu_rows(emu)
    inside = np.full(len(P), np.nan)
    for f in folds:
        with open(out_dir / f'fold{f["fold"]}_posterior.pkl', 'rb') as fh:
            post = pickle.load(fh)
        keep = post['weights'] > 1e-4 * post['weights'].max()
        pts, w = post['points'][keep], post['weights'][keep] / post['weights'][keep].sum()
        s = post['sizes']
        mu = np.vstack([emu.predict_factor(t)[0] for t in pts])          # (n_points, n_rows)
        trunc = getattr(emu, 'truncation_variance', lambda: 0.0)()
        nu = s['site_t_nu']
        sel = np.flatnonzero(P.massif.to_numpy() == f['massif'])
        for i in sel:
            r = row_of[(P.entity.iat[i], int(P.k.iat[i]))]
            var = (P.sig2.iat[i] + trunc + s['s_model_variance'] + s['s_hole_variance'] + s['s_hole_year_variance']
                   + s['s_site_variance'] + (s['s_white_variance'] if P.depth.iat[i] <= w_max_depth else 0))
            if nu is None:
                dist = lambda q: stats.norm.cdf(q, loc=mu[:, r], scale=np.sqrt(var))  # noqa: E731
            else:
                dist = lambda q: stats.t.cdf(q, df=nu, loc=mu[:, r], scale=np.sqrt(var * (nu - 2) / nu))  # noqa: E731
            cdf = float(w @ dist(P.obs.iat[i]))
            inside[i] = (1 - level) / 2 <= cdf <= 1 - (1 - level) / 2
    return inside


def bootstrap(P, cols, rng, n=N_BOOT):
    """Resampled M1 differences cols[0] - cols[1], overall and in the zone, whole rule sites resampled."""
    sites = np.array(sorted(P.site.unique()))
    by = {s: np.flatnonzero(P.site.to_numpy() == s) for s in sites}
    zone = P.band.between(*ft.ZONE).to_numpy()
    d_all, d_zone = [], []
    for _ in range(n):
        draw = rng.choice(sites, len(sites), replace=True)
        idx = np.concatenate([by[s] for s in draw])
        seen = {}
        copy = []
        for s in draw:
            seen[s] = seen.get(s, -1) + 1
            copy.append(np.full(len(by[s]), seen[s]))
        copy = np.concatenate(copy)
        Q = P.iloc[idx].assign(rgi=[f'{r}#{c}' for r, c in zip(P.rgi.to_numpy()[idx], copy)])
        z = zone[idx]
        d_all.append(ft.m1(Q, cols[0]) - ft.m1(Q, cols[1]))
        d_zone.append(ft.m1(Q[z], cols[0]) - ft.m1(Q[z], cols[1]))
    return np.array(d_all), np.array(d_zone)


def cmd_score(pl, F, out, head, fit_dir):
    glaciers = pl._calibration_glaciers()
    E = ft.scored_entities(glaciers, ts.load_borehole_table(), pl.config.fit_exclusions)
    site_of = pl.site_of(glaciers)
    P = frame(E, glaciers, site_of)
    B = json.loads((out / 'folds.json').read_text())
    sc.same_training(B['training'], pl, 'folds.json')
    folds = B['folds']
    fold_of = {f['massif']: f for f in folds}
    if sorted(fold_of) != massifs_with_points(P):
        raise SystemExit('the folds do not match the massifs of the scored points')
    S3 = pd.read_csv(fit_dir / 'site_checks.csv')
    gr = grid_runs(pl.config)
    for fr, run in gr.items():
        add_run(P, ft.col('G', fr), run)
    add_run(P, 'start', start_run(pl.config))
    P['grid'], chosen = ft.lomo(P, {fr: ft.col('G', fr) for fr in GRID})
    fruns = {f['fold']: fold_run(f, pl.config) for f in folds}
    have_runs = all(fruns.values())
    P['step5'] = np.nan
    for f in folds:
        if fruns[f['fold']]:
            col = f'fold{f["fold"]}'
            add_run(P, col, fruns[f['fold']])
            sel = P.massif == f['massif']
            P.loc[sel, 'step5'] = P.loc[sel, col]
    emu = pl.emulator
    row_of = emu_rows(emu)
    rows = np.array([row_of[(e, k)] for e, k in zip(P.entity, P.k)])
    P['step5_emu'] = np.nan
    for f in folds:
        sel = (P.massif == f['massif']).to_numpy()
        P.loc[sel, 'step5_emu'] = emu.predict_factor(np.asarray(f['theta']))[0][rows[sel]]

    # simple rules, each fitted on the other massifs
    air = free_air_means(glaciers)
    P['air'] = P.entity.map(air)
    for c in ('elev_line', 'air_line', 'nearest5'):
        P[c] = np.nan
    ent_mean = P.groupby('entity').agg(obs=('obs', 'mean'), band=('band', 'first'), massif=('massif', 'first'))
    obs_cls = classes_of(P, 'obs', E)
    parts = dict(pd.read_csv(HOLE_PARTS, dtype={'hole': str})[['hole', 'part']].itertuples(index=False))
    epart = entity_part(P, parts)
    cls_rules = {'elev_class': {}, 'suter': {}}
    fitted = {}
    for m in massifs_with_points(P):
        tr, te = P.massif != m, P.massif == m
        b = wls_line(P.band[tr].to_numpy(float), P.obs[tr].to_numpy(), P.rgi[tr].to_numpy())
        P.loc[te, 'elev_line'] = b[0] + b[1] * P.band[te]
        a = wls_line(P.air[tr].to_numpy(float), P.obs[tr].to_numpy(), P.rgi[tr].to_numpy())
        P.loc[te, 'air_line'] = a[0] + a[1] * P.air[te]
        em_tr = ent_mean[ent_mean.massif != m]
        for e in ent_mean.index[ent_mean.massif == m]:
            d = (em_tr.band - ent_mean.band[e]).abs().sort_values(kind='stable')
            P.loc[P.entity == e, 'nearest5'] = em_tr.obs[d.index[:N_NEAREST]].mean()
        names_tr = list(em_tr.index)
        z1, z2 = thresholds(em_tr.band.to_numpy(float), np.array([obs_cls[n]['hole'] for n in names_tr]))
        for e in ent_mean.index[ent_mean.massif == m]:
            bnd = ent_mean.band[e]
            c = 'temperate' if bnd < z1 else 'polythermal' if bnd < z2 else 'cold'
            cls_rules['elev_class'][e] = {t: c for t in ft.TOLS}
            lim = SUTER.get(epart[e], SUTER['other'])
            cls_rules['suter'][e] = {t: ('cold' if bnd >= lim else 'temperate') for t in ft.TOLS}
        fitted[m] = dict(elev_line=b.tolist(), air_line=a.tolist(), elev_class=[z1, z2])

    names = sorted(P.entity.unique())
    methods = (['step5'] if have_runs else []) + ['step5_emu', 'grid', 'start', 'elev_line', 'air_line', 'nearest5']
    bad = {c: sorted(P.entity[~np.isfinite(P[c].to_numpy(float))].unique())
           for c in methods + [ft.col('G', fr) for fr in GRID] if not np.isfinite(P[c].to_numpy(float)).all()}
    if bad:
        raise SystemExit(f'predictions missing: {bad}')
    res = dict(commit=head, training=B['training'], inputs={p.name: sc.sha256(p) for p in (HOLE_PARTS, FREEAIR)},
               points=len(P), entities=len(names), massifs=massifs_with_points(P), grid_chosen=chosen,
               fold_runs=fruns, simple_rule_fits=fitted, methods={})
    for c in methods:
        res['methods'][c] = dict(**run_measures(P, c), classes=class_measures(obs_cls, classes_of(P, c, E), names,
                                                                              site_of))
    for c, mc in cls_rules.items():
        res['methods'][c] = dict(classes=class_measures(obs_cls, mc, names, site_of))

    hp = dict(pd.read_csv(HOLE_PARTS, dtype={'hole': str})[['hole', 'part']].itertuples(index=False))
    fb = {}
    fb['grid'], FBg = frozen_beds(E, chosen, {fr: gr[fr][:-len('_oCP_w0')] for fr in GRID}, hp)
    s0 = start_run(pl.config)[:-len('_oCP_w0')]
    fb['start'], FBs = frozen_beds(E, {m: 'start' for m in res['massifs']}, {'start': s0}, hp)
    if have_runs:
        fb['step5'], FB5 = frozen_beds(E, {f['massif']: f['fold'] for f in folds},
                                       {f['fold']: fruns[f['fold']][:-len('_oCP_w0')] for f in folds}, hp)
        FB5.to_csv(out / 'frozen_bed_step5.csv', index=False)
    for c, v in fb.items():
        res['methods'][c]['frozen_beds'] = v

    P['inside90'] = coverage(P, folds, emu, out, pl.config.s_white_max_depth)
    res['coverage90'] = float(np.nanmean(P['inside90']))
    rng = np.random.default_rng(SEED)
    main_col = 'step5' if have_runs else 'step5_emu'
    d_all, d_zone = bootstrap(P, (main_col, 'grid'), rng)
    p90_all, p90_zone = float(np.quantile(d_all, 0.9)), float(np.quantile(d_zone, 0.9))
    adopt = dict(compared=main_col, diff_M1=res['methods'][main_col]['M1'] - res['methods']['grid']['M1'],
                 diff_zone_M1=res['methods'][main_col]['zone_M1'] - res['methods']['grid']['zone_M1'],
                 p90_diff_M1=p90_all, p90_diff_zone_M1=p90_zone,
                 ci90_diff_M1=np.quantile(d_all, [0.05, 0.95]).tolist(),
                 ci90_diff_zone_M1=np.quantile(d_zone, [0.05, 0.95]).tolist(),
                 coverage90=res['coverage90'])
    adopt['checks'] = {'M1': p90_all <= ADOPT_MARGIN, 'zone M1': p90_zone <= ADOPT_MARGIN,
                       'coverage': COVER[0] <= res['coverage90'] <= COVER[1]}
    adopt['adopt'] = bool(all(adopt['checks'].values())) and have_runs
    if not have_runs:
        adopt['note'] = 'provisional: the fold runs (s5fK) are missing, step5_emu stands in; not a decision'
    res['adoption'] = adopt
    decisive = S3[S3.check3_decisive]
    adopt['check3_decisive_sites'] = decisive.site.tolist()
    th = np.array([f['theta'] for f in folds])
    names_p = list(F['grid']['names'])
    idx = [['refreeze_frac', 'insul_scale', 'advection_scale'].index(n) for n in names_p]
    spread = th[:, idx].std(axis=0, ddof=1)
    own = np.array([sc.result_sd(F)[n] for n in names_p])
    site = np.array([max([0.0] + [max(abs(r[f'shift_{n}']), r[f'sd_without_{n}']) for _, r in decisive.iterrows()])
                     if f'shift_{n}' in S3 else 0.0 for n in names_p])
    res['map_uncertainty'] = dict(names=names_p, result_sd=own.tolist(), fold_sd=spread.tolist(),
                                  decisive_site_reach=site.tolist(), used=np.maximum.reduce([own, spread, site]).tolist(),
                                  reported_single_massif_spread=((len(folds) - 1) * spread).tolist())
    P.to_csv(out / 'benchmark_points.csv', index=False, float_format='%.3f')
    (out / 'benchmark.json').write_text(json.dumps(res, indent=1, default=str))
    for c, v in res['methods'].items():
        print(f'{c:10s}', {k: round(x, 3) for k, x in v.items() if isinstance(x, float)},
              v['classes']['counts'], f"sites {v['classes']['sites_right']}/{v['classes']['sites']}",
              v.get('frozen_beds', ''))
    print('coverage90', round(res['coverage90'], 3), '| adoption', adopt)
    print('map uncertainty', res['map_uncertainty'])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('cmd', choices=['folds', 'score'])
    ap.add_argument('config', nargs='?', default=str(sc.CONFIG))
    ap.add_argument('--fit', default=None)
    ap.add_argument('--out', default=None, help='default: <output_dir>/benchmark')
    ap.add_argument('--stand-in', action='store_true', help='fb1 runs for the grid, start and fold runs (code test)')
    ap.add_argument('--test', action='store_true', help='a code test: skips the committed-tree check')
    a = ap.parse_args()
    global CHECK_RUNS
    if a.stand_in:
        RUNS.update(STAND_IN)
        CHECK_RUNS = False
    head = sc.rules_commit(a.config, test=a.test or a.stand_in)
    pl = sc.pipeline(a.config)
    fit_dir = Path(a.fit) if a.fit else pl.config.output_path / 'fit'
    F = sc.load_fit(fit_dir, pl)
    out = Path(a.out) if a.out else pl.config.output_path / 'benchmark'
    out.mkdir(parents=True, exist_ok=True)
    if a.cmd == 'folds':
        cmd_folds(pl, F, out, head)
    else:
        cmd_score(pl, F, out, head, fit_dir)


if __name__ == '__main__':
    main()
