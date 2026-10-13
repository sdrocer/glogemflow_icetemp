"""The fair run set of week 2: the forcing, decided by the rules fixed before the runs.

Pre-registration: development_history/plans/2026-10-13_fair_test_preregistration.md (committed with
this script before the first run). Runs: full_set_runs.py with --fref/--fins, catchment
CentralEurope_fair, sealed 01161 excluded; tags <tag>_<variant>_r<f_ref>_i<f_ins>_a<f_adv>[_nodeep]_w0.

Decisive, out of sample (f_ref chosen leaving one massif out by M1):
  M1       RMSE below 15 m, each outline with >= 3 points weighted equally (zone: >= 3 zone points)
  classes  correctly classified entities per measured class (cold, polythermal, temperate),
           robust changes only (same at tolerances 0.05 C, the hole's own, 0.5 C)
  warming  week-1 repeat pairs on the new outlines: mean absolute site error and ratio of means
Decision A: FRH adopted unless (1) M1 worse by > 0.3 C overall or in the zone, (2) net robust losses
  of >= 2 polythermal, >= 2 temperate or >= 4 entities overall, or >= 1 of the zone's polythermal and
  temperate entities or >= 2 zone entities, (3) mean absolute site error worse with the 90 % interval
  above 0. Robust if unchanged with each massif removed and at f_ins 0.785 and 1.25; if not robust,
  FRH is central and H a full sensitivity (decision of 2026-10-13).

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/forcing_fair_test.py
       [out_dir] [--tag fair] [--allow-missing]
"""
import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from icetemp.calibration import thermal_structure as ts
from icetemp.calibration.data import DataHandler, fit_exclusions
from icetemp.calibration.validation import STEP5A_GLACIER_IDS

REPO = Path(__file__).resolve().parents[2]
BAYESCAL = REPO / 'data' / 'bayescal'
FIRNICE = 'training_runs/monthly/CentralEurope/PAST/firnice_temperature'
WEEK1_TIME = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week1/time_test/results.csv')
CROSSWALK = REPO / 'data' / 'glenglat_borehole_rgi7.csv'

# ── fixed before the runs ──────────────────────────────────────────────────────────────────
DECIDE = ['H', 'FRH']
REPORTED = ['HR', 'FRSH']
FREFS = [0.10, 0.25, 0.40, 0.55, 0.70, 0.85, 1.00]
FINS0, FADV0 = 1.017, 1.0
FINS_CHECK = [0.785, 1.25]
NODEEP = [('H', 0.55), ('FRH', 0.55)]
DEEP = 15.0
ZONE, ZONE_WIDE = (3005, 3495), (2950, 3550)    # entity band centre, m a.s.l.
MIN_POINTS = 3                                  # outlines with fewer points stay out of M1
MARGIN_M1 = 0.3
LOSS_OVERALL = dict(polythermal=2, temperate=2, total=4)
LOSS_ZONE = dict(warm=1, total=2)               # warm = polythermal or temperate
TOLS = ('0.05', 'hole', '0.5')                  # robustness bracket
TOL_PER_M = 0.0002                              # C per m depth, minimum tolerance
SEALED = {'01161'}
EXPECT = dict(entities=105, points=800, outlines=19, massifs=7,
              classes=dict(cold=83, polythermal=16, temperate=6), zone=dict(cold=13, polythermal=6, temperate=1))
MASSIFS = {
    'Mont Blanc': ['00773', '00774', '00777', '00786'],
    'Monte Rosa': ['01195', '01205', '01225', '01230'],
    'Mischabel-Saas': ['01112', '01127', '01150', '01161'],
    'Bernese Alps': ['01265', '01565', '02596', '02611'],
    'Valais west and Diablerets': ['01039', '01309'],
    'Titlis': ['02669'],
    'Ecrins': ['00245'],
    'Bernina': ['02103', '02216'],
    'Ortles': ['02821'],
    'Oetztal': ['03116'],
}
CLASSES = ('cold', 'polythermal', 'temperate')
N_BOOT, SEED = 4000, 20261013
# ───────────────────────────────────────────────────────────────────────────────────────────

MASSIF_OF = {g: m for m, gs in MASSIFS.items() for g in gs}


def run_name(tag, variant, fref, fins=FINS0, fadv=FADV0, nodeep=False):
    return f'{tag}_{variant}_r{fref:.2f}_i{fins:.3f}_a{fadv:.2f}' + ('_nodeep' if nodeep else '') + '_w0'


def col(variant, fref, fins=FINS0):
    return f'{variant}@{fref:.2f}' + ('' if fins == FINS0 else f'@{fins:.3f}')


# ── scored set ───────────────────────────────────────────────────────────────────────────
def scored_entities(ents, bh_table):
    """entity -> dict(points [(depth, T_obs)], rgi, band, massif, tolerance), points below 15 m."""
    excl = fit_exclusions(ents, ['colle_gnifetti_below_40m'])
    out = {}
    for g in ents:
        if not g.glacier_id or g.glacier_id in SEALED:
            continue
        idx = [i for i, d in enumerate(g.depths) if d > DEEP and (g.glacier_name, float(d)) not in excl]
        if idx:
            if g.glacier_id not in MASSIF_OF:
                raise SystemExit(f'{g.glacier_id} ({g.glacier_name}) is in no massif')
            out[g.glacier_name] = dict(idx=idx, depth=g.depths[idx], obs=g.T_obs[idx], rgi=g.glacier_id,
                                       band=g.elevation_band, massif=MASSIF_OF[g.glacier_id],
                                       tol=ts.borehole_tolerance(g.borehole_id, bh_table))
    return out


def classify_entity(T, depth, tol_mode, hole_tol):
    if tol_mode == 'hole':
        tols = np.maximum(hole_tol, TOL_PER_M * depth)
    else:
        tols = np.maximum(float(tol_mode), TOL_PER_M * depth)
    return ts.classify_profile([ts.classify_point(t, d, tl) for t, d, tl in zip(T, depth, tols)])


def points_frame(E, outputs):
    """Long table: one row per scored point, the model of every run as a column; rows finite in
    all decision runs."""
    rows = []
    for name, e in E.items():
        vals = {}
        for k, out in outputs.items():
            v = out.get(name)
            vals[k] = None if v is None else np.asarray(v[1], dtype=float)
        for j, i in enumerate(e['idx']):
            rows.append(dict(entity=name, rgi=e['rgi'], massif=e['massif'], band=e['band'],
                             depth=float(e['depth'][j]), obs=float(e['obs'][j]),
                             **{k: (np.nan if vals[k] is None or i >= len(vals[k]) else float(vals[k][i]))
                                for k in vals}))
    return pd.DataFrame(rows)


# ── measures ─────────────────────────────────────────────────────────────────────────────
def m1(df, c, min_points=MIN_POINTS, unit='rgi'):
    """RMSE with each outline (unit) weighted equally; outlines with fewer than min_points left out."""
    if df.empty:
        return np.nan
    g = ((df[c] - df.obs) ** 2).groupby(df[unit])
    mse = g.mean()[g.size() >= min_points]
    return float(np.sqrt(mse.mean())) if len(mse) else np.nan


def entity_classes(df, c, E):
    """entity -> {tol_mode: class} for column c (the model, or 'obs')."""
    out = {}
    for name, s in df.groupby('entity'):
        out[name] = {t: classify_entity(s[c].to_numpy(), s.depth.to_numpy(), t, E[name]['tol']) for t in TOLS}
    return out


def class_counts(obs_cls, mod_cls, names):
    """Correctly classified entities per measured class (classes at the hole's tolerance)."""
    out = {c: [0, 0] for c in CLASSES}
    for n in names:
        o = obs_cls[n]['hole']
        out[o][1] += 1
        out[o][0] += int(mod_cls[n]['hole'] == o)
    return out


def robust_changes(obs_cls, a_cls, b_cls, names):
    """Net robust changes from a to b: entities whose correctness changes the same way at every
    tolerance of the bracket. Returns per measured class (gains, losses)."""
    out = {c: [0, 0] for c in CLASSES}
    for n in names:
        ok_a = [a_cls[n][t] == obs_cls[n][t] for t in TOLS]
        ok_b = [b_cls[n][t] == obs_cls[n][t] for t in TOLS]
        c = obs_cls[n]['hole']
        if all(ok_b) and not any(ok_a):
            out[c][0] += 1
        elif all(ok_a) and not any(ok_b):
            out[c][1] += 1
    return out


def m2(obs_cls, mod_cls, names):
    cc = class_counts(obs_cls, mod_cls, names)
    rec = [cc[c][0] / cc[c][1] for c in CLASSES if cc[c][1]]
    return float(np.mean(rec)) if rec else np.nan


def lomo(P, cols_by_fref, massifs=None):
    """Out-of-sample column: each massif scored at the f_ref with the lowest M1 on the others
    (ties: lowest M1 on all massifs). Returns (predictions, {massif: f_ref})."""
    massifs = massifs or sorted(P.massif.unique())
    pred = pd.Series(np.nan, index=P.index)
    chosen = {}
    allm = {fr: m1(P, c) for fr, c in cols_by_fref.items()}
    base = P.massif.str.split('#').str[0]      # resampled copies of one massif leave together
    for m in massifs:
        tr = P[base != m.split('#')[0]]
        sc = {fr: m1(tr, c) for fr, c in cols_by_fref.items()}
        best = min(sc, key=lambda fr: (round(sc[fr], 9), allm[fr]))
        chosen[m] = best
        sel = P.massif == m
        pred[sel] = P.loc[sel, cols_by_fref[best]]
    return pred, chosen


# ── warming over time ────────────────────────────────────────────────────────────────────
def warming_pairs(run_names):
    """Week-1 repeat pairs, outlines from the crosswalk, 01161 and extra sources out; the modelled
    change of every run at the earlier hole's outline and band."""
    import time_test as tt
    res = pd.read_csv(WEEK1_TIME)
    pairs = pd.read_csv(WEEK1_TIME.with_name('pairs.csv'))
    xw = pd.read_csv(CROSSWALK, dtype={'rgi7': str}).set_index('borehole_id').rgi7
    pairs = pairs.merge(res[['pair_id', 'dT_meas']], on='pair_id')
    pairs = pairs[~pairs.extra_source.astype(bool)].copy()
    pairs['outline1'] = pairs.hole1.map(xw)
    pairs['outline2'] = pairs.hole2.map(xw)
    pairs = pairs[~pairs.outline1.isin(SEALED) & ~pairs.outline2.isin(SEALED)].copy()
    bh, prof, meas = tt.load_glenglat()
    p = tt.eligible_profiles(prof, meas)
    months = {(int(r.borehole_id), int(r.id)): r.months for r in p.itertuples()}
    mfs = {r: tt.ModelFiles(r) for r in run_names}
    rows = []
    for r in pairs.itertuples():
        e1 = [months[(r.hole1, int(k))] for k in str(r.profiles1).split()]
        e2 = [months[(r.hole2, int(k))] for k in str(r.profiles2).split()]
        use = None
        for gid, band in ((r.outline1, r.band1), (r.outline2, r.band2)):
            if isinstance(gid, str):
                near, path = mfs[run_names[0]].file_for(gid, band)
                if path is not None and abs(near - band) <= tt.MAX_BAND_GAP:
                    use = (gid, band)
                    break
        row = dict(pair_id=r.pair_id, site=r.site, depth=r.depth, dT_meas=r.dT_meas,
                   outline=None if use is None else use[0],
                   massif=None if use is None else MASSIF_OF.get(use[0]))
        if use is not None:
            for run in run_names:
                _, path = mfs[run].file_for(*use)
                row[run] = (tt.model_at(mfs[run], path, r.depth, e2) - tt.model_at(mfs[run], path, r.depth, e1)
                            if path is not None else np.nan)
        rows.append(row)
    return pd.DataFrame(rows)


def warming_measures(W, cmap, rng=None):
    """Site means of measured and modelled change; cmap: massif -> run column of the variant.
    Returns mean absolute site error and ratio of means (and their site-resampled 90 % ranges)."""
    W = W[W.massif.notna()].copy()
    W['mod'] = [r[cmap[r.massif]] if r.massif in cmap else np.nan for _, r in W.iterrows()]
    W = W[np.isfinite(W['mod']) & np.isfinite(W.dT_meas)]
    S = W.groupby('site').agg(meas=('dT_meas', 'mean'), mod=('mod', 'mean'))
    out = dict(n_sites=len(S), mase=float(np.mean(np.abs(S['mod'] - S.meas))),
               ratio=float(S['mod'].mean() / S.meas.mean()) if S.meas.mean() != 0 else np.nan,
               sites=S.round(3).to_dict(orient='index'))
    return out, S


# ── decision ─────────────────────────────────────────────────────────────────────────────
def decide_A(res_H, res_F, chg_all, chg_zone, warm_diff_ci):
    """Condition by condition; True = FRH adopted."""
    c1 = (res_F['M1'] - res_H['M1'] > MARGIN_M1) or (res_F['zone_M1'] - res_H['zone_M1'] > MARGIN_M1)
    net = {c: chg_all[c][1] - chg_all[c][0] for c in CLASSES}
    tot = sum(net.values())
    znet = {c: chg_zone[c][1] - chg_zone[c][0] for c in CLASSES}
    c2 = (net['polythermal'] >= LOSS_OVERALL['polythermal'] or net['temperate'] >= LOSS_OVERALL['temperate']
          or tot >= LOSS_OVERALL['total'] or znet['polythermal'] + znet['temperate'] >= LOSS_ZONE['warm']
          or sum(znet.values()) >= LOSS_ZONE['total'])
    c3 = warm_diff_ci is not None and warm_diff_ci[0] > 0
    return dict(fails_M1=bool(c1), fails_classes=bool(c2), fails_warming=bool(c3),
                net_losses=net, zone_net_losses=znet, adopt_FRH=not (c1 or c2 or c3))


def evaluate(P, E, obs_cls, W, fins=FINS0, drop_massif=None):
    """Out-of-sample measures of H and FRH and the decision A (optionally one massif removed)."""
    Q = P if drop_massif is None else P[P.massif != drop_massif]
    zone = Q.band.between(*ZONE)
    names = sorted(Q.entity.unique())
    znames = sorted(Q[zone].entity.unique())
    out = {}
    for v in DECIDE:
        pred, chosen = lomo(Q, {fr: col(v, fr, fins) for fr in FREFS})
        Qv = Q.assign(oos=pred)
        cls = entity_classes(Qv, 'oos', E)
        cmap = {m: run_name(TAG, v, fr, fins) for m, fr in chosen.items()}
        wm, S = warming_measures(W, cmap) if W is not None else ({}, None)
        out[v] = dict(chosen=chosen, M1=m1(Qv, 'oos'), zone_M1=m1(Qv[zone], 'oos'), cls=cls,
                      counts=class_counts(obs_cls, cls, names), zone_counts=class_counts(obs_cls, cls, znames),
                      M2=m2(obs_cls, cls, names), warming=wm, site_table=S)
    chg = robust_changes(obs_cls, out['H']['cls'], out['FRH']['cls'], names)
    zchg = robust_changes(obs_cls, out['H']['cls'], out['FRH']['cls'], znames)
    ci = None
    if W is not None and out['H']['site_table'] is not None:
        ci = warming_diff_ci(out['H']['site_table'], out['FRH']['site_table'])
    dec = decide_A(out['H'], out['FRH'], chg, zchg, ci)
    return out, dec, chg, zchg, ci


def warming_diff_ci(SH, SF):
    """90 % interval of MASE(FRH) - MASE(H) from resampling sites."""
    sites = sorted(set(SH.index) & set(SF.index))
    if not sites:
        return None
    rng = np.random.default_rng(SEED)
    eH = np.abs(SH.loc[sites, 'mod'] - SH.loc[sites, 'meas']).to_numpy()
    eF = np.abs(SF.loc[sites, 'mod'] - SF.loc[sites, 'meas']).to_numpy()
    d = [np.mean(eF[k] - eH[k]) for k in (rng.integers(0, len(sites), len(sites)) for _ in range(N_BOOT))]
    return [float(np.percentile(d, 5)), float(np.percentile(d, 95))]


def bootstrap(P, E, obs_cls, n=N_BOOT):
    """Paired massif resampling: drawn copies labelled separately, f_ref re-chosen in each draw.
    Returns 90 % intervals of the M1 differences and the share of draws in which conditions 1 and 2
    of A give the same outcome as on the full set (condition 3 has its own site-resampled interval)."""
    rng = np.random.default_rng(SEED)
    ms = sorted(P.massif.unique())
    cols = {v: {fr: col(v, fr) for fr in FREFS} for v in DECIDE}
    CL = {c: entity_classes(P, c, E) for v in DECIDE for c in cols[v].values()}
    ent_m = P.groupby('entity').massif.first()
    zone_e = set(P[P.band.between(*ZONE)].entity)
    by_m = {m: P[P.massif == m] for m in ms}
    full = None
    d1, dz, outcome = [], [], []
    for b in range(n + 1):
        pick = ms if b == 0 else list(rng.choice(ms, len(ms)))
        Q = pd.concat([by_m[m].assign(massif=f'{m}#{j}', rgi=by_m[m].rgi + f'#{j}', entity=by_m[m].entity + f'#{j}')
                       for j, m in enumerate(pick)], ignore_index=True)
        z = Q.band.between(*ZONE)
        r, picked_cls = {}, {}
        for v in DECIDE:
            Qr = Q.assign(base=Q.entity.str.split('#').str[0])
            pred, chosen = lomo(Qr, cols[v])
            Qv = Qr.assign(oos=pred)
            r[v] = (m1(Qv, 'oos'), m1(Qv[z], 'oos'))
            picked_cls[v] = {e: CL[cols[v][chosen[Qr.loc[Qr.entity == e, 'massif'].iloc[0]]]][e.split('#')[0]]
                             for e in Qr.entity.unique()}
        names = sorted(picked_cls['H'])
        o = {e: obs_cls[e.split('#')[0]] for e in names}
        chg = robust_changes(o, picked_cls['H'], picked_cls['FRH'], names)
        zchg = robust_changes(o, picked_cls['H'], picked_cls['FRH'], [e for e in names if e.split('#')[0] in zone_e])
        dec = decide_A(dict(M1=r['H'][0], zone_M1=r['H'][1]), dict(M1=r['FRH'][0], zone_M1=r['FRH'][1]), chg, zchg, None)
        if b == 0:
            full = dec['adopt_FRH']
            continue
        d1.append(r['FRH'][0] - r['H'][0])
        dz.append(r['FRH'][1] - r['H'][1])
        outcome.append(dec['adopt_FRH'] == full)
    q = lambda x: [float(np.nanpercentile(x, p)) for p in (5, 50, 95)]
    return dict(dM1=q(d1), dM1_zone=q(dz), share_same_outcome_conditions_1_2=float(np.mean(outcome)),
                outcome_full_conditions_1_2=full)


# ── reported extras ──────────────────────────────────────────────────────────────────────
def oos_column(P, v, fins=FINS0):
    pred, chosen = lomo(P, {fr: col(v, fr, fins) for fr in FREFS})
    return P.assign(oos=pred), chosen


def layer_order(s, c, E):
    """Polythermal profile structure by depth: 'cold over temperate', 'temperate over cold' or 'mixed'."""
    s = s.sort_values('depth')
    lab = [ts.classify_point(t, d, max(E[s.entity.iloc[0]]['tol'], TOL_PER_M * d))
           for t, d in zip(s[c], s.depth)]
    if len(set(lab)) < 2:
        return lab[0] if lab else None
    first_t = lab.index('temperate')
    first_c = lab.index('cold')
    tail = lab[max(first_t, first_c):]
    if first_c < first_t and all(x == 'temperate' for x in tail):
        return 'cold over temperate'
    if first_t < first_c and all(x == 'cold' for x in tail):
        return 'temperate over cold'
    return 'mixed'


def frozen_bed_table(run_cols, chosen_by_v, E):
    """Holes that reach the bed, of known depth, deepest reading at most 5 m above the bottom: the
    measured deepest reading and the model at that depth and in those months (each variant at the
    f_ref chosen for the hole's massif)."""
    import time_test as tt
    bh, prof, meas = tt.load_glenglat()
    xw = pd.read_csv(CROSSWALK, dtype={'rgi7': str}).set_index('borehole_id').rgi7
    p = tt.eligible_profiles(prof, meas)
    holes = bh[(bh.to_bed.astype(str).str.lower() == 'true') & bh.depth.notna()]
    rows = []
    mfs = {}
    for h in holes.itertuples():
        gid = xw.get(h.id)
        if gid is None or gid in SEALED or gid not in MASSIF_OF:
            continue
        ph = p[p.borehole_id == h.id]
        if ph.empty:
            continue
        last = ph.sort_values('t_mid').iloc[-1]
        mm = meas[(meas.borehole_id == h.id) & (meas.profile_id == last.id)]
        if mm.empty:
            continue
        dmax = float(mm.depth.max())
        if h.depth - dmax > 5.0:
            continue
        T = float(mm.loc[mm.depth.idxmax(), 'temperature'])
        tol = max(ts.DEFAULT_TOLERANCE if pd.isna(h.temperature_uncertainty) else float(h.temperature_uncertainty),
                  TOL_PER_M * dmax)
        pmp = ts.pressure_melting_point(dmax)
        row = dict(hole=h.id, glacier=h.glacier_name, rgi=gid, massif=MASSIF_OF[gid], depth=dmax, T_meas=T,
                   frozen_meas=T < pmp - tol, deep_cold=T < -2.0)
        for v, chosen in chosen_by_v.items():
            fr = chosen.get(MASSIF_OF[gid])
            if fr is None:
                continue
            rn = run_cols[(v, fr)]
            mfs.setdefault(rn, tt.ModelFiles(rn))
            near, path = mfs[rn].file_for(gid, tt.band_of(h.elevation))
            if path is None or abs(near - tt.band_of(h.elevation)) > tt.MAX_BAND_GAP:
                continue
            Tm = tt.model_at(mfs[rn], path, dmax, [last.months])
            row[f'T_{v}'] = Tm
            row[f'frozen_{v}'] = bool(np.isfinite(Tm) and Tm < pmp - tol)
        rows.append(row)
    return pd.DataFrame(rows)


def simple_rules(P, E, obs_cls):
    """Elevation rules in the same folds: a straight line of temperature on band elevation (M1) and
    two elevation thresholds for the classes (cold above, temperate below, polythermal between)."""
    pred = pd.Series(np.nan, index=P.index)
    cls = {}
    ent = P.groupby('entity').agg(band=('band', 'first'), massif=('massif', 'first'))
    bands = np.unique(ent.band)
    for m in sorted(P.massif.unique()):
        tr, te = P[P.massif != m], P.massif == m
        b, a = np.polyfit(tr.band, tr.obs, 1)
        pred[te] = np.minimum(a + b * P.loc[te, 'band'], 0.0)
        etr = ent[ent.massif != m]
        best, best_n = None, -1
        for lo in bands:
            for hi in bands[bands >= lo]:
                lab = np.where(etr.band > hi, 'cold', np.where(etr.band < lo, 'temperate', 'polythermal'))
                n = sum(lab[i] == obs_cls[e]['hole'] for i, e in enumerate(etr.index))
                if n > best_n:
                    best, best_n = (lo, hi), n
        for e in ent[ent.massif == m].index:
            bb = ent.loc[e, 'band']
            c = 'cold' if bb > best[1] else ('temperate' if bb < best[0] else 'polythermal')
            cls[e] = {t: c for t in TOLS}
    Pv = P.assign(oos=pred)
    names = sorted(P.entity.unique())
    return dict(M1=m1(Pv, 'oos'), counts={c: f'{a} of {b}' for c, (a, b) in class_counts(obs_cls, cls, names).items()},
                M2=m2(obs_cls, cls, names))


def extras(P, E, obs_cls, ents, outputs, W, chosen_by_v):
    """Reported, not decisive."""
    out = {}
    names = sorted(P.entity.unique())
    step5a = P.rgi.isin(STEP5A_GLACIER_IDS)
    for v in DECIDE:
        Pv, _ = oos_column(P, v)
        cls = entity_classes(Pv, 'oos', E)
        r = dict(M1_all_outlines=m1(Pv, 'oos', min_points=1),
                 pooled_rmse=float(np.sqrt(np.mean((Pv.oos - Pv.obs) ** 2))), pooled_bias=float(np.mean(Pv.oos - Pv.obs)),
                 zone_wide_M1=m1(Pv[Pv.band.between(*ZONE_WIDE)], 'oos'),
                 zone_wide_counts={c: f'{a} of {b}' for c, (a, b) in class_counts(
                     obs_cls, cls, sorted(Pv[Pv.band.between(*ZONE_WIDE)].entity.unique())).items()},
                 step5a_M1=m1(Pv[step5a], 'oos'), others_M1=m1(Pv[~step5a], 'oos'))
        warm = lambda k: k in ('polythermal', 'temperate')
        tw = [(warm(obs_cls[n]['hole']), warm(cls[n]['hole'])) for n in names]
        wf = np.mean([m for o, m in tw if o]) if any(o for o, _ in tw) else np.nan
        ck = np.mean([not m for o, m in tw if not o]) if any(not o for o, _ in tw) else np.nan
        r['two_class_skill'] = float(0.5 * (wf + ck))
        em = Pv.groupby('entity').agg(obs=('obs', 'mean'), mod=('oos', 'mean'))
        r['rank_agreement'] = float(stats.spearmanr(em['obs'], em['mod']).correlation)
        poly = [n for n in names if obs_cls[n]['hole'] == 'polythermal']
        r['polythermal_layer_order'] = {n: dict(measured=layer_order(Pv[Pv.entity == n], 'obs', E),
                                                model=layer_order(Pv[Pv.entity == n], 'oos', E)) for n in poly}
        if W is not None:
            cmap = {m: run_name(TAG, v, fr) for m, fr in chosen_by_v[v].items()}
            Wv = W[W.massif.notna()].copy()
            Wv['mod'] = [rr[cmap[rr.massif]] if rr.massif in cmap else np.nan for _, rr in Wv.iterrows()]
            Wv = Wv[np.isfinite(Wv['mod']) & np.isfinite(Wv.dT_meas)]
            x, y = Wv.dT_meas.to_numpy(), Wv['mod'].to_numpy()
            g = Wv.site == 'Grenzgletscher 4250 m'
            r['warming_slope'] = float(np.sum(x * y) / np.sum(x * x)) if len(x) else np.nan
            r['warming_slope_without_grenz4250'] = (float(np.sum(x[~g] * y[~g]) / np.sum(x[~g] ** 2))
                                                    if (~g).any() else np.nan)
        out[v] = r
    # Colle Gnifetti same-site gap (01230 minus 01225, 10-40 m, bands >= 4400), dry firn (>= 4400 m,
    # 10-20 m) and wet/dry ranking at 10-20 m, each variant at its most often chosen f_ref
    for v in DECIDE:
        fr = max(set(chosen_by_v[v].values()), key=list(chosen_by_v[v].values()).count)
        o = outputs[col(v, fr)]
        cg, dry, wd = {'01225': [], '01230': []}, [], []
        for g in ents:
            if g.glacier_id in SEALED or g.glacier_name not in o or o[g.glacier_name] is None:
                continue
            mod = np.asarray(o[g.glacier_name][1], float)[:len(g.depths)]
            d, t = g.depths, g.T_obs
            if g.glacier_id in cg and g.elevation_band >= 4400:
                k = (d >= 10) & (d <= 40) & np.isfinite(mod)
                cg[g.glacier_id] += list(mod[k] - t[k])
            k = (d >= 10) & (d <= 20) & np.isfinite(mod)
            if k.any():
                if g.elevation_band >= 4400:
                    dry += list(mod[k] - t[k])
                wd.append((float(np.mean(t[k])) > -1.0, float(np.mean(mod[k]))))
        out[v]['f_ref_for_site_checks'] = fr
        out[v]['cg_gap_01230_minus_01225'] = (float(np.mean(cg['01230']) - np.mean(cg['01225']))
                                             if cg['01230'] and cg['01225'] else np.nan)
        out[v]['dry_firn_bias_4400'] = float(np.mean(dry)) if dry else np.nan
        lab = np.array([a for a, _ in wd]); val = np.array([b for _, b in wd])
        out[v]['wet_dry_auc_10_20m'] = (float(stats.mannwhitneyu(val[lab], val[~lab]).statistic / (lab.sum() * (~lab).sum()))
                                        if lab.any() and (~lab).any() else np.nan)
    out['simple_rules'] = simple_rules(P, E, obs_cls)
    return out


def bed_types(run, gids):
    """Glacier type per glacier in the last year from the bed log: which bands are frozen."""
    out = {}
    for g in gids:
        f = BAYESCAL / run / FIRNICE / f'firnice_bed_{g}.dat'
        if not f.exists():
            continue
        b = pd.read_csv(f, sep=r'\s+', comment='#', header=None,
                        names=['elev', 'year', 'z', 'T', 'Tmin', 'Tmax', 'pmp', 'zb', 'Tb'])
        last = b[b.year == b.year.max()].sort_values('elev')
        frozen = (last['T'] < last.pmp - 0.1).to_numpy()
        if len(frozen) == 0:
            continue
        if frozen.all():
            t = 'cold-based'
        elif not frozen.any():
            t = 'temperate-based'
        else:
            k = max(1, len(frozen) // 4)
            lo, hi = frozen[:k].mean(), frozen[-k:].mean()
            t = ('frozen front, temperate above' if lo > 0.5 >= hi else
                 'frozen upper part, temperate below' if hi > 0.5 >= lo else 'mixed')
        out[g] = dict(type=t, frozen_share=float(frozen.mean()), n_bands=int(len(frozen)), year=int(b.year.max()))
    return out


def main():
    global TAG
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('out_dir', nargs='?', default='/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week2/fair_test')
    ap.add_argument('--tag', default='fair')
    ap.add_argument('--allow-missing', action='store_true', help='score what exists (smoke checks only)')
    ap.add_argument('--no-boot', action='store_true')
    args = ap.parse_args()
    TAG = args.tag
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    dh = DataHandler(region='CentralEurope')
    dh.load()
    bh_table = ts.load_borehole_table()
    E = scored_entities(dh.calibration_glaciers, bh_table)

    wanted = {col(v, fr): run_name(TAG, v, fr) for v in DECIDE + REPORTED for fr in FREFS}
    wanted |= {col(v, fr, fi): run_name(TAG, v, fr, fi) for v in DECIDE for fr in FREFS for fi in FINS_CHECK}
    outputs, missing = {}, []
    for c, rn in wanted.items():
        f = BAYESCAL / rn / 'training_runs' / 'design_0000.output.pkl'
        if f.exists():
            outputs[c] = pickle.load(open(f, 'rb'))
        else:
            missing.append(rn)
    if missing and not args.allow_missing:
        raise SystemExit(f'{len(missing)} runs missing, e.g. {missing[:3]}')

    P = points_frame(E, outputs)
    data_counts = dict(entities=P.entity.nunique(), points=len(P), outlines=P.rgi.nunique(),
                       massifs=P.massif.nunique())
    for k, v in data_counts.items():
        if v != EXPECT[k]:
            raise SystemExit(f'scored set differs from the pre-registration: {k} {v} != {EXPECT[k]}')
    decision_cols = [c for c in outputs if c.split('@')[0] in DECIDE]
    finite = np.isfinite(P[decision_cols]).all(axis=1) if decision_cols else np.ones(len(P), bool)
    n_dropped = int((~finite).sum())
    P = P[finite].reset_index(drop=True)
    obs_cls = entity_classes(P, 'obs', E)
    counts = dict(entities=P.entity.nunique(), points=len(P), outlines=P.rgi.nunique(), massifs=P.massif.nunique())
    obs_c = pd.Series({n: c['hole'] for n, c in obs_cls.items()})
    zone_names = sorted(P[P.band.between(*ZONE)].entity.unique())
    counts['classes'] = obs_c.value_counts().to_dict()
    counts['zone'] = obs_c[zone_names].value_counts().to_dict()
    for k in ('classes', 'zone'):
        if n_dropped == 0 and counts[k] != EXPECT[k]:
            raise SystemExit(f'measured classes differ from the pre-registration: {k} {counts[k]} != {EXPECT[k]}')
    res = dict(counts=counts, expected=EXPECT, n_dropped_nonfinite=n_dropped, missing_runs=missing)

    W = None
    if not missing:
        W = warming_pairs(sorted(set(wanted.values())))
        W.to_csv(out / 'warming_pairs.csv', index=False)
    main_out, dec, chg, zchg, ci = evaluate(P, E, obs_cls, W)
    res['A'] = dec
    res['A']['warming_mase_diff_ci90'] = ci
    for v in DECIDE:
        o = main_out[v]
        res[v] = dict(chosen_fref=o['chosen'], M1=o['M1'], zone_M1=o['zone_M1'], M2=o['M2'],
                      counts={c: f'{a} of {b}' for c, (a, b) in o['counts'].items()},
                      zone_counts={c: f'{a} of {b}' for c, (a, b) in o['zone_counts'].items()},
                      warming={k: o['warming'].get(k) for k in ('n_sites', 'mase', 'ratio')},
                      at_bound=sum(fr in (FREFS[0], FREFS[-1]) for fr in o['chosen'].values()))
    res['robust_changes'] = dict(overall=chg, zone=zchg)

    # robustness: each massif removed, and the f_ins check
    rob = {}
    for m in sorted(P.massif.unique()):
        _, d, *_ = evaluate(P, E, obs_cls, W, drop_massif=m)
        rob[f'without {m}'] = d['adopt_FRH']
    for fi in FINS_CHECK:
        if all(col(v, fr, fi) in P for v in DECIDE for fr in FREFS):
            _, d, *_ = evaluate(P, E, obs_cls, W, fins=fi)
            rob[f'f_ins {fi}'] = d['adopt_FRH']
    res['robustness'] = rob
    robust = all(v == dec['adopt_FRH'] for v in rob.values())
    res['A']['robust'] = robust
    res['A']['outcome'] = ('FRH adopted' if dec['adopt_FRH'] else 'H stays') if robust else \
        'not robust: FRH central, H as a full sensitivity (decision of 2026-10-13)'
    if not args.no_boot:
        res['bootstrap'] = bootstrap(P, E, obs_cls)

    # reported: the other variants at their own out-of-sample f_ref
    for v in REPORTED:
        if all(col(v, fr) in P for fr in FREFS):
            pred, chosen = lomo(P, {fr: col(v, fr) for fr in FREFS})
            Pv = P.assign(oos=pred)
            cls = entity_classes(Pv, 'oos', E)
            res[v] = dict(chosen_fref=chosen, M1=m1(Pv, 'oos'), zone_M1=m1(Pv[Pv.band.between(*ZONE)], 'oos'),
                          M2=m2(obs_cls, cls, sorted(P.entity.unique())),
                          counts={c: f'{a} of {b}' for c, (a, b) in
                                  class_counts(obs_cls, cls, sorted(P.entity.unique())).items()})
    res['flags'] = dict(fref_at_bound={v: res[v]['at_bound'] for v in DECIDE},
                        stop_sign={v: res[v]['warming'].get('ratio') for v in DECIDE})
    gids = sorted(P.rgi.unique())
    res['bed_types'] = {v: bed_types(run_name(TAG, v, 0.55, nodeep=True), gids) for v, _ in NODEEP}
    if not missing:
        chosen_by_v = {v: main_out[v]['chosen'] for v in DECIDE}
        res['reported'] = extras(P, E, obs_cls, dh.calibration_glaciers, outputs, W, chosen_by_v)
        FB = frozen_bed_table({(v, fr): run_name(TAG, v, fr) for v in DECIDE for fr in FREFS}, chosen_by_v, E)
        FB.to_csv(out / 'frozen_bed.csv', index=False)
        res['frozen_bed'] = {v: dict(
            frozen_found=f"{int((FB.frozen_meas & FB[f'frozen_{v}']).sum())} of {int(FB.frozen_meas.sum())}",
            frozen_deep_cold_found=f"{int((FB.frozen_meas & FB.deep_cold & FB[f'frozen_{v}']).sum())} of "
                                   f"{int((FB.frozen_meas & FB.deep_cold).sum())}",
            temperate_called_frozen=int((~FB.frozen_meas & FB[f'frozen_{v}']).sum()))
            for v in DECIDE if f'frozen_{v}' in FB}
    P.to_csv(out / 'points.csv', index=False)
    (out / 'fair_test.json').write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps({k: res[k] for k in ('counts', 'A', 'robustness')}, indent=1, default=str))


TAG = 'fair'

if __name__ == '__main__':
    main()
