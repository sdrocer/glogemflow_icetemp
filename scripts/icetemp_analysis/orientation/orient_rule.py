"""Stage 1 of the orientation test: rule B (Step 5 plan, revised and approved 21 Oct 2026, step 1b point d), fixed and
committed before any of its runs was read.

Runs, tag o1 (week3/launch_orient_stage1.sh; GloGEM branch orientation-factor): FRHS (firnice_orient 'n') and FRHS
with the orientation factor in three strengths, C (central), W (weaker, c 20), S (stronger, sunny spells)
(write_orient_tables.py), at f_ref 0.03, 0.10, 0.25, 0.40, 0.55, 0.70, 0.85, 1.00 (C also at 1.30, diagnostic),
f_ins 1.017, f_adv 1.0; CentralEurope_fair without 01161. Depth bins pooled per borehole and year (1135 points below
15 m; rescore_hole_year.py). The sunshine of a band is the mean June-August sunshine ratio of its holes (week 1,
orientation_investigation/evidence/bands.csv).

Rule B. The orientation factor is adopted if, for C and for W,
 (1) the shade error inside sites is at least halved, |b_V| <= 0.5 |b_FRHS|. b is the OLS slope of the band bias
     (bands at or above 3000 m; each variant at its own best common f_ref, linear between its runs, 0.03-1.00) on the
     band's sunshine ratio, with site fixed effects and elevation, in C per +0.1 of sunshine. It must hold on all
     sites, and with the Taconnaz and Lys sites left out of both the f_ref choice and the scoring;
 (2) there is no harm overall, as for S (stage 2): M1 (out of sample, f_ref chosen on the grid leaving one massif
     out) not up by more than 0.2 C, overall and in the zone; M1 with the outlines 00774 and 01205 left out of the
     choice and of the scoring not up by more than 0.05 C; no net loss of sites with the right regime, overall and in
     the zone (firn_batch_rule.site_changes).
Single bands or beds that get worse do not veto; they are listed with their change.
Reported, not deciding: S against (1) and (2); the 90 % interval of b (wild cluster bootstrap over sites, 9999
draws); b per massif and with each site left out; the within-site RMS (band biases about their site mean, sites with
at least 2 high bands); the 11 sunny/shaded pairs at 5 sites (one vote per site); M1 above 3800 m; rule 1a as
approved (the in-sample best f_ref of Mont Blanc and of Monte Rosa within 0.2; Taconnaz Glacier cools at f_ref 0.85
and Monte Rosa does not cool by more than 0.2 C); Taconnaz Glacier 3685-3785 m; regime classes; frozen beds; the
time test; every band at depth.

Changes after the max-effort review of 21 Oct, made before any run was read: the within-site measures ((1), the
interval, the within-site RMS, the pairs) leave out entities whose own F differs by more than 0.1 from the F applied
to their model band (write_orient_tables.py applied.csv: Glacier de Taconnaz@4295, both Mont Blanc summit holes),
so the Taconnaz 4285/4295 pair drops; sunny and shaded are told by the applied F; a pair counts only if the model
moves at least 0.2 C closer or farther (thermistor precision), else it is a tie; the two Colle Gnifetti pairs
(01225 against 01230, different climate cells and the 30 mm floor) are flagged.

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/orientation/orient_rule.py [OUT_DIR]
           [--stand-in]    (--stand-in: every variant read from the stage-2 FRHS runs, a test of the scorer only)
           [--frhs-tag T --tag T2 --suffix S]   (the re-check after the burial test: FRHS from tag T, the
           orientation runs from tag T2, both with the run-name suffix S, e.g. b1 o2 _bM; same rule)
"""
import argparse
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import firn_batch_rule as fbr  # noqa: E402
import forcing_fair_test as ft  # noqa: E402
from icetemp.calibration import thermal_structure as ts  # noqa: E402
from icetemp.calibration.data import DataHandler, site_keys  # noqa: E402

W3 = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3')
BANDS = W3 / 'orientation_investigation/evidence/bands.csv'
TABLES = W3 / 'orient_tables'
TAG = 'o1'
VARS = ('FRHS', 'C', 'W', 'S')
DECIDING = ('C', 'W')
GRID = (0.03, 0.10, 0.25, 0.40, 0.55, 0.70, 0.85, 1.00)
EXTRA = {'C': (1.30,)}
FINE = np.round(np.arange(0.03, 1.0001, 0.01), 2)
HIGH = 3000
HALVE, HARMS, HARMS_STRICT = 0.5, 0.2, 0.05
GENERATING_SITES = ('Taconnaz', 'Lysgletscher')
GENERATING_OUTLINES = ('00774', '01205')
PAIRS = [('Grenzgletscher@4405', 'Grenzgletscher@4335'),
         ('Grenzgletscher@4305~01195', 'Grenzgletscher@4295'),
         ('Lysgletscher@4035', 'Lysgletscher@4085'), ('Lysgletscher@4115', 'Lysgletscher@4085'),
         ('Glacier de Taconnaz@4285', 'Glacier de Taconnaz@4215'), ('Glacier de Taconnaz@4285', 'Glacier de Taconnaz@4295'),
         ('Glacier des Bossons@4285', 'Glacier des Bossons@4275'), ('Glacier des Bossons@4285', 'Glacier des Bossons@4295'),
         ('Grenzgletscher@4465~01225', 'Grenzgletscher@4465~01230'), ('Grenzgletscher@4445~01225', 'Grenzgletscher@4445~01230'),
         ('Jungfraufirn@3535', 'Jungfraufirn@3635')]
GROUP_C = ('Taconnaz Glacier@3685', 'Taconnaz Glacier@3785')
PAIR_MIN = 0.2
FLAGGED_PAIRS = {('Grenzgletscher@4465~01225', 'Grenzgletscher@4465~01230'),
                 ('Grenzgletscher@4445~01225', 'Grenzgletscher@4445~01230')}
WEBB = np.array([-np.sqrt(1.5), -1.0, -np.sqrt(0.5), np.sqrt(0.5), 1.0, np.sqrt(1.5)])


NAMES = dict(frhs_tag=TAG, tag=TAG, suffix='')


def run_of(v, fr, stand_in=False):
    if stand_in:
        return ft.run_name('s2', 'FRHS', fr)
    if v == 'FRHS':
        return ft.run_name(NAMES['frhs_tag'], 'FRHS', fr).replace('_w0', f"{NAMES['suffix']}_w0")
    return ft.run_name(NAMES['tag'], 'FRHS', fr).replace('_w0', f"_o{v}{NAMES['suffix']}_w0")


def grid_of(v, stand_in=False):
    g = tuple(fr for fr in GRID if not (stand_in and fr < 0.10))
    return g + (() if stand_in else EXTRA.get(v, ()))


# ── statistics ───────────────────────────────────────────────────────────────────────────
def ols(y, X):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return beta, y - X @ beta


def wild_cluster(y, X, cl, j, B=9999, seed=1):
    """Restricted wild cluster bootstrap (Webb weights) for beta_j = 0: (beta_j, p, 90 % interval)."""
    rng = np.random.default_rng(seed)
    beta, res = ols(y, X)
    keep = [k for k in range(X.shape[1]) if k != j]
    b0, r0 = ols(y, X[:, keep])
    fit0 = X[:, keep] @ b0
    groups = np.unique(cl)
    gi = np.searchsorted(groups, cl)
    XtX = np.linalg.pinv(X.T @ X)

    def tstat(yy):
        bb, rr = ols(yy, X)
        meat = np.zeros((X.shape[1], X.shape[1]))
        for g in range(len(groups)):
            s = X[gi == g].T @ rr[gi == g]
            meat += np.outer(s, s)
        V = XtX @ meat @ XtX
        return bb[j] / np.sqrt(V[j, j]) if V[j, j] > 0 else np.nan
    t_obs = tstat(y)
    tb, bu = np.empty(B), np.empty(B)
    for k in range(B):
        w = rng.choice(WEBB, len(groups))[gi]
        tb[k] = tstat(fit0 + r0 * w)
        bu[k] = ols(X @ beta + res * w, X)[0][j]
    return float(beta[j]), float(np.mean(np.abs(tb) >= abs(t_obs))), (float(np.quantile(bu, .05)), float(np.quantile(bu, .95)))


def design(bb):
    return np.column_stack([bb.ratio_jja - 0.95, (bb.band - 4000) / 100.] +
                           [(bb.site == g).to_numpy(float) for g in sorted(bb.site.unique())])


def slope(bb):
    return float(ols(bb.b.to_numpy(), design(bb))[0][0] * 0.1)


def within_rmse(bb):
    q = bb[bb.groupby('site').entity.transform('size') >= 2]
    dev = q.b - q.groupby('site').b.transform('mean')
    return float(np.sqrt((dev ** 2).mean())), int(q.site.nunique())


# ── the model at any f_ref, linear between a variant's runs ───────────────────────────────
def at_f(P, v, f, grid):
    g = np.array(grid)
    k = int(np.clip(np.searchsorted(g, f) - 1, 0, len(g) - 2))
    w = (f - g[k]) / (g[k + 1] - g[k])
    return (1 - w) * P[ft.col(v, g[k])] + w * P[ft.col(v, g[k + 1])]


def best_common(P, v, grid, sel):
    sc = {f: ft.m1(P[sel].assign(m=at_f(P[sel], v, f, grid)), 'm') for f in FINE}
    return min(sc, key=sc.get)


def band_bias(P, v, grid, sel):
    f = best_common(P, v, grid, sel)
    Q = P[sel].assign(mod=at_f(P[sel], v, f, grid))
    bb = Q.assign(b=Q['mod'] - Q.obs).groupby(['site', 'massif', 'entity', 'band'], as_index=False).agg(
        b=('b', 'mean'), mod=('mod', 'mean'), obs=('obs', 'mean'))
    return f, bb.merge(SUN, on='entity', how='left')


def pairs(bb0, bb1, Fe):
    t0, t1 = bb0.set_index('entity'), bb1.set_index('entity')
    out = []
    for a, b in PAIRS:
        if a not in t0.index or b not in t0.index:
            continue
        sun, sha = (a, b) if Fe.get(a, 1.0) >= Fe.get(b, 1.0) else (b, a)
        dm = t0.obs[sun] - t0.obs[sha]
        d0, d1 = t0['mod'][sun] - t0['mod'][sha], t1['mod'][sun] - t1['mod'][sha]
        e0, e1 = abs(d0 - dm), abs(d1 - dm)
        out.append(dict(sunny=sun, shaded=sha, site=t0.site[sun], dT_meas=dm, dT_FRHS=d0, dT_V=d1,
                        helps=1.0 if e1 < e0 - PAIR_MIN else 0.0 if e1 > e0 + PAIR_MIN else 0.5,
                        flagged=(a, b) in FLAGGED_PAIRS))
    T = pd.DataFrame(out)
    votes = T.groupby('site').helps.mean()
    return T, dict(sites_help=int((votes > .5).sum()), sites_hurt=int((votes < .5).sum()),
                   sites_tie=int((votes == .5).sum()), pairs_help=int((T.helps == 1).sum()),
                   pairs_hurt=int((T.helps == 0).sum()), pairs=len(T))


# ── measures of one variant ──────────────────────────────────────────────────────────────
def lomo_measures(P, E, obs_cls, v, grid):
    cols = {fr: ft.col(v, fr) for fr in grid if fr <= 1.0}
    pred, chosen = ft.lomo(P, cols)
    Q = P.assign(oos=pred)
    zone = Q.band.between(*ft.ZONE)
    keep = ~Q.rgi.isin(GENERATING_OUTLINES)
    pk, _ = ft.lomo(P[keep], cols)
    cls = ft.entity_classes(Q, 'oos', E)
    return dict(chosen=chosen, M1=ft.m1(Q, 'oos'), zone_M1=ft.m1(Q[zone], 'oos'),
                above_M1=ft.m1(Q[Q.band >= 3800], 'oos'), strict_M1=ft.m1(P[keep].assign(oos=pk), 'oos'),
                cls=cls, pred=pred)


def within_measures(P, v, grid, B):
    high = ((P.band >= HIGH) & ~P.entity.isin(EXCLUDED)).to_numpy()
    out_gen = ~P.site.isin(GENERATING_SITES).to_numpy()
    f_all, bb = band_bias(P, v, grid, high)
    f_wo, bb_wo = band_bias(P, v, grid, high & out_gen)
    b, p, ci = wild_cluster(bb.b.to_numpy(), design(bb), bb.site.to_numpy(), 0, B=B)
    rec = dict(f_common=f_all, slope=b * 0.1, slope_p=p, slope_90=(ci[0] * 0.1, ci[1] * 0.1),
               f_common_wo=f_wo, slope_wo=slope(bb_wo), within_rmse=within_rmse(bb)[0],
               within_sites=within_rmse(bb)[1], within_rmse_wo=within_rmse(bb_wo)[0],
               slope_massif={m: slope(bb[bb.massif == m]) for m in ('Mont Blanc', 'Monte Rosa')},
               slope_other=slope(bb[~bb.massif.isin(['Mont Blanc', 'Monte Rosa'])]),
               slope_site_dropped={s: slope(bb[bb.site != s]) for s in sorted(bb.site.unique())},
               within_rmse_site={s: within_rmse(bb[bb.site == s])[0] for s in sorted(bb.site.unique())
                                 if (bb.site == s).sum() >= 2})
    return rec, bb


def rule_1a(P, v, grid):
    g = tuple(fr for fr in grid if fr <= 1.0)
    fine = FINE[(FINE >= g[0]) & (FINE <= g[-1])]
    best = {}
    for m in ('Mont Blanc', 'Monte Rosa'):
        sel = P.massif == m
        sc = {f: ft.m1(P[sel].assign(x=at_f(P[sel], v, f, g)), 'x') for f in fine}
        best[m] = min(sc, key=sc.get)
    tac = P.entity.str.startswith('Taconnaz Glacier') & (P.band >= HIGH)
    mr = (P.massif == 'Monte Rosa') & (P.band >= HIGH)
    c85 = ft.col(v, 0.85)
    return dict(f_best=best, gap=abs(best['Mont Blanc'] - best['Monte Rosa']),
                taconnaz_mean_085=float(P.loc[tac, c85].mean()), monte_rosa_mean_085=float(P.loc[mr, c85].mean()),
                f_best_monte_rosa_to_130=(rule_1a_130(P, v) if v in EXTRA and ft.col(v, EXTRA[v][0]) in P else None))


def rule_1a_130(P, v):
    g = GRID + EXTRA[v]
    sel = P.massif == 'Monte Rosa'
    fine = np.round(np.arange(0.03, 1.3001, 0.01), 2)
    sc = {f: ft.m1(P[sel].assign(x=at_f(P[sel], v, f, g)), 'x') for f in fine}
    return min(sc, key=sc.get)


def frozen_beds(E, runs_by, chosen):
    FB = ft.frozen_bed_table(runs_by, chosen, E)
    out = {}
    for v in chosen:
        F = FB[np.isfinite(FB[f'T_{v}'].astype(float))]
        out[v] = dict(found=int((F.frozen_meas & F[f'frozen_{v}']).sum()), of=int(F.frozen_meas.sum()),
                      temperate_called_frozen=int((~F.frozen_meas & F[f'frozen_{v}']).sum()))
    return FB, out


SUN = None
EXCLUDED = ()


def main():
    global SUN, EXCLUDED
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('out_dir', nargs='?', default=str(W3 / 'orient_stage1'))
    ap.add_argument('--stand-in', action='store_true')
    ap.add_argument('--boot', type=int, default=9999)
    ap.add_argument('--frhs-tag', default=TAG)
    ap.add_argument('--tag', default=TAG)
    ap.add_argument('--suffix', default='')
    a = ap.parse_args()
    NAMES.update(frhs_tag=a.frhs_tag, tag=a.tag, suffix=a.suffix)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    dh = DataHandler(region='CentralEurope', pool_within='hole_year')
    dh.load()
    E = ft.scored_entities(dh.calibration_glaciers, ts.load_borehole_table())
    sites = site_keys([g for g in dh.calibration_glaciers if g.glacier_id not in ft.SEALED])
    runs = {ft.col(v, fr): run_of(v, fr, a.stand_in) for v in VARS for fr in grid_of(v, a.stand_in)}
    missing = [r for r in set(runs.values()) if not ft.output_path(r, 'hole_year').exists()]
    if missing:
        raise SystemExit(f'{len(missing)} hole-year outputs missing, e.g. {sorted(missing)[:3]}')
    cache = {r: pickle.load(open(ft.output_path(r, 'hole_year'), 'rb')) for r in set(runs.values())}
    P = ft.points_frame(E, {c: cache[r] for c, r in runs.items()})
    counts = dict(entities=P.entity.nunique(), points=len(P), outlines=P.rgi.nunique(), massifs=P.massif.nunique())
    expect = dict(entities=ft.EXPECT['entities'], points=1135, outlines=ft.EXPECT['outlines'],
                  massifs=ft.EXPECT['massifs'])
    if counts != expect:
        raise SystemExit(f'scored set differs: {counts} != {expect}')
    finite = np.isfinite(P[list(runs)]).all(axis=1)
    if not finite.all():
        raise SystemExit(f'{int((~finite).sum())} points without a finite model value in some run')
    P['site'] = P.entity.map(sites)
    B = pd.read_csv(BANDS, dtype={'rgi': str})
    ref = B.set_index('entity').site
    bad = {e: (ref.get(e), s) for e, s in P.groupby('entity').site.first().items() if e in ref.index and ref[e] != s}
    if bad:
        raise SystemExit(f'the site table differs from the one the investigation used: {bad}')
    SUN = B[['entity', 'ratio_jja']]
    if SUN.set_index('entity').ratio_jja.reindex(P[P.band >= HIGH].entity.unique()).isna().any():
        raise SystemExit('a high band has no sunshine ratio')
    A = pd.read_csv(TABLES / 'applied.csv', dtype={'rgi': str})
    if set(A.entity) != set(P.entity):
        raise SystemExit('applied.csv and the scored entities differ')
    EXCLUDED = tuple(A[A.exclude_within_site].entity)
    F0 = {v: dict(zip(A.entity, A[f'F_band_{v}'])) for v in ('C', 'W', 'S')}
    F_own = dict(zip(A.entity, A.F_C))

    obs_cls = ft.entity_classes(P, 'obs', E)
    zone = P.band.between(*ft.ZONE)
    names, znames = sorted(P.entity.unique()), sorted(P[zone].entity.unique())
    L, Wm, BB = {}, {}, {}
    for v in VARS:
        L[v] = lomo_measures(P, E, obs_cls, v, grid_of(v, a.stand_in))
        Wm[v], BB[v] = within_measures(P, v, grid_of(v, a.stand_in), a.boot)
        print(v, {k: L[v][k] for k in ('M1', 'zone_M1', 'above_M1', 'strict_M1')},
              {k: Wm[v][k] for k in ('f_common', 'slope', 'slope_wo', 'within_rmse')}, flush=True)

    res = dict(counts=counts, stand_in=a.stand_in, zone=ft.ZONE, grid=GRID, rule='B', variants={}, runs=dict(NAMES),
               excluded_within_site=list(EXCLUDED))
    PT = []
    for v in VARS:
        r = dict(lomo={k: L[v][k] for k in ('chosen', 'M1', 'zone_M1', 'above_M1', 'strict_M1')},
                 within={k: Wm[v][k] for k in Wm[v]}, rule_1a=rule_1a(P, v, grid_of(v, a.stand_in)),
                 counts={c: f'{x} of {n}' for c, (x, n) in ft.class_counts(obs_cls, L[v]['cls'], names).items()},
                 zone_counts={c: f'{x} of {n}' for c, (x, n) in ft.class_counts(obs_cls, L[v]['cls'], znames).items()})
        if v != 'FRHS':
            net, per = fbr.site_changes(obs_cls, L['FRHS']['cls'], L[v]['cls'], names)
            znet, zper = fbr.site_changes(obs_cls, L['FRHS']['cls'], L[v]['cls'], znames)
            T, votes = pairs(BB['FRHS'], BB[v], F0[v])
            PT.append(T.assign(variant=v))
            c1 = dict(all=abs(Wm[v]['slope']) <= HALVE * abs(Wm['FRHS']['slope']),
                      without_generating_sites=abs(Wm[v]['slope_wo']) <= HALVE * abs(Wm['FRHS']['slope_wo']))
            dM = {k: L[v][k] - L['FRHS'][k] for k in ('M1', 'zone_M1', 'strict_M1', 'above_M1')}
            c2 = dict(M1=dM['M1'] <= HARMS, zone_M1=dM['zone_M1'] <= HARMS, strict_M1=dM['strict_M1'] <= HARMS_STRICT,
                      sites=net >= 0, zone_sites=znet >= 0)
            r.update(condition_1=c1, condition_2=c2, dM1=dM, net_sites=net, net_zone_sites=znet, pairs=votes,
                     holds=bool(all(c1.values()) and all(c2.values())),
                     sites_changed={k: d for k, d in per.items() if d['site'] != 'none'},
                     zone_sites_changed={k: d for k, d in zper.items() if d['site'] != 'none'})
        res['variants'][v] = r
    res['adopt'] = bool(all(res['variants'][v]['holds'] for v in DECIDING))

    # frozen beds, time test, every band at depth (reported)
    runs_by = {(v, fr): run_of(v, fr, a.stand_in) for v in VARS for fr in grid_of(v, a.stand_in) if fr <= 1.0}
    FB, fb = frozen_beds(E, runs_by, {v: L[v]['chosen'] for v in VARS})
    FB.to_csv(out / 'frozen_bed.csv', index=False)
    W = ft.warming_pairs(sorted(set(runs_by.values())))
    for v in VARS:
        res['variants'][v]['frozen_beds'] = fb[v]
        cmap = {m: run_of(v, fr, a.stand_in) for m, fr in L[v]['chosen'].items()}
        res['variants'][v]['time_test'] = {k: ft.warming_measures(W, cmap)[0].get(k) for k in ('n_sites', 'mase', 'ratio')}
    Q = P[['entity', 'site', 'rgi', 'massif', 'band', 'depth', 'obs']].copy()
    for v in VARS:
        Q[v] = L[v]['pred']
    Q.to_csv(out / 'points.csv', index=False, float_format='%.3f')
    rows = []
    for e, d in Q.groupby('entity', sort=False):
        rr = dict(band=e, site=d.site.iloc[0], rgi=d.rgi.iloc[0], elev=d.band.iloc[0], n=len(d), obs=d.obs.mean(),
                  F_C=F0['C'][e], F_C_own=F_own[e], excluded_within_site=e in EXCLUDED)
        for v in VARS:
            err = d[v] - d.obs
            rr[f'rmse_{v}'], rr[f'bias_{v}'] = float(np.sqrt((err ** 2).mean())), float(err.mean())
        rows.append(rr)
    T = pd.DataFrame(rows)
    T.to_csv(out / 'bands.csv', index=False, float_format='%.3f')
    for v in DECIDING + ('S',):
        d = T[f'rmse_{v}'] - T.rmse_FRHS
        res['variants'][v]['bands_worse_by_0.2'] = T.loc[d > 0.2, 'band'].tolist()
        res['variants'][v]['bands_better_by_0.2'] = T.loc[d < -0.2, 'band'].tolist()
        res['variants'][v]['group_c_bias'] = {e: float(T.loc[T.band == e, f'bias_{v}'].iloc[0]) for e in GROUP_C}
    res['variants']['FRHS']['group_c_bias'] = {e: float(T.loc[T.band == e, 'bias_FRHS'].iloc[0]) for e in GROUP_C}
    pd.concat(PT).to_csv(out / 'pairs.csv', index=False, float_format='%.3f')
    (out / 'orient_stage1.json').write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps(dict(adopt=res['adopt'], **{v: {k: res['variants'][v].get(k) for k in
                     ('holds', 'condition_1', 'condition_2', 'dM1', 'pairs')} for v in VARS if v != 'FRHS'}),
                     indent=1, default=str))


if __name__ == '__main__':
    main()
