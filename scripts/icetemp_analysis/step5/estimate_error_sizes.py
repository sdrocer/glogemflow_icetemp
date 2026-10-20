"""Step 3 of the plan of 20 Oct: the sizes of the error terms of the Step-5 likelihood, estimated once from the stage-2
FRHS residuals and fixed before any campaign run.

Residuals r = FRHS - reading at every calibration point (all depths, one point per borehole and year, the Colle Gnifetti
points below 40 m left out as in the fit), each massif at the f_ref it chose in stage 2 (0.85 where a massif had no
scored points), f_ins 1.017. Model: r ~ N(mu, D + C) with D the readings' own variances and C =
residual_covariance.model_covariance (band s2 over ell in depth, borehole sh2, borehole-year sy2, site ss2). Restricted
maximum likelihood (REML: the free mean mu, here +0.6 C, is not mistaken for site scatter), from several starts; profile
intervals (2 log-likelihood units); a simulation check (refits of data simulated on the real layout from the estimate);
the fit on points at and above / below 15 m; and the sensitivity arms of the plan:
  A0 the estimate; A1 no site term; A2 the site term at the top of its interval; A3 the other end of the ell/sh2 ridge
  (ell at the far end of its interval, the others refit); A4 estimated without the two Grenz sites; REF today's
  3.5 C^2 over 15 m and no shared terms.

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/step5/estimate_error_sizes.py [OUT_DIR]
           [--n-sim 200]
"""
import argparse
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import forcing_fair_test as ft  # noqa: E402
from icetemp.calibration.data import DataHandler, fit_exclusions, site_keys  # noqa: E402
from icetemp.calibration.residual_covariance import model_covariance  # noqa: E402

W3 = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3')
NAMES = ('s2', 'ell', 'sh2', 'sy2', 'ss2')
LO = np.log([0.01, 1.0, 1e-4, 1e-4, 1e-4])
HI = np.log([30.0, 500.0, 30.0, 30.0, 30.0])
FREF_UNSCORED = 0.85
GRENZ = ('Grenzgletscher tongue', 'Monte Rosa summit saddles')


def residuals():
    chosen = json.load(open(W3 / 'stage2/stage2_hole_year.json'))['summary']['FRHS']['chosen_fref']
    dh = DataHandler(region='CentralEurope', pool_within='hole_year')
    dh.load()
    gl = [g for g in dh.calibration_glaciers if g.glacier_id and g.glacier_id not in ft.SEALED]
    excl = fit_exclusions(gl, ['colle_gnifetti_below_40m'])
    sites = site_keys(gl)
    cache, rows = {}, []
    for g in gl:
        fr = chosen.get(ft.MASSIF_OF.get(g.glacier_id), FREF_UNSCORED)
        run = ft.run_name('s2', 'FRHS', fr)
        if run not in cache:
            cache[run] = pickle.load(open(ft.output_path(run, 'hole_year'), 'rb'))
        o = cache[run].get(g.glacier_name)
        if o is None:
            continue
        mod = np.asarray(o[1], float)[:len(g.depths)]
        for i, d in enumerate(g.depths):
            if (g.glacier_name, float(d)) in excl or not np.isfinite(mod[i]):
                continue
            rows.append(dict(entity=g.glacier_name, site=sites[g.glacier_name], hole=str(g.point_hole[i]),
                             year=int(g.point_year[i]), depth=float(d), sig2=float(g.sigma[i]) ** 2,
                             r=float(mod[i] - g.T_obs[i])))
    return pd.DataFrame(rows)


class Reml:
    def __init__(self, R):
        self.blocks = [dict(e=b.entity.to_numpy(), d=b.depth.to_numpy(), h=list(b.hole), y=b.year.to_numpy(),
                            s=list(b.site), sig2=b.sig2.to_numpy(), r=b.r.to_numpy())
                       for _, b in R.groupby('site', sort=False)]
        self.n = len(R)

    def terms(self, p, r_override=None):
        s2, ell, sh2, sy2, ss2 = np.exp(p)
        logdet = a = b = c = 0.0
        for k, B in enumerate(self.blocks):
            r = B['r'] if r_override is None else r_override[k]
            S = np.diag(B['sig2']) + model_covariance(B['e'], B['d'], s2, ell, hole=B['h'], year=B['y'], site=B['s'],
                                                      sh2=sh2, sy2=sy2, ss2=ss2)
            L = np.linalg.cholesky(S)
            x1 = np.linalg.solve(L, np.ones(len(r)))
            xr = np.linalg.solve(L, r)
            logdet += 2 * np.log(np.diag(L)).sum()
            a += x1 @ x1
            b += x1 @ xr
            c += xr @ xr
        return logdet, a, b, c

    def nll(self, p, r_override=None):
        logdet, a, b, c = self.terms(p, r_override)
        return 0.5 * (logdet + np.log(a) + c - b * b / a)

    def mean(self, p):
        _, a, b, _ = self.terms(p)
        return b / a

    def fit(self, starts, fixed=None, r_override=None):
        fixed = fixed or {}
        free = [k for k in range(5) if k not in fixed]

        def f(q):
            p = np.empty(5)
            p[free] = q
            for k, v in fixed.items():
                p[k] = v
            return self.nll(p, r_override)
        best = None
        for s in starts:
            q0 = np.asarray(s, float)[free]
            res = minimize(f, q0, method='L-BFGS-B', bounds=list(zip(LO[free], HI[free])))
            if best is None or res.fun < best.fun:
                best = res
        p = np.empty(5)
        p[free] = best.x
        for k, v in fixed.items():
            p[k] = v
        return p, float(best.fun)


def starts():
    out = []
    for s2 in (0.5, 2.0):
        for ell in (5.0, 30.0):
            for sh2 in (0.3, 2.0):
                out.append(np.log([s2, ell, sh2, 0.1, 0.5]))
    return out


def profile(model, p_hat, nll_hat, k, half_width=2.5, n_grid=41):
    grid = np.clip(p_hat[k] + np.linspace(-half_width, half_width, n_grid), LO[k], HI[k])
    vals = []
    for v in grid:
        _, f = model.fit([p_hat], fixed={k: v})
        vals.append(f - nll_hat)
    vals = np.array(vals)
    inside = grid[vals <= 2.0]
    return dict(grid=np.exp(grid).round(4).tolist(), dnll=vals.round(3).tolist(),
                lo=float(np.exp(inside.min())) if len(inside) else None,
                hi=float(np.exp(inside.max())) if len(inside) else None)


def as_dict(p):
    return {n: float(v) for n, v in zip(NAMES, np.exp(p))}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('out_dir', nargs='?', default=str(W3 / 'step5_error_sizes'))
    ap.add_argument('--n-sim', type=int, default=200)
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    R = residuals()
    R.to_csv(out / 'residuals.csv', index=False, float_format='%.4f')
    print(f'{len(R)} residuals, {R.entity.nunique()} entities, {R.hole.nunique()} holes, {R.site.nunique()} sites; '
          f'mean {R.r.mean():+.2f} C, sd {R.r.std():.2f} C', flush=True)
    model = Reml(R)
    p0, f0 = model.fit(starts())
    res = dict(n=len(R), n_entities=int(R.entity.nunique()), n_holes=int(R.hole.nunique()),
               n_sites=int(R.site.nunique()), A0=as_dict(p0), mean_A0=float(model.mean(p0)), nll_A0=f0)
    print('A0', res['A0'], f'mean {res["mean_A0"]:+.2f}', flush=True)
    res['profiles'] = {NAMES[k]: profile(model, p0, f0, k) for k in range(5)}
    for n, pr in res['profiles'].items():
        print(f'profile {n}: {pr["lo"]} - {pr["hi"]}', flush=True)
    # arms
    k_ss, k_ell = NAMES.index('ss2'), NAMES.index('ell')
    p1, f1 = model.fit([p0], fixed={k_ss: LO[k_ss]})
    res['A1'] = dict(as_dict(p1), ss2=0.0, dnll=f1 - f0)
    ss_hi = res['profiles']['ss2']['hi']
    p2, f2 = model.fit([p0], fixed={k_ss: np.log(ss_hi)})
    res['A2'] = dict(as_dict(p2), dnll=f2 - f0)
    lo, hi = res['profiles']['ell']['lo'], res['profiles']['ell']['hi']
    far = hi if abs(np.log(hi) - p0[k_ell]) > abs(np.log(lo) - p0[k_ell]) else lo
    p3, f3 = model.fit([p0], fixed={k_ell: np.log(far)})
    res['A3'] = dict(as_dict(p3), dnll=f3 - f0)
    no_g = Reml(R[~R.site.isin(GRENZ)])
    p4, _ = no_g.fit(starts())
    res['A4'] = as_dict(p4)
    res['REF'] = dict(s2=3.5, ell=15.0, sh2=0.0, sy2=0.0, ss2=0.0,
                      dnll=float(model.nll(np.log([3.5, 15.0, 1e-10, 1e-10, 1e-10])) - f0))
    # fits by depth
    for name, sub in (('at_or_above_15m', R[R.depth <= 15]), ('below_15m', R[R.depth > 15])):
        pp, _ = Reml(sub).fit(starts())
        res[f'fit_{name}'] = dict(as_dict(pp), n=len(sub))
    # simulation check
    rng = np.random.default_rng(0)
    s2, ell, sh2, sy2, ss2 = np.exp(p0)
    mu = res['mean_A0']
    chol = [np.linalg.cholesky(np.diag(B['sig2']) + model_covariance(B['e'], B['d'], s2, ell, hole=B['h'], year=B['y'],
                                                                     site=B['s'], sh2=sh2, sy2=sy2, ss2=ss2))
            for B in model.blocks]
    sims = []
    for i in range(a.n_sim):
        r_sim = [mu + L @ rng.standard_normal(L.shape[0]) for L in chol]
        ps, _ = model.fit([p0], r_override=r_sim)
        sims.append(np.exp(ps))
        if (i + 1) % 20 == 0:
            print(f'simulation {i + 1}/{a.n_sim}', flush=True)
    sims = np.array(sims)
    res['simulation'] = {n: dict(p05=float(np.percentile(sims[:, k], 5)), p50=float(np.median(sims[:, k])),
                                 p95=float(np.percentile(sims[:, k], 95))) for k, n in enumerate(NAMES)}
    res['simulation_corr'] = np.corrcoef(np.log(sims).T).round(2).tolist()
    (out / 'error_sizes.json').write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({k: v for k, v in res.items() if k != 'profiles'}, indent=1, default=float))


if __name__ == '__main__':
    main()
