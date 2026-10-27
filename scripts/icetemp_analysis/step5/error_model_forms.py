"""Step 3.1 of the Step 5 plan (20 Oct): the error model's form, chosen by evidence before any fit.

Residuals r = model - reading at every calibration point (all depths, one point per borehole and year, the Colle
Gnifetti points below 40 m left out as in the fit) of one stand-in run near the Step 5 optimum. Each candidate form is
fitted by maximum likelihood (one free mean, several starts) and scored by AIC (BIC reported). Its checks, on its own
fit: the whitened misfit r' S^-1 r / n overall; the leave-one-out z^2 mean at and above / below 15 m; per site the
whitened misfit against the 99 % range of chi^2 with the site's number of points.

Candidates (D = the readings' own variances, always in):
  C0  today: band term exponential in depth (s2, ell; same entity), borehole sh2, borehole and year sy2, site ss2
  C1  C0 + an error of each reading on its own (w2)
  C2  C0 + w2 at and above 15 m only
  C3  C1 with the band term split at 15 m (s2a, ella at and above; s2b, ellb below; no link across)
  C4  C1 with the borehole-and-year term replaced by a borehole term that fades with the years between visits
      (st2, tau)
  C5  C1 with their own w2 and s2 at the ice-core saddles (w2c, s2c; the Monte Rosa summit saddles of the site table)
  C6  C1 with room for rare large misses: a multivariate t per site (nu; S stays the covariance)
  C6b C2 with the same t (the heavy tails on the best Gaussian form)
Sites: 'rule' = holes linked when within 1 km and 150 m of each other (single linkage); 'table' = data.site_keys.

Writes OUT/residuals.csv, OUT/forms.json (estimates, AIC/BIC, checks) and OUT/sites.csv (per site and form).

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/step5/error_model_forms.py OUT
           --run fs1_FRHS_r0.85_i1.017_a1.00_fa_fs1.5_oCP_w0 [--sites rule|table] [--forms C0 C1 ...]
"""
import argparse
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import gammaln
from scipy.stats import chi2

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import forcing_fair_test as ft  # noqa: E402
from icetemp.calibration import thermal_structure as ts  # noqa: E402
from icetemp.calibration.data import DataHandler, fit_exclusions, site_keys  # noqa: E402

SHALLOW = 15.0
SITE_KM, SITE_DZ = 1.0, 150.0
SADDLES = 'Monte Rosa summit saddles'
V_LO, V_HI = 1e-4, 30.0
L_LO, L_HI = 0.3, 1000.0
PARAMS = {    # name -> (start values, lower, upper); variances in C^2, lengths in m, tau in years
    's2': ((0.5, 2.0), V_LO, V_HI), 'ell': ((5.0, 30.0), L_LO, L_HI), 'sh2': ((0.3, 2.0), V_LO, V_HI),
    'sy2': ((0.1,), V_LO, V_HI), 'ss2': ((0.5,), V_LO, V_HI), 'w2': ((0.2,), V_LO, V_HI),
    's2a': ((0.5,), V_LO, V_HI), 'ella': ((2.0,), L_LO, L_HI), 's2b': ((1.0,), V_LO, V_HI),
    'ellb': ((50.0,), L_LO, L_HI), 'st2': ((0.3,), V_LO, V_HI), 'tau': ((10.0,), 0.5, 200.0),
    'w2c': ((0.2,), V_LO, V_HI), 's2c': ((1.0,), V_LO, V_HI), 'nu': ((10.0,), 2.1, 200.0)}
FORMS = {'C0': ('s2', 'ell', 'sh2', 'sy2', 'ss2'),
         'C1': ('s2', 'ell', 'sh2', 'sy2', 'ss2', 'w2'),
         'C2': ('s2', 'ell', 'sh2', 'sy2', 'ss2', 'w2'),
         'C3': ('s2a', 'ella', 's2b', 'ellb', 'sh2', 'sy2', 'ss2', 'w2'),
         'C4': ('s2', 'ell', 'sh2', 'st2', 'tau', 'ss2', 'w2'),
         'C5': ('s2', 'ell', 'sh2', 'sy2', 'ss2', 'w2', 'w2c', 's2c'),
         'C6': ('s2', 'ell', 'sh2', 'sy2', 'ss2', 'w2', 'nu'),
         'C6b': ('s2', 'ell', 'sh2', 'sy2', 'ss2', 'w2', 'nu')}


def rule_sites(holes, bh):
    """hole -> site: holes within SITE_KM and SITE_DZ of each other linked (single linkage)."""
    t = bh.set_index(bh.id.astype(str)).loc[holes]
    lat, lon, z = np.radians(t.latitude.to_numpy()), np.radians(t.longitude.to_numpy()), t.elevation.to_numpy()
    km = 6371.0 * 2 * np.arcsin(np.sqrt(np.sin((lat[:, None] - lat) / 2) ** 2 +
                                        np.cos(lat[:, None]) * np.cos(lat) * np.sin((lon[:, None] - lon) / 2) ** 2))
    link = (km < SITE_KM) & (np.abs(z[:, None] - z) < SITE_DZ)
    lab = np.arange(len(holes))
    changed = True
    while changed:
        new = np.array([lab[link[i]].min() for i in range(len(holes))])
        changed = (new != lab).any()
        lab = new
    names = {}
    for k in np.unique(lab):
        members = t.iloc[np.flatnonzero(lab == k)]
        names[k] = f'{members.glacier_name.mode().iloc[0]} {members.elevation.min():.0f}-{members.elevation.max():.0f}'
    return {h: names[k] for h, k in zip(holes, lab)}


def residuals(run, sites):
    dh = DataHandler(region='CentralEurope', pool_within='hole_year')
    dh.load()
    gl = [g for g in dh.calibration_glaciers if g.glacier_id and g.glacier_id not in ft.SEALED]
    excl = fit_exclusions(gl, ['colle_gnifetti_below_40m'])
    table = site_keys(gl)
    out = pickle.load(open(ft.output_path(run, 'hole_year'), 'rb'))
    rows = []
    for g in gl:
        o = out.get(g.glacier_name)
        if o is None:
            continue
        mod = np.asarray(o[1], float)[:len(g.depths)]
        for i, d in enumerate(g.depths):
            if (g.glacier_name, float(d)) in excl or not np.isfinite(mod[i]):
                continue
            rows.append(dict(entity=g.glacier_name, massif=ft.MASSIF_OF.get(g.glacier_id), band=g.elevation_band,
                             table_site=table[g.glacier_name], hole=str(g.point_hole[i]), year=int(g.point_year[i]),
                             depth=float(d), sig2=float(g.sigma[i]) ** 2, r=float(mod[i] - g.T_obs[i])))
    R = pd.DataFrame(rows)
    hs = rule_sites(sorted(R.hole.unique()), ts.load_borehole_table())
    R['rule_site'] = R.hole.map(hs)
    R['site'] = R.rule_site if sites == 'rule' else R.table_site
    R['saddle'] = R.table_site == SADDLES
    return R


def band_term(e, d, s2, ell, mask=None):
    same = e[:, None] == e[None, :]
    if mask is not None:
        same &= mask[:, None] & mask[None, :]
    return np.where(same, s2 * np.exp(-np.abs(d[:, None] - d[None, :]) / ell), 0.0)


def covariance(form, p, B):
    """S of one site block: the readings' own variances plus the form's terms."""
    e, d, h, y, s, sad = B['e'], B['d'], B['h'], B['y'], B['s'], B['sad']
    S = np.diag(B['sig2'])
    same_h = h[:, None] == h[None, :]
    if form == 'C3':
        S += band_term(e, d, p['s2a'], p['ella'], d <= SHALLOW) + band_term(e, d, p['s2b'], p['ellb'], d > SHALLOW)
    elif form == 'C5':
        S += band_term(e, d, p['s2'], p['ell'], ~sad) + band_term(e, d, p['s2c'], p['ell'], sad)
    else:
        S += band_term(e, d, p['s2'], p['ell'])
    S += p['sh2'] * same_h
    if form == 'C4':
        S += np.where(same_h, p['st2'] * np.exp(-np.abs(y[:, None] - y[None, :]) / p['tau']), 0.0)
    else:
        S += p['sy2'] * (same_h & (y[:, None] == y[None, :]))
    S += p['ss2'] * (s[:, None] == s[None, :])
    if form in ('C2', 'C6b'):
        S += np.diag(p['w2'] * (d <= SHALLOW))
    elif form == 'C5':
        S += np.diag(np.where(sad, p['w2c'], p['w2']))
    elif form != 'C0':
        S += np.diag(np.full(len(d), p['w2']))
    return S


class Model:
    def __init__(self, R, form):
        self.form, self.names = form, FORMS[form]
        self.t = form.startswith('C6')
        self.blocks = [dict(e=b.entity.to_numpy(object), d=b.depth.to_numpy(), h=b.hole.to_numpy(object),
                            y=b.year.to_numpy(), s=b.site.to_numpy(object), sad=b.saddle.to_numpy(bool),
                            sig2=b.sig2.to_numpy(), r=b.r.to_numpy(), site=k) for k, b in R.groupby('site', sort=False)]
        self.n = len(R)
        self.lo = np.log([PARAMS[k][1] for k in self.names])
        self.hi = np.log([PARAMS[k][2] for k in self.names])

    def params(self, q):
        return dict(zip(self.names, np.exp(q[:len(self.names)])))

    def chol(self, p):
        return [np.linalg.cholesky(covariance(self.form, p, B)) for B in self.blocks]

    def gls_mean(self, Ls):
        a = b = 0.0
        for L, B in zip(Ls, self.blocks):
            x1, xr = np.linalg.solve(L, np.ones(len(B['r']))), np.linalg.solve(L, B['r'])
            a, b = a + x1 @ x1, b + x1 @ xr
        return b / a

    def nll(self, q):
        p = self.params(q)
        try:
            Ls = self.chol(p)
        except np.linalg.LinAlgError:
            return 1e12
        mu = q[-1] if self.t else self.gls_mean(Ls)
        f = 0.5 * self.n * np.log(2 * np.pi)
        for L, B in zip(Ls, self.blocks):
            x = np.linalg.solve(L, B['r'] - mu)
            logdet = 2 * np.log(np.diag(L)).sum()
            if not self.t:
                f += 0.5 * (logdet + x @ x)
            else:     # scale matrix S (nu - 2) / nu, so that S is the covariance
                nu, m = p['nu'], len(x)
                c = (nu - 2) / nu
                f += -(gammaln((nu + m) / 2) - gammaln(nu / 2) - 0.5 * m * np.log(nu * np.pi) -
                       0.5 * (logdet + m * np.log(c)) - (nu + m) / 2 * np.log1p(x @ x / (c * nu))) - \
                    0.5 * m * np.log(2 * np.pi)
        return f

    def fit(self):
        grids = np.meshgrid(*[np.log(PARAMS[k][0]) for k in self.names], indexing='ij')
        starts = np.stack([g.ravel() for g in grids], axis=1)
        if self.t:
            starts = np.column_stack([starts, np.zeros(len(starts))])
        bounds = list(zip(self.lo, self.hi)) + ([(-10.0, 10.0)] if self.t else [])
        best = None
        for q0 in starts:
            res = minimize(self.nll, q0, method='L-BFGS-B', bounds=bounds)
            if best is None or res.fun < best.fun:
                best = res
        return best.x, float(best.fun)

    def checks(self, q):
        p = self.params(q)
        Ls = self.chol(p)
        mu = q[-1] if self.t else self.gls_mean(Ls)
        z2, dep, sites = [], [], []
        for L, B in zip(Ls, self.blocks):
            rr = B['r'] - mu
            Q = np.linalg.inv(L @ L.T)
            z = Q @ rr / np.sqrt(np.diag(Q))
            z2.append(z ** 2)
            dep.append(B['d'])
            maha = float(rr @ Q @ rr)
            m = len(rr)
            sites.append(dict(site=B['site'], n=m, misfit_per_n=maha / m,
                              inside_99=bool(chi2.ppf(0.005, m) <= maha <= chi2.ppf(0.995, m))))
        z2, dep = np.concatenate(z2), np.concatenate(dep)
        return dict(mean=float(mu), misfit_per_n=sum(s['misfit_per_n'] * s['n'] for s in sites) / self.n,
                    loo_z2_shallow=float(z2[dep <= SHALLOW].mean()), loo_z2_deep=float(z2[dep > SHALLOW].mean()),
                    sites_outside_99=[s['site'] for s in sites if not s['inside_99']]), sites


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('out')
    ap.add_argument('--run', required=True)
    ap.add_argument('--sites', default='rule', choices=['rule', 'table'])
    ap.add_argument('--forms', nargs='+', default=sorted(FORMS))
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    R = residuals(a.run, a.sites)
    R.to_csv(out / 'residuals.csv', index=False, float_format='%.4f')
    print(f'{a.run}: {len(R)} residuals, {R.entity.nunique()} entities, {R.hole.nunique()} holes, '
          f'{R.site.nunique()} sites ({a.sites}); mean {R.r.mean():+.2f} C, sd {R.r.std():.2f} C', flush=True)
    res, site_rows = dict(run=a.run, sites=a.sites, n=len(R), n_sites=int(R.site.nunique()), forms={}), []
    for form in a.forms:
        M = Model(R, form)
        q, f = M.fit()
        k = len(q)
        chk, srows = M.checks(q)
        res['forms'][form] = dict(params=M.params(q), k=k, nll=f, aic=2 * f + 2 * k, bic=2 * f + k * np.log(M.n),
                                  **chk)
        site_rows += [dict(form=form, **s) for s in srows]
        print(f'{form}: AIC {2 * f + 2 * k:.1f}, misfit/n {chk["misfit_per_n"]:.2f}, shallow '
              f'{chk["loo_z2_shallow"]:.2f}, deep {chk["loo_z2_deep"]:.2f}, sites outside 99 % '
              f'{len(chk["sites_outside_99"])}, {({n: round(v, 3) for n, v in M.params(q).items()})}', flush=True)
    best = min(res['forms'], key=lambda f: res['forms'][f]['aic'])
    for f, v in res['forms'].items():
        v['d_aic'] = v['aic'] - res['forms'][best]['aic']
    (out / 'forms.json').write_text(json.dumps(res, indent=1, default=float))
    pd.DataFrame(site_rows).to_csv(out / 'sites.csv', index=False, float_format='%.3f')


if __name__ == '__main__':
    main()
