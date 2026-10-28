"""Sizes of the depth target's error model, form C6b of the rules plan (27 Oct), by maximum likelihood on one
calibrator's residuals at one theta: the band term (s2, ell), the borehole (sh2), borehole-and-year (sy2) and site
(ss2) terms, each reading's own error at or above the calibrator's s_white_max_depth (w2), and the degrees of freedom
of the per-site t (nu). As in the likelihood: mean zero, the readings' own variances and the emulator's truncation
variance fixed, the emulator mean only. Also the form's own checks of the rules plan (rule 1 and rule 3, checks 1-2).

Usage:
    sizes, info = fit_sizes(cal, theta)      # sizes keyed like the calibrator's fields; info: nll, sizes at a bound
    cal2 = pipeline.make_calibrator(**sizes)
"""
import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln
from scipy.stats import chi2, f as f_dist, norm, t as t_dist

from .residual_covariance import model_covariance

NAMES = ('s2', 'ell', 'sh2', 'sy2', 'ss2', 'w2', 'nu')
FIELD = dict(s2='s_model_variance', ell='s_model_length_m', sh2='s_hole_variance', sy2='s_hole_year_variance',
             ss2='s_site_variance', w2='s_white_variance', nu='site_t_nu')
V_BOUNDS, L_BOUNDS, NU_BOUNDS = (1e-4, 30.0), (0.3, 1000.0), (2.1, 200.0)
BOUNDS = dict(s2=V_BOUNDS, ell=L_BOUNDS, sh2=V_BOUNDS, sy2=V_BOUNDS, ss2=V_BOUNDS, w2=V_BOUNDS, nu=NU_BOUNDS)
STARTS = dict(s2=(0.5, 2.0), ell=(5.0, 30.0), sh2=(0.3, 2.0), sy2=(0.1,), ss2=(0.5,), w2=(0.2,), nu=(10.0,))
OFFSET_P = 1 / 300      # rule 3, check 2: a site's offset is unusual at this two-sided chance under the form's law
SHALLOW = 15.0          # the checks' split into shallow (at and above) and deep readings [m]
N_LAW, LAW_SEED = 20000, 20261028


class SiteBlocks:
    """The calibrator's rows grouped by site, with each row's fixed variance (readings and emulator truncation)."""

    def __init__(self, cal):
        if cal.target != 'depth' or not getattr(cal, '_rows_ok', False):
            raise ValueError("the error model's sizes need a calibrator with target='depth'")
        trunc = getattr(cal.emulator, 'truncation_variance', lambda: 0.0)()
        self.w_max_depth = cal.s_white_max_depth
        site = np.asarray(cal._row_site, dtype=object)
        self.sites = list(dict.fromkeys(cal._row_site))
        self.idx = [np.flatnonzero(site == s) for s in self.sites]
        self.entity, self.depth = cal._row_entity, cal._row_depth
        self.hole = np.asarray(cal._row_hole, dtype=object)
        self.year = np.asarray(cal._row_year, dtype=object)
        self.site = site
        self.D = cal._row_sigma_obs2 + trunc
        self.n = len(self.D)

    def covariance(self, p, i):
        """S of the site block with rows i: the readings' variances plus the model error at sizes p."""
        C = model_covariance(self.entity[i], self.depth[i], p['s2'], p['ell'], hole=self.hole[i], year=self.year[i],
                             site=self.site[i], sh2=p['sh2'], sy2=p['sy2'], ss2=p['ss2'], w2=p['w2'],
                             w_max_depth=self.w_max_depth)
        C[np.diag_indices(len(i))] += self.D[i]
        return C


def residual(cal, theta):
    """obs - emulator mean at theta, in the calibrator's row order (the likelihood's residual)."""
    mean, _ = cal.emulator.predict_factor(np.asarray(theta, dtype=float))
    return cal._row_obs - mean[cal._row_emu_idx]


def neg_log_likelihood(blocks, r, p, gaussian=False):
    """-log likelihood of the residuals r under sizes p: a multivariate t per site whose covariance is the
    block's S (scale S (nu - 2) / nu), or a Gaussian with that S."""
    f = 0.0
    for i in blocks.idx:
        try:
            L = np.linalg.cholesky(blocks.covariance(p, i))
        except np.linalg.LinAlgError:
            return np.inf
        x = np.linalg.solve(L, r[i])
        q, logdet, m = float(x @ x), 2.0 * float(np.log(np.diag(L)).sum()), len(i)
        if gaussian:
            f += 0.5 * (q + logdet + m * np.log(2 * np.pi))
        else:
            nu = p['nu']
            f -= (gammaln((nu + m) / 2) - gammaln(nu / 2) - 0.5 * m * np.log((nu - 2) * np.pi) - 0.5 * logdet
                  - 0.5 * (nu + m) * np.log1p(q / (nu - 2)))
    return f


def fit_sizes(cal, theta, gaussian=False, r=None, start=None):
    """Maximum-likelihood sizes at theta (r: residuals instead, e.g. simulated), on the log scale from a grid of
    starts plus `start` (sizes keyed by the calibrator's field names, e.g. the sizes in use). Returns (sizes keyed by
    the calibrator's field names, info: nll, the sizes at a bound, whether the best start converged)."""
    blocks = SiteBlocks(cal)
    r = residual(cal, theta) if r is None else np.asarray(r, dtype=float)
    if not np.isfinite(r).all():
        raise ValueError(f'{int((~np.isfinite(r)).sum())} residuals are not finite')
    names = [n for n in NAMES if not (gaussian and n == 'nu')]
    lo, hi = np.log([BOUNDS[n][0] for n in names]), np.log([BOUNDS[n][1] for n in names])

    def nll(q):
        v = neg_log_likelihood(blocks, r, dict(zip(names, np.exp(q))), gaussian)
        return v if np.isfinite(v) else 1e12

    grids = np.meshgrid(*[np.log(STARTS[n]) for n in names], indexing='ij')
    starts = np.stack([g.ravel() for g in grids], axis=1)
    if start is not None:
        q = np.log([start[FIELD[n]] for n in names])
        starts = np.vstack([np.clip(q, lo, hi), starts])
    best = None
    for q0 in starts:
        res = minimize(nll, q0, method='L-BFGS-B', bounds=list(zip(lo, hi)))
        if best is None or res.fun < best.fun:
            best = res
    if not best.fun < 1e12:
        raise RuntimeError('no start gave a finite likelihood')
    out = {FIELD[n]: float(v) for n, v in zip(names, np.exp(best.x))}
    if gaussian:
        out['site_t_nu'] = None
    at_bound = [n for n, q, a, b in zip(names, best.x, lo, hi) if min(q - a, b - q) < 1e-6]
    return out, dict(nll=float(best.fun), at_bound=at_bound, converged=bool(best.success))


def sizes_of(cal):
    """The calibrator's current sizes, keyed by NAMES."""
    return {n: getattr(cal, FIELD[n]) for n in NAMES}


def misfit_law_range(m, nu, level=0.99, n=N_LAW, seed=LAW_SEED):
    """Range of sum_k q_k / sum_k m_k under the form's own law, q_k = r_k^T S_k^-1 r_k of the site blocks of sizes m:
    c m_k F(m_k, nu) per site for the t (c = (nu - 2) / nu), chi^2(m_k) for the Gaussian (nu None); by simulation."""
    rng = np.random.default_rng(seed)
    m = np.asarray(m, dtype=int)
    if nu is None:
        q = sum(rng.chisquare(k, n) for k in m)
    else:
        c = (nu - 2) / nu
        q = sum(c * k * rng.f(k, nu, n) for k in m)
    x = q / m.sum()
    return float(np.quantile(x, (1 - level) / 2)), float(np.quantile(x, 1 - (1 - level) / 2))


def checks(cal, theta):
    """The form's own checks at theta with the calibrator's sizes (rules plan, rule 1 and rule 3 checks 1-2).

    overall: r^T S^-1 r / n (about 1 if the form is right) and its 99 % range under the form's own law
    (misfit_law_range); loo_z2 at and above / below SHALLOW: the leave-one-out z^2 = (S^-1 r)_i^2 / (S^-1)_ii, mean
    about 1; per site: q / (c m) against F(m, nu) for the t (chi^2(m) / m for the Gaussian), inside its 99 % range;
    the site's offset, its GLS mean miss over its standard error, unusual below OFFSET_P two-sided under the form's
    own law (t with nu for the t)."""
    blocks = SiteBlocks(cal)
    r = residual(cal, theta)
    p = sizes_of(cal)
    nu = p['nu']
    z2 = np.empty(blocks.n)
    sites, q_all = [], 0.0
    for s, i in zip(blocks.sites, blocks.idx):
        Q = np.linalg.inv(blocks.covariance(p, i))
        Qr = Q @ r[i]
        z2[i] = Qr ** 2 / np.diag(Q)
        q, m = float(r[i] @ Qr), len(i)
        q_all += q
        one = np.ones(m)
        info = float(one @ Q @ one)
        offset, se = float(one @ Qr) / info, info ** -0.5
        if nu is None:
            stat, lo, hi = q / m, chi2.ppf(0.005, m) / m, chi2.ppf(0.995, m) / m
            p_off = 2 * norm.sf(abs(offset) / se)
        else:
            c = (nu - 2) / nu
            stat, lo, hi = q / (c * m), f_dist.ppf(0.005, m, nu), f_dist.ppf(0.995, m, nu)
            p_off = 2 * t_dist.sf(abs(offset) / se / np.sqrt(c), df=nu)
        sites.append(dict(site=s, n=m, misfit_per_n=q / m, law_stat=stat, law_lo=lo, law_hi=hi,
                          inside_99=bool(lo <= stat <= hi), offset=offset, offset_se=se, offset_z=offset / se,
                          offset_p=float(p_off), offset_unusual=bool(p_off < OFFSET_P)))
    shallow = blocks.depth <= SHALLOW
    lo, hi = misfit_law_range([len(i) for i in blocks.idx], nu)
    return dict(misfit_per_n=q_all / blocks.n, misfit_law_lo=lo, misfit_law_hi=hi,
                misfit_inside_99=bool(lo <= q_all / blocks.n <= hi),
                loo_z2_shallow=float(z2[shallow].mean()) if shallow.any() else None,
                loo_z2_deep=float(z2[~shallow].mean()) if (~shallow).any() else None, n=blocks.n,
                sites_outside_99=[x['site'] for x in sites if not x['inside_99']],
                sites_unusual_offset=[x['site'] for x in sites if x['offset_unusual']]), sites
