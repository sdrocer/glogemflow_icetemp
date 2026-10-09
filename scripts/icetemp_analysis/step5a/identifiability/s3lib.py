"""Stage 3a library: targets with their error covariance, two sensitivity methods (quadratic response
surface, SVD + PCGPwM emulator), Fisher summaries and best points. Analysis only.

Units: parameters in prior-sd units (f_ref 0.2887, f_ins 0.3259, f_adv 0.2444; campaign f_ref 0.2598),
temperatures in C. Residual convention for fits: r = y - G(theta). Fisher F = J^T Sigma^-1 J + I.
"""
import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.linalg import cho_factor, cho_solve, solve_triangular
from scipy.optimize import minimize

import target_analysis as ta
from icetemp.calibration import emulator as EM
from icetemp.calibration.priors import Priors

warnings.filterwarnings('ignore')

# The package Emulator assembles its matrix from glacier dicts; here each training output is already
# a row vector over the target's rows. Patched in this process only (no package file changes).
EM.Emulator.assemble_matrix = lambda self, outs: np.column_stack([np.asarray(o, float) for o in outs])

OUT = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/identifiability')
NAMES = ('f_ref', 'f_ins', 'f_adv')
SD = ta.prior_sd()                 # 0.2887, 0.3259, 0.2444
SD_CAMPAIGN = ta.prior_sd(True)    # 0.2598 for f_ref
PRIORS = Priors()
H_FD = 0.1                         # emulator finite-difference half step [prior sd]
TARGETS = ('T1', 'T2', 'T3r')


# ── data ──────────────────────────────────────────────────────────────────────────────────────
def load_all():
    obs = ta.load_obs()
    mods = {e: ta.load_model(e) for e in ta.ENSEMBLES}
    common = ta.common_rows(list(mods.values()), obs)
    return obs, mods, common


def build_target(name, obs, mod, common, ds_run=ta.BEST_RUN):
    """Target with its emulated rows and error covariance Sigma0 = D + S (+ month-to-month term and
    date spread for T3r / literal T3). T1: rows are entities; `base` holds the T2 rows the model is
    emulated on and A the pooling matrix (entities x common T2 rows)."""
    if name == 'T3lit':
        t = ta.target('T3', obs, mod, common)
        t.name = 'T3lit'
        S0 = ta.target_error_cov(ta.target('T3', obs, mod, common), obs, seasonal=ta.seasonal_cov)
        return SimpleNamespace(name=name, rows=t.rows, y=t.y, Y=t.Y, S0=S0, A=None, base=None)
    t = ta.target(name, obs, mod, common)
    if name == 'T1':
        c2 = common[0]
        S0 = ta.target_error_cov(t, obs)
        A = t.A[:, c2]
        t2 = ta.target('T2', obs, mod, common)
        return SimpleNamespace(name=name, rows=t.rows, y=t.y, Y=t.Y, S0=S0, A=A, base=t2)
    if name == 'T2':
        return SimpleNamespace(name=name, rows=t.rows, y=t.y, Y=t.Y, S0=ta.target_error_cov(t, obs), A=None, base=None)
    if name == 'T3r':
        S0 = (ta.target_error_cov(t, obs, seasonal=ta.seasonal_cov)
              + np.diag(ta.date_spread_var(t.rows, obs, mod['Y3'][ds_run])))
        return SimpleNamespace(name=name, rows=t.rows, y=t.y, Y=t.Y, S0=S0, A=None, base=None)
    raise ValueError(name)


def emu_rows(tg):
    """(model values the sensitivity methods are fitted on, pooling matrix or None)."""
    return (tg.base.Y, tg.A) if tg.A is not None else (tg.Y, None)


# ── method (i): quadratic response surface ─────────────────────────────────────────────────────
def quad_basis(X):
    x0, x1, x2 = X.T
    return np.column_stack([np.ones(len(X)), x0, x1, x2, x0 ** 2, x1 ** 2, x2 ** 2, x0 * x1, x0 * x2, x1 * x2])


class Surface:
    """Quadratic in (f_ref, f_ins, f_adv) per row, least squares over the given design points."""

    def __init__(self, theta, Y, sd=SD):
        self.sd = np.asarray(sd, float)
        self.c = theta.mean(axis=0)
        M = quad_basis((theta - self.c) / self.sd)
        self.coef, *_ = np.linalg.lstsq(M, Y, rcond=None)
        fit = M @ self.coef
        ss_res = ((Y - fit) ** 2).sum(axis=0)
        ss_tot = ((Y - Y.mean(axis=0)) ** 2).sum(axis=0)
        with np.errstate(invalid='ignore', divide='ignore'):
            self.r2 = 1.0 - ss_res / ss_tot
        Hm = M @ np.linalg.pinv(M)
        self.loo = (Y - fit) / (1.0 - np.diag(Hm))[:, None]       # PRESS residuals (design x rows)
        self.spread = Y.std(axis=0)
        self.ss_res, self.ss_tot = ss_res, ss_tot

    def predict(self, theta):
        return quad_basis((np.atleast_2d(theta) - self.c) / self.sd) @ self.coef

    def jac(self, theta):
        """(rows x 3): dG / d theta_k in C per prior sd (self.sd)."""
        x = (np.asarray(theta, float) - self.c) / self.sd
        b = self.coef
        g0 = b[1] + 2 * b[4] * x[0] + b[7] * x[1] + b[8] * x[2]
        g1 = b[2] + 2 * b[5] * x[1] + b[7] * x[0] + b[9] * x[2]
        g2 = b[3] + 2 * b[6] * x[2] + b[8] * x[0] + b[9] * x[1]
        return np.column_stack([g0, g1, g2])


# ── method (ii): the package emulator (SVD + surmise PCGPwM) ─────────────────────────────────
class Emu:
    """icetemp.calibration.emulator.Emulator on the target's rows, with the pipeline's
    calibrate_variance (LOO inflation curve); E(theta) = L L^T + truncation variance I, as
    calibrator.residual_and_covariance builds it."""

    def __init__(self, theta, Y, ev=0.999, calibrate=True, seed=None):
        self.emu = EM.Emulator(glaciers=[], explained_variance=ev)
        self.emu.fit(theta, list(Y))
        if calibrate:
            self.emu.calibrate_variance(theta, list(Y))
        self.trunc = self.emu.truncation_variance()
        self.basis = self.emu.basis
        self.p = self.emu.p
        self.n = Y.shape[1]

    def mean(self, theta):
        return self.emu.predict(np.atleast_2d(theta))[0]

    def coeff(self, theta):
        """(mean (n_theta x rows), coefficient variance incl. inflation (n_theta x p))."""
        theta = np.atleast_2d(theta)
        cm, cv = self.emu.backend.predict(self.emu._slice(theta))
        mean = self.emu.mean_[None, :] + cm @ self.basis.T
        curve = getattr(self.emu, '_variance_inflation_curve', None)
        infl = curve.predict(self.emu._nn_distance(theta)) if curve is not None else np.full(len(theta), self.emu.variance_inflation_)
        return mean, np.maximum(cv * infl[:, None], 0.0)

    def jac(self, theta, h=H_FD, sd=SD):
        theta = np.asarray(theta, float)
        J = np.zeros((self.n, 3))
        for k in range(3):
            tp, tm = theta.copy(), theta.copy()
            tp[k] += h * sd[k]
            tm[k] -= h * sd[k]
            mp, mm = self.mean(np.vstack([tp, tm]))
            J[:, k] = (mp - mm) / (2 * h)
        return J

    def E(self, theta):
        _, cv = self.coeff(theta)
        L = self.basis * np.sqrt(cv[0])[None, :]
        return L @ L.T + self.trunc * np.eye(self.n)


def loo_errors(theta, Y, ev=0.999):
    """Emulator leave-one-out prediction errors (design x rows): fit without i, predict at theta_i."""
    out = np.zeros_like(Y)
    for i in range(len(theta)):
        m = np.arange(len(theta)) != i
        e = Emu(theta[m], Y[m], ev=ev, calibrate=False)
        out[i] = e.mean(theta[i])[0] - Y[i]
    return out


# ── Fisher summaries ──────────────────────────────────────────────────────────────────────────
def fisher_data(J, Sigma):
    c = cho_factor(Sigma, lower=True)
    return J.T @ cho_solve(c, J)


def summarize(Fd_full, free=(0, 1)):
    """Posterior summary from the data information (3x3, prior-sd units) for the free parameters:
    spread left (posterior sd / prior sd), correlations, condition number of the posterior precision,
    weakest direction and its posterior sd; data-only correlation f_ref-f_ins."""
    free = list(free)
    Fd = Fd_full[np.ix_(free, free)]
    Fp = Fd + np.eye(len(free))
    C = np.linalg.inv(Fp)
    sdv = np.sqrt(np.diag(C))
    corr = C / np.outer(sdv, sdv)
    ev, evec = np.linalg.eigh(Fp)
    w = evec[:, 0] * np.sign(evec[np.argmax(np.abs(evec[:, 0])), 0])
    evd = np.linalg.eigvalsh(Fd)
    try:
        Cd = np.linalg.inv(Fd)
        corr_d = Cd[0, 1] / np.sqrt(Cd[0, 0] * Cd[1, 1]) if Cd[0, 0] > 0 and Cd[1, 1] > 0 else np.nan
    except np.linalg.LinAlgError:
        corr_d = np.nan
    out = dict(corr_ref_ins=float(corr[0, 1]), corr_data_only=float(corr_d),
               cond_post=float(ev.max() / ev.min()), cond_data=float(evd.max() / max(evd.min(), 1e-12)),
               weak_sd=float(1 / np.sqrt(ev[0])), weak_dir=' '.join(f'{v:+.2f}' for v in w))
    for j, k in enumerate(free):
        out[f'spread_{NAMES[k]}'] = float(sdv[j])
    if len(free) == 3:
        out['corr_ref_adv'] = float(corr[0, 2])
        out['corr_ins_adv'] = float(corr[1, 2])
    out['_w'] = w
    return out


def summarize_J(J, Sigma, free=(0, 1)):
    return summarize(fisher_data(J, Sigma), free)


def told_apart(s):
    return (s['spread_f_ins'] <= 0.5) and (abs(s['corr_ref_ins']) <= 0.8)



# ── design coverage (the calibrator's gate) ───────────────────────────────────────────────────
class Gate:
    """BayesianCalibrator._within_design_coverage on a design: inside the design's box for the free
    parameters and within 2x the median nearest-neighbour spacing (standardised by the design's own
    sd, as Emulator.calibrate_variance sets max_trusted_distance_) of some design point."""

    def __init__(self, theta_design, mult=2.0):
        self.T = np.asarray(theta_design, float)
        self.std = self.T.std(axis=0)
        self.std[self.std == 0] = 1.0
        Z = self.T / self.std
        pair = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
        np.fill_diagonal(pair, np.inf)
        self.spacing = float(np.median(pair.min(axis=1)))
        self.trust = mult * self.spacing
        self.lo, self.hi = self.T.min(axis=0), self.T.max(axis=0)

    def nn(self, theta):
        theta = np.atleast_2d(theta)
        d = (self.T[None, :, :] - theta[:, None, :]) / self.std
        return np.sqrt((d ** 2).sum(-1)).min(axis=1)

    def ok(self, theta, free=(0, 1)):
        theta = np.atleast_2d(theta)
        f = list(free)
        inbox = np.all((theta[:, f] >= self.lo[f] - 1e-12) & (theta[:, f] <= self.hi[f] + 1e-12), axis=1)
        return inbox & (self.nn(theta) <= self.trust)

# ── objective, best points ───────────────────────────────────────────────────────────────────
class Whitener:
    def __init__(self, Sigma):
        self.Lc = np.linalg.cholesky(Sigma)
        self.logdet = 2.0 * np.sum(np.log(np.diag(self.Lc)))

    def chi2(self, R):
        R = np.atleast_2d(R)
        Z = solve_triangular(self.Lc, R.T, lower=True)
        return np.sum(Z ** 2, axis=0)

    def solve(self, B):
        return cho_solve((self.Lc, True), B)


def log_prior(theta):
    return np.array([PRIORS.logpdf(t) for t in np.atleast_2d(theta)])


class SurfaceObjective:
    """-log posterior with G = surface (pooled with A for T1), Sigma = Sigma0 (+ extra fixed term)."""

    def __init__(self, tg, surf, extra=None):
        self.tg, self.surf = tg, surf
        self.W = Whitener(tg.S0 if extra is None else tg.S0 + extra)

    def G(self, theta):
        P = self.surf.predict(theta)
        return P @ self.tg.A.T if self.tg.A is not None else P

    def nlp(self, theta):
        R = self.tg.y[None, :] - self.G(theta)
        return 0.5 * self.W.chi2(R) + 0.5 * self.W.logdet - log_prior(theta)


class EmuObjective:
    """-log posterior with G = emulator mean, Sigma(theta) = Sigma0 + extra + E(theta) (Woodbury on
    the rank-p factor; truncation variance folded into the fixed part). T1: pooled with A."""

    def __init__(self, tg, emu, extra=None, use_E=True):
        self.tg, self.emu, self.use_E = tg, emu, use_E
        A = tg.A
        n = len(tg.y)
        S = tg.S0.copy() if extra is None else tg.S0 + extra
        if use_E:
            S = S + emu.trunc * (A @ A.T if A is not None else np.eye(n))
        self.S = S
        self.W = Whitener(S)
        self.B = A @ emu.basis if A is not None else emu.basis           # rows x p
        self.P0 = self.W.solve(self.B)
        self.G0 = self.B.T @ self.P0

    def G(self, theta):
        mean, _ = self.emu.coeff(theta)
        return mean @ self.tg.A.T if self.tg.A is not None else mean

    def Sigma(self, theta):
        if not self.use_E:
            return self.S
        _, cv = self.emu.coeff(theta)
        Lb = self.B * np.sqrt(cv[0])[None, :]
        return self.S + Lb @ Lb.T

    def nlp(self, theta):
        theta = np.atleast_2d(theta)
        mean, cv = self.emu.coeff(theta)
        G = mean @ self.tg.A.T if self.tg.A is not None else mean
        R = self.tg.y[None, :] - G
        chi = self.W.chi2(R)
        logdet = np.full(len(theta), self.W.logdet)
        if self.use_E:
            U = R @ self.P0                                   # n_theta x p
            for i in range(len(theta)):
                s = np.sqrt(cv[i])
                Mm = np.eye(len(s)) + s[:, None] * self.G0 * s[None, :]
                us = s * U[i]
                chi[i] -= us @ np.linalg.solve(Mm, us)
                logdet[i] += np.linalg.slogdet(Mm)[1]
        return 0.5 * chi + 0.5 * logdet - log_prior(theta)


def design_box(theta_keep):
    return theta_keep.min(axis=0), theta_keep.max(axis=0)


def best_point(obj, gate, free=(0, 1), adv_fixed=1.0, grid_step=0.01, starts3=None):
    """MAP of obj inside the gate (f_adv = adv_fixed unless free). 2p: grid (step 0.01) over the
    design box, then a bounded Nelder-Mead polish; 3p: Nelder-Mead from every design point (and the
    extra starts). Returns (theta, -log posterior, grid (2p) or None)."""
    free = list(free)
    lo, hi = gate.lo[free], gate.hi[free]

    def full(x):
        t = np.array([0.0, 0.0, adv_fixed])
        t[free] = x
        return t

    def f(x):
        t = full(x)
        if not gate.ok(t, free)[0]:
            return 1e12
        return float(obj.nlp(t)[0])

    grid = None
    if len(free) == 2:
        g0 = np.arange(gate.lo[0], gate.hi[0] + 1e-9, grid_step)
        g1 = np.arange(gate.lo[1], gate.hi[1] + 1e-9, grid_step)
        GG = np.array([[a, b, adv_fixed] for a in g0 for b in g1])
        okm = gate.ok(GG, free)
        vals = np.full(len(GG), np.inf)
        idx = np.flatnonzero(okm)
        for i in range(0, len(idx), 1500):
            j = idx[i:i + 1500]
            vals[j] = obj.nlp(GG[j])
        grid = (g0, g1, vals.reshape(len(g0), len(g1)))
        starts = [GG[np.argmin(vals), :2]]
    else:
        starts = [t for t in gate.T] + ([] if starts3 is None else list(starts3))
    best = None
    for s0 in starts:
        x0 = np.asarray(s0, float)
        x0 = x0[free] if len(x0) == 3 else x0
        r = minimize(f, np.clip(x0, lo, hi), method='Nelder-Mead', bounds=list(zip(lo, hi)),
                     options={'maxiter': 4000, 'xatol': 1e-5, 'fatol': 1e-6})
        if best is None or r.fun < best.fun:
            best = r
    return full(best.x), float(best.fun), grid


def grid_posterior(grid):
    """Posterior mean, sd and f_ref-f_ins correlation on the 2p grid of -log posterior values (cells
    inside the gate; uniform cells)."""
    g0, g1, V = grid
    P = np.where(np.isfinite(V), np.exp(-(V - np.nanmin(np.where(np.isfinite(V), V, np.nan)))), 0.0)
    P /= P.sum()
    A0, A1 = np.meshgrid(g0, g1, indexing='ij')
    m0, m1 = (P * A0).sum(), (P * A1).sum()
    v0, v1 = (P * (A0 - m0) ** 2).sum(), (P * (A1 - m1) ** 2).sum()
    c01 = (P * (A0 - m0) * (A1 - m1)).sum()
    return dict(mean_f_ref=float(m0), mean_f_ins=float(m1), spread_f_ref=float(np.sqrt(v0) / SD[0]),
                spread_f_ins=float(np.sqrt(v1) / SD[1]), corr=float(c01 / np.sqrt(v0 * v1)))


def chi2_runs(tg, idx):
    """GLS misfit r^T Sigma0^-1 r of the real runs (rows of tg.Y) and their RMSE."""
    W = Whitener(tg.S0)
    R = tg.Y[idx] - tg.y[None, :]
    return W.chi2(R), np.sqrt(np.mean(R ** 2, axis=1))


# ── where the information comes from ─────────────────────────────────────────────────────────
def info_groups(J, Sigma, labels, free=(0, 1), order=None):
    """Per group of rows: rows, additive share of the data information (f_ref, f_ins diagonal and
    along the weakest posterior direction; row contributions J_i (Sigma^-1 J)_i summed over the
    group), and the posterior summary from the group's rows alone and from all rows without it."""
    labels = np.asarray(labels, dtype=object)
    f = list(free)
    Jf = J[:, f]
    SJ = cho_solve(cho_factor(Sigma, lower=True), Jf)
    Fd = Jf.T @ SJ
    ev, evec = np.linalg.eigh(Fd + np.eye(len(f)))
    w = evec[:, 0]
    out = []
    groups = order if order is not None else sorted(set(labels))
    for g in groups:
        m = labels == g
        if not m.any():
            continue
        c = Jf[m].T @ SJ[m]
        c = 0.5 * (c + c.T)
        rec = dict(group=g, rows=int(m.sum()), share_f_ref=float(c[0, 0] / Fd[0, 0]),
                   share_f_ins=float(c[1, 1] / Fd[1, 1]), share_weak=float(w @ c @ w / (w @ Fd @ w)))
        if len(f) == 3:
            rec['share_f_adv'] = float(c[2, 2] / Fd[2, 2])
        a = summarize_J(J[m], Sigma[np.ix_(m, m)], free)
        rec.update({f'alone_{k}': v for k, v in a.items() if k in ('corr_ref_ins', 'spread_f_ref', 'spread_f_ins', 'spread_f_adv')})
        if (~m).any():
            b = summarize_J(J[~m], Sigma[np.ix_(~m, ~m)], free)
            rec.update({f'without_{k}': v for k, v in b.items() if k in ('corr_ref_ins', 'spread_f_ref', 'spread_f_ins', 'spread_f_adv')})
        out.append(rec)
    return out


# ── per-glacier offset delta (sensitivity only) ───────────────────────────────────────────────
def delta_cov(obs, t1, y1_pred):
    """K_delta as BayesianCalibrator.fit_discrepancy builds it: Discrepancy GP (Matern 3/2 on
    lat/lon/elevation + nugget) fitted to the entity residuals (obs - model) at a representative
    theta; the pure-signal kernel at the entities. Returns (K (entities x entities), kernel text)."""
    from icetemp.calibration.discrepancy import Discrepancy, latlon_to_xyz
    ents = t1.rows.e.to_numpy(int)
    gl = [obs['glaciers'][e] for e in ents]
    lat = np.array([g.latitude for g in gl]); lon = np.array([g.longitude for g in gl])
    elev = np.array([g.elevation for g in gl])
    resid = t1.y - y1_pred
    # pooled observation variance per entity (independent per point), as compute_glacier_residuals
    sig = t1.base.rows.sigma.to_numpy(float)        # common T2 rows, the columns of t1.A
    s2 = []
    for row in t1.A:
        m = row > 0
        s2.append(np.sum(row[m] ** 2 * sig[m] ** 2))
    disc = Discrepancy(nugget_floor={'T_residual': float(np.median(s2))})
    gp = disc.fit('T_residual', lat, lon, resid, elevations=elev)
    k = gp.kernel_.k1 if hasattr(gp.kernel_, 'k1') else gp.kernel_
    return k(latlon_to_xyz(lat, lon, elev)), str(gp.kernel_)


def rows_to_entities(rows, t1):
    """Z (rows x entities of t1): 1 where the row belongs to the entity."""
    ents = t1.rows.e.to_numpy(int)
    return (rows['e'].to_numpy(int)[:, None] == ents[None, :]).astype(float)
