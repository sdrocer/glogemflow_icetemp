"""Stage 3c scratch library: response surfaces and emulators over the alps_s2 design, the GLS
misfit of a target (or a subset of its rows) on a (f_ref, f_ins) grid with f_adv = 1.0, and
posterior summaries in prior-sd units.

Methods (both trained on the 29 runs with f_ref >= 0.1):
  surface   full quadratic in (f_ref, f_ins, f_adv) per row, parameters in prior-sd units; the
            whole prior box f_ref 0.1-1.0, f_ins 0.3-2.0 (above f_ins 1.79 it extrapolates)
  emulator  SVD of the rows keeping 99.9 % of the variance + surmise PCGPwM (the pipeline's
            backend), its mean, restricted to the pipeline's trust region (design box, 2 x the
            median nearest-neighbour spacing).
  emulator_E  as emulator, plus its uncertainty E = basis diag(coefficient variance x a flat
            leave-one-run-out inflation) basis^T + truncation variance, evaluated at the emulator
            optimum and held fixed (a theta-dependent E pulls the optimum into unsampled corners,
            the pathology the pipeline's trust gate documents). Used for reference only: with 29
            runs and none near the optimum, E there is large and coherent.
Both predictions are capped at the pressure-melting point of the row's depth (the model's clamp).

Misfit: chi2(theta) = r^T Sigma^-1 r, r = model - observed; Sigma = D + S (T2) or
D + S + seasonal + date spread (T3r, Stage 2), plus E for emulator_E. Rows may be one-sided
(melting point): r = min(model - min(obs, PMP), 0).
Posterior on the grid: exp(-chi2/2) x prior, f_ref U[0.1, 1] (campaign prior), f_ins truncated
normal (1, 0.35) on [0.3, 2]; f_adv = 1.0. Spreads are posterior sd / prior sd (f_ref 0.2598,
f_ins 0.3259); "lin" values are the linearised (Fisher + unit prior precision) summary at the MAP.
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '4')
import sys
import warnings
from pathlib import Path
from types import SimpleNamespace

REPO = Path('/home/jabeer/projects/glogemflow_development/glogemflow_icetemp')
for _p in (REPO / 'scripts' / 'icetemp_analysis', REPO / 'src'):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve, solve_triangular

import target_analysis as ta
from icetemp.calibration.priors import Priors
from icetemp.calibration.thermal_structure import pressure_melting_point

OUT = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/checks')
PSD = ta.prior_sd(campaign=True)            # f_ref U[0.1, 1] 0.2598, f_ins 0.3259, f_adv 0.2444
PRIOR = Priors(bounds_override={'refreeze_frac': (0.1, 1.0)})
FREF = np.round(np.linspace(0.10, 1.00, 91), 6)
FINS = np.round(np.linspace(0.30, 2.00, 86), 6)
F_ADV = 1.0


def grid(fref=FREF, fins=FINS):
    return np.array([(a, b, F_ADV) for a in fref for b in fins])


def log_prior(TH):
    """Log prior density of (f_ref, f_ins) at f_adv fixed (f_ref uniform: constant)."""
    return PRIOR.insul_scale.logpdf(TH[:, 1])


# ── data ──────────────────────────────────────────────────────────────────────────────────────
def load():
    obs = ta.load_obs()
    mod = ta.load_model('alps_s2')
    keep = mod['design'].keep.to_numpy()
    theta = mod['theta']
    return SimpleNamespace(obs=obs, mod=mod, keep=keep, theta=theta,
                           t2=ta.target('T2', obs, mod), t3=ta.target('T3r', obs, mod))


def pmp(depth):
    return pressure_melting_point(np.asarray(depth, dtype=float))


def sigma_T2(rows):
    return ta.error_cov(rows)


def sigma_T3r(rows, obs, y3_d13):
    """D + S + seasonal (Stage 2) + date spread (design 13) for T3r rows."""
    C = ta.error_cov(rows, seasonal=ta.seasonal_cov)
    return C + np.diag(ta.date_spread_var(rows, obs, y3_d13))


# ── surrogates ────────────────────────────────────────────────────────────────────────────────
class Surface:
    """Full quadratic response surface per row (least squares over the training runs)."""

    def __init__(self, theta, Y, cap=None):
        self.c = theta.mean(axis=0)
        M = self._basis(theta)
        self.coef, *_ = np.linalg.lstsq(M, Y, rcond=None)
        fit = M @ self.coef
        H = M @ np.linalg.pinv(M)
        self.loo = (Y - fit) / (1 - np.diag(H))[:, None]          # PRESS residuals (runs x rows)
        self.loo_rmse = np.sqrt((self.loo ** 2).mean(axis=0))
        ss_tot = ((Y - Y.mean(axis=0)) ** 2).sum(axis=0)
        with np.errstate(invalid='ignore', divide='ignore'):
            self.r2 = 1 - ((Y - fit) ** 2).sum(axis=0) / ss_tot
        self.cap = cap

    def _basis(self, th):
        x = (np.atleast_2d(th) - self.c) / PSD
        return np.column_stack([np.ones(len(x)), x[:, 0], x[:, 1], x[:, 2], x[:, 0] ** 2, x[:, 1] ** 2,
                                x[:, 2] ** 2, x[:, 0] * x[:, 1], x[:, 0] * x[:, 2], x[:, 1] * x[:, 2]])

    def predict(self, th, cap=True):
        out = self._basis(th) @ self.coef
        if cap and self.cap is not None:
            out = np.minimum(out, self.cap[None, :])
        return out

    def subset(self, m):
        s = object.__new__(Surface)
        s.c, s.coef, s.loo, s.loo_rmse, s.r2 = self.c, self.coef[:, m], self.loo[:, m], self.loo_rmse[m], self.r2[m]
        s.cap = None if self.cap is None else self.cap[m]
        return s


class Emu:
    """SVD (explained variance ev) + surmise PCGPwM on the coefficients, as emulator.Emulator."""

    def __init__(self, theta, Y, ev=0.999, cap=None):
        from icetemp.calibration.emulator import SurmisePCGPwM
        F = Y.T                                                  # rows x runs
        self.mu = F.mean(axis=1)
        U, S, _ = np.linalg.svd(F - self.mu[:, None], full_matrices=False)
        cv = np.cumsum(S ** 2) / np.sum(S ** 2)
        self.p = int(min(np.searchsorted(cv, ev) + 1, U.shape[1]))
        self.basis = U[:, :self.p]
        self.S = S
        C = self.basis.T @ (F - self.mu[:, None])
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            self.be = SurmisePCGPwM().fit(theta, C)
        self.trunc = float((S[self.p:] ** 2).sum() / (F.shape[1] * F.shape[0])) if self.p < len(S) else 0.0
        self.cap = cap
        self.infl = 1.0

    def predict(self, th, cap=True):
        """(mean (n_th x rows), coefficient variance (n_th x p) x inflation)."""
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            cm, cvar = self.be.predict(np.atleast_2d(th))
        mean = self.mu[None, :] + cm @ self.basis.T
        if cap and self.cap is not None:
            mean = np.minimum(mean, self.cap[None, :])
        return mean, np.maximum(cvar, 0.0) * self.infl


# ── misfit on a grid ──────────────────────────────────────────────────────────────────────────
def hinge(R, y, one_sided, cap):
    """Residuals with one-sided rows: min(model - min(obs, PMP), 0)."""
    if one_sided is None or not np.any(one_sided):
        return R
    R = R.copy()
    yb = np.minimum(y[one_sided], cap[one_sided])
    G = R[:, one_sided] + y[None, one_sided]
    R[:, one_sided] = np.minimum(G - yb[None, :], 0.0)
    return R


def chi2_fixed(R, Sigma):
    """r^T Sigma^-1 r for each row of R (n_theta x n)."""
    L = np.linalg.cholesky(Sigma)
    W = solve_triangular(L, R.T, lower=True)
    return (W ** 2).sum(axis=0), float(2 * np.log(np.diag(L)).sum())


def chi2_emu(R, Sigma0, basis, cvar, trunc):
    """r^T (Sigma0 + trunc I + B diag(v) B^T)^-1 r and its log det, per theta (Woodbury)."""
    S0 = Sigma0 + trunc * np.eye(len(Sigma0))
    L = np.linalg.cholesky(S0)
    WR = solve_triangular(L, R.T, lower=True)                  # n x nth
    WB = solve_triangular(L, basis, lower=True)                # n x p
    q0 = (WR ** 2).sum(axis=0)
    BtW_r = WB.T @ WR                                          # p x nth
    G = WB.T @ WB                                              # p x p
    ld0 = float(2 * np.log(np.diag(L)).sum())
    out = np.empty(R.shape[0])
    ld = np.empty(R.shape[0])
    for i in range(R.shape[0]):
        v = cvar[i]
        sv = np.sqrt(v)
        K = np.eye(len(v)) + sv[:, None] * G * sv[None, :]
        c = cho_factor(K, lower=True)
        u = sv * BtW_r[:, i]
        out[i] = q0[i] - u @ cho_solve(c, u)
        ld[i] = ld0 + 2 * np.log(np.diag(c[0])).sum()
    return out, ld


def trust_mask_fn(theta_train, mult=2.0):
    """The pipeline's emulator guard (calibrator log_posterior + Emulator.calibrate_variance):
    theta inside the training design's box for the free parameters (f_ref, f_ins) and within
    mult x the median nearest-neighbour spacing (per-dimension standardised) of a design point."""
    sd = theta_train.std(axis=0)
    Z = theta_train / sd
    pair = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(pair, np.inf)
    maxd = mult * float(np.median(pair.min(axis=1)))
    lo, hi = theta_train.min(axis=0), theta_train.max(axis=0)

    def f(TH):
        nn = np.sqrt((((TH / sd)[:, None, :] - Z[None, :, :]) ** 2).sum(-1)).min(axis=1)
        inbox = (TH[:, 0] >= lo[0]) & (TH[:, 0] <= hi[0]) & (TH[:, 1] >= lo[1]) & (TH[:, 1] <= hi[1])
        return inbox & (nn <= maxd)
    return f


def _predict(sur, TH, kind):
    return sur.predict(TH) if kind == 'surface' else sur.predict(TH)[0]


def log_post(sur, y, Sigma, TH, one_sided=None, kind='surface', trust=None):
    """(log posterior, chi2) on the grid TH; Sigma fixed (an emulator's E is plugged in)."""
    R = hinge(_predict(sur, TH, kind) - y[None, :], y, one_sided, sur.cap)
    c2, _ = chi2_fixed(R, Sigma)
    lp = -0.5 * c2 + log_prior(TH)
    if trust is not None:
        lp = np.where(trust(TH), lp, -np.inf)
    return lp, c2


def summarize_grid(TH, logpost):
    w = np.exp(logpost - logpost.max())
    w /= w.sum()
    m = w @ TH[:, :2]
    d = TH[:, :2] - m
    C = (w[:, None] * d).T @ d
    sd = np.sqrt(np.diag(C))
    k = int(np.argmax(logpost))
    return dict(map_f_ref=float(TH[k, 0]), map_f_ins=float(TH[k, 1]), mean_f_ref=float(m[0]), mean_f_ins=float(m[1]),
                sd_f_ref=float(sd[0]), sd_f_ins=float(sd[1]), corr=float(C[0, 1] / (sd[0] * sd[1])),
                spread_f_ref=float(sd[0] / PSD[0]), spread_f_ins=float(sd[1] / PSD[1]), cov=C)


def fit(sur, y, Sigma, one_sided=None, kind='surface', trust=None, refine=True):
    """Grid posterior (coarse grid, then a fine grid over mean +- 6 sd). Returns a dict."""
    TH = grid()
    lp, c2 = log_post(sur, y, Sigma, TH, one_sided, kind, trust)
    s = summarize_grid(TH, lp)
    if refine:
        lo_r, hi_r = max(0.1, s['mean_f_ref'] - 6 * s['sd_f_ref']), min(1.0, s['mean_f_ref'] + 6 * s['sd_f_ref'])
        lo_i, hi_i = max(0.3, s['mean_f_ins'] - 6 * s['sd_f_ins']), min(2.0, s['mean_f_ins'] + 6 * s['sd_f_ins'])
        TH2 = grid(np.linspace(lo_r, hi_r, 61), np.linspace(lo_i, hi_i, 61))
        lp2, c22 = log_post(sur, y, Sigma, TH2, one_sided, kind, trust)
        if np.isfinite(lp2).sum() > 50:
            TH, lp, c2 = TH2, lp2, c22
            s = summarize_grid(TH, lp)
    k = int(np.argmax(lp))
    s['chi2_map'] = float(c2[k])
    s['n'] = int(len(y))
    s['ratio_map'] = float(c2[k] / len(y))
    s['grid'] = (TH, lp)
    return s


def linearised(sur, y, Sigma, th0, one_sided=None, kind='surface', h=0.05):
    """Fisher information of (f_ref, f_ins) at th0 (prior-sd units), Gaussian prior with unit
    precision: corr, spread left. Central differences of the surrogate (capped)."""
    J = np.zeros((len(y), 2))
    for k in range(2):
        tp, tm = np.array(th0, float), np.array(th0, float)
        tp[k] += h * PSD[k]
        tm[k] -= h * PSD[k]
        if kind == 'surface':
            gp, gm = sur.predict(tp)[0], sur.predict(tm)[0]
        else:
            gp, gm = sur.predict(tp)[0][0], sur.predict(tm)[0][0]
        J[:, k] = (gp - gm) / (2 * h)
    if one_sided is not None and np.any(one_sided):
        g0 = sur.predict(np.array(th0, float))[0] if kind == 'surface' else sur.predict(np.array(th0, float))[0][0]
        yb = np.minimum(y, sur.cap) if sur.cap is not None else y
        inactive = one_sided & (g0 >= yb)
        J[inactive] = 0.0
    c = cho_factor(Sigma, lower=True)
    F = J.T @ cho_solve(c, J)
    C = np.linalg.inv(F + np.eye(2))
    sd = np.sqrt(np.diag(C))
    return dict(corr=float(C[0, 1] / (sd[0] * sd[1])), spread_f_ref=float(sd[0]), spread_f_ins=float(sd[1]), F=F, J=J)


def rmse(a, b):
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


def shift_sd(s_a, s_b, ref=None):
    """Optimum shift b - a in prior-sd units, and in units of ref's posterior sd per parameter;
    Mahalanobis distance with the summed posterior covariances (independent data sets)."""
    da = np.array([s_b['mean_f_ref'] - s_a['mean_f_ref'], s_b['mean_f_ins'] - s_a['mean_f_ins']])
    ref = s_a if ref is None else ref
    C = s_a['cov'] + s_b['cov']
    return dict(d_f_ref_prior_sd=float(da[0] / PSD[0]), d_f_ins_prior_sd=float(da[1] / PSD[1]),
                d_f_ref_post_sd=float(da[0] / np.sqrt(ref['cov'][0, 0])),
                d_f_ins_post_sd=float(da[1] / np.sqrt(ref['cov'][1, 1])),
                mahal_sum=float(np.sqrt(da @ np.linalg.solve(C, da))))


# ── context: both targets with their surrogates and error covariances ─────────────────────────
def context(targets=('T2', 'T3r'), Y_override=None):
    """dict name -> SimpleNamespace(rows, y, Y, Sigma, sur, emu, cap). Y_override: {name: Y}
    (design points x rows) to replace the model values (e.g. depth-shifted rows)."""
    import json
    infl = json.load(open(OUT / 'surrogates.json'))
    d = load()
    out = {}
    for name in targets:
        t = d.t2 if name == 'T2' else d.t3
        Y = t.Y if Y_override is None or name not in Y_override else Y_override[name]
        rows = t.rows.copy()
        rows['grp'] = ta.group_of(rows)
        cap = pmp(rows.depth.to_numpy())
        Sigma = sigma_T2(rows) if name == 'T2' else sigma_T3r(rows, d.obs, d.mod['Y3'][ta.BEST_RUN])
        sur = Surface(d.theta[d.keep], Y[d.keep], cap)
        emu = Emu(d.theta[d.keep], Y[d.keep], 0.999, cap)
        emu.infl = infl[name]['inflation']
        out[name] = SimpleNamespace(name=name, rows=rows, y=t.y.copy(), Y=Y, Sigma=Sigma, sur=sur, emu=emu, cap=cap,
                                    d=d)
    return out


class EmuSub:
    """An emulator restricted to a subset of rows (same GP, rows of the basis)."""

    def __init__(self, emu, m):
        self.e, self.m = emu, np.asarray(m)
        self.basis = emu.basis[self.m]
        self.trunc = emu.trunc
        self.cap = None if emu.cap is None else emu.cap[self.m]

    def predict(self, th, cap=True):
        mean, cvar = self.e.predict(th, cap=False)
        mean = mean[:, self.m]
        if cap and self.cap is not None:
            mean = np.minimum(mean, self.cap[None, :])
        return mean, cvar


def emu_cov(sur, th):
    """E at theta for an (EmuSub) emulator: basis diag(coefficient variance x inflation) basis^T
    plus the truncation variance on the diagonal."""
    _, cvar = sur.predict(np.asarray(th, float))
    return (sur.basis * cvar[0][None, :]) @ sur.basis.T + sur.trunc * np.eye(len(sur.basis))


def fit_rows(ctx, m, one_sided=None, methods=('surface', 'emulator'), Sigma=None, y=None):
    """Fit (f_ref, f_ins) at f_adv = 1 on the rows m (bool mask) of a target context.
    surface: quadratic surface over the prior box. emulator: SVD-GP mean inside the pipeline's
    trust region. emulator_E: as emulator, with the emulator's E at the emulator optimum plugged
    in (second pass). Each: grid posterior plus the linearised summary at the MAP."""
    m = np.asarray(m, bool)
    S = ctx.Sigma[np.ix_(m, m)] if Sigma is None else Sigma
    yy = ctx.y[m] if y is None else y
    os_ = None if one_sided is None else np.asarray(one_sided, bool)[m]
    res = {}
    trust = trust_mask_fn(ctx.d.theta[ctx.d.keep])
    for kind in methods:
        if kind == 'surface':
            sur = ctx.sur.subset(m)
            s = fit(sur, yy, S, os_, kind)
            Sl = S
        else:
            sur = EmuSub(ctx.emu, np.flatnonzero(m))
            s = fit(sur, yy, S, os_, 'emulator', trust)
            Sl = S
            if kind == 'emulator_E':
                E = emu_cov(sur, [s['map_f_ref'], s['map_f_ins'], F_ADV])
                Sl = S + E
                s1 = s
                s = fit(sur, yy, Sl, os_, 'emulator', trust)
                s['pass1_map'] = (s1['map_f_ref'], s1['map_f_ins'])
                s['E_median_var'] = float(np.median(np.diag(E)))
                s['E_top_eig'] = float(np.linalg.eigvalsh(E)[-1])
        th0 = np.array([s['map_f_ref'], s['map_f_ins'], F_ADV])
        lin = linearised(sur, yy, Sl, th0, os_, 'surface' if kind == 'surface' else 'emulator')
        s['lin_corr'], s['lin_spread_f_ref'], s['lin_spread_f_ins'] = lin['corr'], lin['spread_f_ref'], lin['spread_f_ins']
        s['J'] = lin['J']
        res[kind] = s
    return res


def fmt(s):
    return (f"MAP f_ref {s['map_f_ref']:.3f} f_ins {s['map_f_ins']:.3f} | mean {s['mean_f_ref']:.3f}+-{s['sd_f_ref']:.3f} "
            f"{s['mean_f_ins']:.3f}+-{s['sd_f_ins']:.3f} | corr {s['corr']:+.2f} spread {s['spread_f_ref']:.2f}/{s['spread_f_ins']:.2f} "
            f"| lin corr {s['lin_corr']:+.2f} spread {s['lin_spread_f_ref']:.2f}/{s['lin_spread_f_ins']:.2f} | chi2/n {s['ratio_map']:.2f}")


def predict_at(ctx, s, kind='surface'):
    th = np.array([s['map_f_ref'], s['map_f_ins'], F_ADV])
    return ctx.sur.predict(th)[0] if kind == 'surface' else ctx.emu.predict(th)[0][0]


# ── reporting helpers ─────────────────────────────────────────────────────────────────────────
def strip(res):
    """Fit results without the grid and Jacobian (for pickling/JSON)."""
    return {k: {kk: (vv.tolist() if isinstance(vv, np.ndarray) else vv) for kk, vv in v.items()
                if kk not in ('grid', 'J')} for k, v in res.items()}


def row(label, target, res, base=None, extra=None):
    """One flat record per method: optimum, spread, corr, shift vs base (prior sd, base post sd)."""
    out = []
    for kind, s in res.items():
        rec = dict(variant=label, target=target, method=kind, n=s['n'], map_f_ref=s['map_f_ref'], map_f_ins=s['map_f_ins'],
                   mean_f_ref=s['mean_f_ref'], mean_f_ins=s['mean_f_ins'], sd_f_ref=s['sd_f_ref'], sd_f_ins=s['sd_f_ins'],
                   corr=s['corr'], spread_f_ref=s['spread_f_ref'], spread_f_ins=s['spread_f_ins'],
                   lin_corr=s['lin_corr'], lin_spread_f_ref=s['lin_spread_f_ref'], lin_spread_f_ins=s['lin_spread_f_ins'],
                   chi2_n=s['ratio_map'])
        if base is not None and kind in base:
            b = base[kind]
            rec['d_f_ref_prior_sd'] = (s['mean_f_ref'] - b['mean_f_ref']) / PSD[0]
            rec['d_f_ins_prior_sd'] = (s['mean_f_ins'] - b['mean_f_ins']) / PSD[1]
            rec['d_f_ref_base_sd'] = (s['mean_f_ref'] - b['mean_f_ref']) / b['sd_f_ref']
            rec['d_f_ins_base_sd'] = (s['mean_f_ins'] - b['mean_f_ins']) / b['sd_f_ins']
            da = np.array([s['mean_f_ref'] - b['mean_f_ref'], s['mean_f_ins'] - b['mean_f_ins']])
            rec['mahal_base'] = float(np.sqrt(da @ np.linalg.solve(np.asarray(b['cov']), da)))
        if extra:
            rec.update(extra.get(kind, {}) if isinstance(extra, dict) and kind in extra else {})
        out.append(rec)
    return out


def table(recs, cols=None):
    df = pd.DataFrame(recs)
    cols = cols or ['variant', 'target', 'method', 'n', 'map_f_ref', 'map_f_ins', 'mean_f_ref', 'mean_f_ins', 'corr',
                    'spread_f_ref', 'spread_f_ins', 'lin_corr', 'lin_spread_f_ins', 'd_f_ref_prior_sd', 'd_f_ins_prior_sd',
                    'mahal_base', 'chi2_n']
    cols = [c for c in cols if c in df]
    return df[cols]


def fit_fins_given_fref(ctx, m, f_ref, one_sided=None, kind='surface', Sigma=None, y=None):
    """1-D posterior of f_ins with f_ref fixed (f_adv = 1): mean, sd, MAP; also the data-only
    optimum (no prior). kind 'surface' or 'emulator' (mean, trust region)."""
    m = np.asarray(m, bool)
    S = ctx.Sigma[np.ix_(m, m)] if Sigma is None else Sigma
    yy = ctx.y[m] if y is None else y
    os_ = None if one_sided is None else np.asarray(one_sided, bool)[m]
    fins = np.linspace(0.3, 2.0, 341)
    TH = np.column_stack([np.full(len(fins), f_ref), fins, np.full(len(fins), F_ADV)])
    if kind == 'surface':
        sur = ctx.sur.subset(m); trust = None
    else:
        sur = EmuSub(ctx.emu, np.flatnonzero(m)); trust = trust_mask_fn(ctx.d.theta[ctx.d.keep])
    lp, c2 = log_post(sur, yy, S, TH, os_, 'surface' if kind == 'surface' else 'emulator', trust)
    w = np.exp(lp - lp.max()); w /= w.sum()
    mu = float(w @ fins); sd = float(np.sqrt(w @ (fins - mu) ** 2))
    ok = np.isfinite(lp)
    k_data = int(np.argmin(np.where(ok, c2, np.inf)))
    return dict(mean=mu, sd=sd, map=float(fins[int(np.argmax(lp))]), data_opt=float(fins[k_data]),
                dchi2_range=float(np.nanmax(np.where(ok, c2, np.nan)) - c2[k_data]), n=int(m.sum()),
                at_bound=bool(fins[k_data] in (fins[0], fins[-1]) or not ok[min(k_data + 1, len(fins) - 1)]))


def rmse_block(ctx2, s, masks2, kind='surface'):
    """Unweighted RMSE over T2 common points (and named subsets) at a fit's MAP, from the T2
    surface (or emulator mean): {name: rmse}."""
    th = np.array([s['map_f_ref'], s['map_f_ins'], F_ADV])
    G = ctx2.sur.predict(th)[0] if kind == 'surface' else ctx2.emu.predict(th)[0][0]
    return {k: rmse(G[m], ctx2.y[m]) for k, m in masks2.items() if np.any(m)}


def run_variants(ctx, variants, base, masks2=None, methods=('surface', 'emulator')):
    """variants: [(label, target, mask, one_sided or None)]; base: baseline fits {target: {method: s}}.
    Returns (records, fits). Records carry RMSE over T2 points at the variant's MAP (masks2)."""
    recs, fits = [], {}
    masks2 = masks2 or {'T2 all': np.ones(len(ctx['T2'].y), bool)}
    for label, tname, m, os_ in variants:
        c = ctx[tname]
        r = fit_rows(c, m, os_, methods)
        fits[(label, tname)] = strip(r)
        rr = row(label, tname, r, base[tname])
        for rec in rr:
            rm = rmse_block(ctx['T2'], r[rec['method']], masks2, rec['method'])
            rec.update({f'rmse {k}': v for k, v in rm.items()})
        recs += rr
    return recs, fits
