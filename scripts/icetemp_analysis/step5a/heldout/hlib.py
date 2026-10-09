"""Stage 3b helpers: targets with their error covariance, surrogates (response surface, SVD+GP
emulator), per-entity GLS misfits on a parameter grid, fits with the prior as a penalty.

Conventions
  theta = (f_ref, f_ins) with f_adv = 1.0; prior: f_ref uniform (bounds only), f_ins truncated
  normal mean 1 sd 0.35 -> penalty ((f_ins - 1) / 0.35)^2 in the objective (= -2 log prior).
  Objective J(theta) = r^T Sigma^-1 r + penalty, r = model - obs (sign irrelevant).
  Fit box = design range of the 29 runs with f_ref >= 0.1 (no extrapolation).
  Prior-sd units: 0.2887 (f_ref, U[0,1]), 0.3259 (f_ins).
"""
import sys
import time
import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve, solve_triangular

REPO = Path('/home/jabeer/projects/glogemflow_development/glogemflow_icetemp')
sys.path.insert(0, str(REPO / 'scripts' / 'icetemp_analysis'))
import target_analysis as ta  # noqa: E402

OUT = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/heldout')
SD3 = ta.prior_sd()                    # 0.2887, 0.3259, 0.2444
SD = SD3[:2]
INS_MEAN, INS_SD = 1.0, 0.35           # f_ins prior (untruncated mean, sd)
F_ADV = 1.0
P_NAT = np.diag([0.0, 1.0 / INS_SD ** 2])   # prior curvature (natural units): uniform f_ref, normal f_ins


def penalty(th):
    th = np.atleast_2d(th)
    return ((th[:, 1] - INS_MEAN) / INS_SD) ** 2


# ── data ──────────────────────────────────────────────────────────────────────────────────────
def load_all():
    """obs, s2 model, joint common rows (771 points), design (32 runs)."""
    obs = ta.load_obs()
    s2 = ta.load_model('alps_s2')
    r7 = ta.load_model('alps_r7')
    common = ta.common_rows([s2, r7], obs)
    design = ta.load_design('alps_s2')
    return obs, s2, common, design


def build_targets(obs, mod, common, date_run=ta.BEST_RUN):
    """T1, T2, T3 (= T3r) with Sigma = D + S (+ seasonal + date spread for T3)."""
    t1 = ta.target('T1', obs, mod, common)
    t2 = ta.target('T2', obs, mod, common)
    t3 = ta.target('T3r', obs, mod, common)
    t1.Sigma = ta.target_error_cov(t1, obs)
    t2.Sigma = ta.target_error_cov(t2, obs)
    t3.Sigma = (ta.target_error_cov(t3, obs, seasonal=ta.seasonal_cov)
                + np.diag(ta.date_spread_var(t3.rows, obs, mod['Y3'][date_run])))
    # entity label per row: T1 rows are entities
    t1.ent = t1.rows['e'].to_numpy(int)
    t2.ent = t2.rows['e'].to_numpy(int)
    t3.ent = t3.rows['e'].to_numpy(int)
    for t in (t1, t2, t3):
        t.system = t.rows['system_short'].to_numpy(object)
        t.n = len(t.y)
    # T2 point index of each T3 row (position in t2.rows)
    pos = {p: i for i, p in enumerate(t2.rows['t2'].to_numpy(int))}
    t3.pt = np.array([pos[p] for p in t3.rows['t2'].to_numpy(int)])
    t2.pt = np.arange(t2.n)
    t1.pt = None
    # T1 pooling matrix restricted to the T2 rows (A over common T2 rows)
    t1.A2 = t1.A[:, common[0]]
    t1.name, t2.name, t3.name = 'T1', 'T2', 'T3'
    return {'T1': t1, 'T2': t2, 'T3': t3}


# ── per-entity GLS structures ─────────────────────────────────────────────────────────────────
class Blocks:
    """Cholesky factor per entity block of a block-diagonal Sigma (rows grouped by `ent`), and
    per-entity quadratic forms of many residual vectors at once. Optional level split:
    GLS entity mean (level) and, with `pts`, GLS point means (for T3 rows)."""

    def __init__(self, Sigma, ent, pts=None, check_block=True):
        self.ent = np.asarray(ent)
        self.ents = np.unique(self.ent)
        self.idx = [np.flatnonzero(self.ent == e) for e in self.ents]
        if check_block:
            mask = self.ent[:, None] != self.ent[None, :]
            off = np.abs(Sigma[mask]).max() if mask.any() else 0.0
            if off > 1e-12:
                raise ValueError(f'Sigma not block-diagonal by entity (max off-block {off})')
        self.L = [np.linalg.cholesky(Sigma[np.ix_(i, i)]) for i in self.idx]
        # level projector pieces: u = L^-1 1
        self.u = [solve_triangular(L, np.ones(len(i)), lower=True) for L, i in zip(self.L, self.idx)]
        self.pts = None
        if pts is not None:
            self.pts = []
            for L, i in zip(self.L, self.idx):
                p = np.asarray(pts)[i]
                up, inv = np.unique(p, return_inverse=True)
                Z = np.zeros((len(i), len(up)))
                Z[np.arange(len(i)), inv] = 1.0
                wZ = solve_triangular(L, Z, lower=True)
                Q, _ = np.linalg.qr(wZ)          # orthonormal basis of the whitened point space
                self.pts.append(Q)

    def parts(self, R):
        """R: (m, n) residual vectors. Returns dict of (m, n_entities) arrays: total, level (GLS
        entity mean), and if pts: points (GLS point means incl. level)."""
        R = np.atleast_2d(R)
        m = R.shape[0]
        tot = np.zeros((m, len(self.idx)))
        lev = np.zeros((m, len(self.idx)))
        pts = np.zeros((m, len(self.idx))) if self.pts is not None else None
        for k, (L, i, u) in enumerate(zip(self.L, self.idx, self.u)):
            W = solve_triangular(L, R[:, i].T, lower=True)       # (n_e, m)
            tot[:, k] = (W ** 2).sum(0)
            lev[:, k] = (u @ W) ** 2 / (u @ u)
            if pts is not None:
                Q = self.pts[k]
                pts[:, k] = ((Q.T @ W) ** 2).sum(0)
        out = dict(total=tot, level=lev)
        if pts is not None:
            out['points'] = pts
        return out


def target_blocks(t):
    return Blocks(t.Sigma, t.ent, pts=getattr(t, 'pt', None) if t.name == 'T3' else None)


# ── surrogates ────────────────────────────────────────────────────────────────────────────────
def poly_basis(theta, form, center):
    """Polynomial basis in standardised parameters x = (theta - center) / prior sd (3 columns).
    forms: 'Q3' full quadratic in 3; 'C2Q' cubic in (f_ref, f_ins) + f_adv, f_adv^2, f_adv x f_ref,
    f_adv x f_ins; 'C2L' cubic in (f_ref, f_ins) + f_adv; 'Q2L' quadratic in (f_ref, f_ins) + f_adv
    and its two interactions."""
    x = (np.atleast_2d(theta) - center) / SD3
    a, b, c = x[:, 0], x[:, 1], x[:, 2]
    one = np.ones(len(x))
    cols = [one, a, b, c]
    if form == 'Q3':
        cols += [a * a, b * b, c * c, a * b, a * c, b * c]
    elif form == 'Q2L':
        cols += [a * a, b * b, a * b, a * c, b * c]
    elif form == 'C2L':
        cols += [a * a, b * b, a * b, a ** 3, b ** 3, a * a * b, a * b * b]
    elif form == 'C2Q':
        cols += [a * a, b * b, a * b, a ** 3, b ** 3, a * a * b, a * b * b, c * c, a * c, b * c]
    elif form == 'L3':
        pass
    else:
        raise ValueError(form)
    return np.column_stack(cols)


class RespSurf:
    """OLS polynomial per row, fitted on the given runs; exact LOO residuals via the hat matrix."""

    def __init__(self, form='Q3'):
        self.form = form

    def fit(self, theta3, Y):
        self.center = theta3.mean(0)
        M = poly_basis(theta3, self.form, self.center)
        self.coef, *_ = np.linalg.lstsq(M, Y, rcond=None)
        H = M @ np.linalg.pinv(M)
        fit = M @ self.coef
        self.loo = (Y - fit) / (1 - np.diag(H))[:, None]       # held-out minus predicted (sign: Y - pred)
        self.loo = -self.loo                                     # predicted - true
        return self

    def predict(self, theta2):
        th = np.atleast_2d(theta2)
        th3 = np.column_stack([th[:, 0], th[:, 1], np.full(len(th), F_ADV)]) if th.shape[1] == 2 else th
        return poly_basis(th3, self.form, self.center) @ self.coef

    def jac(self, theta2, h=1e-4):
        """d model / d theta (natural units), (n_rows, 2), central differences."""
        th = np.asarray(theta2, float)
        J = np.zeros((self.coef.shape[1], 2))
        for k in range(2):
            tp, tm = th.copy(), th.copy()
            tp[k] += h
            tm[k] -= h
            J[:, k] = (self.predict(tp)[0] - self.predict(tm)[0]) / (2 * h)
        return J


class SvdGp:
    """SVD (explained variance ev) + surmise PCGPwM on the coefficients, 3 inputs."""

    def __init__(self, ev=0.999):
        self.ev = ev

    def fit(self, theta3, Y):
        from icetemp.calibration.emulator import SurmisePCGPwM
        F = Y.T
        self.mu = F.mean(axis=1)
        U, S, _ = np.linalg.svd(F - self.mu[:, None], full_matrices=False)
        cv = np.cumsum(S ** 2) / np.sum(S ** 2)
        self.p = int(np.searchsorted(cv, self.ev) + 1)
        self.basis = U[:, :self.p]
        C = self.basis.T @ (F - self.mu[:, None])
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            self.be = SurmisePCGPwM().fit(theta3, C)
        self.S = S
        return self

    def predict(self, theta2, return_var=False, chunk=4000):
        th = np.atleast_2d(theta2)
        th3 = np.column_stack([th[:, 0], th[:, 1], np.full(len(th), F_ADV)]) if th.shape[1] == 2 else th
        outs, vs = [], []
        for s in range(0, len(th3), chunk):
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                cm, cv = self.be.predict(th3[s:s + chunk])
            outs.append(self.mu[None, :] + cm @ self.basis.T)
            vs.append(cv)
        mean = np.vstack(outs)
        if return_var:
            return mean, np.vstack(vs)
        return mean

    def jac(self, theta2, h=0.01):
        th = np.asarray(theta2, float)
        J = np.zeros((len(self.mu), 2))
        for k in range(2):
            tp, tm = th.copy(), th.copy()
            tp[k] += h * SD[k]
            tm[k] -= h * SD[k]
            J[:, k] = (self.predict(tp)[0] - self.predict(tm)[0]) / (2 * h * SD[k])
        return J


def emu_loo(theta3, Y, ev=0.999):
    """Leave-one-run-out predictions of the SVD+GP emulator: (pred - true) per run."""
    out = np.zeros_like(Y)
    for i in range(len(Y)):
        m = np.arange(len(Y)) != i
        e = SvdGp(ev).fit(theta3[m], Y[m])
        out[i] = e.predict(theta3[i:i + 1])[0] - Y[i]
    return out


# ── grid ──────────────────────────────────────────────────────────────────────────────────────
def make_grid(lo, hi, step=0.005):
    fr = np.arange(lo[0], hi[0] + 1e-9, step)
    fi = np.arange(lo[1], hi[1] + 1e-9, step)
    G = np.array([(a, b) for a in fr for b in fi])
    return G, fr, fi


def timer():
    t0 = time.time()
    return lambda: f'{time.time() - t0:.0f} s'


def fisher_parts(B, J):
    """Per-entity Fisher contributions of model Jacobian J (n_rows x 2) under Blocks B:
    dict total / level / points -> (n_entities, 2, 2)."""
    out = dict(total=[], level=[])
    if B.pts is not None:
        out['points'] = []
    for k, (L, i, u) in enumerate(zip(B.L, B.idx, B.u)):
        W = solve_triangular(L, J[i], lower=True)        # (n_e, 2)
        out['total'].append(W.T @ W)
        a = (u @ W) / np.sqrt(u @ u)
        out['level'].append(np.outer(a, a))
        if B.pts is not None:
            Q = B.pts[k]
            P = Q.T @ W
            out['points'].append(P.T @ P)
    return {k: np.array(v) for k, v in out.items()}


def block_part(M, ent):
    """M restricted to its within-entity blocks (zero across entities)."""
    ent = np.asarray(ent)
    return np.where(ent[:, None] == ent[None, :], M, 0.0)


# ── fits on the precomputed grid (h03_grid.pkl) ───────────────────────────────────────────────
METHODS = (('RS', 'E'), ('RS', 'noE'), ('EMU', 'E'), ('EMU', 'noE'), ('RSQ3', 'E'), ('DISC', 'noE'))
METHOD_LABEL = {('RS', 'E'): 'response surface (main)', ('RS', 'noE'): 'response surface, no E',
                ('EMU', 'E'): 'emulator', ('EMU', 'noE'): 'emulator, no E',
                ('RSQ3', 'E'): 'quadratic surface (prototype form)', ('DISC', 'noE'): 'best of the 29 real runs'}


def load_grid():
    import pickle
    with open(OUT / 'h03_grid.pkl', 'rb') as f:
        return pickle.load(f)


def objective_arrays(D, m, ev, target, part='total'):
    """(C per entity on the candidate set, penalty per candidate, held-out SSE per entity, candidate
    thetas). part: 'total', 'level' (GLS entity mean), 'shape' (total - level), 'points' (GLS point
    means), 'depthshape' (points - level), 'm2m' (total - points)."""
    if m == 'DISC':
        src, TH, sse = D['discrete']['noE'], D['discrete']['theta'][:, :2], D['discrete']['sse']
    else:
        src, TH, sse = D['methods'][m][ev], D['G'], D['methods'][m]['sse']
    p = src[target]
    C = {'total': lambda: p['total'], 'level': lambda: p['level'],
         'shape': lambda: p['total'] - p['level'], 'points': lambda: p['points'],
         'depthshape': lambda: p['points'] - p['level'], 'm2m': lambda: p['total'] - p['points']}[part]()
    return C, penalty(TH), sse, TH


def argmin_w(C, pen, w):
    J = C @ w + pen
    k = int(np.argmin(J))
    return k, J


def on_edge(th, lo=(0.10, 0.30), hi=(0.98, 1.79), tol=1e-9):
    th = np.asarray(th)
    return bool((np.abs(th - lo) < tol).any() or (np.abs(th - hi) < tol).any())


class Setup:
    """Data, targets, surrogates and E per method, for subset fits on the grid and Gauss-Newton
    posterior covariances."""

    def __init__(self, methods=('RS', 'EMU', 'RSQ3')):
        self.obs, self.s2, self.common, self.design = load_all()
        self.T = build_targets(self.obs, self.s2, self.common)
        self.D = load_grid()
        self.G = self.D['G']
        self.pen = penalty(self.G)
        keep = self.design['keep'].to_numpy()
        th3 = self.design[['f_ref', 'f_ins', 'f_adv']].to_numpy()
        self.sur = {}
        for m in methods:
            if m == 'EMU':
                self.sur[m] = {k: SvdGp(0.999).fit(th3[keep], self.T[k].Y[keep]) for k in ('T2', 'T3')}
            else:
                form = 'C2Q' if m == 'RS' else 'Q3'
                self.sur[m] = {k: RespSurf(form).fit(th3[keep], self.T[k].Y[keep]) for k in ('T2', 'T3')}
        self._pred = {}

    def E(self, m, target, ev='E'):
        """Emulator term: within-entity blocks of the leave-one-run-out error covariance of the
        method's surrogate (as h03_grid.py); T1: its pooled diagonal."""
        n = self.T[target].n
        if ev == 'noE':
            return np.zeros((n, n))
        if not hasattr(self, '_E'):
            self._E = {}
        if m not in self._E:
            if m == 'EMU':
                import pickle
                loo = pickle.load(open(OUT / 'h01_loo.pkl', 'rb'))['loos']
                L2, L3 = loo['T2']['EMU999'], loo['T3']['EMU999']
            else:
                L2, L3 = self.sur[m]['T2'].loo, self.sur[m]['T3'].loo
            E2f, E3f = L2.T @ L2 / len(L2), L3.T @ L3 / len(L3)
            t1, t2, t3 = self.T['T1'], self.T['T2'], self.T['T3']
            self._E[m] = dict(T1=np.diag(np.diag(t1.A2 @ E2f @ t1.A2.T)), T2=block_part(E2f, t2.ent),
                              T3=block_part(E3f, t3.ent))
        return self._E[m][target]

    def pred_grid(self, m, target):
        key = (m, target)
        if key not in self._pred:
            self._pred = {k: v for k, v in self._pred.items() if k[0] == m}   # keep memory bounded
            self._pred[key] = self.sur[m][target].predict(self.G)
        return self._pred[key]

    def jac(self, m, target, th):
        return self.sur[m][target].jac(np.asarray(th, float))

    # ── level (T1-like) on a subset of T2 rows ──
    def level_t1(self, m, mask2=None, ev='E', th=None):
        """T1-like level over the T2 rows in mask2: per-entity chi2 on the grid (or Fisher at th).
        Returns (C (n_grid x n_ent_sub) or F_e (n_ent_sub x 2 x 2), entity ids)."""
        t2 = self.T['T2']
        mask2 = np.ones(t2.n, bool) if mask2 is None else np.asarray(mask2, bool)
        A, ents = ta.t1_pool(t2.rows, mask2)
        S1 = A @ (t2.Sigma + self.E(m, 'T2', ev)) @ A.T
        v = np.diag(S1)
        if th is None:
            M = self.pred_grid(m, 'T2')
            R = M @ A.T - A @ t2.y
            return R ** 2 / v[None, :], ents
        J = A @ self.jac(m, 'T2', th)
        return np.einsum('ei,ej->eij', J, J) / v[:, None, None], ents

    # ── GLS parts of T2 or T3 on a subset of rows ──
    def parts_sub(self, m, target, mask=None, ev='E', th=None):
        """Per-entity GLS parts (total, level, points) over the rows in mask, on the grid, or the
        per-entity Fisher parts at th. Returns (dict, entity ids)."""
        t = self.T[target]
        mask = np.ones(t.n, bool) if mask is None else np.asarray(mask, bool)
        Sig = (t.Sigma + self.E(m, target, ev))[np.ix_(mask, mask)]
        B = Blocks(Sig, t.ent[mask], pts=t.pt[mask] if target == 'T3' else None)
        if th is None:
            M = self.pred_grid(m, target)[:, mask]
            return B.parts(M - t.y[mask]), B.ents
        J = self.jac(m, target, th)[mask]
        return fisher_parts(B, J), B.ents


def part_of(p, part):
    return {'total': lambda: p['total'], 'level': lambda: p['level'], 'shape': lambda: p['total'] - p['level'],
            'points': lambda: p['points'], 'depthshape': lambda: p['points'] - p['level'],
            'm2m': lambda: p['total'] - p['points']}[part]()


def post_cov(F, prior='exact'):
    """Posterior covariance (natural units) from Fisher F (2x2) and the prior curvature:
    'exact' = uniform f_ref (0) and normal f_ins (1/0.35^2); 'unit' = 1 per prior sd for both."""
    P = P_NAT if prior == 'exact' else np.diag(1.0 / SD ** 2)
    return np.linalg.inv(F + P)


def sd_units(C):
    return C / np.outer(SD, SD)


def mahal(th_a, C_a, th_b, C_b):
    """Distance between two optima in combined sd: sqrt(d^T (C_a + C_b)^-1 d), prior-sd units."""
    d = (np.asarray(th_a) - np.asarray(th_b)) / SD
    Cs = sd_units(C_a) + sd_units(C_b)
    return float(np.sqrt(d @ np.linalg.solve(Cs, d)))


def grid_moments(J, G):
    """Posterior mean and covariance of exp(-J/2) on the grid (natural units)."""
    w = np.exp(-0.5 * (J - J.min()))
    w /= w.sum()
    mu = w @ G
    d = G - mu
    return mu, (d * w[:, None]).T @ d


def summarize_post(C):
    Cs = sd_units(C)
    s = np.sqrt(np.diag(Cs))
    return dict(sd_ref=float(s[0]), sd_ins=float(s[1]), corr=float(Cs[0, 1] / s[0] / s[1]))
