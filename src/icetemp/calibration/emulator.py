"""
Emulator: SVD dimensionality reduction + surmise PCGPwM Gaussian-process emulator.

Two explicit stages, matching the task spec's separate "Dimensionality Reduction" and
"Emulation" deliverables (rather than relying only on PCGPwM's own internal PCA, which would
work but wouldn't expose an interpretable basis or match the requested modular structure):

  1. SVD: stack each training run's flattened depth-temperature profile (evaluated at the
     concatenated (glacier, depth) locations from ALL calibration glaciers -- fixed across
     runs, since firnice_glenglat_lookup pins GloGEM's profile output at each glenglat
     borehole's own observed depths) into a matrix F (n_x, m); reduce to p dominant left
     singular vectors (the "basis") by an explained-variance threshold.
  2. Emulation: fit surmise.emulation.emulator(method='PCGPwM') mapping theta (m, 3) to the
     p-dimensional SVD *coefficients* C (p, m) -- NOT the raw high-dimensional profile -- so
     the GP emulates a small, decorrelated target. predict() reconstructs the full profile as
     mean + basis @ coefficients.

G itself (the expensive simulator being emulated) is the real transient IDL GloGEM model, run
via runner.GloGEMRunner -- this module only consumes its (theta, profile) training pairs.
"""

from dataclasses import dataclass, field

import numpy as np


def _build_x_index(glaciers):
    """Fixed list of (glacier_name, depth) pairs spanning every calibration glacier's own
    observed depths -- the shared 'x' locations every training run's profile must align to
    (guaranteed by firnice_glenglat_lookup placing GloGEM's profile output at these same
    borehole depths on every run, regardless of theta).

    Keyed by glacier_name, NOT glacier_id: multiple GlacierCalibrationData objects can share
    one GloGEM glacier_id (two glenglat boreholes on the same RGI polygon), but glacier_name is
    always unique -- DataHandler groups measurements by glacier_name in the first place. See
    runner.GloGEMRunner.parse_training_output, which disambiguates same-id output files by
    elevation and returns its dict keyed the same way, by glacier_name."""
    index = []
    for g in glaciers:
        for d in g.depths:
            index.append((g.glacier_name, float(d)))
    return index


class EmulatorBackend:
    """Minimal interface an emulator backend must satisfy, so a fallback (e.g. independent
    per-coefficient sklearn GPs) can substitute if surmise is unavailable or fails to build
    for a given numpy/scipy combination (see environment.yaml's numpy<2.2/scipy<1.15 pin note
    -- surmise hard-caps both)."""

    def fit(self, theta, C):
        raise NotImplementedError

    def predict(self, theta):
        """Return (mean, var), each shape (n_theta, p)."""
        raise NotImplementedError


class SurmisePCGPwM(EmulatorBackend):
    def __init__(self, **surmise_kwargs):
        self.surmise_kwargs = surmise_kwargs
        self._emu = None
        self._p = None

    def fit(self, theta, C):
        from surmise.emulation import emulator
        p, m = C.shape
        self._p = p
        x = np.arange(p, dtype=float).reshape(-1, 1)  # 'x' = coefficient index
        self._emu = emulator(x=x, theta=theta, f=C, method='PCGPwM', **self.surmise_kwargs)
        return self

    def predict(self, theta):
        pred = self._emu.predict(x=np.arange(self._p, dtype=float).reshape(-1, 1), theta=theta)
        mean = np.asarray(pred.mean()).T    # surmise returns (p, n_theta) -> (n_theta, p)
        var = np.asarray(pred.var()).T
        return mean, var


class IndependentGPFallback(EmulatorBackend):
    """One independent sklearn GaussianProcessRegressor per SVD coefficient. Simpler and
    weaker than PCGPwM (no cross-coefficient correlation modelling) but has no numpy/scipy
    version constraints -- used only if surmise is unavailable."""

    def __init__(self, **gp_kwargs):
        self.gp_kwargs = gp_kwargs
        self._gps = []

    def fit(self, theta, C):
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
        p, m = C.shape
        self._gps = []
        for k in range(p):
            kernel = ConstantKernel(1.0, (1e-3, 1e3)) * Matern(length_scale=1.0, nu=2.5) \
                + WhiteKernel(noise_level=1e-3)
            gp = GaussianProcessRegressor(kernel=kernel, normalize_y=True,
                                           n_restarts_optimizer=5, **self.gp_kwargs)
            gp.fit(theta, C[k, :])
            self._gps.append(gp)
        return self

    def predict(self, theta):
        means, stds = [], []
        for gp in self._gps:
            m_k, s_k = gp.predict(theta, return_std=True)
            means.append(m_k)
            stds.append(s_k)
        mean = np.column_stack(means)
        var = np.column_stack(stds) ** 2
        return mean, var


@dataclass
class Emulator:
    glaciers: list
    explained_variance: float = 0.99
    backend: EmulatorBackend = None

    x_index: list = field(init=False, default=None)
    basis: np.ndarray = field(init=False, default=None)   # (n_x, p)
    mean_: np.ndarray = field(init=False, default=None)    # (n_x,)
    singular_values_: np.ndarray = field(init=False, default=None)
    p: int = field(init=False, default=None)
    theta_train_: np.ndarray = field(init=False, default=None)

    def __post_init__(self):
        self.x_index = _build_x_index(self.glaciers)
        if self.backend is None:
            try:
                import surmise  # noqa: F401
                self.backend = SurmisePCGPwM()
            except ImportError:
                self.backend = IndependentGPFallback()

    @property
    def n_x(self):
        return len(self.x_index)

    def assemble_matrix(self, training_output):
        """Build F (n_x, m) from a list of per-run outputs (one dict glacier_name -> (depths,
        T_values) per training run, e.g. from GloGEMRunner.parse_training_output, in the same
        order as the theta design matrix). Rows follow self.x_index; a run missing a
        glacier's output leaves those rows as NaN (dropped before SVD, with a warning)."""
        m = len(training_output)
        F = np.full((self.n_x, m), np.nan)
        for j, run in enumerate(training_output):
            for i, (gname, depth) in enumerate(self.x_index):
                entry = run.get(gname)
                if entry is None:
                    continue
                depths, T_values = entry
                depths = np.asarray(depths, dtype=float)
                k = np.argmin(np.abs(depths - depth))
                if abs(depths[k] - depth) < 1e-6:
                    F[i, j] = T_values[k]
        return F

    def fit_svd(self, F):
        """Explicit SVD dimensionality reduction (the 'Dimensionality Reduction' deliverable):
        reduce F's rows (n_x depth-temperature locations) to p dominant basis vectors chosen
        by cumulative explained-variance. Drops any (x) rows that are NaN in every training
        run and any (theta) columns that are NaN anywhere, so a partially-failed design matrix
        (e.g. some IDL runs crashed) can still be used."""
        row_ok = ~np.all(np.isnan(F), axis=1)
        col_ok = ~np.any(np.isnan(F[row_ok]), axis=0)
        F_use = F[np.ix_(row_ok, col_ok)]
        self._row_mask = row_ok
        self._col_mask = col_ok

        mean_ = F_use.mean(axis=1)
        Fc = F_use - mean_[:, None]
        U, S, _ = np.linalg.svd(Fc, full_matrices=False)
        cumvar = np.cumsum(S ** 2) / np.sum(S ** 2)
        p = int(np.searchsorted(cumvar, self.explained_variance) + 1)
        p = min(p, U.shape[1])

        self.mean_ = mean_
        self.basis = U[:, :p]
        self.singular_values_ = S
        self.p = p
        return self.basis, mean_, S

    def fit(self, design_theta, training_output):
        """Full fit: SVD (fit_svd) then emulate theta -> coefficients (backend.fit)."""
        F = self.assemble_matrix(training_output)
        basis, mean_, S = self.fit_svd(F)
        F_use = F[np.ix_(self._row_mask, self._col_mask)]
        theta_use = np.asarray(design_theta)[self._col_mask]
        C = basis.T @ (F_use - mean_[:, None])   # (p, m_used) SVD coefficients

        self.theta_train_ = theta_use
        self.backend.fit(theta_use, C)
        return self

    def predict(self, theta):
        """Predict the full profile (at self.x_index locations, restricted to rows kept by
        fit_svd) for one or more theta points. Returns (mean, var), each (n_theta, n_x_used).
        """
        theta = np.atleast_2d(theta)
        coeff_mean, coeff_var = self.backend.predict(theta)   # (n_theta, p) each
        mean = self.mean_[None, :] + coeff_mean @ self.basis.T
        # Var of a linear combination of independent coefficient predictions:
        # Var(basis @ c) = basis^2 @ Var(c) (ignoring cross-coefficient covariance, which
        # PCGPwM's own predict() does not expose either).
        var = (self.basis ** 2) @ coeff_var.T
        return mean, var.T

    def x_index_used(self):
        return [xi for xi, keep in zip(self.x_index, self._row_mask) if keep]

    def explained_variance_ratio(self):
        return np.cumsum(self.singular_values_ ** 2) / np.sum(self.singular_values_ ** 2)
