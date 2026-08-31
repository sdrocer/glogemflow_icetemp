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


BASAL = 'basal'  # sentinel x_index depth tag, see _build_x_index and assemble_matrix


def _build_x_index(glaciers):
    """Fixed list of (glacier_name, depth) pairs spanning every calibration glacier's own
    observed depths -- the shared 'x' locations every training run's profile must align to
    (guaranteed by firnice_glenglat_lookup placing GloGEM's profile output at these same
    borehole depths on every run, regardless of theta) -- plus one BASAL sentinel row per
    glacier (see assemble_matrix: filled from the model's own deepest RESOLVED grid node,
    only when that's genuinely deeper information than any Track-1 row already covers).

    Keyed by glacier_name, NOT glacier_id: multiple GlacierCalibrationData objects can share
    one GloGEM glacier_id (two glenglat boreholes on the same RGI polygon), but glacier_name is
    always unique -- DataHandler groups measurements by glacier_name in the first place. See
    runner.GloGEMRunner.parse_training_output, which disambiguates same-id output files by
    elevation and returns its dict keyed the same way, by glacier_name."""
    index = []
    for g in glaciers:
        for d in g.depths:
            index.append((g.glacier_name, float(d)))
        index.append((g.glacier_name, BASAL))
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
    # Which columns of theta = (perm_frac, dT_scale, z0) the GP actually sees. Callers still
    # pass and receive FULL 3-vectors -- the slicing happens internally here, so calibrator /
    # validation / writeback need no changes.
    #
    # WHY THIS EXISTS: from 2026-08-18 the calibration holds z0 fixed (it provably cannot reach
    # the real model's physics -- see BayesianCalibrator.fixed_params), which makes its design
    # column CONSTANT. A constant column is degenerate for the GP: confirmed live that
    # surmise PCGPwM then dies with `LinAlgError: Eigenvalues did not converge`, while the
    # identical training data with a varying z0 column fits fine. Adding tiny jitter to lift the
    # degeneracy "works" but is a TRAP: _nn_distance standardizes per dimension by that
    # dimension's own std, so jitter of sd ~0.05 m turns a 0.1 m z0 difference into ~1.8
    # standardized units and lets pure noise dominate the trust gate. Excluding the fixed
    # dimension outright is the only option that is correct in both the GP fit and the distance
    # metric. Leaving z0 free-but-varying would be wrong too -- points far apart in z0 but
    # identical in (perm_frac, dT_scale) are functionally the SAME run, and counting that
    # separation as real distance would inflate every nn_distance.
    active_params: tuple = (0, 1, 2)

    x_index: list = field(init=False, default=None)
    basis: np.ndarray = field(init=False, default=None)   # (n_x, p)
    mean_: np.ndarray = field(init=False, default=None)    # (n_x,)
    singular_values_: np.ndarray = field(init=False, default=None)
    p: int = field(init=False, default=None)
    theta_train_: np.ndarray = field(init=False, default=None)
    variance_inflation_: float = field(init=False, default=1.0)   # see calibrate_variance
    loo_z_std_: float = field(init=False, default=None)            # diagnostic only

    def __post_init__(self):
        self.x_index = _build_x_index(self.glaciers)
        self._glacier_max_obs = {
            g.glacier_name: float(g.depths.max()) for g in self.glaciers if len(g.depths)
        }
        self._variance_inflation_curve = None  # set by calibrate_variance(); see predict()
        if self.backend is None:
            try:
                import surmise  # noqa: F401
                self.backend = SurmisePCGPwM()
            except ImportError:
                self.backend = IndependentGPFallback()

    @property
    def n_x(self):
        return len(self.x_index)

    @staticmethod
    def _local_tolerance(depths, k):
        """Half the local grid spacing around depths[k] (GloGEM's own fit_dzstep = [1,5,20] m
        fine/medium/coarse bands, settings.pro:488-489) -- how close an observed depth must be
        to this grid node to accept it as that node's real-world counterpart. Not a fixed
        constant: derived from the grid actually returned, so it stays correct regardless of
        which region/grid definition produced it."""
        n = len(depths)
        if n < 2:
            return 1e-6
        if k == 0:
            return (depths[1] - depths[0]) / 2.0
        if k == n - 1:
            return (depths[k] - depths[k - 1]) / 2.0
        return min(depths[k] - depths[k - 1], depths[k + 1] - depths[k]) / 2.0

    def assemble_matrix(self, training_output):
        """Build F (n_x, m) from a list of per-run outputs (one dict glacier_name -> (depths,
        T_values) per training run, e.g. from GloGEMRunner.parse_training_output, in the same
        order as the theta design matrix). Rows follow self.x_index; a run missing a
        glacier's output leaves those rows as NaN (dropped before SVD, with a warning).

        Track-1 (depth) rows: matched to the model's nearest OWN grid depth, accepted within
        _local_tolerance -- NOT an exact match. GloGEM's output grid (1,2,...,9,14,19,...,259 m)
        essentially never lands exactly on an arbitrary glenglat measurement depth, so an exact-
        match requirement (the previous behaviour) silently discarded ~88% of otherwise-valid,
        within-resolved-depth comparisons purely from grid quantization, not genuine data loss.

        Track-2 (BASAL) rows: filled from the model's own deepest RESOLVED grid node for that
        glacier (the last non-NaN entry -- see firnice_temperature_model.pro's tl_fit[..,tt-2],
        a real geothermal-flux-driven conduction result, not a numerical artifact of
        truncation). Only filled when the glacier's deepest OBSERVATION exceeds that resolved
        node -- i.e. when this is genuinely new information Track-1 can't already provide.
        Otherwise left NaN, so fit_svd's existing all-NaN-row drop removes the redundant row
        automatically (same mechanism that already excludes glaciers with no training output
        at all, e.g. failed runs)."""
        m = len(training_output)
        F = np.full((self.n_x, m), np.nan)
        for j, run in enumerate(training_output):
            for i, (gname, depth) in enumerate(self.x_index):
                entry = run.get(gname)
                if entry is None:
                    continue
                depths, T_values = entry
                depths = np.asarray(depths, dtype=float)
                T_values = np.asarray(T_values, dtype=float)
                if depth == BASAL:
                    valid = np.flatnonzero(~np.isnan(T_values))
                    if valid.size == 0:
                        continue
                    deepest_k = valid.max()
                    obs_max = self._glacier_max_obs.get(gname)
                    if obs_max is not None and obs_max <= depths[deepest_k]:
                        continue  # Track-1 already covers this glacier's deepest observation
                    F[i, j] = T_values[deepest_k]
                else:
                    k = np.argmin(np.abs(depths - depth))
                    if abs(depths[k] - depth) <= self._local_tolerance(depths, k):
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
        theta_use = self._slice(np.asarray(design_theta)[self._col_mask])
        C = basis.T @ (F_use - mean_[:, None])   # (p, m_used) SVD coefficients

        self.theta_train_ = theta_use
        self.backend.fit(theta_use, C)
        return self

    def _slice(self, theta):
        """Restrict a theta array (n, 3) to the columns this emulator was trained on."""
        theta = np.atleast_2d(np.asarray(theta, dtype=float))
        if tuple(self.active_params) == (0, 1, 2):
            return theta
        return theta[:, list(self.active_params)]

    def _nn_distance(self, theta):
        """Standardized nearest-neighbor distance, per query theta, to self.theta_train_ (the
        emulator's actual training design) -- per-dimension standardized by theta_train_'s own
        std so perm_frac/dT_scale/z0's very different natural scales contribute comparably.
        A property of theta-space location alone, shared by every output row/glacier -- see
        calibrate_variance's docstring for why this replaced raw predicted variance as the
        covariate. Returns shape (n_theta,)."""
        theta = self._slice(np.atleast_2d(theta))
        theta_std = self.theta_train_.std(axis=0)
        theta_std[theta_std == 0] = 1.0
        diffs = (self.theta_train_[None, :, :] - theta[:, None, :]) / theta_std
        return np.sqrt(np.sum(diffs ** 2, axis=2)).min(axis=1)

    def predict(self, theta):
        """Predict the full profile (at self.x_index locations, restricted to rows kept by
        fit_svd) for one or more theta points. Returns (mean, var), each (n_theta, n_x_used).

        var is corrected by calibrate_variance's fitted curve if available (preferred --
        location-dependent on theta-space distance, see that method's docstring), else by the
        flat variance_inflation_ scalar (1.0 = no-op until either has been run, so an
        un-calibrated emulator behaves exactly as before). getattr with a default guards
        against unpickling an emulator saved before _variance_inflation_curve existed.
        """
        theta = np.atleast_2d(theta)
        coeff_mean, coeff_var = self.backend.predict(self._slice(theta))   # (n_theta, p) each
        mean = self.mean_[None, :] + coeff_mean @ self.basis.T
        # Var of a linear combination of independent coefficient predictions:
        # Var(basis @ c) = basis^2 @ Var(c) (ignoring cross-coefficient covariance, which
        # PCGPwM's own predict() does not expose either).
        var = ((self.basis ** 2) @ coeff_var.T).T

        curve = getattr(self, '_variance_inflation_curve', None)
        if curve is not None:
            nn_dist = self._nn_distance(theta)             # (n_theta,) -- one distance per
            inflation = curve.predict(nn_dist)              # theta, applied to ALL its rows:
            var = var * inflation[:, None]                  # extrapolation risk is a property
        else:                                               # of theta-space location, not of
            var = var * self.variance_inflation_             # which glacier/depth row it is.
        return mean, var

    def predict_factor(self, theta):
        """Predict, returning the emulator covariance as a FACTOR rather than a diagonal.

        Returns (mean, L) with mean (n_x_used,) and L (n_x_used, p), such that the emulator's
        predictive covariance is exactly `L @ L.T`. Single theta only.

        WHY THIS EXISTS. `predict` returns `var = (basis**2) @ coeff_var`, which is the DIAGONAL
        of `basis @ diag(coeff_var) @ basis.T` -- a rank-p (p = 6 here) matrix over ~695 rows.
        Using that diagonal in a likelihood asserts ~695 INDEPENDENT model errors where the
        emulator only ever claims p uncertain coefficients, and it is ANTI-conservative in the
        directions that matter: measured on campaign 7, the diagonal understates the variance
        along the leading basis direction -- a coherent whole-profile warm/cold shift, i.e. the
        direction dT_scale moves -- by 233x (true 201.6, diagonal 0.860). Since sigma_emu is
        86-98% of the diagonal error budget, this is the dominant term, not a refinement.
        Compounding it, only 498 of campaign 7's 695 rows carry DISTINCT predictions (one
        entity's 41 observations collapse to 17), so the diagonal asserts independence between
        rows whose predictive correlation is exactly 1.

        The trace is preserved by construction: sum(diag(L @ L.T)) == sum(var from predict()).
        """
        theta = np.atleast_2d(np.asarray(theta, dtype=float))
        if theta.shape[0] != 1:
            raise ValueError(f'predict_factor takes ONE theta, got {theta.shape[0]}')
        coeff_mean, coeff_var = self.backend.predict(self._slice(theta))
        mean = self.mean_ + (coeff_mean @ self.basis.T)[0]

        curve = getattr(self, '_variance_inflation_curve', None)
        if curve is not None:
            inflation = float(curve.predict(self._nn_distance(theta))[0])
        else:
            inflation = float(self.variance_inflation_)
        # basis is (n_x_used, p); scaling its columns by sqrt(coeff_var * inflation) gives a
        # factor whose outer product is the inflated rank-p covariance.
        L = self.basis * np.sqrt(np.maximum(coeff_var[0] * inflation, 0.0))[None, :]
        return mean, L

    def truncation_variance(self):
        """Per-row variance DISCARDED by the SVD truncation, as a scalar.

        `fit_svd` keeps components up to `explained_variance` (0.99), so ~1% of the training
        signal is thrown away. That part is genuinely near-independent across rows -- it is what
        the retained basis cannot represent -- and so is the one non-arbitrary DIAGONAL term the
        emulator contributes. Without it, replacing the diagonal with the rank-p factor would
        leave the likelihood claiming exactly zero uncertainty in the 689 directions orthogonal
        to the basis, which is worse than the diagonal it replaces.

        Estimated from the singular-value spectrum and spread uniformly over rows; this is an
        approximation (the true truncation error is not uniform), deliberately on the small side.
        """
        S = np.asarray(self.singular_values_, dtype=float)
        p = self.basis.shape[1]
        if p >= len(S):
            return 0.0
        n_rows = self.basis.shape[0]
        n_train = max(getattr(self, 'theta_train_', np.zeros((1, 1))).shape[0], 1)
        return float((S[p:] ** 2).sum() / (n_train * n_rows))

    def x_index_used(self):
        return [xi for xi, keep in zip(self.x_index, self._row_mask) if keep]

    def explained_variance_ratio(self):
        return np.cumsum(self.singular_values_ ** 2) / np.sum(self.singular_values_ ** 2)

    def calibrate_variance(self, design_theta, training_output, n_folds=None, seed=0,
                            trust_multiplier=5.0, trust_spacing_multiplier=2.0):
        """Leave-one-out cross-validate THIS emulator's own predictive variance against real
        held-out design points, and fit a LOCATION-DEPENDENT variance-inflation curve
        (_variance_inflation_curve) correcting for measured overconfidence -- predict() applies
        it automatically from then on, so every downstream consumer (BayesianCalibrator's
        likelihood, MCMC, LOO validation) is corrected transparently, with no changes needed on
        their end. Also sets max_trusted_distance_ (see below) -- BayesianCalibrator's
        _within_design_coverage gate uses it to hard-exclude candidate theta too far from real
        training data, rather than relying on the inflation curve alone.

        Confirmed empirically (full 100-fold LOO-CV against real CentralEurope production
        data): standardized residuals z=(true-pred)/pred_std had std(z)~2.1-2.4 overall (1.0 =
        perfectly calibrated -- some run-to-run variation from PCGPwM's own fitting
        stochasticity), with severe location-dependence -- std(z)~1.2-1.4 near the design's
        center/middle, ~3.8 near its edges. The emulator DOES correctly sense it's less
        reliable near the edges (predicted std grows ~2.3x center-to-edge), but true error
        grows even faster (~2.9x), so the stated uncertainty doesn't keep pace -- exactly the
        gap a KO calibration likelihood (which divides the residual penalty by predicted
        variance) can exploit by treating a poor, edge-of-design fit as "forgivably uncertain"
        rather than genuinely untrustworthy. A single GLOBAL inflation factor was tried first
        and confirmed insufficient at the worst corners (dT_scale still pinned at its design
        ceiling after applying it) -- the true relationship needs to scale with how far into
        sparse territory a given prediction actually is.

        This fits a MONOTONIC (isotonic) regression of E[z^2] on a standardized theta-space
        NEAREST-NEIGHBOR distance (see _nn_distance) -- NOT the emulator's own raw predicted
        variance. Raw variance was tried first and rejected: pooled directly across all output
        rows/glaciers, it produced an almost completely flat curve, because different rows have
        wildly different NATURAL variance scales (a stable deep-ice temperature vs. a variable
        shallow firn one) -- "raw_var=0.5" doesn't mean the same thing for every row, breaking
        the assumption that its magnitude alone is a comparable cross-row proxy for
        extrapolation risk. Theta-space distance has no such problem: it is a property of the
        INPUT location alone, identical for every row/glacier a given theta predicts, and
        directly reproduces the pattern independently confirmed via a coarser theta-distance
        tercile split (std(z)~1.2 near the design's center, ~3.8 near its edges). Isotonic (not
        a fixed parametric curve) because the only assumption physically justified here is that
        overconfidence gets no BETTER as distance grows -- y_min=1.0 additionally forbids the
        curve from ever DEFLATING variance below the GP's own raw estimate.

        WHY a hard exclusion distance exists at all, beyond the inflation curve above: confirmed
        empirically (against a real corner the calibration kept converging to) that no amount of
        HONEST variance inflation can prevent the KO likelihood from preferring a poorly-fitting,
        far-from-training-data theta over a well-fitting, well-covered one -- log-posterior at
        that corner peaked around a 300-1000x inflation multiplier, far beyond what leave-one-out
        cross-validation itself supports there (~25x). Pushing the curve past what the data
        justifies just to defeat one specific outcome would not be a calibration, it would be
        reverse-engineering a number. The correct response to "we have essentially no real
        information this far from training data" is to refuse to draw any conclusion there at
        all, not to assign it a large-but-still-finite variance and let the likelihood decide.
        max_trusted_distance_ is set to the smallest theta-space distance at which the fitted
        curve's required inflation first exceeds `trust_multiplier` times its own baseline value
        (the inflation right at distance 0, which by construction -- the curve is fit
        increasing=True -- is its global minimum) -- i.e. where the curve's own PAVA-pooled
        breakpoints show a genuine, data-determined jump in required inflation, not an
        arbitrarily chosen cutoff.

        FALLBACK WHEN NO JUMP IS FOUND: if the curve never exceeds `trust_multiplier * baseline`
        anywhere in the tested range (`beyond` empty), an earlier version of this method trusted
        theta all the way out to `dist.max()` -- the single farthest distance any LOO fold
        happened to test. Confirmed to fail in practice: the theta-expansion campaign (145-point
        design, CentralEurope 25-glacier) produced a completely FLAT curve (7.36x at every
        tested distance from 0 to 3.0 -- LOO-CV found no detectable distance/error relationship
        at all, plausibly because that design's added points clustered tightly in one region
        rather than spreading distance-diversity evenly across the LOO folds), so `dist.max()`
        silently set max_trusted_distance_ to 1.549 -- not because the data supported trusting
        that far, but because no jump was found to stop it sooner. MCMC then found and settled
        on a theta at standardized distance 1.169 from the nearest real training point (zero
        points anywhere near perm_frac>0.8 with dT_scale in [3.5,4.3]), with R-hat=1.0 and huge
        ESS (the sampler converged cleanly onto this unsupported corner), and LOO RMSE got
        WORSE than the prior campaign (2.649 -> 4.561 degC) despite the emulator's own
        diagnostics (LOO z-std, SVD components) all looking better -- i.e. exactly the
        exploit-a-technically-admitted-but-unsupported-corner failure this hard gate exists to
        prevent, just relocated by the flat-curve fallback rather than blocked by it. A single
        lucky-far LOO fold (or, as here, a total absence of signal) should not single-handedly
        set the whole trust radius. Fixed to the median of tested distances when no jump is
        found -- "trust extends as far as most (not all) of our tested evidence supports," a
        deliberately conservative middle value rather than the single most extreme point
        sampled. When a real jump IS found (the normal case, e.g. the 25-glacier campaign's
        clean 2.2x-then-172x step at ~0.85), that data-determined breakpoint is used exactly as
        before -- this fallback change only bites when the curve gives no signal to work with.

        Expensive (refits the full emulator once per held-out point, ~6-7s each for this
        project's design) -- call ONCE per campaign, after fit(), not inside any per-theta loop.
        n_folds subsamples for a faster (less precise) estimate; None = full leave-one-out.
        """
        from sklearn.isotonic import IsotonicRegression

        design_theta = np.asarray(design_theta)
        n = len(design_theta)
        fold_idx = (range(n) if n_folds is None
                    else np.random.default_rng(seed).choice(n, size=n_folds, replace=False))
        design_active = self._slice(design_theta)
        theta_std = design_active.std(axis=0)
        theta_std[theta_std == 0] = 1.0

        dist_chunks, z2_chunks = [], []
        for i in fold_idx:
            train_idx = [j for j in range(n) if j != i]
            # active_params MUST be propagated: a fold emulator built with the default
            # (0,1,2) would hit the degenerate-constant-column LinAlgError this class exists
            # to avoid whenever a parameter is held fixed.
            fold_emu = Emulator(glaciers=self.glaciers, explained_variance=self.explained_variance,
                                 active_params=self.active_params)
            fold_emu.fit(design_theta[train_idx], [training_output[j] for j in train_idx])

            true_col = fold_emu.assemble_matrix([training_output[i]])[:, 0]
            true_used = true_col[fold_emu._row_mask]

            mean_pred, var_pred = fold_emu.predict(design_theta[i])  # variance_inflation_/curve
            mean_pred, var_pred = mean_pred[0], var_pred[0]           # is a no-op on a fresh fold_emu

            # nearest-neighbor standardized distance from this held-out theta to the REMAINING
            # (fold's own) training design -- one scalar per fold, shared by every row it
            # contributes below (see _nn_distance / this method's docstring for why distance,
            # not raw variance, is the covariate).
            diffs = (design_active[train_idx] - design_active[i]) / theta_std
            nn_dist = np.sqrt(np.sum(diffs ** 2, axis=1)).min()

            # floor var_pred before dividing -- a near-zero (but not exactly zero) predicted
            # variance blows z^2 up to nonsense (seen directly: a single ~1e-20 value turned
            # loo_z_std_ into ~1e13 and corrupted an earlier version of this fit), the same
            # failure mode the flat-scalar version's std_pred = sqrt(max(var_pred, 1e-12))
            # guarded against.
            valid = ~np.isnan(true_used)
            var_floored = np.maximum(var_pred[valid], 1e-6)
            dist_chunks.append(np.full(int(valid.sum()), nn_dist))
            z2_chunks.append(((true_used[valid] - mean_pred[valid]) ** 2) / var_floored)

        dist = np.concatenate(dist_chunks)
        z2 = np.concatenate(z2_chunks)

        self.loo_z_std_ = float(np.sqrt(np.mean(z2)))       # kept for reporting/backward compat
        self.variance_inflation_ = self.loo_z_std_ ** 2      # flat fallback, see predict()

        curve = IsotonicRegression(y_min=1.0, out_of_bounds='clip', increasing=True)
        curve.fit(dist, z2)
        # DEGENERACY CHECK. A curve whose thresholds are all equal is a CONSTANT: it returns the
        # same inflation at every theta-space distance, so the location-dependent extrapolation
        # guard documented above is silently inert and only the hard max_trusted_distance_ gate is
        # doing any work. Campaign 7 shipped exactly this -- predict() returned
        # 8.718175525142334 at every distance from 0.0 to 0.29 -- and nothing anywhere reported
        # it. Warn loudly rather than raise: a degenerate curve is a usable (if blunt) fallback,
        # but it must never pass unnoticed again.
        if np.ptp(np.asarray(curve.y_thresholds_, dtype=float)) <= 0:
            import warnings
            warnings.warn(
                f'variance-inflation curve is DEGENERATE (constant '
                f'{float(curve.y_thresholds_[0]):.4f} at every theta distance): the '
                f'location-dependent extrapolation guard is inert and only max_trusted_distance_ '
                f'is active. See calibrate_variance\'s FALLBACK docstring.',
                RuntimeWarning, stacklevel=2)
            self.inflation_curve_degenerate_ = True
        else:
            self.inflation_curve_degenerate_ = False
        self._variance_inflation_curve = curve

        # ── trust radius: set from the DESIGN'S OWN GEOMETRY, not from the isotonic curve ──
        # The curve's breakpoints have no principled relationship to how far apart the design
        # points actually are, and deriving the trust radius from them failed in BOTH directions
        # on real campaigns:
        #   campaign 3 (145 pts, 3D): breakpoint gave 1.549 = 6.1x the design's median
        #     point-to-point spacing. MCMC converged (R-hat 1.00) onto a corner 1.169 away with
        #     no real support, and LOO RMSE nearly doubled. TOO LOOSE.
        #   campaign 5 (150 pts, 2D): breakpoint gave 0.157 = 1.00x median spacing -- i.e. INSIDE
        #     the design's own typical spacing, so 35.3% of the design box was refused and the
        #     accepted region shattered into NINE disconnected islands. emcee cannot cross an
        #     -inf barrier, so walkers were trapped in separate islands and reported a spurious
        #     "bimodal" posterior with R-hat 1.59. TOO TIGHT.
        #   campaign 4 (245 pts, 3D) also sat at 1.00x median spacing and was equally shattered;
        #     it only appeared healthy (R-hat 1.00) because its optimum happened to fall just
        #     INSIDE an island rather than just outside. That was luck, not correctness.
        #
        # The gate's job is to refuse genuine EXTRAPOLATION, not to refuse interpolation between
        # neighbouring design points -- interpolating there is precisely what a GP is for. So the
        # radius is now scaled to the design's own natural length scale: the median distance from
        # each design point to its nearest neighbour. A query may be up to
        # `trust_spacing_multiplier` times that far from the design before being refused.
        #
        # 2.0 validated against all three campaigns above (grid scan of the full design box):
        #   campaign 5: 2.1% of box refused, accepted region is a SINGLE connected component
        #               holding 100% of the accepted area (vs 9 islands before)
        #   campaign 3: threshold 0.506 vs the bad corner at 1.169 -> still REJECTED, 2.3x margin
        #   campaign 4: threshold 0.330 vs its posterior at 0.157 -> accepted, comfortably
        # 1.5x also connects the space but leaves 10.3% refused; 1.0x still shatters (9 islands).
        Z = self.theta_train_ / np.where(self.theta_train_.std(axis=0) == 0, 1.0,
                                          self.theta_train_.std(axis=0))
        pair = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
        np.fill_diagonal(pair, np.inf)
        self.design_spacing_ = float(np.median(pair.min(axis=1)))
        self.max_trusted_distance_ = float(trust_spacing_multiplier * self.design_spacing_)

        # kept purely as a diagnostic -- what the old curve-derived rule WOULD have said, so a
        # future campaign can still see whether the curve detected anything unusual.
        baseline = float(curve.predict([0.0])[0])
        beyond = curve.X_thresholds_[curve.y_thresholds_ > trust_multiplier * baseline]
        self.curve_trust_distance_ = float(beyond.min()) if len(beyond) else None
        return curve
