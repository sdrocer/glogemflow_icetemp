"""
BayesianCalibrator: Kennedy-O'Hagan Bayesian calibration, y = G(x;theta) + delta(x) + epsilon.

  - G(x;theta): the Emulator (SVD + surmise PCGPwM), standing in for the expensive real IDL
    GloGEM model. Its per-point predictive variance is Sigma_emu(theta).
  - delta(x): a spatial GP over ONE weighted-mean temperature residual per calibration
    glacier (discrepancy.Discrepancy, reused from the already-validated notebook 05 design).
    Its kernel hyperparameters (length-scale, signal variance) are fit ONCE, via sklearn's own
    marginal-likelihood optimizer, on the POOLED GLOBAL residual field at a representative
    theta -- decoupled from the theta-sampling loop below, exactly as calibration_scheme_
    prompt.md specifies ("to avoid repeating the overfitting failure" seen when hyperparameters
    were fit on a ~6-glacier regional subset). delta is analytically marginalised (a zero-mean
    GP prior integrates out to exactly the multivariate-normal marginal used in log_likelihood
    below) rather than sampled -- MCMC only explores theta.
  - epsilon: glenglat's own reported temperature_uncertainty (DataHandler's `sigma`, already
    inflated for equilibrium=='estimated' profiles).

Likelihood, evaluated at each theta proposal (n = number of calibration glaciers used):
    r(theta) ~ N(0, K_delta + diag(sigma_emu(theta)^2 + sigma_obs^2))
where r_i(theta) is glacier i's weighted-mean (G_emu(theta) - y_obs) residual and K_delta is
the FIXED spatial covariance matrix (pure Matern*Constant signal, no nugget -- the nugget's
role is superseded here by the explicit, theta-dependent sigma_emu^2 + sigma_obs^2 term, so it
is deliberately excluded from K_delta to avoid double-counting noise).

Sampling: emcee (gradient-free affine-invariant ensemble MCMC) -- the emulator backend
(PCGPwM or the sklearn fallback) is not differentiable, so a gradient-based sampler (NUTS)
would require re-implementing the emulator in an autodiff framework; emcee has no such
requirement (locked decision, see the calibration plan's sampler comparison).
"""

from dataclasses import dataclass, field

import numpy as np

from .emulator import BASAL

# 'basal' rows (see emulator.BASAL) compare the model's deepest RESOLVED grid node against a
# glacier's deepest AVAILABLE glenglat observation -- these are rarely at the same physical
# depth (e.g. Grenzgletscher: 41 m modeled vs 348 m observed; see the calibration write-up on
# GloGEM's per-band ice-thickness ceiling). This is a genuine near-bed thermal-regime test
# (geothermal flux + insulation from ice thickness), not a same-depth comparison, so it needs a
# much looser observation-noise floor than an actual borehole reading at that exact depth.
BASAL_SIGMA_INFLATION = 5.0

# ── model-error term S_model ────────────────────────────────────────────────────────────
# The single most consequential defect found in the 2026-08-20 adversarial review: fit_discrepancy
# DISCARDS the discrepancy GP's own fitted WhiteKernel (4.77 degC^2) on the stated ground that
# sigma_obs^2 + sigma_emu^2 supersedes it. It does not. glenglat's sigma_obs is THERMISTOR
# PRECISION (median 0.10 degC^2), while the real within-entity residual scatter is ~2.05 degC (sd)
# with a THETA-INDEPENDENT floor of 2.25-3.81 degC^2 measured on the real 150-point IDL training
# output -- no reachable theta can absorb it, because it is model discrepancy, not noise.
# Consequence, Mahalanobis/n (should be ~1 for an honestly calibrated likelihood):
#     62.5 at theta0 | 13.0 at campaign 7's own posterior mode | 76.5 at (0.176, 1.380)
# i.e. the likelihood was over-sharp by one to two orders of magnitude, which is why campaign 7's
# credible intervals came out 6-28x too narrow, why dT_scale pinned at its 5.0 bound with sd
# 4.1e-4, and why the high-recall region of theta-space was excluded by THOUSANDS of log-units.
# Restoring an explicit model-error term returns Mahalanobis/n to 0.95-1.18.
#
# AMPLITUDE. Three independent estimates converge: the discarded GP nugget (4.77), a
# profile-likelihood ML fit at fixed L=15 (3.6, with a 2-log-unit interval [3.3, 4.2]), and the
# theta-independent within-entity floor on real IDL output (2.25-3.81). 3.5 is the headline value;
# vary it via the dataclass field for the sensitivity arms {2.3, 3.5, 5.0}.
#
# FIXED, NOT FITTED. Measurement says fitting would be safe (s2 is identified to +-25% from 4343
# within-entity pairs, and the theta MAP is unchanged over s2 in [1.5, 3.5]), but fixing keeps the
# MCMC 2-D and removes any argument about a free variance parameter absorbing signal. Judgement
# call on a measured-to-be-benign question.
S_MODEL_VARIANCE = 3.5      # degC^2
#
# CORRELATION LENGTH, in DEPTH LAG. Exponential (Matern nu=1/2). 15 m, justified by the MEASURED
# within-entity residual variogram range (14.0 m at interior theta, 14.1 at the prior mean, 17.2
# at campaign 7's posterior) -- NOT by thermal_structure.SEASONAL_DEPTH_M, which is a depth
# THRESHOLD (`deep = depths > seasonal_depth_m`) rather than a lag scale; the two share units and
# nothing else. Exponential rather than a smoother kernel because the measured variogram is
# strongly non-smooth at the origin AND because smooth kernels drive the contrast variance of
# near-coincident depths to zero quadratically -- with a minimum within-entity depth gap of 0.13 m
# (data._pool_by_depth_bin stores each bin's MEAN depth, not its centre), that let a single pair
# carry 62-74% of the total quadratic form.
S_MODEL_LENGTH_M = 15.0


def _weighted_mean_and_var(values, weights, per_value_var):
    """Weighted mean and its variance, for a weighted average of independent quantities with
    per-value variance `per_value_var` (used to pool per-depth residuals/obs-error/emulator-
    variance onto one value per glacier, using the same depth weights as DataHandler)."""
    w = np.asarray(weights, dtype=float)
    v = np.asarray(values, dtype=float)
    pv = np.asarray(per_value_var, dtype=float)
    wsum = w.sum()
    mean = np.sum(w * v) / wsum
    var = np.sum((w ** 2) * pv) / (wsum ** 2)
    return mean, var


def _gaussian_logpdf(residual, Sigma):
    """log N(residual | 0, Sigma) via a Cholesky factorisation.

    Uses scipy's cho_factor/cho_solve rather than np.linalg.cholesky followed by two
    np.linalg.solve calls. The latter was the original implementation and is correct but
    wasteful: np.linalg.solve runs a general LU factorisation (O(n^3)) on each of the two
    TRIANGULAR factors, throwing away the triangularity the Cholesky just established.
    cho_solve dispatches to LAPACK's potrs, which is two O(n^2) triangular back-substitutions.
    Irrelevant at the mean-based likelihood's n=107; it is not irrelevant at the depth-resolved
    likelihood's n=695, and this runs on every MCMC evaluation (32 walkers x 30000 steps).

    Returns -inf on a non-positive-definite Sigma, so a bad proposal is rejected by the sampler
    rather than raising.
    """
    from scipy.linalg import cho_factor, cho_solve

    # cho_factor raises ValueError -- NOT LinAlgError -- on a NaN/inf entry, so the except
    # clause below does not catch it and a single bad Sigma would abort a 30000-step chain
    # instead of the sampler simply rejecting that proposal. A NaN is reachable in practice:
    # sigma_emu2 comes from the emulator, whose predictive variance can degenerate. Guard
    # explicitly rather than passing check_finite=False, which would let NaNs through into the
    # factorisation and return a silently meaningless log-likelihood.
    if not np.isfinite(Sigma).all():
        return -np.inf
    try:
        c, low = cho_factor(Sigma, lower=True)
    except np.linalg.LinAlgError:
        return -np.inf
    alpha = cho_solve((c, low), residual)
    quad = float(residual @ alpha)
    # only the lower triangle of `c` is meaningful (the upper holds untouched input values),
    # but the DIAGONAL is always the Cholesky factor's own -- safe to read directly.
    logdet = 2.0 * float(np.sum(np.log(np.diag(c))))
    n = len(residual)
    return -0.5 * (quad + logdet + n * np.log(2 * np.pi))

@dataclass
class BayesianCalibrator:
    emulator: object       # emulator.Emulator, already fit
    glaciers: list          # calibration glaciers, must match emulator.glaciers (order-free)
    priors: object           # priors.Priors
    track: str = 'all'      # 'all' | 'depth' (Track-1 only) | 'basal' (Track-2 only)
    s_model_variance: float = S_MODEL_VARIANCE   # degC^2; 0.0 restores pre-2026-08-20 behaviour
    s_model_length_m: float = S_MODEL_LENGTH_M   # depth-lag correlation length [m]
    use_elevation: bool = True         # pass elevation as a 3rd discrepancy coordinate
    use_emulator_factor: bool = True   # carry the emulator's rank-p covariance FACTOR instead of
    # its diagonal (see Emulator.predict_factor). False restores pre-2026-08-20 behaviour.
    fixed_params: dict = None   # e.g. {'z0': 15.0}; None/{} = calibrate all three. Parameters
    # named here are HELD FIXED and excluded from the search space; the rest are calibrated.
    #
    # SUPERSEDES the earlier `fixed_perm_frac` field, which fixed perm_frac at 1.0 on the
    # grounds that perm_frac and dT_scale "enter ONLY as a product (ins = perm_frac * ICE_FRAC
    # * dT_scale * dT_firn_band) and so are structurally non-identifiable". That justification
    # was WRONG, and was disproven on 2026-08-18: it describes the PYTHON ANALYTICAL SURROGATE
    # (physics.cp_model_single:70-74), not the real IDL forward model the emulator is trained
    # on. In the real model, perm_frac appears at exactly ONE line --
    # firnice_temperature_model.pro:89, `z_perm_b = firnice_perm_depth * firnice_perm_frac_b`,
    # scaling the Herron-Langway meltwater percolation depth -- and does NOT appear in the
    # surface boundary condition at all (lines 163-168 use only dT_scale and ICE_FRAC). So
    # perm_frac and dT_scale are INDEPENDENT physical mechanisms there (percolation depth vs
    # surface insulation amplitude), not a confounded product. Verified empirically on the real
    # 245-run training set: spearman(perm_frac, output) = +0.618 (p=3e-27), i.e. strongly
    # influential and separately identifiable. Fixing it discarded a real parameter.
    #
    # Conversely `z0` is the one that genuinely cannot be calibrated as the model currently
    # stands: it appears NOWHERE in firnice_temperature_model.pro, and although
    # initialise_firnicetemp_spinup.pro:189 reads firnice_z0_firn_b to build the initial C&P
    # profile, glogem.pro applies the per-glacier/design-point override AFTER that (line 420
    # builds the profile; lines 422/425/428 then write the array, with nothing recomputing the
    # profile). So a design point's z0 can never reach the physics. Verified empirically:
    # spearman(z0, output) = -0.031, p=0.63 -- indistinguishable from zero, against a validated
    # control (dT_scale 0.203 vs 4.989 changes output by 11.35 degC, and all 245 output vectors
    # are distinct). This is why z0's posterior never narrowed in any of the four campaigns.
    #
    # See emulator.BASAL. theta should still be FIT from track='depth' only: fitting against a
    # POOLED depth+basal residual let a handful of basal rows, comparing the model's own
    # thickness-capped "fake bed" against a much deeper real borehole, drag the fit to a prior
    # bound (confirmed empirically: even heavily downweighting basal rows via
    # BASAL_SIGMA_INFLATION just flips it to the OPPOSITE bound rather than settling at an
    # interior value -- these rows were never well-posed to inform profile curvature at all).
    # track='basal' is for the non-gating basal DIAGNOSTIC only (see pipeline.py/validation.py):
    # evaluate a Track-1-fitted theta's own residual against basal rows, never fit theta from it.

    def __post_init__(self):
        if self.track not in ('all', 'depth', 'basal'):
            raise ValueError(f"track must be 'all', 'depth', or 'basal', got {self.track!r}")
        # Per-dimension bounds of the emulator's ACTUAL training design (not the nominal prior
        # range) -- see log_posterior's design-coverage gate. dT_scale's prior is a truncated
        # Gaussian(mean=1,std=1) on [0.2,5.0]; a 100-point LHS draw from that naturally puts
        # almost no samples near the 4-sigma tail, so the real design covers only ~[0.22,4.01]
        # despite the nominal upper bound of 5.0. Confirmed the emulator's predictive variance
        # inflates ~100-200x approaching that unsampled edge (from ~0.02-0.13 near the design's
        # center to ~9+ at dT_scale=5), and the KO likelihood -- dividing the residual penalty
        # by that variance -- rewards this as "forgivably uncertain" fit rather than correctly
        # distrusting an unsampled, extrapolated region; MCMC/MAP exploit it by walking straight
        # for the edge regardless of whether the underlying physics is actually ambiguous there.
        self._design_lo = self.emulator.theta_train_.min(axis=0)
        self._design_hi = self.emulator.theta_train_.max(axis=0)
        # Keyed by glacier_name (always unique -- DataHandler groups measurements by name),
        # NOT glacier_id: multiple GlacierCalibrationData objects can share one GloGEM
        # glacier_id (two glenglat boreholes on the same RGI polygon). glacier_id is only
        # needed downstream, for the IDL-facing residual file (see writeback.py).
        self._glacier_by_name = {g.glacier_name: g for g in self.glaciers}
        self._x_used = self.emulator.x_index_used()
        # group x_index_used rows by glacier_name, preserving each row's position in the
        # emulator's output vector (needed to slice mean/var predictions per glacier). Track
        # filtering happens HERE, by skipping rows rather than pre-slicing self._x_used --
        # `row` must stay an index into the emulator's own full, unfiltered prediction vector
        # (self.emulator.predict()'s output order), which compute_glacier_residuals indexes
        # into directly via `rows` below; renumbering a pre-filtered list would silently
        # misalign predictions with the wrong rows.
        self._rows_by_glacier = {}
        for row, (gname, depth) in enumerate(self._x_used):
            is_basal_row = (depth == BASAL)
            if self.track == 'depth' and is_basal_row:
                continue
            if self.track == 'basal' and not is_basal_row:
                continue
            self._rows_by_glacier.setdefault(gname, []).append(row)

        self._calib_glacier_names = [
            gname for gname in self._rows_by_glacier if gname in self._glacier_by_name
        ]
        self._lat = np.array([self._glacier_by_name[n].latitude for n in self._calib_glacier_names])
        self._lon = np.array([self._glacier_by_name[n].longitude for n in self._calib_glacier_names])
        # ELEVATION as a third discrepancy coordinate. discrepancy.py has supported it since the
        # elevation split (anisotropic/ARD Matern, elevation in metres against a km sphere
        # embedding), and campaign 7's config calls it "Required, not cosmetic" -- but it was
        # never wired in here, so every campaign through 7 fit K_delta on lat/lon ALONE. After
        # the split up to 24 entities of one glacier share essentially identical lat/lon (median
        # inter-entity separation 0.171 km), so on lat/lon alone the fitted length scale pins at
        # its 10 km LOWER BOUND and co-located entities correlate ~0.991.
        # NOT A KNOWN IMPROVEMENT -- treat as an experiment with a check. Measured on campaign 7,
        # adding elevation raises the GP marginal likelihood by +10.7 but lands an elevation ARD
        # length scale of 1530 m (larger than most within-glacier spans), barely moves the fitted
        # nugget (4.77 -> 4.25), and DROPS the effective rank from 8.5 to 3.8. If it makes things
        # worse once S_model is in place, the honest reading is that the elevation split needs an
        # along-flow coordinate this data cannot supply -- a finding, not a failure.
        self._elev = np.array([self._glacier_by_name[n].elevation for n in self._calib_glacier_names])

        # per-glacier observation index (which g.depths/g.T_obs entry each x_used row maps
        # to), depth weights, and observation variance, aligned to _x_used row order. Cached
        # once here rather than recomputed inside compute_glacier_residuals -- that method
        # runs on every MCMC likelihood evaluation (tens of thousands of times per run), and
        # the mapping itself doesn't depend on theta.
        self._obs_idx = {}
        self._weights = {}
        self._sigma_obs2 = {}
        for gname in self._calib_glacier_names:
            g = self._glacier_by_name[gname]
            rows = self._rows_by_glacier[gname]
            depths_used = [self._x_used[r][1] for r in rows]
            is_basal = np.array([d == BASAL for d in depths_used])
            idx = np.array([
                np.argmax(g.depths) if b else np.argmin(np.abs(g.depths - d))
                for d, b in zip(depths_used, is_basal)
            ])
            self._obs_idx[gname] = idx
            self._weights[gname] = g.weights[idx]
            sigma2 = g.sigma[idx] ** 2
            sigma2[is_basal] *= BASAL_SIGMA_INFLATION ** 2
            self._sigma_obs2[gname] = sigma2

        self._build_row_structures()

        self.discrepancy = None
        self._K_delta = None  # cached pure-signal covariance at the calibration locations
        self._Sigma_base = None  # cached theta-INDEPENDENT part of Sigma (see fit_discrepancy)

    # ── row-level structures (flat across entities) ────────────────────────────────────
    def _build_row_structures(self):
        """Flatten every entity's used rows into one vector, and build the pooling matrix A.

        A is (n_entities, n_rows) with A[k] holding entity k's normalised depth weights, so
        `A @ row_vector` reproduces compute_glacier_residuals' weighted mean EXACTLY. This lets
        the row-level covariance terms (model error, the emulator's rank-p factor) be built once
        at row level and pooled to entity level as `A M A.T`, instead of being approximated by
        _weighted_mean_and_var -- which pools them as INDEPENDENT and so shrinks the emulator
        contribution by ~1/n_rows (measured correct/coded ratio: median 3.4, p95 11.7).

        Skipped for track='basal': emulator.BASAL is the STRING 'basal' sitting in the depth
        column, so |depth_a - depth_b| is undefined there. That track is a non-gating diagnostic
        (see the fixed_params docstring) and keeps the original diagonal treatment.
        """
        self._rows_ok = False
        if self.track == 'basal':
            return
        emu_idx, ent_idx, depth, sig2, obs = [], [], [], [], []
        for k, gname in enumerate(self._calib_glacier_names):
            g = self._glacier_by_name[gname]
            rows = self._rows_by_glacier[gname]
            idx = self._obs_idx[gname]
            depths_used = [self._x_used[r][1] for r in rows]
            if any(d == BASAL for d in depths_used):
                return   # mixed track='all' -- fall back to the original diagonal path
            emu_idx.extend(rows)
            ent_idx.extend([k] * len(rows))
            depth.extend(float(d) for d in depths_used)
            sig2.extend(self._sigma_obs2[gname])
            obs.extend(g.T_obs[idx])
        if not emu_idx:
            return

        self._row_emu_idx = np.asarray(emu_idx, dtype=int)
        self._row_entity = np.asarray(ent_idx, dtype=int)
        self._row_depth = np.asarray(depth, dtype=float)
        self._row_sigma_obs2 = np.asarray(sig2, dtype=float)
        self._row_obs = np.asarray(obs, dtype=float)

        n_ent, n_row = len(self._calib_glacier_names), len(self._row_emu_idx)
        A = np.zeros((n_ent, n_row))
        start = 0
        for k, gname in enumerate(self._calib_glacier_names):
            w = np.asarray(self._weights[gname], dtype=float)
            A[k, start:start + len(w)] = w / w.sum()
            start += len(w)
        self._A = A
        self._rows_ok = True

    def _row_model_covariance(self):
        """S_model: within-entity exponential covariance in DEPTH LAG, zero across entities."""
        if self.s_model_variance <= 0:
            return np.zeros((len(self._row_depth),) * 2)
        d = self._row_depth
        same = self._row_entity[:, None] == self._row_entity[None, :]
        return np.where(same, self.s_model_variance
                         * np.exp(-np.abs(d[:, None] - d[None, :]) / self.s_model_length_m), 0.0)

    # ── residuals at a given theta ──────────────────────────────────────────────────────
    def compute_glacier_residuals(self, theta):
        """Weighted-mean (obs - emulator prediction) residual per calibration glacier, plus
        the pooled observation and emulator predictive variance. Returns arrays ordered as
        self._calib_glacier_names."""
        mean_pred, var_pred = self.emulator.predict(theta)
        mean_pred, var_pred = mean_pred[0], var_pred[0]

        residual = np.zeros(len(self._calib_glacier_names))
        sigma_obs2 = np.zeros_like(residual)
        sigma_emu2 = np.zeros_like(residual)
        for k, gname in enumerate(self._calib_glacier_names):
            g = self._glacier_by_name[gname]
            rows = self._rows_by_glacier[gname]
            idx = self._obs_idx[gname]
            obs_vals = g.T_obs[idx]
            pred_vals = mean_pred[rows]
            pred_var = var_pred[rows]
            w = self._weights[gname]

            r_i, _ = _weighted_mean_and_var(obs_vals - pred_vals, w, np.zeros_like(w))
            residual[k] = r_i
            _, sigma_obs2[k] = _weighted_mean_and_var(obs_vals, w, self._sigma_obs2[gname])
            _, sigma_emu2[k] = _weighted_mean_and_var(pred_vals, w, pred_var)

        return residual, sigma_obs2, sigma_emu2

    # ── discrepancy hyperparameter pre-fit (once, outside the MCMC loop) ───────────────
    def fit_discrepancy(self, theta0=None):
        """Fit delta(x)'s spatial-covariance hyperparameters once, via sklearn's marginal-
        likelihood optimizer, on residuals at a representative theta (prior mean by default).
        Caches the pure-signal (Constant*Matern, no WhiteKernel) covariance matrix at the
        calibration-glacier locations for reuse in every log_likelihood() call."""
        from .discrepancy import Discrepancy

        if theta0 is None:
            theta0 = np.array([
                self.priors.perm_frac.mean(), self.priors.dT_scale.mean(), self.priors.z0.mean(),
            ])
        residual0, sigma_obs2_0, sigma_emu2_0 = self.compute_glacier_residuals(theta0)

        # 'T_residual' is a temperature-space discrepancy, not one of Discrepancy's default
        # per-parameter nugget floors -- use the median pooled observation-error variance
        # across calibration glaciers as the noise floor instead (same role: stop the GP from
        # over-interpolating below the level its own data can actually resolve).
        nugget_floor = {'T_residual': float(np.median(sigma_obs2_0))}
        self.discrepancy = Discrepancy(nugget_floor=nugget_floor)
        gp = self.discrepancy.fit('T_residual', self._lat, self._lon, residual0,
                                   elevations=self._elev if self.use_elevation else None)
        # pure spatial-signal kernel (drop the WhiteKernel term -- superseded by the explicit,
        # theta-dependent sigma_obs^2 + sigma_emu^2 in log_likelihood, avoiding double-counting)
        signal_kernel = gp.kernel_.k1 if hasattr(gp.kernel_, 'k1') else gp.kernel_
        from .discrepancy import latlon_to_xyz
        X = latlon_to_xyz(self._lat, self._lon,
                           self._elev if self.use_elevation else None)
        self._K_delta = signal_kernel(X)
        self._build_sigma_base()
        return self.discrepancy

    def _build_sigma_base(self):
        """Cache the theta-INDEPENDENT part of Sigma, at ENTITY level.

        Sigma(theta) = K_delta + A [ diag(sigma_obs^2 + sigma_trunc^2) + S_model ] A.T
                                + (A L_emu(theta)) (A L_emu(theta)).T
        Only the last term moves with theta, and it is a rank-p (p=6) update, so caching the rest
        turns each likelihood evaluation into one 107x6 matmul plus a Cholesky. Measured on the
        n=695 row-level version this is 20.2 -> 6.55 ms/eval, i.e. 5.4 h -> 1.75 h of MCMC.
        """
        if not getattr(self, '_rows_ok', False):
            self._Sigma_base = None
            return
        trunc = getattr(self.emulator, 'truncation_variance', lambda: 0.0)()
        D = np.diag(self._row_sigma_obs2 + trunc)
        M = D + self._row_model_covariance()
        self._Sigma_base = self._K_delta + self._A @ M @ self._A.T

    # ── KO log-posterior ─────────────────────────────────────────────────────────────
    def residual_and_covariance(self, theta):
        """(residual, Sigma) at ENTITY level -- the exact pair log_likelihood scores.

        Exposed so the calibration diagnostic (mahalanobis_per_dof) scores the SAME quantities the
        sampler does, rather than a reconstruction that could drift out of step with it.
        """
        if self._K_delta is None:
            raise RuntimeError('call fit_discrepancy() before log_likelihood()/run_mcmc()')
        if getattr(self, '_rows_ok', False) and self.use_emulator_factor \
                and self._Sigma_base is not None:
            mean, L = self.emulator.predict_factor(theta)
            residual = self._A @ (self._row_obs - mean[self._row_emu_idx])
            M = self._A @ L[self._row_emu_idx]          # (n_entities, p) -- the ONLY theta-
            Sigma = self._Sigma_base + M @ M.T           # dependent part, a rank-p update
            return residual, Sigma
        # legacy path: diagonal emulator variance, no model-error term (pre-2026-08-20).
        residual, sigma_obs2, sigma_emu2 = self.compute_glacier_residuals(theta)
        return residual, self._K_delta + np.diag(sigma_obs2 + sigma_emu2)

    def log_likelihood(self, theta):
        residual, Sigma = self.residual_and_covariance(theta)
        return _gaussian_logpdf(residual, Sigma)

    def mahalanobis_per_dof(self, theta):
        """r.T Sigma^-1 r / n -- the likelihood's own calibration check.

        Should be ~1 if the covariance model is honest: it says "the residuals are as big as the
        uncertainty budget claims". Much greater than 1 means the likelihood is OVER-SHARP, i.e.
        it thinks it knows the answer far better than it does, and will exclude good theta with
        false confidence. Measured on campaign 7 BEFORE the S_model fix: 62.5 at theta0, 13.0 at
        its own posterior mode, 76.5 at (0.176, 1.380). After: 0.95-1.18.

        adoption gates on this -- see adoption.MAHALANOBIS_BOUNDS.
        """
        from scipy.linalg import cho_factor, cho_solve
        residual, Sigma = self.residual_and_covariance(theta)
        if not np.isfinite(Sigma).all():
            return float('nan')
        try:
            c, low = cho_factor(Sigma, lower=True)
        except np.linalg.LinAlgError:
            return float('nan')
        return float(residual @ cho_solve((c, low), residual) / len(residual))

    def log_posterior(self, theta):
        lp = self.priors.logpdf(theta)
        if not np.isfinite(lp):
            return -np.inf
        if not self._within_design_coverage(theta):
            return -np.inf
        ll = self.log_likelihood(theta)
        if not np.isfinite(ll):
            return -np.inf
        return lp + ll

    def _expand_theta(self, theta_free):
        """theta_free is the vector actually being searched/sampled (the FREE parameters, in
        _FULL_NAMES order). Returns the full (perm_frac, dT_scale, z0) 3-tuple every likelihood
        method (log_posterior, compute_glacier_residuals, ...) expects -- those methods are
        UNCHANGED by which parameters are fixed, only the search space is."""
        if not self.fixed_params:
            return np.asarray(theta_free, dtype=float)
        out, it = [], iter(np.asarray(theta_free, dtype=float))
        for name in self._FULL_NAMES:
            out.append(float(self.fixed_params[name]) if name in self.fixed_params else next(it))
        return np.array(out)

    def _free_param_names(self):
        fixed = self.fixed_params or {}
        return tuple(n for n in self._FULL_NAMES if n not in fixed)

    _FULL_NAMES = ('perm_frac', 'dT_scale', 'z0')

    def _within_design_coverage(self, theta):
        """Reject (-inf, via log_posterior) any theta outside the emulator's ACTUAL training
        design coverage -- see __post_init__'s _design_lo/_design_hi docstring. Only checked
        for FREE parameters (self._free_param_names()): a deliberately fixed parameter is a
        physical choice (see that field's docstring), not something this guard needs to police,
        and may legitimately sit fractionally outside the observed design range (e.g. 1.0 vs an
        observed max of ~0.999 from a finite LHS sample).

        Also rejects theta whose NEAREST-NEIGHBOR distance to the real training design exceeds
        emulator.max_trusted_distance_ (see Emulator.calibrate_variance) -- the per-dimension
        box check above admits "corners" that sit inside every dimension's individual range but
        are still far from any real training combination (e.g. dT_scale near its max together
        with z0 near its min, when the real design never actually explored that combination
        jointly). Confirmed empirically that no amount of honest variance inflation prevents the
        KO likelihood from preferring exactly such a corner over a well-covered, better-fitting
        point -- this is a hard exclusion, not a softer reweighting, because the right response
        to "no real information here" is to refuse a claim, not widen the uncertainty until the
        likelihood math is satisfied. getattr guards against an emulator that hasn't had
        calibrate_variance() run (no distance limit applied, matching prior behaviour)."""
        for name in self._free_param_names():
            i = self._FULL_NAMES.index(name)
            if not (self._design_lo[i] <= theta[i] <= self._design_hi[i]):
                return False
        max_dist = getattr(self.emulator, 'max_trusted_distance_', None)
        if max_dist is not None:
            nn_dist = self.emulator._nn_distance(np.asarray(theta).reshape(1, -1))[0]
            if nn_dist > max_dist:
                return False
        return True

    def _effective_bounds(self):
        """Prior bounds intersected with the emulator's actual design coverage, for each FREE
        parameter -- used to seed find_map's LHS starts and run_mcmc's walkers inside the
        region log_posterior will actually accept, rather than wasting starts/walkers in
        territory that's just going to return -inf."""
        prior_bounds = self.priors.bounds()
        eff = {}
        for name in self._free_param_names():
            i = self._FULL_NAMES.index(name)
            b_lo, b_hi = prior_bounds[name]
            eff[name] = (max(b_lo, self._design_lo[i]), min(b_hi, self._design_hi[i]))
        return eff

    # ── MAP point estimate ──────────────────────────────────────────────────────────
    def find_map(self, n_starts=24, seed=0):
        """Multi-start Nelder-Mead MAP search. A SINGLE start (even from the prior mean) is
        unreliable for this likelihood: an LHS sweep against real production data found local
        optima differing by tens of thousands of log-posterior units depending on the start,
        with perm_frac/dT_scale racing to their prior bounds and z0 landing essentially
        anywhere in [5,200] once they do. LHS-sample diverse starts across the EFFECTIVE bounds
        (prior bounds intersected with the emulator's actual training design coverage -- see
        _effective_bounds/_within_design_coverage; z0 sampled log-uniformly, matching its
        prior) and keep the best converged optimum -- used both to seed run_mcmc()'s walkers
        and by validation.Validator's per-fold MAP mode. Always returns a full (perm_frac,
        dT_scale, z0) 3-tuple, even when some parameters are fixed (see _expand_theta) -- every
        downstream consumer expects the full theta."""
        from scipy.optimize import minimize
        from scipy.stats import qmc

        eff_bounds = self._effective_bounds()
        free_names = self._free_param_names()
        ndim = len(free_names)
        lo = np.array([eff_bounds[n][0] if n != 'z0' else np.log(eff_bounds[n][0]) for n in free_names])
        hi = np.array([eff_bounds[n][1] if n != 'z0' else np.log(eff_bounds[n][1]) for n in free_names])
        # z0's starts are drawn LOG-uniformly to match its loguniform prior, so its column (and
        # only its column) is exponentiated back. Guarded because z0 is not necessarily a free
        # parameter any more: since 2026-08-18 the default is to hold z0 FIXED (it cannot reach
        # the real model's physics -- see fixed_params), in which case it is absent from
        # free_names entirely and an unguarded .index('z0') raises
        # "ValueError: tuple.index(x): x not in tuple". That is exactly what killed campaign 5's
        # step 4 after its 150 IDL runs had already completed.
        if 'z0' in free_names:
            z0_col = free_names.index('z0')
        else:
            z0_col = None
        u = qmc.LatinHypercube(d=ndim, seed=seed).random(n_starts)
        starts = lo + u * (hi - lo)
        if z0_col is not None:
            starts[:, z0_col] = np.exp(starts[:, z0_col])
        nm_bounds = [eff_bounds[n] for n in free_names]  # explicit optimizer-level bounds too,
        # not just the log_posterior -inf gate -- Nelder-Mead's simplex can otherwise wander
        # (and occasionally get stuck) once several vertices land on a -inf plateau.

        best_x, best_lp = None, -np.inf
        for theta0 in starts:
            res = minimize(lambda t: -self.log_posterior(self._expand_theta(t)), theta0,
                            method='Nelder-Mead', bounds=nm_bounds,
                            options={'maxiter': 3000, 'xatol': 1e-6, 'fatol': 1e-6})
            lp = -res.fun
            if res.success and lp > best_lp:
                best_x, best_lp = res.x, lp
        if best_x is None:
            # every start failed -- fall back to the prior mean of each FREE parameter
            prior_by_name = self.priors.as_dict()
            return self._expand_theta([prior_by_name[n].mean() for n in self._free_param_names()])
        return self._expand_theta(best_x)

    def grid_search(self, n_per_dim=41):
        """Exhaustive coarse grid over the FREE parameters, returning (theta, log_posterior).

        An independent check that the MCMC actually found the global mode. It is not a substitute
        for sampling -- it gives a point, not a distribution -- but it cannot get stuck, which is
        exactly the failure the sampler is prone to here.

        WHY THIS IS NEEDED. Campaign 7 reported R-hat 1.00 and ESS 15533 while its posterior mode
        sat 223 log-units BELOW a gate-accepted grid point. R-hat measures whether the walkers
        agree with EACH OTHER, not whether they found the best region -- and run_mcmc starts all
        of them in a 1e-2*prior_std ball around a single find_map result, so a shared bad start
        produces confident agreement on the wrong answer. Independent searches also disagreed
        about where the old likelihood's MAP was (three grids gave three different dT_scale
        locations), which is itself a symptom of the over-sharp covariance this class now fixes.

        Cost at the default: 41^2 = 1681 evaluations at ~1.2 ms = ~2 s. Negligible next to MCMC.
        """
        eff = self._effective_bounds()
        free = self._free_param_names()
        axes = [np.linspace(*eff[n], n_per_dim) for n in free]
        best_theta, best_lp = None, -np.inf
        for point in np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1).reshape(-1, len(free)):
            theta = self._expand_theta(point)
            lp = self.log_posterior(theta)
            if lp > best_lp:
                best_theta, best_lp = theta, lp
        return best_theta, best_lp

    # ── MCMC ────────────────────────────────────────────────────────────────────────
    def run_mcmc(self, n_walkers=32, n_steps=2000, burn_in=500, thin=5, seed=42,
                  progress=False):
        """Sample the KO posterior over theta = (perm_frac, dT_scale, z0) with emcee.
        Returns (sampler, flat_samples) -- flat_samples has burn-in removed and is thinned.

        Walkers start in a small ball around a MAP point estimate, NOT scattered across the
        full prior. Confirmed empirically this project's first production run (walkers drawn
        directly from Priors.rvs, spanning z0's full 40x log-uniform [5,200] range) had r_hat
        2.0-2.2 after 5000 steps, and running 6x longer (30000 steps) did not fix it (r_hat
        1.7-2.9, if anything worse for 2 of 3 parameters) -- consistent with poorly-conditioned
        affine-invariant stretch moves from an initial ensemble geometry unrelated to the
        target's local covariance, not simply an under-sampled chain. Small-ball
        initialization around a good starting point is emcee's own standard guidance.
        """
        import emcee

        if self._K_delta is None:
            self.fit_discrepancy()

        rng = np.random.default_rng(seed)

        map_theta = self.find_map(seed=seed)  # always a full 3-tuple, see _expand_theta
        free_names = self._free_param_names()
        ndim = len(free_names)
        prior_stds = {'perm_frac': self.priors.perm_frac.std(), 'dT_scale': self.priors.dT_scale.std(),
                      'z0': self.priors.z0.std()}
        prior_std = np.array([prior_stds[n] for n in free_names])
        full_names = ('perm_frac', 'dT_scale', 'z0')
        map_free = np.array([map_theta[full_names.index(n)] for n in free_names])

        # Perturb by a small fraction of each parameter's OWN prior std, not a fixed
        # absolute jitter -- z0 (std ~tens) and perm_frac (std ~0.1-0.3) are on completely
        # different scales, so one constant ball radius would be wildly wrong for one of them.
        p0 = map_free + 1e-2 * prior_std * rng.standard_normal((n_walkers, ndim))
        eff_bounds = self._effective_bounds()  # prior bounds ∩ emulator's actual design
        # coverage (see _within_design_coverage) -- walkers must start where log_posterior can
        # return a finite value, not just inside the nominal prior range.
        for i, name in enumerate(free_names):
            lo, hi = eff_bounds[name]
            p0[:, i] = np.clip(p0[:, i], lo, hi)

        def log_posterior_free(theta_free):
            return self.log_posterior(self._expand_theta(theta_free))

        # Box-bound clipping alone is not enough once the nearest-neighbor-distance gate
        # (_within_design_coverage) is active: that region can have an irregular shape whose
        # true boundary sits well inside the per-dimension box, so a walker can clip to a
        # "valid-looking" box position that is still outside real design coverage. Confirmed
        # directly against real production data: with the MAP landing essentially exactly ON
        # the trust boundary (a genuine constrained-optimum result -- see this method's
        # docstring -- not a bug), 19/32 walkers spawned at log_posterior=-inf under box-clip
        # alone, which corrupted emcee's own move-acceptance arithmetic (RuntimeWarning:
        # invalid value in scalar subtract, from -inf minus -inf) and produced unusable mixing
        # (tau ~1650 steps, R-hat 1.24 after 30000 steps). Shrink any invalid walker's offset
        # from map_free (guaranteed finite -- it's find_map's own result) by half, repeatedly,
        # until it lands somewhere valid.
        for w in range(n_walkers):
            shrink = 1.0
            while not np.isfinite(log_posterior_free(p0[w])) and shrink > 1e-6:
                shrink *= 0.5
                p0[w] = map_free + shrink * (p0[w] - map_free)
            if not np.isfinite(log_posterior_free(p0[w])):
                p0[w] = map_free.copy()  # last resort: the exact MAP

        sampler = emcee.EnsembleSampler(n_walkers, ndim, log_posterior_free)
        sampler.run_mcmc(p0, n_steps, progress=progress)

        flat_free = sampler.get_chain(discard=burn_in, thin=thin, flat=True)
        if not self.fixed_params:
            flat_samples = flat_free
        else:
            # re-expand to the full (perm_frac, dT_scale, z0) shape downstream consumers expect
            # (posterior_samples.csv, PARAM_NAMES-keyed reads, ResidualWriter, ...) -- each
            # fixed parameter's column is simply constant, never sampled.
            flat_samples = np.column_stack([
                (np.full(len(flat_free), float(self.fixed_params[n])) if n in self.fixed_params
                 else flat_free[:, self._free_param_names().index(n)])
                for n in self._FULL_NAMES
            ])
        return sampler, flat_samples

    def convergence_diagnostics(self, sampler, burn_in=500, thin=5):
        """R-hat (split-chain) and effective sample size via arviz, plus emcee's own
        integrated-autocorrelation-time estimate. Only over the parameters actually sampled
        (see _free_param_names) -- when a parameter is fixed, `sampler` itself is a lower-dim emcee
        sampler (run_mcmc's ndim), so there is no perm_frac column to report R-hat on."""
        import arviz as az

        chain = sampler.get_chain(discard=burn_in, thin=thin)  # (n_steps, n_walkers, ndim)
        idata = az.from_dict(posterior={
            name: chain[:, :, i].T  # arviz wants (chain, draw)
            for i, name in enumerate(self._free_param_names())
        })
        summary = az.summary(idata)
        try:
            tau = sampler.get_autocorr_time(discard=burn_in, quiet=True)
        except Exception:
            tau = None
        return summary, tau
