"""
Discrepancy: spatial Gaussian Process model-discrepancy term delta(x) for the KO calibration.

y = G(x; theta) + delta(x) + epsilon

delta(x) captures "missing physics" -- systematic spatial bias the 3 global parameters alone
cannot fix (regional accumulation/aspect/valley effects, etc). This is a faithful port of the
already-drafted, reviewed scheme from notebooks/05_firnicetemp_calibration.ipynb ("Gaussian
Process residual correction" cell), the direct fix for the diagnosed nearest-neighbour
failure (undamped correction copied up to 650+ km): a zero-mean GP with a distance-decaying
kernel automatically shrinks to "trust the global model" far from data and blends multiple
nearby calibrated glaciers instead of copying the single nearest one.

Trained on residuals at the GLOBAL calibration set (not just one region) -- fitting kernel
hyperparameters on ~6 regional points is exactly the failure mode that killed the CE-specific
regression (calibration_scheme_prompt.md, Tier-3a).
"""

from dataclasses import dataclass, field

import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel

EARTH_R_KM = 6371.0088

# Matern length-scale search range [km]: notebook 05's own choice, wide enough to let
# marginal-likelihood optimisation find anything from "very local" to "near-continental".
LENGTH_SCALE_INIT_KM = 80.0
# Lower bound dropped 10.0 -> 0.02 km on 2026-08-20. After the elevation split the MEDIAN
# inter-entity separation is 0.171 km, so a 10 km floor is two orders of magnitude coarser than
# the structure the split exists to expose: sklearn's optimiser pinned the fit AT the lower bound
# on campaign 7 (ConvergenceWarning, "close to the specified lower bound 10.0") and co-located
# entities came out correlated ~0.991 -- i.e. the GP was forced to treat one glacier's 24 bands
# as a single location no matter what the data said. 0.02 km (20 m) is below the closest real
# separation, so the bound no longer binds and the marginal likelihood is free to choose.
LENGTH_SCALE_BOUNDS_KM = (0.02, 500.0)


# Elevation length-scale search range [m], used only when latlon_to_xyz is given elevations.
# Wide enough to span "discrepancy changes within one glacier's elevation span" (~50 m) to
# "effectively elevation-independent" (~3000 m, i.e. larger than any single glacier's range).
ELEV_LENGTH_SCALE_INIT_M = 300.0
# Upper bound raised 3000 -> 12000 m on 2026-08-20: the horizontal dimensions were pinning at
# their own 500 km upper bound (ConvergenceWarning on dimensions 1 and 2), which is the ARD
# kernel's way of saying "this coordinate carries no information" -- a legitimate answer it must
# be able to express without hitting a wall, otherwise the pinned value is an artefact of the
# bound rather than a finding.
ELEV_LENGTH_SCALE_BOUNDS_M = (20.0, 12000.0)


def latlon_to_xyz(lat, lon, elevation=None):
    """Embed (lat, lon) on the unit sphere (scaled by Earth radius) so a standard isotropic
    kernel behaves like a great-circle-distance kernel.

    If `elevation` (m a.s.l.) is given, it is appended as a FOURTH coordinate and the caller
    must use an ANISOTROPIC kernel (see Discrepancy._fit_one) -- the first three columns are in
    KILOMETRES while elevation is in METRES, so an isotropic length scale of 10-500 km would
    make a 1000 m elevation difference effectively invisible. Giving elevation its own length
    scale is what makes the extra coordinate do anything at all.

    WHY ELEVATION IS NEEDED (2026-08-19): calibration entities are being split from one per
    glacier to one per (glacier, 10 m elevation band) -- 24 -> 125 entities for CentralEurope --
    so that deep boreholes on a glacier's tongue stop being matched against the model column at
    that glacier's summit (Grenzgletscher pooled boreholes from 2015 m to 4485 m into a single
    entity assigned 4450 m, so its 348 m tongue borehole was compared against a 35.9 m model
    column). After that split, up to 28 entities of the same glacier sit at essentially
    IDENTICAL lat/lon; on lat/lon alone the kernel would treat them as one location, giving a
    singular or near-singular covariance and collapsing exactly the along-flow structure the
    split exists to expose.
    """
    lat_r, lon_r = np.radians(lat), np.radians(lon)
    x = EARTH_R_KM * np.cos(lat_r) * np.cos(lon_r)
    y = EARTH_R_KM * np.cos(lat_r) * np.sin(lon_r)
    z = EARTH_R_KM * np.sin(lat_r)
    cols = [x, y, z]
    if elevation is not None:
        cols.append(np.asarray(elevation, dtype=float))
    return np.column_stack(cols)


@dataclass
class Discrepancy:
    """General-purpose spatial residual GP -- one fit per named scalar-per-glacier field.

    Used twice in the full pipeline, both times via the same zero-mean, distance-decaying
    machinery:
      - calibrator.BayesianCalibrator.fit_discrepancy(): ONE fit, on a temperature-space
        residual ('T_residual'), whose kernel this analytically marginalises into the KO
        likelihood -- the delta(x) term proper.
      - writeback.ResidualWriter: one fit PER calibration parameter (perm_frac/dT_scale/z0),
        on (posterior-calibrated value - Tier-2 transfer-model prediction), to fan the KO
        result out from the ~65 calibrated glaciers to all 217 glenglat glaciers (and beyond)
        -- the direct, reviewed replacement for the undamped k-NN residual copy.

    nugget_floor: dict field_name -> minimum WhiteKernel noise variance, so the GP cannot
    over-interpolate below the level its own data can resolve and keeps some shrinkage even
    exactly at a training location. For the per-parameter writeback fits, this is set from
    each parameter's Tier-1 grid-search step (variance ~ step^2/12); for the calibrator's
    'T_residual' fit it is set from the median observation-error variance instead (see
    BayesianCalibrator.fit_discrepancy).
    """

    nugget_floor: dict = field(default_factory=lambda: {
        'perm_frac': (0.2) ** 2 / 12,   # PERM_FRAC_GRID step = 0.2
        'dT_scale': (0.1) ** 2 / 12,    # DT_SCALE_GRID step = 0.1
        'z0': (5.0) ** 2 / 12,          # Z0_GRID step = 5.0
    })
    random_state: int = 42

    def __post_init__(self):
        self._gps = {}
        self._locations = {}

    def _fit_one(self, X_xyz, residuals, nugget_floor):
        # 3 columns -> isotropic great-circle kernel (the original behaviour, unchanged).
        # 4 columns -> ANISOTROPIC (ARD) kernel: the three spatial columns are kilometres on the
        # sphere embedding, the fourth is elevation in metres, so they cannot share a length
        # scale (see latlon_to_xyz). The spatial three are given the same initial value and the
        # same bounds, so they stay near-isotropic unless the data genuinely says otherwise --
        # a mild departure from an exact great-circle metric, acceptable because all calibration
        # glaciers sit in one region where the embedding is close to a local tangent plane.
        n_dim = np.asarray(X_xyz).shape[1]
        if n_dim == 4:
            ls = [LENGTH_SCALE_INIT_KM] * 3 + [ELEV_LENGTH_SCALE_INIT_M]
            ls_bounds = [LENGTH_SCALE_BOUNDS_KM] * 3 + [ELEV_LENGTH_SCALE_BOUNDS_M]
        else:
            ls, ls_bounds = LENGTH_SCALE_INIT_KM, LENGTH_SCALE_BOUNDS_KM
        kernel = (
            ConstantKernel(1.0, (1e-3, 1e3))
            * Matern(length_scale=ls, length_scale_bounds=ls_bounds, nu=1.5)
            + WhiteKernel(noise_level=max(nugget_floor, 1e-6),
                          noise_level_bounds=(max(nugget_floor, 1e-6) * 0.5, nugget_floor * 50 + 10.0))
        )
        gp = GaussianProcessRegressor(kernel=kernel, normalize_y=False, n_restarts_optimizer=10,
                                       random_state=self.random_state)
        gp.fit(X_xyz, residuals)
        return gp

    def fit(self, param_name, latitudes, longitudes, residuals, elevations=None):
        """Fit the spatial residual GP for one parameter ('perm_frac' | 'dT_scale' | 'z0').

        residuals: calibrated_value - global_transfer_model_prediction, one per glacier (see
        writeback.compute_residuals). normalize_y=False keeps the prior mean exactly 0 --
        "no nearby data => trust the global Tier-2 model unmodified".
        """
        X = latlon_to_xyz(np.asarray(latitudes), np.asarray(longitudes), elevations)
        gp = self._fit_one(X, np.asarray(residuals), self.nugget_floor[param_name])
        self._gps[param_name] = gp
        self._locations[param_name] = X
        return gp

    def predict(self, param_name, latitudes, longitudes, elevations=None):
        """Posterior (mean, std) of delta(x) for one parameter at new locations.

        `elevations` MUST be supplied iff the corresponding fit() supplied them -- a GP trained
        on 4 columns cannot predict from 3. Guarded below rather than left to fail deep inside
        sklearn with an opaque shape error."""
        gp = self._gps[param_name]
        n_trained = self._locations[param_name].shape[1]
        if (elevations is None) != (n_trained == 3):
            raise ValueError(
                f"discrepancy '{param_name}' was fit on {n_trained} coordinate columns but "
                f"predict was called {'without' if elevations is None else 'with'} elevations")
        X = latlon_to_xyz(np.asarray(latitudes), np.asarray(longitudes), elevations)
        mean, std = gp.predict(X, return_std=True)
        return mean, std

    def predict_all(self, latitudes, longitudes, elevations=None):
        """Posterior (mean, std) for all three parameters at once. Returns two dicts keyed by
        param name, each mapping to an array of shape (n_locations,)."""
        means, stds = {}, {}
        for name in self._gps:
            m, s = self.predict(name, latitudes, longitudes, elevations)
            means[name], stds[name] = m, s
        return means, stds

    def diagnostics(self, param_name):
        gp = self._gps[param_name]
        matern = gp.kernel_.k1.k2 if hasattr(gp.kernel_, 'k1') else None
        return {
            'kernel': str(gp.kernel_),
            'length_scale_km': getattr(matern, 'length_scale', np.nan),
            'log_marginal_likelihood': gp.log_marginal_likelihood_value_,
            'n_train': self._locations[param_name].shape[0],
        }
