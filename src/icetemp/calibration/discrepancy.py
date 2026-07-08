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
LENGTH_SCALE_BOUNDS_KM = (10.0, 500.0)


def latlon_to_xyz(lat, lon):
    """Embed (lat, lon) on the unit sphere (scaled by Earth radius) so a standard isotropic
    kernel behaves like a great-circle-distance kernel."""
    lat_r, lon_r = np.radians(lat), np.radians(lon)
    x = EARTH_R_KM * np.cos(lat_r) * np.cos(lon_r)
    y = EARTH_R_KM * np.cos(lat_r) * np.sin(lon_r)
    z = EARTH_R_KM * np.sin(lat_r)
    return np.column_stack([x, y, z])


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
        kernel = (
            ConstantKernel(1.0, (1e-3, 1e3))
            * Matern(length_scale=LENGTH_SCALE_INIT_KM, length_scale_bounds=LENGTH_SCALE_BOUNDS_KM, nu=1.5)
            + WhiteKernel(noise_level=max(nugget_floor, 1e-6),
                          noise_level_bounds=(max(nugget_floor, 1e-6) * 0.5, nugget_floor * 50 + 10.0))
        )
        gp = GaussianProcessRegressor(kernel=kernel, normalize_y=False, n_restarts_optimizer=10,
                                       random_state=self.random_state)
        gp.fit(X_xyz, residuals)
        return gp

    def fit(self, param_name, latitudes, longitudes, residuals):
        """Fit the spatial residual GP for one parameter ('perm_frac' | 'dT_scale' | 'z0').

        residuals: calibrated_value - global_transfer_model_prediction, one per glacier (see
        writeback.compute_residuals). normalize_y=False keeps the prior mean exactly 0 --
        "no nearby data => trust the global Tier-2 model unmodified".
        """
        X = latlon_to_xyz(np.asarray(latitudes), np.asarray(longitudes))
        gp = self._fit_one(X, np.asarray(residuals), self.nugget_floor[param_name])
        self._gps[param_name] = gp
        self._locations[param_name] = X
        return gp

    def predict(self, param_name, latitudes, longitudes):
        """Posterior (mean, std) of delta(x) for one parameter at new locations."""
        gp = self._gps[param_name]
        X = latlon_to_xyz(np.asarray(latitudes), np.asarray(longitudes))
        mean, std = gp.predict(X, return_std=True)
        return mean, std

    def predict_all(self, latitudes, longitudes):
        """Posterior (mean, std) for all three parameters at once. Returns two dicts keyed by
        param name, each mapping to an array of shape (n_locations,)."""
        means, stds = {}, {}
        for name in self._gps:
            m, s = self.predict(name, latitudes, longitudes)
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
