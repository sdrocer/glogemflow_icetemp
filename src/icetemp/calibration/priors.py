"""
Priors: prior distributions for the KO calibration parameters (perm_frac, dT_scale, z0).

Per calibration_scheme_prompt.md / the task spec: log-uniform priors for parameters spanning
orders of magnitude (z0), Gaussian priors for parameters with an established literature/
default mean (dT_scale), uniform elsewhere (perm_frac -- no directional prior belief, only
physical bounds).

theta = (perm_frac, dT_scale, z0) are GLOBAL scalars shared across all calibration glaciers in
one GloGEM training run (see runner.write_calibration_override_single) -- matching the task
spec's "Calibrate perm_frac, dT_scale, and z0" as three scalar parameters, not per-glacier
ones. Prior means therefore reference GloGEM's own scalar defaults (settings.pro:167-169:
firnice_perm_frac=1.0, firnice_dT_scale=1.0, firnice_z0_firn=15.0), the best "literature mean"
available without conditioning on per-glacier covariates.
"""

from dataclasses import dataclass

import numpy as np
from scipy import stats

from .physics import PERM_FRAC_BOUNDS, DT_SCALE_BOUNDS, Z0_BOUNDS

PARAM_NAMES = ('perm_frac', 'dT_scale', 'z0')

# Gaussian prior for dT_scale: mean = settings.pro default (1.0, "no correction" baseline);
# std = 1.0 -- wide enough that +/-2 std roughly spans the full [0.2, 5.0] bound, so the prior
# is only mildly informative and lets the KO likelihood dominate once data is available.
DT_SCALE_PRIOR_MEAN = 1.0
DT_SCALE_PRIOR_STD = 1.0


@dataclass
class Priors:
    """Prior distributions for (perm_frac, dT_scale, z0), as scipy.stats frozen distributions.

    Each distribution exposes .pdf/.logpdf/.ppf/.rvs, used respectively by:
      - calibrator.py: log-prior terms in the KO log-posterior (.logpdf)
      - design.py: Latin Hypercube sampling via the inverse-CDF trick (.ppf)
      - plots.py: prior-density curves to compare against the posterior (.pdf)
    """

    perm_frac: object = None
    dT_scale: object = None
    z0: object = None

    def __post_init__(self):
        if self.perm_frac is None:
            lo, hi = PERM_FRAC_BOUNDS
            self.perm_frac = stats.uniform(loc=lo, scale=hi - lo)
        if self.dT_scale is None:
            lo, hi = DT_SCALE_BOUNDS
            a, b = (lo - DT_SCALE_PRIOR_MEAN) / DT_SCALE_PRIOR_STD, (hi - DT_SCALE_PRIOR_MEAN) / DT_SCALE_PRIOR_STD
            self.dT_scale = stats.truncnorm(a, b, loc=DT_SCALE_PRIOR_MEAN, scale=DT_SCALE_PRIOR_STD)
        if self.z0 is None:
            lo, hi = Z0_BOUNDS
            self.z0 = stats.loguniform(lo, hi)

    def as_dict(self):
        return {'perm_frac': self.perm_frac, 'dT_scale': self.dT_scale, 'z0': self.z0}

    def logpdf(self, theta):
        """Joint log-prior density at theta = (perm_frac, dT_scale, z0). Returns -inf outside
        the support of any marginal (keeps emcee's random walk inside physical bounds)."""
        pf, ds, z0 = theta
        lp = (self.perm_frac.logpdf(pf) + self.dT_scale.logpdf(ds) + self.z0.logpdf(z0))
        return lp if np.isfinite(lp) else -np.inf

    def rvs(self, size=1, random_state=None):
        """Draw `size` samples of theta directly from the prior (e.g. emcee walker init)."""
        rng = np.random.default_rng(random_state)
        pf = self.perm_frac.rvs(size=size, random_state=rng)
        ds = self.dT_scale.rvs(size=size, random_state=rng)
        z0 = self.z0.rvs(size=size, random_state=rng)
        return np.column_stack([pf, ds, z0])

    def bounds(self):
        return {'perm_frac': PERM_FRAC_BOUNDS, 'dT_scale': DT_SCALE_BOUNDS, 'z0': Z0_BOUNDS}
