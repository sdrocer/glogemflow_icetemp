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
    # Optional per-parameter (lo, hi) overrides, e.g. {'perm_frac': (0.02, 1.0)}. DIAGNOSTIC USE:
    # campaign 5's posterior pinned against the LOWER bound of both free parameters
    # (perm_frac 0.125 with a 0.1 floor, dT_scale 0.243 with a 0.2 floor, posterior sd 0.006 and
    # 0.014), i.e. the likelihood wants to go below the parameterisation. Widening lets us locate
    # where the unconstrained optimum actually is, which quantifies how far outside the
    # parameterisation the model is asking to be.
    #
    # This is legitimate to sample because the TRAINING path does not clip: the flat override
    # applier apply_firnicetemp_calibration.pro assigns firnice_*_b directly with no re-clip
    # (unlike the _bayes/_knn appliers, which do clip and are NOT used by training runs). And
    # the values in physics.py are described there as "settings.pro /
    # apply_firnicetemp_calibration_knn.pro re-clip bounds" -- operational limits, not hard
    # physics. Sub-floor values are a DIAGNOSTIC to locate the optimum; adopting one for
    # production would need separate physical justification.
    bounds_override: dict = None

    def __post_init__(self):
        b = self.bounds_override or {}
        if self.perm_frac is None:
            lo, hi = b.get('perm_frac', PERM_FRAC_BOUNDS)
            self.perm_frac = stats.uniform(loc=lo, scale=hi - lo)
        if self.dT_scale is None:
            lo, hi = b.get('dT_scale', DT_SCALE_BOUNDS)
            a, bb = (lo - DT_SCALE_PRIOR_MEAN) / DT_SCALE_PRIOR_STD, (hi - DT_SCALE_PRIOR_MEAN) / DT_SCALE_PRIOR_STD
            self.dT_scale = stats.truncnorm(a, bb, loc=DT_SCALE_PRIOR_MEAN, scale=DT_SCALE_PRIOR_STD)
        if self.z0 is None:
            lo, hi = b.get('z0', Z0_BOUNDS)
            self.z0 = stats.loguniform(lo, hi)

    def as_dict(self):
        return {'perm_frac': self.perm_frac, 'dT_scale': self.dT_scale, 'z0': self.z0}

    def bounds(self):
        b = self.bounds_override or {}
        return {'perm_frac': b.get('perm_frac', PERM_FRAC_BOUNDS),
                'dT_scale': b.get('dT_scale', DT_SCALE_BOUNDS),
                'z0': b.get('z0', Z0_BOUNDS)}

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
