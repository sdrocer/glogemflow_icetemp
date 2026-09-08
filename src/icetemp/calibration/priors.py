"""
Priors: prior distributions for the KO calibration parameters (perm_frac, dT_scale,
advection_scale).

theta = (perm_frac, dT_scale, advection_scale) are GLOBAL scalars shared across all
calibration glaciers in one GloGEM training run (see runner.write_calibration_override_single).

z0 (formerly the third calibrated parameter) has been RETIRED from theta: calibrator.py's own
diagnostics found spearman(z0, output) = -0.031, p=0.63 across four campaigns -- indistinguishable
from zero -- because z0 only ever shapes the analytical C&P SPINUP profile
(initialise_firnicetemp_spinup.pro) and never reaches the transient physics afterward. A design
point's z0 therefore could never reach the model output it was being calibrated against; its
posterior never narrowed in any campaign. z0 now stays fixed at its settings.pro default (15.0 m)
for every run.

advection_scale (f_adv) replaces it: a per-band multiplier on the depth-averaged advection
velocity u, applied every substep in firnice_temperature_model.pro (GloGEM/procedures/processing/
firnice_temperature_model.pro), so -- unlike z0 -- it genuinely reaches the transient physics that
produces the output being calibrated against. Prior: Gaussian, mean 1.0 (unscaled baseline
physics, GloGEM/procedures/initialise/settings.pro:173's firnice_adv_scale default), std 0.25,
truncated to [0.3, 1.7] -- informative around the physical baseline (an advection field 3x too
fast or 70% suppressed is not a plausible correction), unlike dT_scale's much wider, only mildly
informative prior below.

perm_frac keeps its uniform prior on physical bounds only (no directional prior belief) and
dT_scale keeps its Gaussian prior centered on the settings.pro default -- both unchanged from
before this swap.
"""

from dataclasses import dataclass

import numpy as np
from scipy import stats

from .physics import PERM_FRAC_BOUNDS, DT_SCALE_BOUNDS

PARAM_NAMES = ('perm_frac', 'dT_scale', 'advection_scale')

# Gaussian prior for dT_scale: mean = settings.pro default (1.0, "no correction" baseline);
# std = 1.0 -- wide enough that +/-2 std roughly spans the full [0.2, 5.0] bound, so the prior
# is only mildly informative and lets the KO likelihood dominate once data is available.
DT_SCALE_PRIOR_MEAN = 1.0
DT_SCALE_PRIOR_STD = 1.0

# Gaussian prior for advection_scale: mean = settings.pro default (1.0, unscaled baseline
# physics); std = 0.25, truncated to [0.3, 1.7] (+/-2.8 std) -- see module docstring.
ADVECTION_SCALE_PRIOR_MEAN = 1.0
ADVECTION_SCALE_PRIOR_STD = 0.25
ADVECTION_SCALE_BOUNDS = (0.3, 1.7)


@dataclass
class Priors:
    """Prior distributions for (perm_frac, dT_scale, advection_scale), as scipy.stats frozen
    distributions.

    Each distribution exposes .pdf/.logpdf/.ppf/.rvs, used respectively by:
      - calibrator.py: log-prior terms in the KO log-posterior (.logpdf)
      - design.py: Latin Hypercube sampling via the inverse-CDF trick (.ppf)
      - plots.py: prior-density curves to compare against the posterior (.pdf)
    """

    perm_frac: object = None
    dT_scale: object = None
    advection_scale: object = None
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
        if self.advection_scale is None:
            lo, hi = b.get('advection_scale', ADVECTION_SCALE_BOUNDS)
            a, bb = ((lo - ADVECTION_SCALE_PRIOR_MEAN) / ADVECTION_SCALE_PRIOR_STD,
                     (hi - ADVECTION_SCALE_PRIOR_MEAN) / ADVECTION_SCALE_PRIOR_STD)
            self.advection_scale = stats.truncnorm(
                a, bb, loc=ADVECTION_SCALE_PRIOR_MEAN, scale=ADVECTION_SCALE_PRIOR_STD)

    def as_dict(self):
        return {'perm_frac': self.perm_frac, 'dT_scale': self.dT_scale,
                'advection_scale': self.advection_scale}

    def bounds(self):
        b = self.bounds_override or {}
        return {'perm_frac': b.get('perm_frac', PERM_FRAC_BOUNDS),
                'dT_scale': b.get('dT_scale', DT_SCALE_BOUNDS),
                'advection_scale': b.get('advection_scale', ADVECTION_SCALE_BOUNDS)}

    def logpdf(self, theta):
        """Joint log-prior density at theta = (perm_frac, dT_scale, advection_scale). Returns
        -inf outside the support of any marginal (keeps emcee's random walk inside physical
        bounds)."""
        pf, ds, adv = theta
        lp = (self.perm_frac.logpdf(pf) + self.dT_scale.logpdf(ds) + self.advection_scale.logpdf(adv))
        return lp if np.isfinite(lp) else -np.inf

    def rvs(self, size=1, random_state=None):
        """Draw `size` samples of theta directly from the prior (e.g. emcee walker init)."""
        rng = np.random.default_rng(random_state)
        pf = self.perm_frac.rvs(size=size, random_state=rng)
        ds = self.dT_scale.rvs(size=size, random_state=rng)
        adv = self.advection_scale.rvs(size=size, random_state=rng)
        return np.column_stack([pf, ds, adv])
