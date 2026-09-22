"""
Priors: prior distributions for the KO calibration parameters (refreeze_frac, insul_scale,
advection_scale).

theta = (refreeze_frac, insul_scale, advection_scale) are GLOBAL scalars shared across all
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
fast or 70% suppressed is not a plausible correction), unlike insul_scale's much wider, only mildly
informative prior below.

refreeze_frac (f_ref) replaces perm_frac, retired because it was dead on bare ice, limited by
water availability rather than by depth, and coarsely quantised by the layer grid. f_ref scales
the melt+rain entering the column directly (fit_water in firnice_temperature_model.pro), so it
acts wherever there is any melt at all. Uniform prior on [0, 1] -- a fraction, no directional
belief.

insul_scale (f_ins) replaces dT_scale, retired with the CART decision tree that supplied its
offset: the tree jumped discontinuously across T_amplitude and elevation thresholds, and the
offset it scaled was an empirical surface-temperature correction rather than a physical process.
f_ins instead multiplies the Calonne (2011) snow/firn conductivity, so winter snow insulation
emerges from the seasonally varying snowpack the model already tracks.
"""

from dataclasses import dataclass

import numpy as np
from scipy import stats

from .physics import REFREEZE_FRAC_BOUNDS, INSUL_SCALE_BOUNDS

PARAM_NAMES = ('refreeze_frac', 'insul_scale', 'advection_scale')

# Gaussian prior for insul_scale: mean = 1.0 (Calonne 2011 as published); std = 0.35, so
# +/-2 std spans roughly [0.3, 1.7] -- the span of the published snow-conductivity relations
# (Sturm 1997 at 0.59, Yen 1981 at 1.09) plus room on either side.
INSUL_SCALE_PRIOR_MEAN = 1.0
INSUL_SCALE_PRIOR_STD = 0.35

# Gaussian prior for advection_scale: mean = settings.pro default (1.0, unscaled baseline
# physics); std = 0.25, truncated to [0.3, 1.7] (+/-2.8 std) -- see module docstring.
ADVECTION_SCALE_PRIOR_MEAN = 1.0
ADVECTION_SCALE_PRIOR_STD = 0.25
ADVECTION_SCALE_BOUNDS = (0.3, 1.7)


@dataclass
class Priors:
    """Prior distributions for (refreeze_frac, insul_scale, advection_scale), as scipy.stats frozen
    distributions.

    Each distribution exposes .pdf/.logpdf/.ppf/.rvs, used respectively by:
      - calibrator.py: log-prior terms in the KO log-posterior (.logpdf)
      - design.py: Latin Hypercube sampling via the inverse-CDF trick (.ppf)
      - plots.py: prior-density curves to compare against the posterior (.pdf)
    """

    refreeze_frac: object = None
    insul_scale: object = None
    advection_scale: object = None
    # Optional per-parameter (lo, hi) overrides, e.g. {'refreeze_frac': (0.02, 1.0)}. DIAGNOSTIC USE:
    # campaign 5's posterior pinned against the LOWER bound of both free parameters
    # (refreeze_frac 0.125 with a 0.1 floor, insul_scale 0.243 with a 0.2 floor, posterior sd 0.006 and
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
        if self.refreeze_frac is None:
            lo, hi = b.get('refreeze_frac', REFREEZE_FRAC_BOUNDS)
            self.refreeze_frac = stats.uniform(loc=lo, scale=hi - lo)
        if self.insul_scale is None:
            lo, hi = b.get('insul_scale', INSUL_SCALE_BOUNDS)
            a, bb = (lo - INSUL_SCALE_PRIOR_MEAN) / INSUL_SCALE_PRIOR_STD, (hi - INSUL_SCALE_PRIOR_MEAN) / INSUL_SCALE_PRIOR_STD
            self.insul_scale = stats.truncnorm(a, bb, loc=INSUL_SCALE_PRIOR_MEAN, scale=INSUL_SCALE_PRIOR_STD)
        if self.advection_scale is None:
            lo, hi = b.get('advection_scale', ADVECTION_SCALE_BOUNDS)
            a, bb = ((lo - ADVECTION_SCALE_PRIOR_MEAN) / ADVECTION_SCALE_PRIOR_STD,
                     (hi - ADVECTION_SCALE_PRIOR_MEAN) / ADVECTION_SCALE_PRIOR_STD)
            self.advection_scale = stats.truncnorm(
                a, bb, loc=ADVECTION_SCALE_PRIOR_MEAN, scale=ADVECTION_SCALE_PRIOR_STD)

    def as_dict(self):
        return {'refreeze_frac': self.refreeze_frac, 'insul_scale': self.insul_scale,
                'advection_scale': self.advection_scale}

    def bounds(self):
        b = self.bounds_override or {}
        return {'refreeze_frac': b.get('refreeze_frac', REFREEZE_FRAC_BOUNDS),
                'insul_scale': b.get('insul_scale', INSUL_SCALE_BOUNDS),
                'advection_scale': b.get('advection_scale', ADVECTION_SCALE_BOUNDS)}

    def logpdf(self, theta):
        """Joint log-prior density at theta = (refreeze_frac, insul_scale, advection_scale). Returns
        -inf outside the support of any marginal (keeps emcee's random walk inside physical
        bounds)."""
        pf, ds, adv = theta
        lp = (self.refreeze_frac.logpdf(pf) + self.insul_scale.logpdf(ds) + self.advection_scale.logpdf(adv))
        return lp if np.isfinite(lp) else -np.inf

    def rvs(self, size=1, random_state=None):
        """Draw `size` samples of theta directly from the prior (e.g. emcee walker init)."""
        rng = np.random.default_rng(random_state)
        pf = self.refreeze_frac.rvs(size=size, random_state=rng)
        ds = self.insul_scale.rvs(size=size, random_state=rng)
        adv = self.advection_scale.rvs(size=size, random_state=rng)
        return np.column_stack([pf, ds, adv])
