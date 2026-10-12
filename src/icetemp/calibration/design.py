"""
DesignSampler: Latin Hypercube design over theta = PARAM_NAMES for training runs.

Samples scipy.stats.qmc.LatinHypercube in the unit cube, then maps each dimension through its
prior's inverse-CDF (.ppf) -- the standard LHS-in-quantile-space trick, which handles the
uniform refreeze_frac prior and the truncated-Gaussian insul_scale / advection_scale priors
through one interface.

Dimensions are taken from Priors.as_dict() rather than named individually, so a parameter swap
needs no change here.
"""

from dataclasses import dataclass, field

import numpy as np
from scipy.stats import qmc

from .priors import Priors, PARAM_NAMES


@dataclass
class DesignSampler:
    priors: Priors = field(default_factory=Priors)
    seed: int = 42
    fixed_params: dict = None   # name -> value, held constant; the LHS spans the other parameters

    def sample(self, n, scramble=True):
        """Return an (n, len(PARAM_NAMES)) array of theta design points."""
        dists = self.priors.as_dict()
        fixed = self.fixed_params or {}
        free = [name for name in PARAM_NAMES if name not in fixed]
        sampler = qmc.LatinHypercube(d=len(free), scramble=scramble, seed=self.seed)
        u = sampler.random(n=n)
        cols = {name: dists[name].ppf(u[:, k]) for k, name in enumerate(free)}
        return np.column_stack([cols[name] if name in cols else np.full(n, float(fixed[name]))
                                for name in PARAM_NAMES])

    def sample_around(self, n, center, half_width, scramble=True):
        """n Latin-hypercube points spread uniformly over center +- half_width (arrays over
        PARAM_NAMES; free parameters, clipped to the prior bounds), fixed parameters constant:
        the second design stage, around a first result."""
        bounds = self.priors.bounds()
        fixed = self.fixed_params or {}
        free = [name for name in PARAM_NAMES if name not in fixed]
        u = qmc.LatinHypercube(d=len(free), scramble=scramble, seed=self.seed).random(n=n)
        cols = {}
        for k, name in enumerate(free):
            i = PARAM_NAMES.index(name)
            lo = max(bounds[name][0], center[i] - half_width[i])
            hi = min(bounds[name][1], center[i] + half_width[i])
            if not lo < hi:
                raise ValueError(f'{name}: empty box {lo}..{hi} around {center[i]}')
            cols[name] = lo + u[:, k] * (hi - lo)
        return np.column_stack([cols[name] if name in cols else np.full(n, float(fixed[name]))
                                for name in PARAM_NAMES])

    def sample_df(self, n, scramble=True):
        """Same as sample(), returned as a DataFrame with PARAM_NAMES columns (convenience for
        persisting/reloading a design, e.g. before a multi-hour training run)."""
        import pandas as pd
        theta = self.sample(n, scramble=scramble)
        return pd.DataFrame(theta, columns=PARAM_NAMES)
