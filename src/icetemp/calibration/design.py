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

    def sample(self, n, scramble=True):
        """Return an (n, len(PARAM_NAMES)) array of theta design points."""
        dists = self.priors.as_dict()
        d = len(PARAM_NAMES)
        sampler = qmc.LatinHypercube(d=d, scramble=scramble, seed=self.seed)
        u = sampler.random(n=n)
        return np.column_stack([dists[name].ppf(u[:, k])
                                for k, name in enumerate(PARAM_NAMES)])

    def sample_df(self, n, scramble=True):
        """Same as sample(), returned as a DataFrame with PARAM_NAMES columns (convenience for
        persisting/reloading a design, e.g. before a multi-hour training run)."""
        import pandas as pd
        theta = self.sample(n, scramble=scramble)
        return pd.DataFrame(theta, columns=PARAM_NAMES)
