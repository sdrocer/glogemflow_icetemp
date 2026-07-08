"""
DesignSampler: Latin Hypercube design over theta = (perm_frac, dT_scale, z0) for training runs.

Samples scipy.stats.qmc.LatinHypercube in the unit cube, then maps each dimension through its
prior's inverse-CDF (.ppf) -- the standard LHS-in-quantile-space trick. This automatically
space-fills z0 in LOG-space (since loguniform.ppf is the inverse of a log-uniform CDF) without
any separate log-handling logic, and equally naturally handles the truncated-Gaussian dT_scale
prior and the uniform perm_frac prior through the same interface.
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
        """Return an (n, 3) array of theta = (perm_frac, dT_scale, z0) design points."""
        sampler = qmc.LatinHypercube(d=3, scramble=scramble, seed=self.seed)
        u = sampler.random(n=n)  # (n, 3) in [0, 1)
        theta = np.column_stack([
            self.priors.perm_frac.ppf(u[:, 0]),
            self.priors.dT_scale.ppf(u[:, 1]),
            self.priors.z0.ppf(u[:, 2]),
        ])
        return theta

    def sample_df(self, n, scramble=True):
        """Same as sample(), returned as a DataFrame with PARAM_NAMES columns (convenience for
        persisting/reloading a design, e.g. before a multi-hour training run)."""
        import pandas as pd
        theta = self.sample(n, scramble=scramble)
        return pd.DataFrame(theta, columns=PARAM_NAMES)
