"""The model-error covariance of the depth target, in one place for the calibrator and for estimating its sizes
(plan of 20 Oct, Steps 2b and 3). Two residuals (model minus reading) share error when they share
  their band (entity):        s2 * exp(-|depth difference| / ell)
  their borehole:             sh2   (the hole's spot in the band)
  their borehole and year:    sy2   (one visit)
  their site:                 ss2   (all bands of one site)
and each reading at or above w_max_depth has an error of its own, w2 (its exact depth and month; rules plan of
27 Oct, form C6b). The readings' own errors (the diagonal) are added by the caller. With sh2 = sy2 = ss2 = w2 = 0
this is the covariance used before 2026-10-20."""
import numpy as np


def model_covariance(entity, depth, s2, ell, hole=None, year=None, site=None, sh2=0.0, sy2=0.0, ss2=0.0, w2=0.0,
                     w_max_depth=None):
    entity = np.asarray(entity)
    d = np.asarray(depth, dtype=float)
    C = np.zeros((len(d), len(d)))
    if s2 > 0:
        C += np.where(entity[:, None] == entity[None, :], s2 * np.exp(-np.abs(d[:, None] - d[None, :]) / ell), 0.0)
    if sh2 > 0 or sy2 > 0:
        if hole is None or year is None or any(h is None for h in hole) or (sy2 > 0 and any(y is None for y in year)):
            raise ValueError('borehole terms need every point\'s hole and year (pool_within="hole_year")')
        h = np.asarray(hole, dtype=object)
        same = h[:, None] == h[None, :]
        C += sh2 * same
        if sy2 > 0:
            y = np.asarray(year)
            C += sy2 * (same & (y[:, None] == y[None, :]))
    if ss2 > 0:
        if site is None:
            raise ValueError('the site term needs every point\'s site')
        s = np.asarray(site, dtype=object)
        C += ss2 * (s[:, None] == s[None, :])
    if w2 > 0:
        C[np.diag_indices(len(d))] += w2 * (d <= w_max_depth if w_max_depth is not None else np.ones(len(d), bool))
    return C
