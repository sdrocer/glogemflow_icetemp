"""
Thermal-structure classification: temperate / cold / polythermal.

Uses the SAME real pressure-melting-point (PMP) physics GloGEM itself applies as a clamp
(GloGEM/procedures/processing/firnice_temperature_model.pro), applied IDENTICALLY to model
predictions and real glenglat observations -- see classify_point's docstring for why that
symmetry is a deliberate decision, not an oversight.

Companion diagnostics (robustness, elevation_spread) exist because two things were confirmed
empirically against the real CentralEurope held-out glacier set, not assumed:
  1. Half of 24 real held-out glaciers have their warmest observed point within 0.3 degC of the
     PMP threshold -- classification is genuinely tolerance-sensitive for about half the data,
     not a rare edge case (see robustness()).
  2. Most calibration glaciers pool multiple real glenglat boreholes under one glacier_name,
     spanning real elevation ranges up to 2470 m (Grenzgletscher: Colle Gnifetti firn saddle
     ~3900-4485m pooled with the lower glacier tongue ~2500-2600m) -- and the real model was
     only ever asked to predict at ONE representative elevation per glacier_name (see
     runner.GloGEMRunner.write_glenglat_lookup), so a pooled classification is not necessarily
     a fair "whole glacier" statement (see elevation_spread()).
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .data import GLENGLAT_DIR

# PMP formula ported from GloGEM/procedures/processing/firnice_temperature_model.pro line 61:
#     pmp_profile[jp] = (fit_dz[1,jp] * 0.9d0 / 10.0d0) * (-0.00742d0)
# fit_dz[1,*] confirmed as depth in metres (prepare_output_firnicetemp.pro line 84: "Depth in
# m"). E.g. -0.067 degC at 100 m depth, -0.334 degC at 500 m -- small, but a real, physically
# derived threshold, not an arbitrary one.
PMP_COEF = 0.09 * (-0.00742)  # degC per metre depth


def pressure_melting_point(depth_m):
    """PMP(z), degC -- the warmest ice can be at depth z (m) without melting."""
    return depth_m * PMP_COEF


# Fallback tolerance for the ~34% of boreholes with no reported glenglat temperature_uncertainty
# (551/835 report one; median 0.10 degC, 75th percentile 0.20 degC). 0.2 degC matches real
# thermistor precision and sits at roughly that 75th percentile -- a deliberately
# conservative-but-not-excessive default, not a guess.
DEFAULT_TOLERANCE = 0.2

# Bracket used by robustness() to test whether a classification is sensitive to the exact
# tolerance chosen. Exists because of finding (1) in the module docstring, not as a generic
# hedge -- with half the real held-out set within 0.3 degC of the threshold, a single tolerance
# number needs a sensitivity check, not blind trust.
ROBUSTNESS_BRACKET = (0.05, 0.5)

# Elevation spread beyond which a pooled glacier-level classification stops being a trustworthy
# "whole glacier" statement (see module docstring finding (2)). A blunt, documented cutoff, not
# a statistical test: >100 m spans several GloGEM 10 m elevation bands many times over, wide
# enough that along-flow structure (e.g. a cold-front/temperate-upglacier pattern) could exist
# within it, undetectable by one pooled point.
ELEVATION_SPREAD_CAUTION_M = 100.0

# Depth of zero annual amplitude: below this, the annual surface-temperature wave has
# effectively damped out. A reading SHALLOWER than this can read "cold" simply because the
# borehole was logged in or shortly after winter, on a glacier that is temperate at every depth
# that actually matters -- a seasonal artifact, not persistent thermal structure.
#
# 15 m is the middle of the ~10-20 m range commonly cited for ice/firn. NOT derived from
# GloGEM's own thermal constants: firnice_temperature_model.pro computes per-layer heat
# capacity (cap_fit, line 77) but the conductivity/diffusivity constants needed to derive a
# project-specific damping depth were not located, so this is an honest literature default,
# deliberately exposed as a parameter rather than hardcoded as if it were GloGEM-derived.
#
# WHY THIS IS A FLAG AND NOT A HARD FILTER: excluding shallow readings outright was tested
# against the real CentralEurope held-out set and rejected on cost -- dropping everything
# shallower than 15/20/30 m leaves 5/9/14 of 24 glaciers with NO usable data at all (of 428
# total points, only 281/243/182 survive). With polythermal recall already resting on just 10
# true-positive glaciers, that much attrition would make any campaign-to-campaign comparison
# too noisy to interpret. Flagging preserves the sample while staying honest about which calls
# are confirmed by persistent deep readings and which could be seasonal.
SEASONAL_DEPTH_M = 15.0


def classify_point(T, depth_m, tolerance):
    """Single depth reading -> 'temperate' or 'cold'.

    Applies the SAME depth-adjusted PMP threshold to model predictions and real observations
    alike. An earlier version of this analysis judged predictions against the C&P analytical
    surrogate's raw 0 degC cap (physics.cp_model_single has no pressure term at all) while
    judging real observations against the true depth-adjusted PMP -- a real asymmetry, however
    small in magnitude, that judged the model slightly more generously than the data at depth.
    Fixed on an explicit decision: model and data are judged by the identical rule, always.
    """
    return 'temperate' if T >= pressure_melting_point(depth_m) - tolerance else 'cold'


def classify_profile(labels):
    """Point labels -> one profile-level call: 'temperate' if every point is temperate, 'cold'
    if every point is cold, 'polythermal' if both appear."""
    u = set(labels)
    if u == {'temperate'}:
        return 'temperate'
    if u == {'cold'}:
        return 'cold'
    return 'polythermal'


def classify(T_values, depths, tolerance):
    """classify_profile over classify_point for each (T, depth) pair -- the usual entry point."""
    labels = [classify_point(t, d, tolerance) for t, d in zip(T_values, depths)]
    return classify_profile(labels)


def robustness(T_values, depths, bracket=ROBUSTNESS_BRACKET):
    """'robust' if the profile classification is UNCHANGED across `bracket` (default 0.05-0.5
    degC), else 'borderline'. Lets a downstream reader distinguish a solid classification from
    one sitting on a knife-edge, rather than reporting every glacier with equal implied
    confidence -- see module docstring finding (1) for why this matters in practice, not just
    in principle."""
    lo, hi = bracket
    return 'robust' if classify(T_values, depths, lo) == classify(T_values, depths, hi) else 'borderline'


def seasonal_confidence(T_values, depths, tolerance, seasonal_depth_m=SEASONAL_DEPTH_M):
    """Is this profile's classification confirmed by readings deep enough to be free of the
    annual surface-temperature wave? Returns one of:

      'confirmed_at_depth' -- classification is UNCHANGED when only readings deeper than
                              seasonal_depth_m are used. The call rests on persistent
                              structure, not a possible winter artifact.
      'shallow_dependent'  -- deep readings exist, but restricting to them CHANGES the
                              classification. The call depends on shallow points that could be
                              seasonal.
      'shallow_only'       -- no readings deeper than seasonal_depth_m at all, so the
                              classification cannot be confirmed either way.

    See SEASONAL_DEPTH_M for why this is a flag rather than a hard depth filter.
    """
    depths = np.asarray(depths, dtype=float)
    T_values = np.asarray(T_values, dtype=float)
    deep = depths > seasonal_depth_m
    if not np.any(deep):
        return 'shallow_only'
    if classify(T_values, depths, tolerance) == classify(T_values[deep], depths[deep], tolerance):
        return 'confirmed_at_depth'
    return 'shallow_dependent'


def load_borehole_table(glenglat_dir=GLENGLAT_DIR):
    """glenglat's borehole table, with `id` cast to str to match data.py's own convention.

    The cast is load-bearing, not cosmetic: borehole.csv parses `id` as int64, while
    GlacierCalibrationData.borehole_id is a STRING (data.py:230 does its own
    `borehole['id'].astype(str)`, and data.py:414 stores `str(bh_row['borehole_id'])`). Without
    the cast here, borehole_tolerance's `==` comparison is int64-vs-str, never matches, and
    silently returns DEFAULT_TOLERANCE for 100% of boreholes -- so the per-borehole
    temperature_uncertainty this module's whole tolerance design rests on would never be used,
    with no error raised. Confirmed live before the fix: all 24 CentralEurope held-out glaciers
    received 0.2 where the intended values span 0.05-0.5, and correcting it flips Taconnaz
    Glacier's observed class (polythermal -> cold, its reported uncertainty being 0.10 not 0.2).
    """
    df = pd.read_csv(Path(glenglat_dir) / 'borehole.csv')
    df['id'] = df['id'].astype(str)
    return df


def borehole_tolerance(borehole_id, borehole_df, default=DEFAULT_TOLERANCE):
    """Real tolerance from glenglat's own temperature_uncertainty column for this borehole,
    falling back to `default` if missing or if borehole_id isn't found.

    NOTE ON GRANULARITY: GlacierCalibrationData pools potentially many real boreholes under one
    glacier_name (see module docstring finding (2)), retaining only ONE representative
    borehole_id for the whole pooled profile -- so this currently returns one tolerance for an
    entire glacier's pooled readings, not a true per-reading value for multi-borehole glaciers.
    This mirrors the SAME single-representative-point limitation elevation_spread() flags for
    elevation; fixing it fully would need per-reading borehole traceability that the pooling
    step (data._pool_by_depth_bin) does not currently retain. Documented here rather than
    silently assumed away -- deliberately not overclaiming precision the data doesn't carry.
    """
    row = borehole_df.loc[borehole_df['id'] == str(borehole_id)]
    if row.empty:
        return default
    val = row.iloc[0]['temperature_uncertainty']
    return default if pd.isna(val) else float(val)


@dataclass
class ElevationSpread:
    """How much real spatial spread a glacier_name's pooled boreholes span, versus the SINGLE
    representative elevation the model was actually asked to predict at (see
    runner.GloGEMRunner.write_glenglat_lookup: one elevation per glacier_name, from a single
    representative borehole). See module docstring finding (2)."""
    glacier_name: str
    n_locations: int
    elevation_range_m: float
    elevation_min: float
    elevation_max: float

    @property
    def caution(self):
        return self.elevation_range_m > ELEVATION_SPREAD_CAUTION_M


def elevation_spread(glacier_name, borehole_df):
    """None if the glacier isn't found or has no elevation data."""
    sub = borehole_df[borehole_df['glacier_name'] == glacier_name]
    locs = sub[['latitude', 'longitude']].drop_duplicates()
    elevs = sub['elevation'].dropna()
    if elevs.empty:
        return None
    return ElevationSpread(
        glacier_name=glacier_name, n_locations=len(locs),
        elevation_range_m=float(elevs.max() - elevs.min()),
        elevation_min=float(elevs.min()), elevation_max=float(elevs.max()),
    )
