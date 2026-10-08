"""Observation matching: the model read at each reading's date (runner.match_model_to_obs)."""
import calendar

import numpy as np
import pandas as pd

from icetemp.calibration.data import _reading_months
from icetemp.calibration.runner import match_model_to_obs

DEPTHS = np.array([1.0, 2.0, 3.0])


def end_of_month_profiles(fn, hyd_years=(1972, 1973)):
    """Profile file layout: hydrological-year label, calendar month, one column per depth,
    each row the value at the END of its month."""
    rows = []
    for hy in hyd_years:
        for m in list(range(10, 13)) + list(range(1, 10)):
            y = hy - 1 if m >= 10 else hy
            t = y + m / 12   # end of month m
            rows.append([hy, m] + [fn(t, z) for z in DEPTHS])
    return pd.DataFrame(rows, columns=['year', 'month'] + list(DEPTHS))


class Point:
    def __init__(self, depths, obs_months):
        self.depths = np.asarray(depths, float)
        self.obs_months = obs_months


def test_linear_in_time_recovered_exactly():
    # a field linear in time is recovered exactly at any day of a month
    df = end_of_month_profiles(lambda t, z: 3.0 * t - z)
    ym, f = 197207, 10 / 31
    g = Point([2.0], [(np.array([ym]), np.array([1.0]), np.array([f]))])
    _, v = match_model_to_obs(df, g)
    t = 1972 + (6 + f) / 12
    assert abs(v[0] - (3.0 * t - 2.0)) < 1e-9


def test_sine_phase_within_a_day():
    period, amp = 1.0, 8.0
    fn = lambda t, z: amp * np.sin(2 * np.pi * t / period)
    df = end_of_month_profiles(fn)
    for ym, day in ((197203, 15), (197207, 10), (197211, 20)):
        y, m = divmod(ym, 100)
        nd = calendar.monthrange(y, m)[1]
        g = Point([1.0], [(np.array([ym]), np.array([1.0]), np.array([day / nd]))])
        _, v = match_model_to_obs(df, g)
        t = y + (m - 1 + day / nd) / 12
        # chord of a sine over one month: error below amp * (2 pi / 12)^2 / 8
        assert abs(v[0] - fn(t, 1.0)) < amp * (2 * np.pi / 12) ** 2 / 8 + 1e-9


def test_last_day_and_first_month_give_month_end():
    df = end_of_month_profiles(lambda t, z: 10 * t + z)
    for ym in (197203, 197110):           # an ordinary month, the run's first month
        g = Point([3.0], [(np.array([ym]), np.array([1.0]), np.array([1.0 if ym == 197203 else 0.3]))])
        _, v = match_model_to_obs(df, g)
        y, m = divmod(ym, 100)
        hy = y + 1 if m >= 10 else y
        expected = df[(df.year == hy) & (df.month == m)][3.0].iloc[0]
        assert abs(v[0] - expected) < 1e-12


def test_old_two_entry_months_unchanged():
    df = end_of_month_profiles(lambda t, z: 5 * t + z)
    g2 = Point([1.0], [(np.array([197205, 197206]), np.array([2, 1]))])
    g3 = Point([1.0], [(np.array([197205, 197206]), np.array([2.0, 1.0]), np.array([1.0, 1.0]))])
    assert np.allclose(match_model_to_obs(df, g2)[1], match_model_to_obs(df, g3)[1])


def test_reading_months_ranges():
    assert _reading_months(197207, 197207, 0.4, 197110, 202509) == ([197207], [1.0], [0.4])
    yms, w, f = _reading_months(197011, 197202, 0.5, 197110, 202509)
    assert yms == [197110, 197111, 197112, 197201, 197202]
    assert abs(sum(w) - 1) < 1e-12 and f == [0.5] * 5
