"""Known answers of the firn ledger reference (scripts/icetemp_analysis/firn_ledger.py)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts' / 'icetemp_analysis'))
import firn_ledger as fl  # noqa: E402


def test_layers_add_and_melt_from_the_top():
    # +0.5, +0.3, then -0.6 removes the 0.3 layer and 0.3 of the 0.5 layer
    out = fl.run([0.5, 0.3, -0.6], 2000)
    assert np.allclose([o[0] for o in out], [0.5, 0.8, 0.2])
    assert [o[1] for o in out] == [1, 2, 1] and out[-1][2] == 2000


def test_mass_conserving_removal_of_thin_layers():
    # three thin layers 0.1 each, then -0.25: exactly 0.25 goes, 0.05 of the oldest remains
    out = fl.run([0.1, 0.1, 0.1, -0.25], 2000)
    assert np.isclose(out[-1][0], 0.05) and out[-1][2] == 2000


def test_ageing_after_20_years():
    # a single layer survives through year 2020 (age 20) and is ice in 2021
    out = fl.run([1.0] + [0.0] * 21, 2000)
    assert out[20][1] == 1 and out[20][2] == 2000
    assert out[21][1] == 0 and out[21][0] == 0.0


def test_steady_state_start():
    out = fl.run([0.4], 2000, b_ref=0.4)
    # 21 layers before the run (1979-1999), +1 in 2000, 1979 ages out: 21 layers of 0.4
    assert np.isclose(out[0][0], 21 * 0.4) and out[0][1] == 21 and out[0][2] == 1980
    assert fl.run([0.0], 2000, b_ref=-0.2)[0] == (0, 0, 0)


def test_depth_on_a_density_profile():
    z = np.array([0.0, 1.0, 2.0])
    d = np.array([500.0, 500.0, 500.0])
    assert np.isclose(fl.depth(0.75, z, d), 1.5)
    assert np.isclose(fl.depth(5.0, z, d), 2.0)   # capped at the profile bottom
