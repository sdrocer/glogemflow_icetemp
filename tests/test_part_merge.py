"""merge_part_outputs: each point from the run of its hole's part (Step 5 plan, step 1c)."""
from types import SimpleNamespace

import numpy as np
import pytest

from icetemp.calibration.runner import GloGEMRunner


def entity(name, holes, depths):
    return SimpleNamespace(glacier_name=name, point_hole=np.array(holes, dtype=object), depths=np.array(depths, float))


def outputs(depths, value):
    d = np.array(list(depths) + [60.0])
    return (d, np.full(len(d), value) + np.arange(len(d)) * 0.01)


def test_points_take_their_part():
    g = entity('A@3005', ['1', '2', '1'], [16.0, 20.0, 30.0])
    runs = {p: {'A@3005': outputs(g.depths, v)} for p, v in (('G', -1.0), ('N', -5.0), ('S', -2.0))}
    m = GloGEMRunner.merge_part_outputs(runs, [g], {'1': 'N', '2': 'S'})
    d, t = m['A@3005']
    assert np.array_equal(d, [16.0, 20.0, 30.0, 60.0])
    assert np.allclose(t, [-5.0, -1.99, -4.98, -4.97])     # last: part of the deepest point's hole (N)


def test_missing_part_gives_none_and_other_parts_may_be_missing():
    g = entity('B@4005', ['3'], [18.0])
    runs = {'G': {'B@4005': None}, 'N': {'B@4005': outputs(g.depths, -7.0)}, 'S': {}}
    assert GloGEMRunner.merge_part_outputs(runs, [g], {'3': 'N'})['B@4005'][1][0] == -7.0
    assert GloGEMRunner.merge_part_outputs(runs, [g], {'3': 'G'})['B@4005'] is None


def test_depth_grids_must_agree_and_holes_are_needed():
    g = entity('C@3505', ['4'], [17.0])
    bad = {'G': {'C@3505': outputs([17.0], -1.0)}, 'N': {'C@3505': outputs([17.5], -1.0)}}
    with pytest.raises(ValueError):
        GloGEMRunner.merge_part_outputs(bad, [g], {'4': 'G'})
    g.point_hole = None
    with pytest.raises(ValueError):
        GloGEMRunner.merge_part_outputs({'G': {}}, [g], {})
