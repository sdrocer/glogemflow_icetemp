"""runner.parse_training_output: which model band an entity is compared with."""
from types import SimpleNamespace as E

import numpy as np

from icetemp.calibration.runner import GloGEMRunner


def write_profile(path, elev, value):
    lines = [f'Point elevation  {elev:4d} masl: Depth in m', 'Year  Month   1   2   3']
    for hy, m in ((1972, 10), (1972, 11)):
        lines.append(f'{hy:4d}{m:4d}' + f'{value:8.3f}' * 3)
    path.write_text('\n'.join(lines) + '\n')


def entity(elev, band):
    return E(glacier_name=f'Test@{band}', glacier_id='00001', elevation=elev, elevation_band=band,
             depths=np.array([2.0]), obs_months=[(np.array([197111]), np.array([1.0]), np.array([1.0]))])


def test_band_edge_prefers_own_band(tmp_path):
    write_profile(tmp_path / 'temp_ID1_00001.dat', 4295, -5.0)   # sorts first
    write_profile(tmp_path / 'temp_ID2_00001.dat', 4305, -7.0)
    out = GloGEMRunner().parse_training_output([entity(4300.0, 4305)], firnice_dir=tmp_path)
    assert abs(out['Test@4305'][1][0] - (-7.0)) < 1e-9


def test_nearest_band_still_wins(tmp_path):
    write_profile(tmp_path / 'temp_ID1_00001.dat', 4295, -5.0)
    write_profile(tmp_path / 'temp_ID2_00001.dat', 4315, -7.0)
    out = GloGEMRunner().parse_training_output([entity(4302.0, 4305)], firnice_dir=tmp_path)
    assert abs(out['Test@4305'][1][0] - (-5.0)) < 1e-9
