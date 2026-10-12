"""Synthetic calibration setup: four entities and an emulator fit on 24 runs of a smooth model."""
from types import SimpleNamespace

import numpy as np
import pytest

from icetemp.calibration.design import DesignSampler
from icetemp.calibration.emulator import Emulator
from icetemp.calibration.priors import Priors


def glacier(name, gid, band, depths, lat, lon, rng):
    depths = np.asarray(depths, dtype=float)
    return SimpleNamespace(glacier_name=name, base_glacier_name=name.split('@')[0], glacier_id=gid,
                           glogem_region='CentralEurope', elevation=float(band), elevation_band=band,
                           latitude=lat, longitude=lon, depths=depths,
                           T_obs=-4.0 + 0.02 * depths + rng.normal(0, 0.5, len(depths)),
                           weights=np.where((depths >= 10) & (depths <= 30), 3.0, 1.0),
                           sigma=np.full(len(depths), 0.3))


def model(theta, g):
    f_ref, f_ins, f_adv = theta
    z = g.depths
    return (-6.0 + 2.5 * f_ref * np.exp(-z / 12.0) + 1.2 * (f_ins - 1.0) * np.exp(-z / 30.0)
            - 0.4 * (f_adv - 1.0) + 0.01 * z)


@pytest.fixture(scope='module')
def setup():
    rng = np.random.default_rng(0)
    gl = [glacier('A@3005', '00001', 3005, [2, 5, 9, 14, 19, 24], 45.9, 7.8, rng),
          glacier('B@3505', '00002', 3505, [1, 3, 6, 12, 18], 46.0, 7.9, rng),
          glacier('C@4455', '01230', 4455, [5, 15, 25, 35, 45, 55], 45.93, 7.88, rng),
          glacier('D@2605', '00004', 2605, [4, 8, 16, 32], 46.8, 10.8, rng)]
    theta = DesignSampler(Priors(), seed=3).sample(24)
    outputs = [{g.glacier_name: (g.depths, model(t, g) + rng.normal(0, 0.02, len(g.depths)))
                for g in gl} for t in theta]
    emu = Emulator(glaciers=gl, explained_variance=0.999).fit(theta, outputs)
    return gl, theta, outputs, emu
