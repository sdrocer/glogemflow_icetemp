"""Data rules of the Step 5 rules plan (27 Oct): sites by the written rule (rule 2), deep old ice below 40 m where the
measured accumulation is 0.6 m w.e. per year or less (rule 6 b, Janosch 28 Oct)."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from icetemp.calibration.data import FIT_EXCLUSIONS, fit_exclusions, low_accumulation_entities, rule_site_keys

HOLES = pd.DataFrame(dict(
    id=[1, 2, 3, 4, 5],
    glacier_name=['A', 'A', 'B', 'A', 'C'],
    latitude=[46.0, 46.006, 46.012, 46.0, 46.3],      # 1-2 and 2-3 about 0.67 km apart, 1-3 1.3 km
    longitude=[7.9, 7.9, 7.9, 7.9, 7.9],
    elevation=[3000.0, 3100.0, 3200.0, 3300.0, 3000.0]))


def ent(name, holes, depths=(5.0, 50.0)):
    return SimpleNamespace(glacier_name=name, point_hole=np.array(holes, dtype=object),
                           depths=np.asarray(depths, float), glogem_region='CentralEurope', glacier_id='00001',
                           elevation_band=4000)


def test_rule_sites_chain_and_split():
    gl = [ent('A@3005', ['1']), ent('A@3105', ['2']), ent('B@3205', ['3']), ent('A@3305', ['4']), ent('C@3005', ['5'])]
    s = rule_site_keys(gl, HOLES)
    assert s['A@3005'] == s['A@3105'] == s['B@3205'] == 'A 3000-3200'     # chained, though 1 and 3 are 1.3 km apart
    assert s['A@3305'] == 'A 3300-3300'                                   # 200 m above hole 3: its own site
    assert s['C@3005'] == 'C 3000-3000'


def test_rule_sites_refuse_split_entities():
    with pytest.raises(ValueError):
        rule_site_keys([ent('A@3005', ['1', '4'])], HOLES)
    with pytest.raises(ValueError):
        rule_site_keys([SimpleNamespace(glacier_name='A@3005', point_hole=None)], HOLES)


def test_low_accumulation_rule(tmp_path):
    t = tmp_path / 'acc.csv'
    pd.DataFrame(dict(entity=['A@3005', 'B@3205', 'C@3005', 'D@1', 'E@1'], acc_central=[0.5, 0.6, np.nan, 0.2, 0.61],
                      use=['used', 'used', 'no measurement', 'borrowed or basin value', 'used'])).to_csv(t, index=False)
    assert low_accumulation_entities(t) == {'A@3005': 0.5, 'B@3205': 0.6}
    rules = dict(FIT_EXCLUSIONS, test=dict(accumulation_table=t, max_accumulation=0.6, min_depth=40.0))
    gl = [ent('A@3005', ['1'], (5.0, 40.0, 41.0, 90.0)), ent('B@3205', ['3'], (50.0,)), ent('C@3005', ['5'], (60.0,))]
    assert fit_exclusions(gl, ['test'], rules) == {('A@3005', 41.0), ('A@3005', 90.0), ('B@3205', 50.0)}


def test_low_accumulation_table_is_complete():
    low = low_accumulation_entities()
    assert all(0 < v <= 0.6 for v in low.values())
    assert low["Grenzgletscher@4465~01225"] == 0.6          # on the line: left out (Janosch, 28 Oct)
    assert 'Glacier des Bossons@4305~00773' in low       # Dome du Gouter, 0.5 m w.e./a measured

