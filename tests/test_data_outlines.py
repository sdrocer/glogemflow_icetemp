"""Loader rules of 2026-10-13: holes on the outline that contains them, entities split by outline,
deep-only profiles admitted, missing covariates filled per hole, and the table of left-out holes."""
import numpy as np
import pandas as pd
import pytest

from icetemp.calibration import data
from icetemp.calibration.data import DataHandler


@pytest.fixture
def mini(tmp_path, monkeypatch):
    """Glacier G: holes 1-3 in band 3005 (1, 2 on outline 00001; 3 lies in 00002 although its GLIMS
    id points to 00001; 2 is 'estimated' and has no derived covariates), hole 4 in band 2505 with
    readings below 30 m only."""
    g = tmp_path / 'glenglat'
    g.mkdir()
    pd.DataFrame(dict(id=[1, 2, 3, 4], glacier_name='G', glims_id='GA', latitude=46.0, longitude=8.0,
                      elevation=[3001.0, 3004.0, 3008.0, 2501.0], date_min='2000-08-01',
                      mass_balance_area='accumulation', temperature_uncertainty=0.1)
                 ).to_csv(g / 'borehole.csv', index=False)
    pd.DataFrame(dict(borehole_id=[1, 2, 3, 4], id=1, date_min='2000-08-01', date_max='2000-08-01',
                      equilibrium=['true', 'estimated', 'true', 'true'])).to_csv(g / 'profile.csv', index=False)
    pd.DataFrame(dict(borehole_id=[1, 1, 1, 2, 2, 3, 3, 4, 4], profile_id=1,
                      depth=[5.0, 15.0, 25.0, 12.0, 20.0, 15.0, 25.0, 40.0, 80.0],
                      temperature=[-3.0, -2.5, -2.0, -1.5, -1.0, -4.0, -3.5, -0.5, -0.2])
                 ).to_csv(g / 'measurement.csv', index=False)
    derived = tmp_path / 'derived.csv'
    pd.DataFrame(dict(borehole_id=[1, 3, 4], profile_id=1, T_maat=-8.0, T_amplitude=20.0, dT_firn=4.4,
                      regime='maritime', glogem_region='CentralEurope')).to_csv(derived, index=False)
    lookup = tmp_path / 'lookup.csv'
    pd.DataFrame(dict(glims_id=['GA', 'GB'], rgi7_official=['x', 'y'], rgi7_glogem=['RGI70-11.00001', 'RGI70-11.00002'],
                      area_km2=1.0, cenlat=46.0, cenlon=8.0)).to_csv(lookup, index=False)
    outlines = tmp_path / 'outlines.csv'
    pd.DataFrame(dict(borehole_id=[1, 2, 3, 4], glacier_name='G', rgi7=['00001', '00001', '00002', '00001'])
                 ).to_csv(outlines, index=False)
    monkeypatch.setattr(data, '_load_wgms_ela', lambda *a, **k: {})
    monkeypatch.setattr(DataHandler, '_era5_supplement', lambda self, row, year: (-7.0, 19.0, 'maritime', 'CentralEurope'))

    def load(**flags):
        dh = DataHandler(glenglat_dir=g, derived_csv=derived, glims_lookup_csv=lookup,
                         borehole_outline_csv=outlines, region='CentralEurope', **flags)
        dh.load()
        return dh
    return load


def test_new_rules(mini):
    dh = mini()
    ents = {e.glacier_name: e for e in dh.calibration_glaciers}
    assert set(ents) == {'G@3005~00001', 'G@3005~00002', 'G@2505'}
    assert ents['G@3005~00002'].glacier_id == '00002'
    # the 'estimated' hole without derived covariates keeps its readings
    assert {12.0, 20.0} <= set(ents['G@3005~00001'].depths)
    assert ents['G@2505'].depths.max() == 80.0
    assert dh.left_out.empty


def test_old_rules_and_left_out_reasons(mini):
    dh = mini(outline_by_position=False, entry_any_below=None, fill_missing_covariates=False)
    ents = {e.glacier_name: e for e in dh.calibration_glaciers}
    assert set(ents) == {'G@3005'}
    assert ents['G@3005'].glacier_id == '00001'
    assert not ({12.0, 20.0} & set(ents['G@3005'].depths))
    reasons = dict(zip(dh.left_out.borehole_id, dh.left_out.reason))
    assert reasons == {'2': 'other region or no covariates', '4': 'no reading at 10-30 m'}


def test_hole_without_outline_is_reported(mini, tmp_path):
    pd.DataFrame(dict(borehole_id=[1, 2, 3, 4], glacier_name='G', rgi7=['00001', '00001', '00002', np.nan])
                 ).to_csv(tmp_path / 'outlines.csv', index=False)
    lk = pd.read_csv(tmp_path / 'lookup.csv')
    lk.loc[0, 'glims_id'] = 'GX'                       # hole 4's GLIMS id no longer matches either
    lk.to_csv(tmp_path / 'lookup.csv', index=False)
    dh = mini()
    reasons = dict(zip(dh.left_out.borehole_id, dh.left_out.reason))
    assert reasons.get('4') == 'in the data but on no RGI7 outline: not run'
