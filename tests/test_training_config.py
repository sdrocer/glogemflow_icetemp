"""The generated GloGEM training config: the thermal forcing switches, defaults = before 2026-10-12."""
import re

from icetemp.calibration.runner import GloGEMRunner


def settings(path):
    text = path.read_text()
    return dict(re.findall(r"^(firnice_\w+)\s*=\s*'([^']*)'", text, flags=re.M)), text


def test_defaults_keep_todays_physics(tmp_path):
    r = GloGEMRunner(run_dir=tmp_path / 'runs', catchments_dir=tmp_path / 'catchments')
    s, text = settings(r.write_training_config())
    assert (s['firnice_forcing_toff'], s['firnice_hl1980'], s['firnice_lh_water'], s['firnice_perm_limit'],
            s['firnice_thermal_lapse']) == ('y', 'n', 'n', 'hl', 'floor')
    assert not re.search(r'\{\w+\}', text)        # every template field filled


def test_switches_reach_the_config(tmp_path):
    r = GloGEMRunner(run_dir=tmp_path / 'runs', catchments_dir=tmp_path / 'catchments', forcing_toff='n',
                     hl1980='y', lh_water='y', perm_limit='profile', catchment='CentralEurope_gate')
    s, text = settings(r.write_training_config())
    assert (s['firnice_forcing_toff'], s['firnice_hl1980'], s['firnice_lh_water'], s['firnice_perm_limit']) == \
        ('n', 'y', 'y', 'profile')
    assert "catchment_selection   = 'CentralEurope_gate'" in text


def test_lapse_and_cold_percolation_reach_the_config(tmp_path):
    r = GloGEMRunner(run_dir=tmp_path / 'runs', catchments_dir=tmp_path / 'catchments', perm_limit='cold',
                     thermal_lapse='reanalysis')
    s, _ = settings(r.write_training_config())
    assert (s['firnice_perm_limit'], s['firnice_thermal_lapse']) == ('cold', 'reanalysis')


def test_firn_ledger_reaches_the_config(tmp_path):
    r = GloGEMRunner(run_dir=tmp_path / 'runs', catchments_dir=tmp_path / 'catchments', perm_limit='firn')
    assert settings(r.write_training_config())[0]['firnice_perm_limit'] == 'firn'


def test_mb_hist_file_reaches_the_config(tmp_path):
    hist = lambda r: re.findall(r"^mb_hist_file\s*=\s*'([^']*)'", r.write_training_config().read_text(), flags=re.M)
    assert hist(GloGEMRunner(run_dir=tmp_path / 'runs', catchments_dir=tmp_path / 'catchments')) == ['']
    assert hist(GloGEMRunner(run_dir=tmp_path / 'runs', catchments_dir=tmp_path / 'catchments',
                             mb_hist_file='/x/h.dat')) == ['/x/h.dat']


def test_snow_refreeze_reaches_the_config(tmp_path):
    r = GloGEMRunner(run_dir=tmp_path / 'runs', catchments_dir=tmp_path / 'catchments')
    assert settings(r.write_training_config())[0]['firnice_snow_refreeze'] == 'n'
    r = GloGEMRunner(run_dir=tmp_path / 'runs', catchments_dir=tmp_path / 'catchments', snow_refreeze='y')
    assert settings(r.write_training_config())[0]['firnice_snow_refreeze'] == 'y'


def test_firn_switches_reach_the_config(tmp_path):
    def read(r):
        s, text = settings(r.write_training_config())
        age = re.search(r'^firnice_firn_age\s*=\s*(\d+)', text, flags=re.M).group(1)
        return age, s['firnice_firn_obs_dir'], s['firnice_firn_obs_mask'], s['firnice_firn_column']
    r = GloGEMRunner(run_dir=tmp_path / 'runs', catchments_dir=tmp_path / 'catchments')
    assert read(r) == ('20', '', 'n', 'fixed')
    r = GloGEMRunner(run_dir=tmp_path / 'runs', catchments_dir=tmp_path / 'catchments', firn_age=30,
                     firn_obs_dir='/x/firnobs', firn_obs_mask='y', firn_column='ledger')
    assert read(r) == ('30', '/x/firnobs', 'y', 'ledger')


def test_mass_balance_refreezing_reaches_the_config(tmp_path):
    def read(r):
        text = r.write_training_config().read_text()
        return tuple(re.search(rf"^{k}\s*=\s*'(\w)'", text, flags=re.M).group(1)
                     for k in ('refreezing_parametrised', 'refreezing_full'))
    assert read(GloGEMRunner(run_dir=tmp_path / 'runs', catchments_dir=tmp_path / 'catchments')) == ('y', 'n')
    assert read(GloGEMRunner(run_dir=tmp_path / 'r2', catchments_dir=tmp_path / 'c2', mb_refreezing='full')) == ('n', 'y')
