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


def test_snow_refreeze_reaches_the_config(tmp_path):
    r = GloGEMRunner(run_dir=tmp_path / 'runs', catchments_dir=tmp_path / 'catchments')
    assert settings(r.write_training_config())[0]['firnice_snow_refreeze'] == 'n'
    r = GloGEMRunner(run_dir=tmp_path / 'runs', catchments_dir=tmp_path / 'catchments', snow_refreeze='y')
    assert settings(r.write_training_config())[0]['firnice_snow_refreeze'] == 'y'
