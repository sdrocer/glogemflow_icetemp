"""Leave-one-system-out folds: data.glacier_systems."""
from types import SimpleNamespace as E

from icetemp.calibration.data import glacier_systems


def test_monte_rosa_is_one_system():
    ents = [E(glacier_name='Grenzgletscher@4445', base_glacier_name='Grenzgletscher', glacier_id='01225'),
            E(glacier_name='Grenzgletscher@4455', base_glacier_name='Grenzgletscher', glacier_id='01230'),
            E(glacier_name='Gornergletscher@2525', base_glacier_name='Gornergletscher', glacier_id='01225'),
            E(glacier_name='Breithornplateau@3815', base_glacier_name='Breithornplateau', glacier_id='01225'),
            E(glacier_name='Hintereisferner@2555', base_glacier_name='Hintereisferner', glacier_id='03116'),
            E(glacier_name='Somewhere@3000', base_glacier_name='Somewhere', glacier_id=None)]
    s = glacier_systems(ents)
    monte_rosa = {s[e.glacier_name] for e in ents[:4]}
    assert len(monte_rosa) == 1
    assert s['Hintereisferner@2555'] == 'Hintereisferner'
    assert s['Somewhere@3000'] == 'Somewhere'
    assert len(set(s.values())) == 3
