"""The fair-test scorer (scripts/icetemp_analysis/forcing_fair_test.py) on cases with a known answer."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location(
    'fft', Path(__file__).resolve().parents[1] / 'scripts/icetemp_analysis/forcing_fair_test.py')
fft = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fft)


def frame():
    """Two massifs, two outlines each, three points per outline; runs A@x are the observation plus
    a known offset, so M1 is that offset."""
    rows = []
    for m, rgis in (('M1', ['a', 'b']), ('M2', ['c', 'd'])):
        for r in rgis:
            for k in range(3):
                rows.append(dict(entity=f'{r}@3000', rgi=r, massif=m, band=3000, depth=20.0 + k, obs=-1.0))
    P = pd.DataFrame(rows)
    for fr, off in ((0.10, 2.0), (0.55, 0.5), (1.00, 1.0)):
        P[fft.col('A', fr)] = P.obs + off
    return P


def test_m1_weights_outlines_and_drops_small_ones():
    P = frame()
    P.loc[P.rgi == 'a', 'A@0.55'] += 3.0             # outline a: error 3.5, others 0.5
    assert np.isclose(fft.m1(P, 'A@0.55'), np.sqrt((3.5 ** 2 + 3 * 0.5 ** 2) / 4))
    small = P[~((P.rgi == 'b') & (P.depth > 20.0))]  # outline b keeps 1 point: left out
    assert np.isclose(fft.m1(small, 'A@0.55'), np.sqrt((3.5 ** 2 + 2 * 0.5 ** 2) / 3))


def test_lomo_picks_the_best_fref_on_the_other_massif():
    P = frame()
    P.loc[P.massif == 'M2', 'A@1.00'] = P.loc[P.massif == 'M2', 'obs']   # 1.00 perfect on M2 only
    pred, chosen = fft.lomo(P, {fr: fft.col('A', fr) for fr in (0.10, 0.55, 1.00)})
    assert chosen == {'M1': 1.00, 'M2': 0.55}        # M1 is scored with M2's best, and vice versa
    assert np.allclose(pred[P.massif == 'M1'], P.loc[P.massif == 'M1', 'obs'] + 1.0)


def test_lomo_leaves_resampled_copies_out_together():
    P = frame()
    P.loc[P.massif == 'M1', 'A@1.00'] = P.loc[P.massif == 'M1', 'obs']   # 1.00 perfect on M1 only
    Q = pd.concat([P[P.massif == 'M1'].assign(massif='M1#0'), P[P.massif == 'M1'].assign(massif='M1#1'),
                   P[P.massif == 'M2'].assign(massif='M2#2')], ignore_index=True)
    _, chosen = fft.lomo(Q, {fr: fft.col('A', fr) for fr in (0.10, 0.55, 1.00)})
    assert chosen['M1#0'] == 0.55                    # its twin M1#1 must not teach it 1.00


def test_robust_changes_need_all_tolerances():
    obs = {'e1': {t: 'polythermal' for t in fft.TOLS}, 'e2': {t: 'cold' for t in fft.TOLS}}
    a = {'e1': {t: 'polythermal' for t in fft.TOLS}, 'e2': {t: 'cold' for t in fft.TOLS}}
    b = {'e1': {t: 'cold' for t in fft.TOLS},                                   # robust loss
         'e2': {'0.05': 'temperate', 'hole': 'cold', '0.5': 'cold'}}            # not robust
    ch = fft.robust_changes(obs, a, b, ['e1', 'e2'])
    assert ch['polythermal'] == [0, 1] and ch['cold'] == [0, 0]


def test_decision_conditions():
    H = dict(M1=2.0, zone_M1=2.0)
    none = {c: [0, 0] for c in fft.CLASSES}
    assert fft.decide_A(H, dict(M1=2.29, zone_M1=2.0), none, none, [-0.1, 0.1])['adopt_FRH']
    assert not fft.decide_A(H, dict(M1=2.31, zone_M1=2.0), none, none, None)['adopt_FRH']
    two_poly = dict(none, polythermal=[0, 2])
    assert not fft.decide_A(H, dict(M1=2.0, zone_M1=2.0), two_poly, none, None)['adopt_FRH']
    zone_warm = dict(none, temperate=[0, 1])
    assert not fft.decide_A(H, dict(M1=2.0, zone_M1=2.0), none, zone_warm, None)['adopt_FRH']
    assert not fft.decide_A(H, dict(M1=2.0, zone_M1=2.0), none, none, [0.01, 0.2])['adopt_FRH']
