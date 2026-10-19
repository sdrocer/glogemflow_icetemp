"""The perfect-firn test of the batch of 19 Oct, site by site (reported beside the rule of firn_batch_rule.py).

At f_ref 0.85 (the value every massif but Monte Rosa and Mont Blanc chose), for the borehole sites of the 9 GLAMOS
glaciers: mean measured temperature, and the bias (model minus measured) and RMSE of FRH, FRHL (ledger fed GloGEM's
balance), FRHLo (fed the observed balance) and FRHLO (also the observed firn mask). Writes perfect_firn_sites.csv.

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/firn/perfect_firn_sites.py [out_dir]
"""
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import forcing_fair_test as ft
from icetemp.calibration import thermal_structure as ts
from icetemp.calibration.data import DataHandler

V = ('FRH', 'FRHL', 'FRHLo', 'FRHLO')
G9 = ['01112', '02103', '02596', '01150', '01039', '01309', '02611', '01225', '02216']


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3/batch_rule')
    dh = DataHandler(region='CentralEurope')
    dh.load()
    E = ft.scored_entities(dh.calibration_glaciers, ts.load_borehole_table())
    outs = {ft.col(v, 0.85): pickle.load(open(ft.BAYESCAL / ft.run_name('b19', v, 0.85) / 'training_runs' /
                                              'design_0000.output.pkl', 'rb')) for v in V}
    P = ft.points_frame(E, outs)
    P['site'] = P.entity.str.split('@').str[0]
    rows = []
    for s, d in P[P.rgi.isin(G9)].groupby('site'):
        r = dict(site=s, rgi=d.rgi.iloc[0], n=len(d), band_lo=d.band.min(), band_hi=d.band.max(), obs=d.obs.mean())
        for v in V:
            e = d[ft.col(v, 0.85)] - d.obs
            r[f'bias_{v}'], r[f'rmse_{v}'] = e.mean(), np.sqrt((e ** 2).mean())
        rows.append(r)
    T = pd.DataFrame(rows)
    T.to_csv(out / 'perfect_firn_sites.csv', index=False, float_format='%.3f')
    print(T.round(2).to_string(index=False))


if __name__ == '__main__':
    main()
