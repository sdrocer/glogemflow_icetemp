import sys; sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2')
import numpy as np, pandas as pd
import target_analysis as ta
import s2lib as L
pd.set_option('display.width', 250, 'display.max_columns', 60, 'display.max_rows', 300)
obs = ta.load_obs(); mod = ta.load_model('alps_s2')
rows, Y, y = L.month_rows(obs, mod)
Rres = Y - y[None, :]
P = L.pairs(rows)
print('pairs', len(P), 'points', P.t2.nunique(), 'same borehole', P.same_bh.sum())
G = L.gamma(P, Rres)
keep = mod['design'].keep.to_numpy()
lagnames = [n for _, _, n in L.LAG_BINS]
def table(P, G, label, balanced=False):
    b = L.binned(P, G, ['zb', 'lagbin'], balanced)
    d13 = b[ta.BEST_RUN].unstack('zb').reindex(index=lagnames, columns=L.BANDS)
    med = b.median(axis=1).unstack('zb').reindex(index=lagnames, columns=L.BANDS)
    med29 = b.loc[:, np.flatnonzero(keep)].median(axis=1).unstack('zb').reindex(index=lagnames, columns=L.BANDS)
    n = P.groupby(['zb', 'lagbin']).size().unstack('zb').reindex(index=lagnames, columns=L.BANDS)
    npt = P.groupby(['zb', 'lagbin']).t2.nunique().unstack('zb').reindex(index=lagnames, columns=L.BANDS)
    out = pd.concat({'d13': d13.round(3), 'med32': med.round(3), 'med29': med29.round(3), 'pairs': n, 'points': npt}, axis=1)
    print(f'\n== {label} ==')
    print(out.to_string())
table(P, G, 'all pairs, pair-weighted')
table(P, G, 'all pairs, point-balanced', balanced=True)
s = P.same_bh.to_numpy()
table(P[s].reset_index(drop=True), G[:, s], 'same borehole, pair-weighted')
table(P[~s].reset_index(drop=True), G[:, ~s], 'different boreholes, pair-weighted')
# who supplies each bin
print('\npairs per glacier system and lag bin:')
print(P.groupby(['system_short', 'lagbin']).size().unstack('lagbin').reindex(columns=lagnames).fillna(0).astype(int).to_string())
print('\nobservation-error part 0.5(sig_i^2+sig_j^2) mean by band:', P.groupby('zb').var_obs.mean().round(4).to_dict())
P.to_pickle('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2/pairs.pkl')
