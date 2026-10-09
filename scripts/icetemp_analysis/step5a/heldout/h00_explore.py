import sys, time, warnings
sys.path.insert(0, 'scripts/icetemp_analysis')
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
import target_analysis as ta
t0 = time.time()
obs = ta.load_obs()
s2 = ta.load_model('alps_s2'); r7 = ta.load_model('alps_r7')
common = ta.common_rows([s2, r7], obs)
print('load', time.time() - t0)
for name in ('T1', 'T2', 'T3r'):
    t = ta.target(name, obs, s2, common)
    print(name, t.y.shape, t.Y.shape, np.isfinite(t.Y).all())
t3 = ta.target('T3r', obs, s2, common)
print(t3.rows.columns.tolist())
print(t3.rows.groupby('e').size().describe())
t2 = ta.target('T2', obs, s2, common)
print(t2.rows.columns.tolist())
print(t2.rows.groupby('system_short').size())
print(t2.rows.groupby('e').size().describe())
d = ta.load_design('alps_s2')
print(d.describe())
print(ta.prior_sd(), ta.prior_sd(True))
print(s2.keys())
