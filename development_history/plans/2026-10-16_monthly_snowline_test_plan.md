# The monthly snow-line test, and the fixes pushed to main and GloGEMflow

*Plan of 2026-10-16, evening, for your sign-off. Terms are in the glossary at the end.*

## Resume here (written 2026-10-16, when you left for the week)

**Where it stands:**
- **Not yet approved.** Your last comment set the daily switch aside ("go ahead in our todos").
  The plan below was updated for that and waits for your sign-off.
- **Nothing in it has been executed yet:** no pushes, no runs, no commits, no memory updates. Plan
  mode allowed only this file.
- **Next week:**
  - **If this conversation is resumed,** everything is here.
  - **In a new session,** memory points here (`project_thermal_forcing_plan.md` names this file).
  - Either way: approve, or comment, and I start with step 0 below.

**Step 0 after sign-off (before the "Order of work"):**
1. **Save what lives only in this session's `/tmp` scratch folder.** `/tmp` is a RAM disk, so a
   reboot wipes it.
   - The step-B check `stepB/compare_stepB.py` exists nowhere else; Appendix B below is a full copy.
   - Commit it as `scripts/icetemp_analysis/snowline/compare_snow.py`, and copy its outputs
     (`stepB/*.csv`, `README.md`, `verdict.txt`) to `…/week2/snowcal/stepB/`.
   - If the folder is gone, the step-B numbers are still in chapter 9 (`c56f030`), and the script
     regenerates the tables.
2. **The memory updates this plan could not make:**
   - `project_thermal_forcing_plan.md`: a new RESUME note pointing to this plan;
   - `project_snowline_port.md`: the step-B verdict, the monthly test, the daily-versus-monthly
     result;
   - a new memory for the daily-versus-monthly comparison and your decision;
   - a new memory for the 360-day thermal year;
   - `MEMORY.md`: its "Ice-temp plan (READ FIRST)" line still says RESUME 2026-10-15 and names
     the old handover note.

**Git state at the time of writing** (all clean unless noted):

| repository / worktree | branch | commit | note |
|---|---|---|---|
| `GloGEM` (main tree) | GloGEMflow, local | `d075c98` | target fix on top; 33 commits ahead of origin, not pushed |
| `GloGEM_wt_snowline` | snowline-calibration | `c7a477c` | the port and its three fixes; local until the big merge |
| `GloGEM_wt_fixtarget` | fix/geodetic-target-period | `a1fb2fa` | pushed (the remote branch exists); to go into main by fast-forward |
| `GloGEM_wt_fixflow` | fix/geodetic-target-period-glogemflow | `5dc500b` | local, no upstream; to be pushed onto origin/GloGEMflow |
| `GloGEM_wt_gate3`, `GloGEM_wt_mainref` | detached | `ce3fc75`, `fb7f9c0` | gate-3 runs; main reference |
| origin/main, origin/GloGEMflow | – | `fb7f9c0`, `dac9fa4` | – |
| `_wt_v3` (AMOC session) | glogemflow-v3 | `550c701` | holds the model-changes page `docs/source/development/model-changes-2026.md` (G31, G32); message that session before editing |
| `glogemflow_icetemp` | master | `c56f030` | 90 ahead of origin, not pushed |

**Uncommitted in `glogemflow_icetemp`:**
- your yaml (line 38) and the deleted notebook pdf: yours, left alone;
- four ch09 SVGs that differ only in date stamps and random ids (regenerated at 16:46): restore
  them with `git checkout`;
- `firn_model/`, untracked since 14 Oct: not touched.

**Runs on disk** (under `/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week2/snowcal/`):
- `boot/`: daily, `c0.4`–`c3.0` calibrations and `f0.4`–`f3.0` forward runs; `boot_result/`;
- `fwdB_none_fixed`, `fwdB_boot7`: monthly forward runs of step B;
- `refit_*`: monthly re-fits;
- `fixtest/`: the reference for section 3's tests;
- `gate3/`, `data_glamos3/`.

## Context

- **Your question** was whether the thermal model should run on the daily GloGEM, which reproduces
  the snow lines better.
- **Your decision**, after the numbers: not now; go ahead with the todos.
- **This plan therefore:**
  - records the daily-versus-monthly comparison (section 1);
  - runs the monthly snow-line test, adjusted to what the comparison showed (section 2);
  - pushes the two fixes to main and GloGEMflow, as you asked (section 3);
  - keeps the snowline branch local until the big merge (section 4).
- **Approving the plan accepts the recommendations** in "Decisions"; comment on any you want
  differently.

| # | your answer | in this plan |
|---|---|---|
| – | daily GloGEM for the thermal model: not now | section 1: a short record for the book |
| 1 | the monthly snow-line test | section 2, adjusted after section 1 |
| 2 | the target fix directly on main, pushed, documented | section 3 |
| 3 | push the fix to GloGEMflow | section 3 |
| 4 | the refreezing fix directly, too | section 3: on main and GloGEMflow |
| 5 | the snowline calibration with the big merge later | section 4: the branch stays local |

## 1. Daily GloGEM for the thermal model: not now (a record for the book)

**What the comparison showed** (existing runs; the 12 glaciers with snow lines, 2000–2019, 4015
readings):
- **May to August:** the monthly model is as good as the daily one.
- **September:** the daily model is better at every c_prec (bias +60 to +104 m, against +238 m
  monthly).
  - The monthly model decides snow or rain from the month's mean temperature, so the cold days'
    snow becomes rain.
  - Snow that falls and melts within one month never shows in a month-end value.
- **The summer's highest snow line**, which decides where firn survives, is the same in both
  models (0 m difference on average).
- **A daily thermal model** would take 10–30× longer per run, needs large work in GloGEM and the
  pipeline, and needs every calibration redone.
- **Kept for later:**
  - a daily mass balance feeding a monthly thermal model;
  - snowfall spread over the month's daily temperatures, as the melt already is.

**For the book:**
- chapter 9 gets a short section with the by-month table and a figure (the snow line by month, both
  models against the images);
- the comparison script is committed as `snowline/daily_vs_monthly.py`.
- This section is also the reason for the adjustment in section 2.

**A side finding, fixed later (decision 6):**
- The monthly thermal model gives every month 30 days (`rf_dt` in `settings.pro`), so its year has
  360 days.
- Heat conduction, ice flow and geothermal heat run 1.4 % short each year, while the meltwater heat
  is complete. That acts like raising f_ref by 1.4 %, far less than one step of the f_ref scan.
- **The fix:** the step length from the month's real length, for the thermal model only (the full
  refreezing model shares `rf_dt`).
  - It is made together with the thermal reruns (C), so that it lands in one rerun set.
  - It is then recorded on the model-changes page.
- The finding goes into memory now.

## 2. The monthly snow-line test, adjusted

**What it is.**
- **Today:** we chose c_prec with the daily model, Alexandra's way, and carried it into the monthly
  model.
- **Step B showed** that this does not work: the monthly snow line behaves differently.
- **The test:** choose c_prec with the monthly model itself, against the same satellite snow lines.

```
 today (the daily route)                          the monthly test
   daily GloGEM, scan c_prec                        monthly GloGEM, c_prec = 0.4, 0.6, ... 3.0
     daily snow line vs satellite                     each: DDF and T_off fitted to the target
     -> best c_prec                                   monthly snow line vs satellite
   monthly GloGEM with that c_prec                    -> best c_prec, bootstrap over years
     -> monthly snow line too high in September       -> the c_prec the thermal model's snow asks for
```

**The adjustment from section 1:**
- **It uses only the readings a monthly model can follow:** the snow line's rise through the summer,
  up to each year's highest reading (per glacier).
- **The September readings after that point are left out:**
  - no monthly c_prec can fix them (in step B, c_prec moved the September bias by only 17 m);
  - left in, they would dominate the fit.
- **This is also the part that matters for the thermal model:** the rise, and how high it gets,
  decide where firn survives.

```
 snow-line
 elevation        ● highest reading of the year
    ▲          ●     ○
    │       ●           ○      ● used: May up to the highest reading
    │    ●                     ○ left out: the readings after it
    │ ●
    └──────────────────────────▶ May ... September
```

**How a monthly model meets a reading on one day.**
- The model gives the snow line at the end of each month.
- A reading in between is compared with the line joining the two month ends (interpolated).
- Using only the month-end value would shift the comparison by about 80 m (step B).

**What it can show.**
- **If it works:** we get a c_prec that fits the model the thermal runs use, checked against the
  snow lines and against GLAMOS.
- **If it does not** (the curves are flat, or ask for the end of the range):
  - that is a finding: c_prec cannot fix the monthly snow cycle;
  - it becomes a stated limitation, and the thermal reruns (C) go ahead on the standard mass
    balance.

**Rules, fixed before the runs:**
- **Kept values:** a glacier keeps its c_prec only if its bootstrap 90 % range reaches neither 0.4
  nor 3.0.
- **The T_off test**, with this morning's flag rule.
- **The adoption rule**, which you decide on after M4:
  - (a) on the kept glaciers, the monthly snow line fits better (weighted RMSE on the readings
    above, interpolated);
  - (b) the winter balance does not move away from GLAMOS (Aletsch; Morteratsch only as a partial
    outline, flagged);
  - (c) C1 breaks none of A's conditions (later).

**Steps (no change to GloGEM's code):**
- **M1 runs.**
  - Fourteen monthly calibrations of the 12 glaciers (2000–2019 target), each with every glacier's
    c_prec fixed at one value (`calibrate_cprec_file`).
  - Each is followed by a monthly forward run, 1999–2025, with `snowline_definition = 'p10'` and
    the output of the step-B forward runs.
- **M2 bootstrap.**
  - `cprec_bootstrap.py` gets a monthly option: the snow line read from `Snowline_mon` and
    interpolated to each reading's date.
  - Readings: Sommer's raw snow lines, 2000–2025, May up to each year's highest reading.
  - 1000 draws over the years.
- **M3:** the monthly re-fit of the 23 glaciers with the kept values (as in A2), and the T_off test.
- **M4:** the snow check. The step-B script is moved into the repository as
  `snowline/compare_snow.py`, so that the rule is committed. Then your decision.
- **Time:** about 1.5 hours of runs on the busy host, four at a time as today.
- **Effort:**
  - runs as background scripts;
  - the bootstrap and the check in the main session, since they decide the recipe;
  - a reviewer only before the thermal run set, as agreed.

## 3. The two fixes, directly on main and GloGEMflow

**The two commits**, both already documented in detail:
- **The target period:** `a1fb2fa`, *fix(calibration): average the geodetic target over all its
  hydrological years, after one spin-up year* (the bug, OGGM's convention, the spin-up evidence, the
  effect).
- **Full refreezing:** `07d694d`, *fix(refreezing): initialise the full refreezing model's state for
  every glacier* (the symptom, the cause in the modular split, why working runs are unchanged).

**On main** (`origin/main` is still `fb7f9c0`, the base of `a1fb2fa`):
1. In a worktree on `main`, fast-forward `main` to `a1fb2fa`.
2. Cherry-pick the refreezing fix with `-x`, which records its origin.
3. Test (below), then `git push origin main`.
4. Delete the remote branch `fix/geodetic-target-period`, whose commit is then on main.

**On the shared GloGEMflow** (`origin/GloGEMflow` is `dac9fa4`):
1. The prepared local branch already holds the target fix on `dac9fa4` (`5dc500b`). The refreezing
   fix is added with `-x`.
2. Test, then push it onto `origin/GloGEMflow`, a fast-forward of two commits.
3. Our 32 other local GloGEMflow commits stay local.

**On our local GloGEMflow:**
- **It already has the target fix** (`d075c98`); the refreezing fix is cherry-picked as well.
- **`snowline-calibration` is rebased on it**; its own copy of the refreezing fix drops out.
- **The same change then has two hashes**, one here and one on the remote. They merge without
  conflict at the big merge.

**Tests before each push** (the repository's own test calibration of Aletsch and Morteratsch,
monthly and daily):
- **Default refreezing unchanged:** with the refreezing fix added, the results equal those with the
  target fix alone (the `fixtest/fixtarget_*` runs of this afternoon).
- **Full refreezing fixed:** a copy with `refreezing_full = 'y'` stops with `TL_RF` undefined on
  `fb7f9c0` and runs after the fix.
- **On GloGEMflow:** the same two checks on the GloGEMflow branch.

**Model-changes page (G31, G32):**
- Their status becomes "on main and GloGEMflow".
- The AMOC session has untested work in `_wt_v3` and asked to be messaged before that file is
  edited, so I send it the two status lines first.

## 4. The snowline calibration (your point 5)

- **It goes into main with the big merge of GloGEMflow**, which then contains the snowline
  calibration.
- **Until then the branch stays local**, rebased on GloGEMflow whenever that moves.
- **Gate 3's two open cases** (Abramov, Golubin) are tested on her own data folder before that
  merge.

## Order of work

1. **M1 runs started in the background** (about 1.5 hours, four at a time).
2. **While they run:**
   - section 3: tests, the pushes to main and GloGEMflow, the remote branch deleted;
   - section 4: the rebase;
   - the AMOC session messaged, then the G31/G32 status lines.
3. **The book and memory:**
   - the chapter 9 section from section 1, with its figure and `daily_vs_monthly.py`;
   - memory: the comparison, the decision, the 360-day finding;
   - commits as each piece lands.
4. **When M1 is done:** M2 → M3 → M4, then your decision on adoption.

## Decisions for you

| # | question | recommendation |
|---|---|---|
| 1 | the thermal model stays on the monthly GloGEM | decided by you: yes |
| 2 | the monthly test fits the snow line's rise, up to each year's highest reading | yes |
| 3 | the monthly test compares the interpolated monthly snow line (not the month-end value) | yes |
| 4 | snow lines of 2000–2025 (all of Sommer's years) | yes |
| 5 | delete the remote branch `fix/geodetic-target-period` after the direct push to main | yes |
| 6 | the 360-day thermal year is fixed together with the thermal reruns (C), not now | yes |

## Verification

- **Section 1:** `daily_vs_monthly.py` reproduces the numbers above.
- **Section 2:**
  - every fixed-value run closes its target;
  - where the minimum lies inside the range, the bootstrap's value agrees with a direct look at the
    curves;
  - M3–M4 as A2 and B today, with the rules above.
- **Section 3:**
  - the tests above, before each push;
  - after the pushes, `git ls-remote` shows `main` and `GloGEMflow` at the new commits;
  - gate 1 of the snowline branch is still byte-identical against the updated local GloGEMflow.

## Glossary

| term | meaning |
|---|---|
| c_prec | precipitation factor: multiplies the reanalysis precipitation |
| DDF, T_off | degree-day (melt) factor; temperature offset, fitted when DDF reaches its bounds |
| daily / monthly GloGEM | the same model run in daily or monthly steps; the thermal model runs on the monthly one |
| snow line, p10 | the elevation below which 10 % of the snow-covered glacier area lies (Sommer's statistic, now also the model's) |
| highest snow line | the highest snow line of a summer, usually late August or September; close to the firn line |
| interpolated | the model's snow line at a reading's date, read between the two month-end values |
| bias, weighted RMSE | mean of model minus image; the root-mean-square difference, weighted by 1/σ |
| σ (sigma) | the stated uncertainty of one satellite snow-line reading |
| bootstrap | resampling the observation years to see how firmly the data fix a value |
| f_ref | the thermal model's refreezing factor (scales the meltwater heat) |
| `rf_dt` | the length of one thermal solver step (10 days in monthly runs: a 30-day month in three steps) |
| A1, A2, B, C, M1–M4 | steps of the plans of 16 Oct: A1 the daily bootstrap, A2 its re-fit, B the snow check, C the thermal reruns, M1–M4 the monthly test |
| A's conditions | the adoption conditions of rule A: M1 on held-out massifs, no robust loss of temperature classes, warming error not worse |
| GLAMOS | Glacier Monitoring Switzerland (measured mass balances) |
| AMOC session | the other Claude session, working on the AMOC study in the worktree `_wt_v3` |
| fast-forward | moving a branch forward to a later commit that already contains it, without a merge |
| cherry-pick (-x) | copying a commit onto another branch; -x writes the original hash into the message |
| `main`, GloGEMflow | the shared main branch, and the shared GloGEMflow branch, of the GloGEM repository |
| RAM disk (`/tmp`) | a folder held in memory: fast, but emptied at every reboot; this session's scratch folder lives there |

---

## Appendix A: `daily_vs_monthly.py` (draft, for section 1)

- **To be committed** as `scripts/icetemp_analysis/snowline/daily_vs_monthly.py`.
- **Self-contained:** it reads the satellite CSVs and the runs on `/scratch_net`, not the `/tmp`
  scratch folder.
- **Checked on 16 Oct** (run read-only, printing only): it reproduces section 1's numbers.
  - September bias +238 against +60 to +104 m;
  - 0 m difference in the highest snow line;
  - 10 of 12 glaciers inside the daily range.
  - One value differs by 1 m from rounding.

```python
"""Daily versus monthly GloGEM snow lines on the same satellite readings (book chapter 9; plan of 16 Oct, evening).

Monthly: fwdB_none_fixed (standard monthly calibration, hydrological years 1999-2025), Snowline_mon interpolated
between month ends to each reading's date. Daily: boot/f<c> (c_prec fixed at each scan value 0.4-3.0, melt fitted
to the same geodetic target, 1999-2019), the p10 snow line of the reading's day.
Readings: Sommer TSLA sla_raw_m, qa >= 0.666, May-September, hydrological years 2000-2019, weights 1/sigma.
Prints (1) the weighted bias and RMSE by month, (2) the summer's highest snow line (July-September).

Usage: python daily_vs_monthly.py [SNOWCAL_DIR] [TSLA_DIR]
"""
import calendar
import glob
import os
import sys

import numpy as np
import pandas as pd

S = sys.argv[1] if len(sys.argv) > 1 else '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week2/snowcal'
T = sys.argv[2] if len(sys.argv) > 2 else \
    '/home/jabeer/projects/glogemflow_development/glogemflow_icetemp/data/snowlines/tsla_csv'
TSLA = ['00773', '00774', '01112', '01127', '01205', '01225', '01230', '02216', '02596', '02611', '02821', '03116']
Y0_MON, NY_MON = 1999, 27   # monthly run: hydrological years 1999-2025
Y0_DAY = 1999               # daily runs: hydrological years 1999-2019, 1999 the spin-up year
YEARS = (2000, 2019)


def monthly(fn):
    out = {}
    for ln in open(fn).readlines()[1:]:
        p = ln.split()
        if p:
            out[p[0]] = np.array(p[1:], float).reshape(NY_MON, 12)
    return out


def month_end(arr, year, month):
    i = (year + 1 if month >= 10 else year) - Y0_MON
    return arr[i, (month - 10) % 12] if 0 <= i < NY_MON else np.nan


def interp(arr, t):
    prev = t.replace(day=1) - pd.Timedelta(days=1)
    a, b = month_end(arr, prev.year, prev.month), month_end(arr, t.year, t.month)
    return a + (b - a) * t.day / calendar.monthrange(t.year, t.month)[1]


def daily(fwd):
    m = {}
    for ln in open(f'{fwd}/daily/CentralEurope/PAST/PAST_original/'
                   'centraleurope_Snowline_day_r1_CentralEurope_snowcal.dat'):
        p = ln.split()
        if p[0] == 'ID' or len(p) < 300:
            continue
        m.setdefault(p[0], {})[int(p[1])] = np.array(p[3:], float)
    return m


def wb(mod, obs, w):
    ok = np.isfinite(mod)
    r, w = mod[ok] - obs[ok], w[ok]
    return np.sum(w * r) / w.sum(), np.sqrt(np.sum(w * r * r) / w.sum())


mon = monthly(f'{S}/fwdB_none_fixed/monthly/CentralEurope/PAST/PAST_original/'
              'centraleurope_Snowline_mon_r1_CentralEurope_fair.dat')
runs = sorted(glob.glob(f'{S}/boot/f*'), key=lambda p: float(os.path.basename(p)[1:]))
cs = [float(os.path.basename(p)[1:]) for p in runs]
day = [daily(p) for p in runs]

obs = []
for gid in TSLA:
    d = pd.read_csv(f'{T}/tsla_RGI70-11.{gid}.csv', parse_dates=['date'])
    d = d[(d.qa >= 0.666) & d.month.between(5, 9) & d.hydro_year.between(*YEARS)].copy()
    d['gid'] = gid
    obs.append(d)
o = pd.concat(obs, ignore_index=True)
o['w'] = 1 / o.sigma_m
o['mon'] = [interp(mon[g], t) for g, t in zip(o.gid, o.date)]
doy = o.date.dt.dayofyear.values.copy()
doy[o.date.dt.is_leap_year.values & (doy > 59)] -= 1          # the runs' 365-day calendar
k = 365 * (o.date.dt.year.values - Y0_DAY + 1) + doy - 274     # day index from 1 Oct of the first run year
o['col'] = k % 365
assert ((Y0_DAY + k // 365) == o.hydro_year).all()
D = []
for c, m in zip(cs, day):
    o[f'd{c:.1f}'] = [m[g][h][j] for g, h, j in zip(o.gid, o.hydro_year, o.col)]
    D.append(f'd{c:.1f}')

print(f'{o.gid.nunique()} glaciers, {len(o)} readings, hydrological years {YEARS[0]}-{YEARS[1]}')
print('month     n  bias_mon  bias_day(min..max)  wrmse_mon  wrmse_day(min..max)')
for mo, s in o.groupby('month'):
    bm, rm = wb(s.mon.values, s.sla_raw_m.values, s.w.values)
    bd, rd = zip(*[wb(s[c].values, s.sla_raw_m.values, s.w.values) for c in D])
    print(f'{mo:5d} {len(s):5d} {bm:8.0f}  {min(bd):7.0f}..{max(bd):5.0f}  {rm:9.0f}  {min(rd):7.0f}..{max(rd):5.0f}')

# the summer's highest snow line, July-September
rows = []
for g in TSLA:
    yrs = range(YEARS[0], YEARS[1] + 1)
    m_max = np.mean([mon[g][y - Y0_MON, 9:12].max() for y in yrs])
    d_max = [np.mean([m[g][y][273:365].max() for y in yrs]) for m in day]
    s = o[(o.gid == g) & (o.month >= 7)]
    n = s.groupby('hydro_year').size()
    grp = s[s.hydro_year.isin(n[n >= 3].index)].groupby('hydro_year')
    rows.append(dict(gid=g, mon_max=m_max, day_lo=min(d_max), day_hi=max(d_max), n_years=grp.ngroups,
                     obs_max=grp.sla_raw_m.max().mean(), mon_on_obs=grp.mon.max().mean(),
                     day_on_obs_lo=min(grp[c].max().mean() for c in D),
                     day_on_obs_hi=max(grp[c].max().mean() for c in D)))
t = pd.DataFrame(rows).set_index('gid')
print(t.round(0).to_string())
inside = int(((t.mon_max >= t.day_lo) & (t.mon_max <= t.day_hi)).sum())
print(f'monthly minus the middle of the daily range: {(t.mon_max - (t.day_lo + t.day_hi) / 2).mean():.0f} m; '
      f'monthly inside the daily range on {inside} of {len(t)} glaciers')
```

## Appendix B: copy of `stepB/compare_stepB.py` (the step-B check)

- **Verbatim** from `/tmp/claude-553859/…/scratchpad/stepB/compare_stepB.py`, the only copy.
- **It writes its tables next to itself** (`OUT`).
- **For M4** it becomes `snowline/compare_snow.py`:
  - the runs as arguments;
  - the readings limited to the rise up to each year's highest reading (section 2).

```python
"""Step B: BEFORE (fwdB_none_fixed) vs AFTER (fwdB_boot7) snow line, winter balance, ELA."""
import calendar
import glob
import os

import numpy as np
import pandas as pd

OUT = os.path.dirname(os.path.abspath(__file__))
RUNS = {'before': '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week2/snowcal/fwdB_none_fixed',
        'after': '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week2/snowcal/fwdB_boot7'}
SUB = 'monthly/CentralEurope/PAST/PAST_original'
DATA = '/home/jabeer/projects/glogemflow_development/glogemflow_icetemp/data/snowlines'
FIX = '/itet-stor/jabeer/glazioarch/GlacioBaseData/MassBalance/glacier'
MH = '/scratch_net/vierzack06_fourth/mhuss/MassBalance_all'
Y0, NY = 1999, 27
CHANGED = ['00773', '01205', '01225', '02216', '02596', '02611', '03116']
TSLA = ['00773', '00774', '01112', '01127', '01205', '01225', '01230', '02216', '02596', '02611', '02821', '03116']
MEAS_W = (1, 2, 4, 5)  # GLAMOS methods with a measured winter balance

tg = pd.read_csv(f'{DATA}/targets.csv')
NAME = {r.rgi7[-5:]: r['name'] for _, r in tg.iterrows()}
MATCH = pd.read_csv(f'{DATA}/reference/glamos_vs_geodetic.csv', dtype={'rgi7': str})


def table(run, var):
    fn = f'{RUNS[run]}/{SUB}/centraleurope_{var}_r1_CentralEurope_fair.dat'
    lines = open(fn).readlines()
    out = {}
    for ln in lines[1:]:
        p = ln.split()
        if p:
            out[p[0]] = np.array(p[1:], float)
    return out


def calib(run):
    fn = glob.glob(f'{RUNS[run]}/monthly/CentralEurope/calibration/*.dat')[0]
    return pd.read_csv(fn, sep=r'\s+', dtype={'ID': str}).set_index('ID')


def model_month(arr, year, month):
    Y = year + 1 if month >= 10 else year
    i = Y - Y0
    return arr[i, (month - 10) % 12] if 0 <= i < NY else np.nan


def model_interp(arr, t):
    nd = calendar.monthrange(t.year, t.month)[1]
    prev = t.replace(day=1) - pd.Timedelta(days=1)
    a, b = model_month(arr, prev.year, prev.month), model_month(arr, t.year, t.month)
    return a + (b - a) * t.day / nd


def wstats(d, mcol):
    ok = np.isfinite(d[mcol]) & np.isfinite(d.sla_raw_m)
    w, e = 1 / d.sigma_m[ok], (d[mcol] - d.sla_raw_m)[ok]
    if ok.sum() == 0:
        return np.nan, np.nan, 0
    return np.sqrt(np.sum(w * e * e) / w.sum()), np.sum(w * e) / w.sum(), int(ok.sum())


# ---------------- 1. snow line ----------------
sl = {r: {k: v.reshape(NY, 12) for k, v in table(r, 'Snowline_mon').items()} for r in RUNS}
obs = []
for gid in TSLA:
    d = pd.read_csv(f'{DATA}/tsla_csv/tsla_RGI70-11.{gid}.csv', parse_dates=['date'])
    d = d[(d.qa >= 0.666) & d.month.between(5, 9) & d.hydro_year.between(Y0, Y0 + NY - 1)].copy()
    d['gid'] = gid
    for r in RUNS:
        d[f'{r}_mon'] = [model_month(sl[r][gid], t.year, t.month) for t in d.date]
        d[f'{r}_int'] = [model_interp(sl[r][gid], t) for t in d.date]
    obs.append(d)
obs = pd.concat(obs, ignore_index=True)
obs['group'] = np.where(obs.gid.isin(CHANGED), 'changed7', 'unchanged5')
obs.to_csv(f'{OUT}/snowline_matched_obs.csv', index=False, float_format='%.1f')


def stats_block(d):
    rec = {}
    for var in ('mon', 'int'):
        for r in RUNS:
            rm, b, n = wstats(d, f'{r}_{var}')
            rec[f'{r}_wrmse_{var}'] = rm
            rec[f'{r}_wbias_{var}'] = b
        rec[f'n'] = n
        rec[f'd_wrmse_{var}'] = rec[f'after_wrmse_{var}'] - rec[f'before_wrmse_{var}']
    return rec


rows = []
for gid, d in obs.groupby('gid'):
    rows.append(dict(rgi7=gid, name=NAME[gid], changed=gid in CHANGED,
                     yrs=f'{d.hydro_year.min()}-{d.hydro_year.max()}', **stats_block(d)))
for grp, d in obs.groupby('group'):
    rows.append(dict(rgi7=f'POOLED_{grp}', name='pooled (obs-weighted 1/sigma)', changed=grp == 'changed7',
                     yrs=f'{d.hydro_year.min()}-{d.hydro_year.max()}', **stats_block(d)))
SLG = pd.DataFrame(rows)
# glacier-mean of per-glacier scores (each glacier weighted equally) as a sensitivity
for grp, ch in (('changed7', True), ('unchanged5', False)):
    s = SLG[(SLG.changed == ch) & ~SLG.rgi7.str.startswith('POOLED')]
    rec = dict(rgi7=f'GLACIERMEAN_{grp}', name='mean of per-glacier scores', changed=ch, yrs='')
    for c in s.columns:
        if c.startswith(('before_', 'after_', 'd_')):
            rec[c] = s[c].mean()
    rec['n'] = s.n.sum()
    rows.append(rec)
SLG = pd.DataFrame(rows)
SLG.to_csv(f'{OUT}/snowline_per_glacier.csv', index=False, float_format='%.1f')

rows = []
for grp, dg in obs.groupby('group'):
    for m, d in dg.groupby('month'):
        rows.append(dict(group=grp, month=m, **stats_block(d)))
SLM = pd.DataFrame(rows)
SLM.to_csv(f'{OUT}/snowline_by_month.csv', index=False, float_format='%.1f')

# per glacier and month (changed only)
rows = []
for (gid, m), d in obs[obs.gid.isin(CHANGED)].groupby(['gid', 'month']):
    rows.append(dict(rgi7=gid, name=NAME[gid], month=m, **stats_block(d)))
pd.DataFrame(rows).to_csv(f'{OUT}/snowline_by_glacier_month_changed7.csv', index=False, float_format='%.1f')

# sensitivity: drop inferred full-snow / snow-free observations
rows = []
for grp, d in obs[obs.flag.isna()].groupby('group'):
    rows.append(dict(group=grp, variant='no_inferred_flags', **stats_block(d)))
pd.DataFrame(rows).to_csv(f'{OUT}/snowline_pooled_noflag.csv', index=False, float_format='%.1f')


# ---------------- 2/3. GLAMOS winter balance, ELA ----------------
def read_fix(name):
    lines = open(f'{FIX}/{name}_fix.dat', encoding='latin-1').readlines()
    h = [s.strip() for s in lines[0].split(';')]
    n, zlo, zhi = int(h[4]), float(h[5]), float(h[6])
    edges = np.linspace(zlo, zhi, n + 1)
    recs = {}
    for l in lines:
        if l.startswith('#') or not l.strip():
            continue
        p = l.split()
        f = lambda s: np.nan if s.lower() == 'nan' else float(s.rstrip('.'))
        b = np.array([f(x) for x in p[12:]])
        bw, nw, ba, na, A = (b[i * n:(i + 1) * n] for i in range(5))
        recs[int(p[4][:4])] = dict(meth=int(p[0]), Bw=f(p[5]) / 1000, Ba=f(p[6]) / 1000, ELA=f(p[7]),
                                   A=f(p[9]), bw=bw / 1000, nw=nw, A_eb=A)
    return edges, recs


def read_mh(names):
    """mhuss bands_<yr>.dat summed over several glaciers: lo, hi, area, fxd_ann, meas_ann, fxd_win, meas_win."""
    sets = []
    for nm in names:
        d = f'{MH}/{nm}/massbalance/results'
        s = {}
        for fn in sorted(glob.glob(f'{d}/bands_*.dat')):
            rows = [[float(x) for x in l.split()[:7]] for l in open(fn).readlines()[1:] if len(l.split()) >= 7]
            s[int(os.path.basename(fn)[6:10])] = np.unique(np.array(rows), axis=0)
        sets.append(s)
    lo = min(a[:, 0].min() for s in sets for a in s.values())
    hi = max(a[:, 1].max() for s in sets for a in s.values())
    a0 = next(iter(sets[0].values()))
    w = a0[0, 1] - a0[0, 0]
    edges = np.arange(lo, hi + w / 2, w)
    n = len(edges) - 1
    recs = {}
    for yr in sorted(set.intersection(*[set(s) for s in sets])):
        A = np.zeros(n); W = np.zeros(n)
        for s in sets:
            for r in s[yr]:
                i = int(round((r[0] - lo) / w))
                if 0 <= i < n and r[2] > 0:
                    A[i] += r[2]; W[i] += r[2] * r[5]
        with np.errstate(all='ignore'):
            bw = np.where(A > 0, W / A, np.nan) / 1000
        recs[yr] = dict(meth=-1, Bw=np.nansum(bw * A) / A.sum(), ELA=np.nan, A=A.sum(), bw=bw, nw=np.zeros(n), A_eb=A)
    return edges, recs


def mband(run, rid):
    rows = [l.split() for l in open(f'{RUNS[run]}/{SUB}/mb_elevation/bwelev_{rid}.dat').readlines()[1:] if l.strip()]
    z = np.array([float(r[0]) for r in rows])
    v = np.array([[float(x) for x in r[1:]] for r in rows])
    v[v <= -98] = np.nan
    return z, v


def on_bins(z, v, edges):
    out = np.full((len(edges) - 1, v.shape[1]), np.nan)
    for i in range(len(edges) - 1):
        m = (z >= edges[i]) & (z < edges[i + 1])
        if m.any():
            with np.errstate(all='ignore'):
                out[i] = np.nanmean(v[m], axis=0)
    return out


WB = {r: table(r, 'Winter_balance_sfc') for r in RUNS}
EL = {r: table(r, 'ELA') for r in RUNS}
AR = {r: table(r, 'Area') for r in RUNS}

# rgi7, label, source, outline ratio, flag, counts for condition (b)
ratio = {r.short: r.area_ratio for _, r in MATCH.iterrows()}
SRC = [('02596', 'fix:aletsch', True), ('02216', 'fix:pers', True),
       ('01112', 'fix:hohsaas', False), ('01309', 'fix:sexrouge', False), ('01039', 'fix:tortin', False),
       ('01150', 'fix:alphubel', False), ('02103', 'fix:corvatsch', False),
       ('02216', 'mh:morteratsch+pers', False), ('02611', 'mh:fiescher', False), ('01225', 'mh:gorner', False)]

bw_rows, ela_rows, band_rows, band_prof = [], [], [], []
for rid, src, cond_b in SRC:
    kind, nm = src.split(':')
    if kind == 'fix':
        edges, recs = read_fix(nm)
        rat = ratio.get(nm, np.nan)
    else:
        edges, recs = read_mh(nm.split('+'))
        rat = np.nan
    yrs_all = [y for y in recs if 2000 <= y <= 2025 and np.isfinite(recs[y]['Bw'])]
    if kind == 'fix':
        yrs_meas = [y for y in yrs_all if recs[y]['meth'] in MEAS_W]
    elif nm == 'morteratsch+pers':
        yrs_meas = [y for y in yrs_all if 2021 <= y <= 2025]  # winter stakes on both (peg_*_w.dat)
    else:
        yrs_meas = []  # no winter stakes in 2000-2025 (fiescher: first winter survey 2026; gorner: none)
    A_src = np.nanmean([recs[y]['A'] for y in yrs_all]) if yrs_all else np.nan
    if kind == 'mh':
        rat = A_src / AR['before'][rid][2003 - Y0 + 0] if rid in AR['before'] else np.nan
    flag = '' if (np.isfinite(rat) and 0.85 <= rat <= 1.15) else f'outline ratio {rat:.2f}'
    if kind == 'mh':
        flag = ('mhuss evaluation, not GLAMOS-published; ' + ('winter stakes 2021-2025' if yrs_meas else
                'NO winter stakes 2000-2025 -> winter is model-inferred, info only') + '; ' + flag)
    for sel, yrs in (('measured', yrs_meas), ('all', yrs_all)):
        if not yrs:
            continue
        o = np.mean([recs[y]['Bw'] for y in yrs])
        rec = dict(rgi7=rid, name=NAME[rid], source=src, changed=rid in CHANGED, years_sel=sel, n_years=len(yrs),
                   yrs=f'{min(yrs)}-{max(yrs)}', area_obs_km2=A_src, outline_ratio=rat, flag=flag.strip('; '),
                   Bw_obs=o)
        for r in RUNS:
            m = np.mean([WB[r][rid][y - Y0] for y in yrs])
            rec[f'Bw_{r}'] = m
            rec[f'absdiff_{r}'] = abs(m - o)
        rec['cond_b_pass'] = rec['absdiff_after'] <= rec['absdiff_before'] + 1e-9
        rec['counts_for_cond_b'] = cond_b and sel == 'measured' and rid in CHANGED
        bw_rows.append(rec)
        # ELA
        ey = [y for y in yrs if np.isfinite(recs[y]['ELA'])]
        if ey:
            o = np.mean([recs[y]['ELA'] for y in ey])
            rec2 = dict(rgi7=rid, name=NAME[rid], source=src, changed=rid in CHANGED, years_sel=sel, n_years=len(ey),
                        yrs=f'{min(ey)}-{max(ey)}', flag=flag.strip('; '), ELA_obs=o)
            for r in RUNS:
                rec2[f'ELA_{r}'] = np.mean([EL[r][rid][y - Y0] for y in ey])
                rec2[f'bias_{r}'] = rec2[f'ELA_{r}'] - o
            ela_rows.append(rec2)
        # bands
        ob = np.array([recs[y]['bw'] for y in yrs]).T
        Aw = np.array([recs[y]['A_eb'] for y in yrs]).T
        nwt = np.array([recs[y]['nw'] for y in yrs]).T
        zc = 0.5 * (edges[1:] + edges[:-1])
        rec3 = dict(rgi7=rid, name=NAME[rid], source=src, changed=rid in CHANGED, years_sel=sel, n_years=len(yrs),
                    yrs=f'{min(yrs)}-{max(yrs)}', bin_m=edges[1] - edges[0], flag=flag.strip('; '))
        for r in RUNS:
            z, v = mband(r, rid)
            mb = on_bins(z, v, edges)[:, [y - Y0 for y in yrs]]
            ok = np.isfinite(ob) & np.isfinite(mb) & (Aw > 0)
            w = np.where(ok, Aw, 0)
            e = np.where(ok, mb - ob, 0)
            rec3[f'n_bandyears'] = int(ok.sum())
            rec3[f'bias_{r}'] = np.sum(w * e) / w.sum()
            rec3[f'rmse_{r}'] = np.sqrt(np.sum(w * e * e) / w.sum())
            okp = ok & (nwt > 0)
            if okp.sum():
                wp = np.where(okp, Aw, 0)
                rec3['n_bandyears_with_pts'] = int(okp.sum())
                rec3[f'bias_ptbins_{r}'] = np.sum(wp * e) / wp.sum()
                rec3[f'rmse_ptbins_{r}'] = np.sqrt(np.sum(wp * e * e) / wp.sum())
            with np.errstate(all='ignore'):
                mprof = np.nanmean(np.where(ok, mb, np.nan), axis=1)
                oprof = np.nanmean(np.where(ok, ob, np.nan), axis=1)
            for i in np.where(np.isfinite(oprof))[0]:
                band_prof.append(dict(rgi7=rid, source=src, years_sel=sel, run=r, z_lo=edges[i], z_hi=edges[i + 1],
                                      A_obs=np.nanmean(Aw[i]), n_winter_pts=int(nwt[i].sum()), n_years=int(ok[i].sum()),
                                      bw_obs=oprof[i], bw_model=mprof[i]))
        band_rows.append(rec3)

BW = pd.DataFrame(bw_rows); BW.to_csv(f'{OUT}/glamos_Bw_glacierwide.csv', index=False, float_format='%.3f')
ELA = pd.DataFrame(ela_rows); ELA.to_csv(f'{OUT}/glamos_ELA.csv', index=False, float_format='%.0f')
BB = pd.DataFrame(band_rows); BB.to_csv(f'{OUT}/glamos_Bw_bands_summary.csv', index=False, float_format='%.3f')
BP = pd.DataFrame(band_prof)
BP = BP.pivot_table(index=['rgi7', 'source', 'years_sel', 'z_lo', 'z_hi', 'A_obs', 'n_winter_pts', 'bw_obs'],
                    columns='run', values='bw_model').reset_index()
BP.to_csv(f'{OUT}/glamos_Bw_band_profiles.csv', index=False, float_format='%.3f')

# ---------------- 4. sanity: unchanged glaciers ----------------
rows = []
VARS = ['Snowline_mon', 'Winter_balance_sfc', 'Annual_Balance_sfc', 'ELA', 'Accumulation_sfc', 'Snowmelt_sfc',
        'Icemelt_sfc', 'Volume', 'Area', 'AAR', 'Precipitation_mon', 'Melt_mon']
for var in VARS:
    a, b = table('before', var), table('after', var)
    for gid in a:
        if gid in CHANGED:
            continue
        d = b[gid] - a[gid]
        i = int(np.nanargmax(np.abs(d))) if np.isfinite(d).any() else 0
        rows.append(dict(var=var, rgi7=gid, name=NAME.get(gid, ''), max_absdiff=np.nanmax(np.abs(d)),
                         mean_absdiff=np.nanmean(np.abs(d)), mean_before=np.nanmean(a[gid]), idx_max=i))
ca, cb = calib('before'), calib('after')
for gid in ca.index:
    if gid in CHANGED:
        continue
    for c in ('DDFsnow', 'Cprec', 'T_off', 'Ba', 'Bw'):
        rows.append(dict(var=f'calib_{c}', rgi7=gid, name=NAME.get(gid, ''), max_absdiff=abs(cb.loc[gid, c] - ca.loc[gid, c]),
                         mean_absdiff=abs(cb.loc[gid, c] - ca.loc[gid, c]), mean_before=ca.loc[gid, c], idx_max=0))
SAN = pd.DataFrame(rows); SAN.to_csv(f'{OUT}/sanity_unchanged16.csv', index=False, float_format='%.4g')
cal = ca[['DDFsnow', 'DDFice', 'Cprec', 'T_off', 'Bw', 'Ba']].join(cb[['DDFsnow', 'DDFice', 'Cprec', 'T_off', 'Bw', 'Ba']],
                                                                    lsuffix='_before', rsuffix='_after')
cal.insert(0, 'name', [NAME.get(g, '') for g in cal.index])
cal.insert(1, 'changed', [g in CHANGED for g in cal.index])
cal.to_csv(f'{OUT}/calibration_before_after.csv', float_format='%.3f')

# ---------------- 5. adoption rule ----------------
pc = SLG[SLG.rgi7 == 'POOLED_changed7'].iloc[0]
cond_a = pc.after_wrmse_mon < pc.before_wrmse_mon
# year-block bootstrap of the pooled wRMSE difference (resample hydrological years, all 7 glaciers together)
rng = np.random.default_rng(0)
d7 = obs[obs.group == 'changed7']
yrs = np.sort(d7.hydro_year.unique())
byy = {y: d7[d7.hydro_year == y] for y in yrs}
boot = []
for _ in range(2000):
    dd = pd.concat([byy[y] for y in rng.choice(yrs, len(yrs))])
    boot.append([wstats(dd, 'after_mon')[0] - wstats(dd, 'before_mon')[0], wstats(dd, 'after_int')[0] - wstats(dd, 'before_int')[0]])
boot = np.array(boot)
BOOT = pd.DataFrame(dict(variant=['model_month', 'model_interp'],
                         d_wrmse=[pc.d_wrmse_mon, pc.d_wrmse_int],
                         p05=np.percentile(boot, 5, axis=0), p50=np.percentile(boot, 50, axis=0),
                         p95=np.percentile(boot, 95, axis=0), frac_improved=(boot < 0).mean(axis=0)))
BOOT.to_csv(f'{OUT}/snowline_pooled_bootstrap_changed7.csv', index=False, float_format='%.2f')
print(BOOT.round(2))
cb_rows = BW[BW.counts_for_cond_b]
cond_b = bool(cb_rows.cond_b_pass.all())
with open(f'{OUT}/verdict.txt', 'w') as f:
    f.write(f'(a) pooled wRMSE changed7 May-Sep (model_month): before {pc.before_wrmse_mon:.1f} m, after {pc.after_wrmse_mon:.1f} m -> {cond_a}\n')
    f.write(f'    sensitivity model_interp: before {pc.before_wrmse_int:.1f}, after {pc.after_wrmse_int:.1f}\n')
    for _, r in cb_rows.iterrows():
        f.write(f'(b) {r.rgi7} {r["name"]} [{r.source}, {r.yrs}, n={r.n_years}] obs {r.Bw_obs:.3f} before {r.Bw_before:.3f} '
                f'after {r.Bw_after:.3f}: |d| {r.absdiff_before:.3f} -> {r.absdiff_after:.3f} pass={r.cond_b_pass} {r.flag}\n')
    f.write(f'(b) all -> {cond_b}\n')

pd.set_option('display.width', 250, 'display.max_columns', 40)
print(SLG.round(1).to_string(index=False))
print(SLM.round(1).to_string(index=False))
print(BW.round(3).to_string(index=False))
print(ELA.round(0).to_string(index=False))
print(BB.round(3).to_string(index=False))
print(SAN.sort_values('max_absdiff', ascending=False).groupby('var').head(2).round(4).to_string(index=False))
print(cal.round(3).to_string())
print(open(f'{OUT}/verdict.txt').read())
```
