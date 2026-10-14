# Snow and firn from data, a probabilistic product, and the road to v1.0

*Plan of 2026-10-14 (afternoon). The plans of 2026-10-12/13 are archived in
`glogemflow_icetemp/development_history/plans/`; today's results are in book chapter 8 (icetemp
`42278cd`) and in `plans/2026-10-14_plan_additions.md` (your two additions, `a6bdbce`).*

## Context

Today's decisions, all by rules fixed before the runs:

| question | outcome |
|---|---|
| forcing (rule A) | offset-free FRH is central; H is carried as a full sensitivity (the result flips without Bernina) |
| meltwater depth limit ('front', rule B) | failed its gate twice, so 'hl' stays. The rule had two stable states: water stops at 4 m from a cold start, passes from a temperate start |
| meltwater refreezing in the winter snow on ice | fixes the false frozen tongues (Bossons 84 → 3 % frozen) and improves M1 and warming, but loses Tête Rousse's measured frozen bed, so it stays off |
| zone target D (3000–3500 m) | not met; the time box is extended (your rule) |

**The common cause found today: too much winter snow, and so too much firn, on small high
glaciers.**
- At Corvatsch and Hohsaas, GLAMOS shows a higher ELA (by 125–150 m), less winter snow (by a factor
  of 1.3–1.9) and negative band balances where the model has firn. Tête Rousse and Altels look the
  same.
- So the meltwater soaks through firn that in reality is ice, down to a close-off depth (53–80 m)
  deeper than the whole column.
- The same excess snow fills the thin ice columns, which is why the snow-refreezing switch breaks
  Tête Rousse.

You added two tools and one deliverable:
- the GLAMOS-driven firn model (Matthias; firn maps and depths for all Swiss glaciers), which you
  will load into the repo;
- snowline calibration: Alexandra's code and paper, which you will give me access to, and the
  Alpine snowline product (ESSD preprint essd-2026-35);
- a probabilistic map: the probability of being polythermal, and where cold ice is likely.

**First diagnostic result (calibrations already run, 4 zone glaciers).** Allowing the
precipitation factor below 1.3 brings the winter snow to GLAMOS, but barely moves the ELA:

| glacier (RGI7) | c_prec (today → 0.8 bound) | winter balance, m w.e. (GLAMOS) | ELA, m (calibration period mean) |
|---|---|---|---|
| Corvatsch (02103) | 1.30 → 0.95 | 1.20 → 0.87 (≈ 0.83, 2013–25) | 3164 → 3175 |
| Hohsaas (01112) | 1.30 → 0.80 | 1.66 → 1.03 (≈ 0.97) | 3205 → 3224 |
| Altels (01265) | 1.30 → 0.80 | 1.55 → 0.94 | 2956 → 3009 |
| Tête Rousse (00777) | 1.30 → 0.80 | 1.63 → 1.00 | 3179 → 3186 |

The calibration with the 1.3 bound, run on these 4 glaciers alone, reproduces today's snow and ELA:
c_prec identical, winter balance within 1 %, ELA within 10 m. DDF and T_off trade off slightly
differently, for example Altels 4.35 / 0.00 against 3.73 / 0.59, because on 4 glaciers the cell
means that seed T_off differ. Why the ELA barely moves:

```
 calibration: one glacier-wide target (geodetic, 2000-2020), two knobs that trade off
   less snow (c_prec down)  ->  needs less melt to hit the target (DDF down)
   snow and melt both shrink  ->  where they balance (the ELA) hardly moves
 a second observation is needed to pin WHERE snow survives: the snow line / firn line
```

So the precipitation bound alone cannot fix the firn. The snow line and the firn data are what
pin it.

## Step 1: finish the diagnostic (0.5–1 day; approved)

**What and why:** does the corrected snow fix the zone, and does Tête Rousse keep its frozen bed?
The answer shows how much the snow matters, before we build anything.

- **icetemp:** a `--calibration-dir` option in `full_set_runs.py`. The runner already supports a
  calibration source (`runner.py: calibration_source_dir`, the catchment renaming in
  `copy_calibration_data`).
- **Runs** on the 4 zone glaciers (`--ids`, own catchment), about 5 min each:
  - FRH at f_ref 0.85 with the calibration of the 1.3 bound ('ref') and of the 0.8 bound ('low');
  - the same two with the snow-refreezing switch (GloGEM worktree `snow-refreeze`);
  - the 4 runs without deepening, for the bed.
- **Compare, period-matched against GLAMOS:**
  - firn line and ELA by year (`firnice_firnline`), and winter snow on 1 May;
  - zone classes and RMSE at the boreholes;
  - Tête Rousse hole 64 and the bed types.
- This is a diagnostic only: the regional calibration is not changed anywhere. No note is sent; you
  will give me access to Alexandra's snowline calibration instead (Step 2).

## Step 2: read the snow and firn tools you bring in, then decide (after you upload them)

**What and why:** you will give me two things:
- Alexandra's snowline-calibration code and her published paper;
- Matthias's firn model, loaded into the repo.

I read and run them first. We then decide together what goes into the model and how.

- **Alexandra's snowline calibration** (paper read: von der Esch et al., in review, GloGEM in the
  Indus basin):
  - **How it works:** c_prec is scanned from 0.8 to 3.0 in steps of 0.2. For each value, DDF_snow
    (1.75–4.5) is set so the geodetic balance 2000–2019 is matched; T_off is used only when DDF hits
    a bound. Of these equally good pairs, it keeps the one whose simulated transient snow line
    best matches the observed one (weighted RMSE, spline minimum).
  - **Its results:** on 12 test glaciers the snow-line error fell from 255 to 113 m and the
    winter-balance error halved (0.42 → 0.21 m w.e.), with the annual balance unchanged.
  - This is exactly the second constraint our Step-1 result says is missing.
  - **For us:**
    - the Alpine transient snow lines (Sommer et al. 2026) replace the HMA product. They cover
      glaciers of at least 1 km²: of ours, Hohsaas, Aletsch, Fiescher, Gorner, Lys, Morteratsch,
      Bossons, Taconnaz and 01230, but not Corvatsch, Altels or Tête Rousse;
    - GloGEM already writes a monthly snow line (commit `5483c21`). Her calibration ran daily, our
      thermal runs are monthly, so whether to calibrate daily and carry the parameters over is
      decided once I have her code;
    - small glaciers (below 1 km²) keep the geodetic calibration in her framework. For our zone
      glaciers, that means a rule for glaciers without snow lines, e.g. a regional relation or
      GLAMOS winter balances where they exist.
  - **Her code arrives by the end of the week** as a separate branch of the GloGEM repo. It was
    written when glogem.pro was still one big script that held all procedures, including the
    snowline calibration.
    - The goal is to port it into its own procedure on our GloGEMflow branch, e.g.
      `procedures/calibration/calibrate_snowline.pro`, next to the existing
      `calibrate_single_glacier.pro`, with a switch.
    - Off, the calibration stays byte-identical.
    - Main session, one reviewer at max effort, and a gate that reproduces her published numbers on
      a test glacier where possible.
  - **Test, after the port:** recalibrate the calibration glaciers that have snow lines
    (diagnostic). Compare the snow, snow line and ELA with GLAMOS, period-matched.
  - **Before that, possible now:** score the current calibration's monthly snow line against the
    Sommer snow lines for those glaciers.
- **Matthias's firn model** (`glogemflow_icetemp/firn_model/firn_change.pro`):
  - **Rules:** a layer per positive annual balance (m w.e.), melt removes the newest layers first,
    layers older than 20 years become ice, and firn counts after two consecutive years. Outputs:
    firn thickness, age, years since firn loss, annual firn area.
  - **Inputs:** the GLAMOS distributed mass-balance grids of
    `/scratch_net/vierzack06_fourth/mhuss/MassBalance_all/`, which are readable locally: 1970–2025
    grids, plus 100 m band balances in `bands_YYYY.dat`.
  - **Coverage of our glaciers:** Corvatsch, Hohsaas, Aletsch, Fiescher, Gorner, Morteratsch, Sex
    Rouge, Tortin and Alphubel. Not Altels, Tête Rousse, Bossons, Taconnaz (France) or Lys
    (Italy). Fee is the outline of the sealed Chessjengletscher (01161), so it stays out of
    development checks.
  - **Run:** with local paths in a copy; your file stays as it is.
  - **Test:** firn extent and thickness by year, averaged to GloGEM's bands, against the thermal
    firn mask (`firnice_firnline`) and the Herron–Langway meltwater depth.
- **Data that help both** (found today):
  - **The Alpine snowline product** (Sommer et al., ESSD preprint essd-2026-35; Zenodo
    10.5281/zenodo.18223929): transient and end-of-summer snow lines, 2000–2025, for the 408 Alpine
    glaciers of at least 1 km², RMSE about 100 m. It covers Hohsaas, Aletsch, Fiescher, Gorner,
    Lys, Bossons and Taconnaz, but not Corvatsch, Altels or Tête Rousse (too small).
  - **Unpublished Sentinel-2 surface-type maps** (snow / firn and old snow / ice), 10 m,
    2015–2021, in `glazioarch/GlacioData/Sentinel2_GlacierSurfaceType_2015-2021/`. They give an
    independent check of a firn mask.
- **Output:** one table per glacier (snow, snow line, firn extent and depth: model against data),
  where the misfit concentrates, and a short options note for our decision.

## Order of work this week

```
 now -> end of week (until Alexandra's branch arrives):
   Step 1  diagnostic: thermal runs of the 4 zone glaciers with the 1.3 and 0.8 calibrations
   Step 2  firn model on our 9 Swiss calibration glaciers -> firn mask and meltwater depth vs model
           snow-line scoring of today's calibration against the Sommer snow lines (Python)
   Step 4  probabilistic product: design note; find Daniel's 2024 plan of action
   (small)  outdated comment at firnice_temperature_model.pro:75-81, f_adv note (no runs going)
 end of week -> : port the snowline calibration into GloGEMflow (gate), recalibrate, then Step 3
```

## Step 3: put snow and firn into the model (after the Step 2 decision; a physics step: its own short plan, one reviewer at max, a pre-registered rule)

The two likely parts, to be confirmed by Step 2:

**3a. The mass balance: the snow calibrated with snow lines (Alexandra's method).**
- Today one target (the geodetic mass balance) sets c_prec, DDF and T_off one after another. Adding
  the observed snow lines lets the calibration choose the pair (snow, melt) that both closes the
  budget and puts the snow line where it is observed.
- For our thermal runs we would use the recalibrated parameters for the calibration glaciers,
  documented, until the mass-balance side decides. That is your call (see the decisions).

**3b. The thermal model: real firn thickness for the meltwater depth (Matthias's rules, inside GloGEM).**
- A firn-layer bookkeeping per band, kept by the mass balance and not by the temperatures, with the
  firn model's rules:

  ```
  each year, per band:  a new firn layer after two years with snow left
                        melt removes the newest layers first
                        layers older than 20 years become ice
  meltwater reaches     min( firn thickness , Herron-Langway close-off )
  ```

- It is checked against the firn model's thickness (Swiss glaciers), the Sentinel-2 surface types
  and the snow lines. It runs wherever GloGEM runs, including French and Italian glaciers and the
  projections.
- This is the state-independent criterion the failed depth rule lacked: taken from the mass balance,
  not from the temperatures it changes. So it has no two stable states.
- After 3a and 3b, the snow-refreezing switch is judged again, on realistic snow.

Each part gets the usual sequence:
- a gate (defaults byte-identical, mass balance unchanged where it should be);
- a run set at 7 f_ref values, f_ref chosen leaving one massif out;
- a rule fixed in advance: A's conditions, the frozen-bed table, the known fronts, and zone target D
  (≥ 2 more zone entities right than H).

## Step 4: the probabilistic product (design now, compute at production)

- **Per glacier, band and year:**
  - the probability of being cold, polythermal or temperate;
  - the probability of a frozen bed;
  - the elevation range where cold ice is likely.
- **Where the spread comes from:**
  - the posterior parameter sets of the Bayesian calibration;
  - the forcing (FRH central, H as structural sensitivity);
  - the measurement tolerance of the classes.
- **Check on held-out massifs:** do the probabilities mean what they say? Measured with a
  reliability diagram and a Brier score against the boreholes. Parameter spread alone usually
  underestimates structural error. If it does here, the probabilities are widened (e.g. by the
  neighbouring-borehole floor of about 1.5–1.9 °C), or the gap is stated.
- **Daniel's 2024 plan of action:** match its "probability of being temperate". I search for it
  first; you may know where it is.
- **Scoring:** the three-class skill and the frozen-bed table are reported as probabilities
  (Brier) beside the hard classes.

## Step 5: v1.0 and production

```
 Steps 1-3 done (time box extended for the zone, your rule)
   -> Alpine ensemble with the adopted model -> T2 Bayesian calibration -> Z2-Z4 measured once
   -> v1.0 freeze review (ultracode panel, you opt in) -> Chessjengletscher (sealed) scored once
   -> cost benchmark -> hindcast 1940-2025 -> projections (GMIP4 GCMs x SSPs x posterior sets,
      H as sensitivity) -> probabilistic maps: regime, frozen bed, cold-ice range, years of change
```

Open items carried along:
- the frozen tongues of Bossons and Taconnaz (until 3a/3b allow judging snow refreezing again);
- the massif split of f_ref without T_off (does it shrink once snow and firn are right?);
- the spin-up on the 2003 geometry (your earlier question; the thickness diagnostic waits for your
  go);
- the outdated comment at `firnice_temperature_model.pro:75-81`, and the finding that f_adv scales
  only horizontal advection (a small commit when no runs are going).

## Earlier to-dos carried over (added 2026-10-14 on your check)

| item | when | how |
|---|---|---|
| **rule C**: does any variant (HR: best M1 2.46 °C; FRSH) beat FRH + 'hl' on A's conditions without being worse in the zone? | this week | Python, on the fair runs |
| **sunshine residual check** (old plan 1d): do FRH's errors below 15 m follow summer sunshine? | this week | Python, with the week-1 exposure tables |
| **FRSH (station lapse rates)**: an adoption rule, fixed before it is scored | with Step 3 | pre-registered like A |
| **model-changes document** for the other developers before the main merge (every switch, fix and finding since June) | running, from this week | main session (writing) |
| **for the mass-balance side, through you**: T_off grows with height (93 % of the high glaciers at both bounds); the RGI7 band aspect coding | your call | text given earlier |
| **housekeeping**: your yaml line 38 (`explained_variance` 0.99 → 0.999, your file); the unmerged worktrees `perc-wetting-front` (closed) and `snow-refreeze` (until it is judged again); the LowLatitudes bal_month caveat | as it comes | small commits, no runs going |

## Effort per task (as agreed)

| step | how |
|---|---|
| 1 runs and comparison | main session; runs as background scripts |
| 2 reading Alexandra's code and paper, Matthias's firn model; the data comparisons | main session for the code (it must be understood exactly); diagnostic agents (Opus, high) for the per-glacier comparisons; I check their numbers |
| 3a/3b design and code | main session at max; one independent reviewer (Opus, max); gate and rule as above |
| 4 design | main session; the reliability check later as a diagnostic agent (Opus, high) |
| 5 v1.0 freeze review | ultracode panel, only after you opt in |

## Decisions for you

| topic | recommendation |
|---|---|
| Step 1 as above (diagnostic, 4 glaciers) | approved |
| Step 2: the firn model (in the repo) and Alexandra's method (paper read; her branch by the end of the week) | firn model and snow-line scoring start now; the port into GloGEMflow and the recalibration once her branch is up |
| for our thermal runs, use a snowline-constrained recalibration of the calibration glaciers until the mass-balance side adopts one | decide after Step 2 |
| firn thickness from the mass balance as the meltwater depth limit (3b) | decide after Step 2 |
| the probabilistic product as the main deliverable beside the hard classes | yes |
| Daniel's 2024 plan of action | if you know where it is, tell me; otherwise I search |

## Verification

- **Step 1:** the 'ref' calibration reproduces today's snow and ELA (done: c_prec identical, winter
  balance within 1 %, ELA within 10 m; DDF and T_off trade off slightly differently). The thermal
  runs with 'ref' are therefore compared with 'low', not with the fair run.
- **Steps 2–3:** every data comparison is period-matched. The gate is defaults byte-identical. The
  rules are fixed and their scorer committed before the run set; a known-answer check precedes
  each scorer.
- **Step 4:** the reliability check is computed on held-out massifs only.
- Book chapter 8 (then 9), memory and commits as each step lands. No config files are committed;
  config.pro stays yours.
