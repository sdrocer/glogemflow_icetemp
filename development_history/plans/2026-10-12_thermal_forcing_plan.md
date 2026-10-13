# The thermal forcing at high elevation: the temperature offset, the melt it drives, and the meltwater depth

*Plan of 2026-10-12. Written so each step can be defended in a paper: every change follows a
published method or removes an inconsistency shown by data the model was not tuned to, and every
change is measured with paired runs. The previous plan (30 Sep) is archived in
`development_history/plans/2026-09-30_snow_comparison_target_plan*.md`.*

## Context

Two GloGEM runs on all 24 calibration glaciers (icetemp `ebaaf43`) showed that below 15 m the
model predicts unseen glaciers about as well as the 7 test glaciers (median 1.7–2.0 °C), except in
cold high firn, where it is far too warm:

| site | measured, 15–40 m | modelled | first clue |
|---|---|---|---|
| Mont Blanc summit, 4785–4805 m (outline 00786) | −14 to −17 °C | 0 °C | mass-balance offset T_off = +8.4 °C |
| Colle Gnifetti, outline 01230 | −13.0 to −13.4 °C | −8.5 to −9.7 °C | T_off = 2.67 °C |
| Colle Gnifetti, outline 01225 (same saddle) | −13.0 to −13.4 °C | −12.1 to −12.9 °C | T_off = 0.67 °C |
| Taconnaz, 3400–3800 m (outline 00774) | −3 to −9 °C | 0 °C down to 110 m | T_off only 0.12 °C |

You asked how best to solve the offset problem, to check the meltwater depth, and to judge
everything as we would have to defend it in a paper. Nothing here has been run or changed yet;
three read-only investigations of the code and data underlie it.

## What we found

### 1. The offset closes each glacier's mass budget; it is not the local air temperature

How GloGEM's calibration sets it (`settings.pro:81-97`, `procedures/calibration/`), following Huss
& Hock (2015):

```
 target: the glacier's geodetic mass balance 2000–2020 (Hugonnet et al. 2021)
 step 1  precipitation factor c_prec        bounded 1.3 – 2.3       ── enough? stop
 step 2  degree-day factor DDF              bounded 1.75 – 4.5      ── enough? stop
 step 3  temperature offset T_off           NO bound, ±1 °C steps   ── until the budget closes
 and before step 1: every glacier starts from the mean T_off of its 0.1° grid cell (previous pass)
```

So a glacier's T_off is (its cell neighbours' average) + (whatever its own budget still needed after
the other two factors hit their bounds). Mont Blanc 00786 = 2.64 (cell) + 5.75 (own). Colle Gnifetti
on 01230 inherits its cell's 2.67 °C from tiny summit glaciers that needed +7 to +10 °C, although
01230's own precipitation factor sits at its upper bound: it needed more snow, not more warmth.
In CentralEurope (4073 glaciers) the median is 0.10 °C, 8 % exceed 2 °C and 38 exceed 5 °C, almost
all positive and on small summit glaciers (several with the regional fill value as target).

Two tests with data the model was never tuned to:

*Same site, two outlines.* The Colle Gnifetti boreholes lie metres apart on outlines 01225 and
01230, which also use different ERA5 cells. The measured firn is the same on both.

| | measured 10–40 m | thermal surface temperature at the saddle (run mean) | modelled firn |
|---|---|---|---|
| 01225 (T_off 0.67) | −13.0 to −13.5 °C | −14.03 with T_off, −14.70 without | −12.1 to −12.9 °C |
| 01230 (T_off 2.67) | −13.0 to −13.4 °C | −12.06 with T_off, −14.73 without | −8.5 to −9.7 °C |

Without the offset the two air temperatures agree to 0.03 °C; with it they are 2 °C apart and the
firn 3.5 °C apart.

*Cold dry firn as a thermometer.* Where there is almost no melt, firn at 10–20 m is close to the mean
annual air temperature (a textbook relation, e.g. Cuffey & Paterson 2010; used in the Alps by
Haeberli & Funk 1991, Suter et al. 2001, Gilbert & Vincent 2013). At the sites with a large offset
(Mont Blanc, Colle Gnifetti on 01230), measured firn minus band air temperature is:

| band air temperature (mass-balance lapse rate) | mean | RMS |
|---|---|---|
| without T_off | +0.25 °C | 0.77 °C |
| with T_off | −5.3 °C | 6.0 °C (firn colder than its own air: impossible without melt) |

With the thermal model's steeper lapse rate the numbers are +1.4 °C without and −4.1 °C with the
offset; firn slightly warmer than its air is expected (winter snow damps the cold), so this does not
settle the lapse rate, only the offset. Where the offset is small (< 0.8 °C) the test cannot tell
the two apart, and all 22 of our outlines have offsets ≥ 0, so negative offsets go untested (a
caveat for the paper). The air temperature was
rebuilt in Python from the ERA5 files GloGEM reads and reproduces GloGEM's logged surface
temperature to 0.005 °C, so this test needs no model run.

*Where it enters the thermal model.* At one line (`glogem.pro:530`, the band air temperature
`tg`), and from there into every atmospheric input: surface temperature, the melt and rain whose
refreezing heats the firn, snow cover, firn mask, the meltwater depth, the initial profile and the
whole spin-up forcing. No setting keeps it out.

### 2. The meltwater depth

```
 surface ─┬─ summer meltwater (melt + rain) × f_ref
          ▼  fills the column from the top: each layer is warmed to 0 °C before water goes deeper
   5 m   ███ 0 °C
  20 m   ███ 0 °C      a layer already at 0 °C passes the water on untouched
  59 m   ███ 0 °C
  99 m   ▓▓▓ ← the limit: the pore close-off depth (830 kg/m³, Herron & Langway 1980), 75–120 m
               on high-accumulation bands — water left over is discarded
```

- **A deviation from the cited formula:** in Herron & Langway (1980) the depth of the first
  densification stage does not depend on the accumulation rate; the code makes it proportional to it
  (`firnice_perm_depth.pro:19`, since `780520a`, no reference). Two unit slips: metres of water are
  converted with the ice density (917 instead of 1000 kg/m³), and the latent heat per metre of water
  uses the ice density (`settings.pro:297-298`), so each metre of meltwater counts 9 % short.
- **Fixing them changes Taconnaz little.** Its column is temperate to the bottom because a lot of
  water arrives every summer (about 0.5 m w.e. or more at 3400–3800 m), fills the column from the top,
  and the limit (pore close-off, 75–120 m, often below the bottom of the column) lets the surplus heat
  reach the deepest cold layer. The thermal model's own density profile reaches 830 kg/m³ at 34–39 m,
  so the latent heat goes deeper than the firn the conduction sees.

## How to solve the offset problem — the options, judged as a reviewer would

| | option | for | against |
|---|---|---|---|
| **A** | thermal model's **energy** inputs from the location's climate (no T_off); its **mass** inputs from the calibrated budget | both tests above; no new parameter; the same site gets the same forcing whatever outline it is on | the melt that refreezes in the firn is no longer the mass balance's melt — reported, not hidden |
| B | a regional offset instead of each glacier's | keeps a bias correction | no evidence of a regional air-temperature bias; the offsets are closure terms that trade off against the other two factors (equifinality, e.g. Rounce et al. 2020) |
| C | cap the offset (e.g. 2 °C) | simple | any bound is arbitrary |
| D | recalibrate the mass balance with a bounded offset | fixes the source | the global calibration is shared with GMIP4 and not ours; we report the finding to that side |
| E | leave out glaciers with large offsets | simple | hides it, and summit glaciers are the cold-firn sites we need |

**Recommendation: A.** The split follows what each input physically is:

```
 ENERGY into the firn — from the climate of the place          MASS — from the calibrated budget
   surface temperature                                            snow on the ice (insulation)
   melt and rain that refreeze (same melt model and DDF)         firn or ice (the firn mask)
   initial profile, spin-up temperatures, HL temperature         accumulation (HL depth, burial)
```

The mass balance itself keeps its calibration, so firn extent (checked against WGMS ELAs in
chapter 6) and glacier-wide mass budgets are unchanged. Paper wording, in short: *"The thermal model
is forced by the downscaled air temperature without the glacier-specific temperature offset of the
mass-balance calibration: the offset closes each glacier's budget after the precipitation and melt
factors reach their bounds, and is not a property of the local climate. At Colle Gnifetti, boreholes
on two outlines whose offsets differ by 2 °C would otherwise be modelled 3.5 °C apart."*

## The plan

```
 Step 0  dry-firn check of the air temperature (Python, no run)
 Step 1  GloGEM: switches, defaults = today; code checks (IDL harnesses)
 Step 2  icetemp: pass the switches into the training config; scoring for N runs
 Step 3  gate runs on 3 glaciers (minutes), then 5 runs on all 24 glaciers in parallel (~30 min)
 Step 4  stop: you see the table and decide what becomes the default
 Step 5  book, model-changes text, memory; then the Alpine ensemble is redone before any campaign
```

Steps 1–2 keep today's behaviour as the default, so nothing changes until you decide in Step 4.

### Step 0 — the air temperature checked against dry firn (Python only)

A script rebuilds GloGEM's band air temperature from the ERA5 files exactly as GloGEM does
(grid point from the cell's glacier centroid, monthly lapse rate, the thermal lapse rate above the
grid elevation) and compares its 10-year mean before each reading with the measured 10–20 m firn
temperature, with and without the offset. Two site selections, written into the script before it
is run: (a) every entity at or above 4400 m (the dry-snow zone: Colle Gnifetti, Mont Blanc, upper
Bossons); (b) every entity colder than −5 °C at 10–20 m, which shows where melt makes cold firn
warmer than its air. Honest note for the paper: (a) is the set we already looked at; (b) is the wider
check. Output: a table and one figure (measured firn against air temperature, with and without the
offset). *Pass:* the rebuild reproduces the logged winter surface temperatures to 0.01 °C.

### Step 1 — GloGEM (branch GloGEMflow): switches, defaults = today

| switch | 'y' / new behaviour | why |
|---|---|---|
| `firnice_forcing_toff = 'n'` | the thermal model's energy inputs without T_off: surface temperature, the melt and rain that refreeze in the firn, the Herron–Langway temperature, the initial profile, the spin-up forcing | the evidence in section 1 |
| `firnice_hl1980 = 'y'` | the pore close-off depth from the closed form of Herron & Langway (1980): stage 1 independent of accumulation, depth integrated exactly | the code deviates from the cited formula |
| `firnice_lh_water = 'y'` | latent heat per metre of **water** (1000 kg/m³) in the thermal model only (the mass balance's full refreezing scheme keeps its own constant) | unit slip; equals f_ref × 1.09 |
| `firnice_perm_limit = 'profile'` | diagnostic only (run P): water stops where the thermal model's own density profile reaches 830 kg/m³ (default 'hl' = today) | measures how much the deep limit matters |

- New include `procedures/processing/firnice_thermal_forcing.pro` computes the thermal surface
  temperature, melt and rain from the offset-free temperature with the same formulas, in the same
  order, as `accumulation.pro` and `meltmodel.pro` (degree-day model '1', monthly; anything else
  stops the run with a message). It writes no mass-balance variable.
- Changed reads: `firnice_temperature_model.pro:196, 208, 237, 254`; `firnice_reference_forcing.pro:77,
  101` (the spin-up forcing; the reference balance, snow and firn mask stay as they are);
  `spinup_firnicetemp_refyear.pro:12-15`; `initialise_firnicetemp_spinup.pro:34-35`;
  `firnice_perm_depth.pro:18-26`. Mass inputs keep the offset (snow cover, firn mask, accumulation
  for burial and the HL depth): they were calibrated together with c_prec and DDF, and removing the
  offset there alone would, for example, bury far too much snow at Colle Gnifetti (c_prec at its
  upper bound there).
- A diagnostic log line per glacier: firn-band melt + rain with and without the offset (how much the
  latent heat changes), and a wet-depth record per year and band (deepest node brought to 0 °C, water
  discarded at the limit), to see what happens at Taconnaz.
- Also: fix the `@`-include trap at `glogem.pro:406-407` (K14 in the GMIP4 log; harmless in our runs).
- *Checks before any run* (IDL harnesses, like the snow-conduction test): the new melt and rain code
  equals the mass balance's bit for bit when the offset is kept; the HL depth equals a numerical
  integration of Herron & Langway's equations to 1 mm; the latent-heat change equals rescaling f_ref
  by 1000/917 to 1e-10 °C. The snow-conduction harness is updated to the new names.
- *Commits:* one per switch, each checked byte-identical with its default before committing.

### Step 2 — icetemp

- The runner and the configuration pass the switches into the training config (defaults =
  today); `full_set_runs.py` takes a variant name, uses a fresh run tag per variant (an old `.done`
  file would skip a run) and records the GloGEM commit and a clean status in the run folder.
- `full_set_score.py` scores any number of runs; `compare_run_files.py` compares every output file
  byte by byte; a test checks the training config.
- No config file is committed: each run folder keeps the configuration it used.

### Step 3 — the runs (prior-mean parameters 0.55 / 1.017 / 1.0, no tuning)

*Gate runs first* (minutes), on a small catchment of their own (00786, 00774, 03116) so the shared
catchment file is untouched: with the defaults every file must equal today's run byte for byte; with
the offset left out, Hintereisferner (T_off exactly 0.00) must stay byte-identical and Mont Blanc must
cool.

*Then five runs on all 24 glaciers in parallel* (vierzack03, about 30 min; IDL licences checked
first), each paired with today's run `ce25_pred_w0`:

| run | offset in thermal forcing | HL as published | latent heat per m water | purpose |
|---|---|---|---|---|
| R | yes | no | no | the new code reproduces today, all 308 files byte-identical |
| F | **no** | no | no | the forcing fix alone |
| H | yes | **yes** | **yes** | the meltwater-depth fixes alone |
| FH | **no** | **yes** | **yes** | both |
| P | yes | yes | yes, plus water stops where the model's own density reaches 830 kg/m³ | diagnostic only: how much of Taconnaz's deep heat is the deep limit |

### Checks, fixed before the runs

| check | passes if |
|---|---|
| code reproduces today (R) | every output file byte-identical; F, H, FH leave the mass-balance files identical |
| F touches only the offset (F) | the three outlines with T_off exactly 0.00 (01039, 02216, Hintereisferner 03116) byte-identical to today |
| same site, two outlines (F, FH) | Colle Gnifetti 01230 minus 01225 at 10–40 m: from about 3.4 °C to below 1 °C |
| no offset left in the result (F, FH) | the correlation of deep bias with T_off across outlines no longer significant (permutation test) |
| nothing else harmed (F) | deep RMSE of the 17 other glaciers without Mont Blanc and Taconnaz worse by at most 0.05 °C |
| H | a bug fix: adopted whatever it does; its effect reported |

Reported in one table: RMSE and bias below 15 m and on all points, for the 7 test glaciers and the 17
others, with and without Mont Blanc and Taconnaz; firn and ice points; Mont Blanc; Taconnaz by band;
the change per glacier against its offset; the interaction (FH − F) − (H − today).

*What I expect, so the result can surprise us:* F brings Mont Blanc from 0 °C to about −15 °C and
Colle Gnifetti on 01230 to about −12.5 °C, and changes little elsewhere (most offsets are below
0.8 °C). Neither F nor H cures Taconnaz (offset 0.12; HL depth changes by under 10 m). That will
leave one question open, for you: where the meltwater's heat should stop (today water reaches 45–120
m through layers the model's own density profile already treats as ice from about 35–49 m; run P
measures this).

### Step 4 — stop: you see the table

You decide which switches become the default for the ice-temperature runs (the settings defaults
stay as they are; GMIP4 runs have the thermal module off).

### Step 5 — afterwards

- Book: a new chapter 8 (the thermal forcing: the offset, the two tests, the meltwater depth), a
  dated correction in chapter 6 (it says the offset "corrects the downscaled air temperature").
- Model-changes text for the GMIP4 session (handed over through you; we do not edit their page), and
  the offset finding for the mass-balance side: unbounded, inherited from the cell mean, up to
  +9.7 °C on summit glaciers with fill-value targets.
- Memory updated as results land.
- Then the Alpine ensemble is redone with the adopted forcing before any calibration campaign: every
  calibration number so far (Step 5a, the trial calibration) used the old forcing, and Colle Gnifetti
  on 01230 was one of the sites that shaped it.

## Decisions for you

| topic | recommendation |
|---|---|
| how to treat the offset | A: energy inputs without it, mass inputs with it |
| Herron–Langway and the latent-heat unit | treat both as bug fixes (switch H) |
| run P (water stops at the model's own close-off depth) | yes, as a diagnostic only; it costs no extra time |
| report the offset finding to the mass-balance side | yes, through you |
| where the meltwater's heat should stop | after the runs: the published self-consistent option (the Herron–Langway density profile used for conduction too) or water retention in firn (new physics) — not now |

## Verification

- Step 0: the rebuild matches GloGEM's logged surface temperature to 0.01 °C.
- Step 1: harnesses (melt/rain identity, HL against numerical integration, latent heat = f_ref
  rescaling); greps show the thermal model no longer reads the old variables.
- Step 3: gate runs and R byte-identical; all 111 entities with profiles; every spin-up converged;
  the checks above.
- Tests: `tests/` all pass, plus the training-config test.

## Housekeeping

- GloGEM edits only on GloGEMflow, small commits, nothing changed while runs are going; no process
  killed; IDL licences checked before each launch.
- No config files committed.
- The GMIP4 session's list of ice-temperature files in GloGEM (move to icetemp, drop or keep) is a
  separate cleanup, later, with you.
