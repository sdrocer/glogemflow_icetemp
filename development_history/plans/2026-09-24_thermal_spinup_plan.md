# A fast thermal spin-up that runs the same physics as the main run

## Context — what problem this solves

The **spin-up** gives the ice-temperature model a sensible starting state before the main run
(1991–2020 in calibration). It matters more than it sounds: the boreholes mostly sample ice tens
of metres deep, and that deep ice is almost entirely spin-up state — 30 years of main run barely
reach it.

Today's spin-up has two problems:

**1. It is slow for thick glaciers, and gets slower as thickness squared.** A thick column settles
towards equilibrium only a little each simulated year. Measured on RGI7 geometry, the settling time
is ≤ 15 yr for Tête Rousse but up to **391 yr** for Grenz/Gorner. So the current 200-cycle spin-up
is still **~0.5 °C away from equilibrium** on the deep glaciers (measured against a 3000-cycle
brute-force run), and a 3000-cycle run of three glaciers takes an hour. That cannot scale to the
~4000 Alpine glaciers, let alone globally — and the thick glaciers, which need the most cycles, are
also the ones with the most bands, so they cost the most per cycle.

**2. It does not run the main run's physics.** All verified in the code:

| | spin-up today | main run |
|---|---|---|
| ice flow (advection) | **off** | on |
| snow, rain | **none** | real snowpack and rain |
| melt | simplified formula on monthly means | full melt model |
| climate | **1990 alone** (see below) | 1991–2020 |
| surface lapse rate | at least 6.5 °C/km | reanalysis |
| season | January → December | starts in October |

The climate one is a silent bug: the setting says 1961–1990, but the calibration runs only load
1990–2020, so the spin-up quietly uses one year — 1990, the 3rd warmest of the 30 (+0.8 °C). The
check that should catch a missing period never fires, because 12 months counts as "found".

**Why this matters for the calibration parameters.** A parameter that does nothing during the
spin-up can only reach deep ice through the 30-year main run. Today `f_adv` does nothing in the
spin-up (advection off), `f_ins` works only through firn density (no snow), and `f_ref` only
through melt (no rain). That plausibly contributes to the three parameters looking alike. The
honest expectation: fixing this should mainly help `f_adv` on large glaciers. `f_ref` and `f_ins`
may still look alike, because both mainly shift the annual mean temperature just below the
seasonal layer.

**Does advection belong in the spin-up? Yes.** Estimated at equilibrium: for Tête Rousse it
changes temperatures by ≤ 0.13 °C, but for Grenz/Gorner the along-glacier flow alone shifts them
by a median ~2.6 °C — larger than any parameter's effect — and today's 30-year main run reaches
only 7 % of that on the thickest band. With the new method it costs essentially nothing, so it goes
in everywhere.

## Decisions (yours)

| topic | decision |
|---|---|
| climate for the spin-up | **mean 1961–1990 climate**, as the flow-model spin-up does |
| lapse rate for the thermal surface | **reanalysis up to the climate grid point (~2400 m); at least 6.5 °C/km above it** — spin-up *and* main run |
| vertical ice motion | burial/emergence = the band's annual **surface mass balance** |
| ice arriving from the band above | mapped by **relative depth** |
| convergence depth | **15 m** |

## The ideas behind the plan, explained

### A. The drift correction — how the spin-up gets fast

Think of the deep ice as a big pot of water cooling down: each year it closes only a small part of
the remaining gap to its final temperature. Brute force simply waits — thousands of years for thick
glaciers.

But after a few simulated years you can *see* where the column is heading and how fast. If you know
how fast something is changing, and how quickly that change fades, you can predict where it will
end up — and jump straight there. Then run the real model a few more years to check. If it still
moves a little, measure again and jump again. Usually two or three jumps are enough.

```
 temperature at 60 m depth
   ^                                             true equilibrium
   |                                   . . . . . . . . . . . . . . . . . . .
   |                        . . ' ' '
   |               . . ' '        brute force: creeps up for ~1000 years
   |         . ' '
   |     . '      jump  ────────────────────────────────►  X  (lands here)
   |   .'
   | .'
   +───────────────────────────────────────────────────────────►  simulated years
     |<─ 5 real ─>|<3>| measure the drift, jump; repeat 2–3×; then verify
```

**Why this is safe:** the drift is measured on the *real* model. The jump is only a prediction.
If the prediction is a bit off, we land a little short or long, the next measurement shows it,
and we correct. We stop only when **the real model no longer drifts** — so the result is, by
construction, the real model's own equilibrium. A simplified copy of the physics can never
decide the answer. That matters in this codebase, where most of today's bugs ran cleanly and
were silently wrong.

**Cost:** about 16–22 simulated years per glacier, whatever its thickness — against 200 today, which
is still ~0.5 °C short on thick glaciers.

### B. The jump calculation, term by term

For each band, the jump δ (a temperature change per layer) is found from

```
   L · δ   +   H · (δ from the band above  −  δ)   =   − D
```

| symbol | what it is | plain words |
|---|---|---|
| **D** | measured drift | how much each layer warmed or cooled over the last simulated year, per unit time: *where is the column still heading, and how fast?* |
| **δ** | the jump | the temperature change we add to each layer |
| **L** | settling of the column | heat conduction between layers plus vertical ice motion — the same coefficients the model uses. It answers: *if I shift the profile by δ, how fast would conduction and burial/emergence move it?* |
| **H** | along-glacier exchange | see below |

The equation reads: *find the jump that cancels the drift* — the slow settling that would take
centuries, done in one step.

Boundary conditions, and the layers left alone:
- **surface: δ = 0.** The climate sets the surface temperature; we never change it.
- **bed: no change in heat flow.** The geothermal flux is fixed, so the jump must not alter it.
- **layers that received meltwater last year are left alone.** Refreezing sets their temperature,
  not slow diffusion, and the real model handles them in the next simulated year.
- **layers at the melting point that are still being pushed warmer are left alone.** They cannot
  get warmer. Ones being cooled are allowed to move.

**Order of solving:** highest band first (nothing flows into it), then each band below, using the
jump just found for the band above. Ice only flows downhill, so one pass is exact.

### C. The term H = u / dx — "how often the ice gets replaced from upglacier"

```
   band i+1 (upglacier)          band i
   ───────────────────           ───────────────────
        ice ──────── moves downhill at speed u ────►
   |<──────────── dx: distance between the two bands ───────────>|
```

- **u** — ice speed in that layer (m/yr). Fastest near the surface, zero at the bed.
- **dx** — along-glacier distance to the band above: 10 m of elevation ÷ tan(slope). For example,
  57 m on a 10° slope.
- **u/dx** — how many times per year the ice in that layer is replaced by ice arriving from the
  band above. At u = 50 m/yr and dx = 57 m that is ~0.9 per year: after about a year, the layer's
  ice has come from upglacier.

It is called a *relaxation rate* because it keeps pulling band i's temperature towards the upglacier
temperature at that speed. The model caps it for numerical stability (at most 0.8 of the distance
per time step). The jump uses the same cap, so it aims at what the model actually does.

### D. Where the refreezing energy inflation comes from

Each temperature node stands for a slice of ice. Around 9 m the layers change from 1 m to 5 m apart:

```
 depth
  8.0 m  ── node 8
                    ┐  refreezing loop's slice for node 9: 1 m (8 → 9 m)
  8.5 m  ┐          │
  9.0 m  │  ── node 9┘
         │   conduction's slice for node 9: halfway to each neighbour
         │   = 8.5 m → 11.5 m = 3 m of ice
 11.5 m  ┘
 14.0 m  ── node 10
```

When meltwater reaches the 9 m node, the refreezing loop asks: *how much water does it take to warm
this slice to 0 °C?* — and computes that for **1 m** of ice. It then sets the node to 0 °C. But for
conduction that node is **3 m** of ice. So water that could warm only 1 m has warmed 3 m:
**three times too much energy.** At 59 m the gaps are 5 m above and 20 m below, a 12.5 m slice
against 5 m: **2.5×**.

**This is a regression I introduced.** Before yesterday's stencil fix (`dac9fa4`), conduction also
treated the node as a 1 m slice, so the two agreed — both wrong in the same way. The fix corrected
conduction but not the refreezing loop. It matters wherever meltwater percolates past 9 m, which is
common in firn — exactly where `f_ref` acts. The fix is to use the same slice in the refreezing
loop (`firnice_temperature_model.pro:176, 192`).

### E. Relative-depth mapping — ice arriving from a thinner band

```
 band i+1 (upglacier), 50 m thick         band i, 100 m thick
 surface ──────────────                   surface ──────────────
  10 m ●  ζ = 0.2  ──────────────────────►  ● 20 m   ζ = 0.2
  25 m ●  ζ = 0.5  ──────────────────────►  ● 50 m   ζ = 0.5
  50 m ▓  bed      ──────────────────────►  ● 100 m  bed
 (today: band i's layers at 60–100 m take band i+1's "below-bed" value — bedrock, not ice)
```

Today band i takes band i+1's temperature at the **same depth**. When band i+1 is thinner — true
in about half the bands — there is no ice at that depth upglacier, and the model copies the bed
value downwards instead. Relative depth (ζ = depth ÷ thickness) instead takes the upglacier
temperature at the **same fraction of the thickness**, because ice flows roughly parallel to the
bed.

**One thing to get right (your call: solve it the same way):** once ice moves along relative-depth
levels, the vertical motion is also measured relative to those levels — both parts of the flow use
one coordinate. Otherwise part of the vertical motion would be counted twice or not at all. The
prototype notebook derives that form and checks it conserves heat on a toy two-band glacier that
thickens downstream.

### F. The lapse rate, above and below the climate grid point

```
 elevation
  4500 m  ●  Colle Gnifetti band      ▲  above the grid point we extrapolate:
          │                           │  use at least 6.5 °C/km
  2400 m  ■  ERA5 climate grid point ─┤  ERA5 lapse rate is fitted over terrain up to here
          │                           ▼  below: the reanalysis lapse rate
  2000 m  ●  glacier tongue
```

ERA5's lapse rate describes the terrain around the grid point, up to ~2400 m. Extending a shallow
near-ground rate 2 km upwards makes summits too warm — that was the documented reason for the
floor (commit `46b2e0c`). The same argument says nothing about bands *below* the grid point, so
there the reanalysis rate stays. The two meet at the grid point, so there is no jump. This applies
to the **thermal surface only**. Melt and mass balance keep the reanalysis rate, which is
consistent with keeping energy and mass separate.

### G. Why 15 m for the convergence test

The convergence test only looks at ice below the seasonal wave. **15 m** is the commonly cited
depth of zero annual amplitude for Alpine glaciers, so it is easy to defend. In ice the annual
swing at 15 m is about 1.4 % of the surface swing (damping depth ~3.5 m); in firn it is below
0.2 %. The test also compares the same month in successive years, so what is left of the seasonal
wave cancels anyway. Defined as *"layers at and below 15 m"*, it starts at the 19 m node on
today's grid and at ~15.7 m on the planned geometric grid — 15 m is always a true lower bound.
**My view: a good choice.**

## The plan

### Phase 0 — fix what is broken, land what is finished
1. **Fix the refreezing slice** (D). It is a regression, so it goes first. Regression tell: without
   percolation past 9 m the output is unchanged; with it, less warming at 9 m and 59 m.
2. **Commit fix 3** (column to borehole depth; verified 482/482 observations inside the column).
3. **Convergence depth = first node at or below 15 m** (G).
4. **Revert the interim "stop after 30 cycles" edits.** Your goal — stop as soon as equilibrium is
   reached — is exactly what the drift correction does. The interim versions measured either
   inaccurate (up to 0.98 °C off) or expensive (running into the 500-cycle cap).
5. **Rerun the 3000-cycle brute-force reference** in the background, so Phase 1 validates against
   today's physics including fix 1.
6. Save to memory: plans explain each step in plain words, with schematics; and the history book
   below is kept up to date as work lands.
7. **Start the model-development history book** (see next section).

### The history book — how the model got where it is
A new folder in the icetemp repo, one Markdown file per chapter, figures next to it:

```
glogemflow_icetemp/development_history/
  README.md                          table of contents + how chapters are written
  figures/                           image files the chapters show inline
  01_starting_point.md               glenglat, CART calibration, perm_frac / dT_scale
  02_parameter_redesign.md           f_ref, f_ins, f_adv; energy and mass kept separate
  03_why_the_parameters_look_alike.md  f_ref / f_ins degeneracy, weak perm_frac, inert snow
  04_depth_window.md                 bed node, implicit solver, stencil, column to borehole depth
  05_thermal_spinup.md               this plan: one-year climate bug, drift correction, physics
```

Each chapter follows the same pattern: *what we believed → what we measured → what was wrong →
what we changed and why → what is still open*, with the numbers and commit hashes, in plain words.
**Figures sit inside the text, right where the concept is explained.** Markdown does this well —
three kinds, chosen per figure:

| kind | used for | how it appears |
|---|---|---|
| **PNG plots** from the notebooks | results: profiles vs boreholes, drift curves, sensitivity | `![caption](figures/drift_vs_bruteforce.png)` — the image shows inline on GitHub, in VS Code's preview and in a PDF |
| **SVG schematics** drawn with matplotlib | concepts: the node slices (D), relative depth (E), the lapse rate (F) | the same way, as clean line drawings that stay sharp when zoomed |
| **Mermaid diagrams** written as text in the file | flow of steps: the spin-up loop, the order the solver works through bands | GitHub and VS Code draw them automatically; no image file needed |

The ASCII sketches in this plan are the drafts; in the book they become proper SVG schematics. The
figures are made by one script per chapter (`development_history/figures/make_ch05.py`) so they
can be redrawn when numbers change.
- **Markdown**, not PDF: it shows directly on GitHub with the figures inline, is easy to edit, and
  diffs cleanly. A PDF of the whole book, figures included, can be made with one `pandoc` command
  whenever you need one.
- Chapters 01–04 are backfilled from the git history, the notebooks and my saved notes. I write
  chapter 05 as each phase of this plan lands, so it is finished when the work is.
- Figures for the book go in `development_history/figures/`. That is an explicit exception to "no figures
  unasked"; figures anywhere else still need your go-ahead.
- Committed in its own commits, separate from code. The book's figures are committed too, because
  the chapters need them to display. Checked: the repo's `.gitignore` only excludes the top-level
  `figures/`, not `development_history/figures/`.

### Phase 1 — build the drift correction (current physics; default unchanged)
1. **Python prototype notebook** (`notebooks/05_spinup_drift_correction.ipynb`): one real column
   (real grid, the model's conduction, seasonal surface forcing, geothermal bed, vertical motion).
   Show that the jumps reach the brute-force equilibrium in ≤ 4 steps. Also derive and test the
   relative-depth vertical-motion form (E).
2. **IDL**, in `spinup_firnicetemp_thermal.pro`, behind `firnice_spinup_solver = 'cycle' | 'drift'`
   (default `'cycle'`): 5 real years, then blocks of 3; measure D, solve for δ (B), jump; stop after
   ≥ 2 jumps once every jump is < 0.01 °C; verify with 5 plain years. After 8 jumps without success,
   fall back to plain cycling and flag `converged = 0`.
3. **Validate** against the brute-force reference: ≤ 0.02 °C everywhere; 50 extra plain years change
   nothing below 15 m by more than 0.01 °C; runtime on the same host.

**→ Check-in with you**, with the numbers, before any physics changes.

### Phase 2 — make the spin-up run the main run's physics
One commit each, so every change can be attributed.
1. **Mean 1961–1990 climate** from the full reanalysis series, which is still in memory as
   `rtemp/rprec`. The grid-cell indices must be saved under their own names right after
   downscaling, because the solver's arrays overwrite `bb`/`cc`. At least 20 years per calendar
   month, otherwise stop loudly. The spin-up ends in September 1990 and the main run continues in
   October 1990, as the flow model does — no year-by-year run-in, which saves 30 simulated years.
2. **Same forcing as the main run**, computed without side effects: the full melt formula
   (sub-monthly variability, degree-day factor by surface type), high-elevation precipitation
   reduction, a snowpack with the January reset, and rain. Check it is identical to the main
   run's first year.
3. **Season** runs October → September, like the main run.
4. **Firn mask:** a thermal-only mask from the 1961–1990 balance, used in the spin-up and the main
   run's first 5 years. The mass-balance mask is untouched.
5. **Lapse rate** (F), in the main-run thermal model as well.
6. **Vertical motion** = surface mass balance, measured by depth, firn-aware. This replaces the
   layer-index profile and the snow-reservoir term.
7. **Relative-depth mapping** of ice arriving from the band above (E).
8. Rename initialise's `tt` (the mean annual temperature array), which the model's layer count
   silently overwrites.

**→ Check-in with you** before switching defaults, with the effect of each step.

### Phase 3 — switch over and redo the science
1. One brute-force run under the new physics (3 glaciers, once) to confirm the drift correction
   still matches.
2. Defaults become `'drift'` and the new physics.
3. Python rejects runs whose spin-up flagged `converged = 0`, fewer than 20 reference years, or a
   failed identity check.
4. Rerun the 32-point, 7-glacier sensitivity ensemble and notebook 02.

## Checks against changes that run cleanly but are wrong

| silent failure | check |
|---|---|
| reference climate short, or from the wrong grid cell | ≥ 20 years per calendar month (hard stop); raw series equals `temp` where they overlap |
| spin-up forcing drifts from the main run | identity check on surface temperature, melt, rain and snow in the first main-run year |
| snowpack grows without limit | same-month snow equal across years; never above twice the annual snowfall |
| jumps settle somewhere that is not an equilibrium | 50 extra plain years move nothing below 15 m by more than 0.01 °C |
| a parameter never reaches deep ice | preflight notebook 01: each parameter ×1.5 must move the spin-up state below 15 m in ≥ 10 % of Grenz bands; under today's physics `f_adv` must move it by exactly 0, which proves the check works |
| spin-up leaks into mass balance | mass-balance output bit-identical across spin-up modes |
| season mismatch | last spin-up month = September 1990, first main-run month = October 1990 |
| the numerical cap, not physics, sets the result | share of capped bands logged; warn above 5 % |

## Noted, not in this plan
- The standalone ice-flow estimate gives unphysical speeds on steep, thick bands (up to
  ~1400 m/yr), held back only by the numerical cap.
- The model year is 360 days while speeds are per 365.25 days (1.4 %). Document it.
- Solving all bands at once instead of one by one is the main lever for all-Alps and global runs. A
  separate performance project.
- Using the 1961–1990 balance for vertical motion is a calibration-window simplification. Long
  projections will need the evolving balance.

## Status of the earlier depth-window plan

| fix | state |
|---|---|
| 1 bed node in output | committed `47d0c30` |
| 2a implicit solver default | committed `daca009` |
| 2b conduction stencil | committed `dac9fa4` — its refreezing side effect is Phase 0.1 |
| grid dependants tied to depth | committed `b64cf35` |
| 3 column to borehole depth | verified, not yet committed → Phase 0.2 |
| 2c geometric grid | not started |
| 4 strain heating fix + sizing | not started |
| 5 re-measure separability | becomes Phase 3.4 |
