# 5. The thermal spin-up

> **In short.** The boreholes mostly sample ice tens of metres deep, and that ice is almost
> entirely spin-up state: 30 years of main run barely reach it. So the spin-up largely decides
> what the calibration sees. In September 2026 we found that it used one year of climate
> instead of thirty, was too slow to reach equilibrium on thick glaciers, and ran different
> physics from the main run. **Main lesson: the spin-up is part of the model, not a
> preliminary — it must run the same physics, and it must actually reach equilibrium.**
>
> *This chapter is written as the work lands. Status (2026-09-24): Phase 0 done; drift correction
> prototyped on one column; IDL implementation next.*

## What we believed

The thermal spin-up (`spinup_firnicetemp_thermal.pro`) cycled a 1961–1990 mean climate until
the deep layers stopped changing (≥ 200 cycles, tolerance 0.05 °C), then ran the reference
period year by year up to the start of the main run. We assumed that:

- it used the 1961–1990 climate the setting names;
- 200 cycles reached equilibrium;
- apart from simplifications that do not matter for the deep ice, it ran the main run's physics.

## What we measured

### It used one year of climate, not thirty

Calibration runs load only 1990–2020. The spin-up asked for 1961–1990, found just 1990 and
averaged that single year. 1990 is the third-warmest year of the thirty, **0.8 °C** warmer than
the 1961–1990 mean. The check for a missing period never fired, because it only tests that each
calendar month occurs at least once — and 12 months of 1990 pass.

### It was far too slow on thick glaciers

A column settles towards equilibrium over a time that grows with the **square of its
thickness**, roughly *τ ≈ 4H² / (π²κ)*:

| glacier (RGI7) | settling time τ |
|---|---|
| Glacier de Tête Rousse (00777) | ≤ 15 yr |
| Grenz/Gorner (01225) | up to 391 yr |

Each cycle only closes a fraction of the remaining gap, so thick columns need thousands of
cycles. Against a 3000-cycle brute-force reference (same physics), the error left by each
stopping rule was:

| stopping rule | cycles (00777 / 01225 / 01230) | remaining error °C (00777 / 01225 / 01230) |
|---|---|---|
| 200 cycles minimum, change < 0.05 °C (default) | 200 / 200 / 200 | 0.001 / **0.452** / **0.534** |
| change per cycle < 0.05 °C | 30 / 112 / 100 | 0.018 / 0.683 / 0.981 |
| change per cycle < 0.01 °C | 30 / 254 / 138 | 0.018 / 0.357 / 0.777 |
| projected remaining drift < 0.05 °C | 57 / 500 (cap) / 500 (cap) | 0.010 / 0.126 / 0.091 |

A small change per cycle does not mean the column is close to equilibrium: a thick column
changes slowly *because* it settles slowly. The projected-drift rule is accurate but hits the
500-cycle cap and doubles the runtime (765 s against 366 s for three glaciers). None of these
rules scales to ~4000 Alpine glaciers, let alone globally.

### It ran different physics from the main run

| | spin-up | main run |
|---|---|---|
| ice flow (advection) | off | on |
| snow, rain | none | real snowpack and rain |
| melt | simplified degree-day formula on monthly means | full melt model |
| climate | 1990 only (see above) | 1991–2020 |
| thermal surface lapse rate | at least 6.5 °C/km | reanalysis |
| season | January → December | starts in October |

This matters for calibration: a parameter that does nothing in the spin-up can only reach the
deep ice through the 30-year main run. With advection off, `f_adv` could not shape the deep
state at all.

### Advection matters at equilibrium, but only for large glaciers

Estimated steady-state effect of advection, and how much of it 30 years of main run reach on
the thickest band:

| glacier (RGI7) | vertical | along-glacier | reached in 30 yr |
|---|---|---|---|
| Tête Rousse (00777) | 0.05 °C | 0.01 °C | 86 % |
| Grenz/Gorner (01225) | 0.64 °C | ~2.6 °C (median) | 7 % |
| Monte Rosa east face (01230)¹ | 0.23 °C | ~1.1 °C | 22 % |

¹ Listed as "Grenzgletscher" in glenglat; the RGI7 outline is the Monte Rosa east face.

On Grenz/Gorner the along-glacier effect alone is larger than any calibration parameter's.

## What we changed and why

*(Filled in as each phase lands; see the plan for the full reasoning.)*

### Phase 0: clearing the ground (2026-09-24)

- The refreezing loop now uses the same slice of ice as conduction (GloGEM `f9024f3`,
  chapter 4).
- Column to the borehole depth committed (GloGEM `6d05a87`).
- Convergence checked from the first node at or below 15 m (GloGEM `59e9835`).
- The interim stopping rules (30-cycle minimum, projected drift) were reverted: the drift
  correction replaces them.

### Phase 1: the drift correction, tested on one column

`notebooks/05_spinup_drift_correction.ipynb` copies one IDL column into Python: the same grid,
time step, boundary conditions, melting-point clamp and upwind advection. Its exact equilibrium
can be computed directly, so every result is checked against the true answer.

![Brute force vs drift correction on one column](figures/ch05_drift_vs_bruteforce.svg)

On seven test columns (cold and warm, 39–259 m, burial, strong emergence, seasonal snow the jump
does not know about), the drift correction reaches the exact equilibrium to ≤ 0.0003 °C in
**11–17 simulated years**. Brute force still has 4.2 °C left after 200 years on the cold 259 m
column, and 0.09 °C after 3000.

**What we learned on the way.** The first version failed on every column whose base reaches the
melting point. There were two separate reasons:

1. *The jump did not know that ice cannot pass 0 °C.* From a cold start it predicted warming of
   up to +95 °C at depth; the clamp cut it off and the next jump swung back, in a loop. Fix: the
   melting point is built into the jump. Layers that would pass it are held at it, and layers
   already there and still being warmed stay there.
2. *The bed behaves differently once it melts.* With a cold bed the geothermal **flux** is fixed.
   Under emergence, a uniform warming of the column can then only escape by diffusing down
   against rising ice, which takes about 9000 years. So the jump read a drift of 0.001 °C/yr as
   "8 °C still to go". Once the bed is at the melting point, the model fixes its
   **temperature** instead and the extra heat melts ice. The jump now switches its bed
   condition to match.

```mermaid
flowchart LR
    A[measured drift D] --> B{bed at the<br/>melting point?}
    B -- no --> C[fixed-flux bed]
    B -- yes --> D[fixed-temperature bed]
    C --> E[solve the jump;<br/>hold layers at 0 °C]
    D --> E
    E --> F[add the jump,<br/>run 3 real years]
```

**Relative depth.** The notebook also derives the burial speed that matches relative-depth
advection: *surface balance × (1 − share of the ice flux passing above that level)*. It equals
the surface balance at the surface and is zero at the bed by construction. The current code
uses *surface balance × (1 − ζ)*, which is only right for plug flow.

The approach, in one diagram:

```mermaid
flowchart TD
    A[Run 5 real years<br/>with the full model] --> B[Measure the drift D:<br/>how fast each layer<br/>is still changing]
    B --> C[Solve for the jump δ that<br/>cancels the drift, band by band<br/>from the top of the glacier down]
    C --> D[Add δ to the profile]
    D --> E[Run 3 real years]
    E --> F{Jump < 0.01 °C<br/>after ≥ 2 jumps?}
    F -- no --> B
    F -- yes --> G[Verify: 5 plain years]
    G --> H[Start main run]
```

Decisions taken (2026-09-24):

| topic | decision | why |
|---|---|---|
| climate | mean 1961–1990, spin-up ends Sep 1990, main run from Oct 1990 | as the flow-model spin-up does |
| thermal lapse rate | reanalysis below the climate grid point (~2400 m), at least 6.5 °C/km above | ERA5's rate is fitted on terrain up to the grid point; extrapolated 2 km up it makes summits too warm |
| vertical ice motion | the band's annual surface mass balance | kinematic condition at the surface |
| ice arriving from upglacier | mapped by relative depth, with vertical motion in the same coordinate | ice flows roughly parallel to the bed; a thinner band above has no ice at the same absolute depth |
| convergence depth | 15 m | the depth of zero annual amplitude on Alpine glaciers |

## What is still open

- The drift correction in IDL, including the along-glacier coupling between bands, validated
  against the 3000-cycle brute-force run `truthref2_w0`.
- Making the spin-up physics match the main run (plan Phase 2), then rerunning the sensitivity
  ensemble (Phase 3).

## Sources

- Plan: [plans/2026-09-24_thermal_spinup_plan.md](plans/2026-09-24_thermal_spinup_plan.md)
- Runs: `data/bayescal/truthref_w0` (3000 cycles, before the refreezing fix), `truthref2_w0`
  (after it), `tol001_w0`, `projdrift_w0`, `longspin_w0`
- Code: `GloGEM/procedures/initialise/spinup_firnicetemp_thermal.pro`,
  `procedures/processing/downscale_gcmdata_monthly.pro`, `procedures/flow/spinup_flowmodel.pro`
