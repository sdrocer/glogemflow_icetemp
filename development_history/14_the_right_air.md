# 14. The right air, and what it reveals

> **In short.** Chapter 13 showed that ERA5's free air matches the mountain stations far better than the 2 m air the
> model used, and that the ice-temperature model gets worse with it.
>
> On 23 Oct Janosch turned the question around:
> - **Better input is kept**, and the model's new errors are read as compensations to find and fix.
> - **A standing rule** follows: input data are judged by their own observations, never by model skill.
> - **The free air becomes the base of the thermal model**, without any surface cooling.
> - **The fixes come from the parameter space first**, before new physics. The candidate for the deep error is the
>   past-climate offset of the spin-up.
>
> *Status (23 Oct, 13:45): the base (`fb1`) and the knob runs (`fd2`) are running, 130 runs on five parts per band.
> Plan: `plans/2026-10-23_freeair_base_plan.md`.*

## What we believed

- **An input is adopted when the model gets better with it.** The free-air test of chapter 13 was built that way: a
  rule against today's model (o4).
- **The 2 m air's errors were a known limitation** to be fixed together with a surface energy balance in later work.

## What changed: the standing rule

- **The rule (Janosch, 23 Oct):** input and forcing data are chosen by how well they match observations of the input
  itself (stations, radiosondes, gauges), never by whether they improve the model's results. If better input makes
  the model worse, the model needs fixing. This holds for every future forcing dataset.
- **Why it matters here.** The cold 2 m air made the model look right for the wrong reasons. Its winter cold stood in
  for:
  - the snow surface's own radiative cooling;
  - the colder past that the deep ice remembers;
  - f_ref pressed against its upper limit.

  With the right air these three show themselves (chapter 13):

| where (free air, no cooling, against the measurements) | error | likely cause |
|---|---|---|
| below 3800 m, every depth | as before (< 0.2 °C change) | — |
| above 3800 m, 15–20 m | +0.5 °C | the snow surface is not colder than the air |
| above 3800 m, 30–200 m (Mont Blanc sites, Fiescher) | +2.3 to +3.8 °C | the deep ice's colder past; advection tuned to the old air |
| f_ref | 0.40–0.70 instead of 0.85–1.0 | it no longer makes up for the cold air |

- **Gilbert's surface cooling stays out** (Janosch, 23 Oct): no empirical relations that add processes.

## Parameters before physics

**Janosch's question:** can the parameter space absorb these errors instead of new physics? It can, under two
conditions:
- **The parameter acts where the error is.** What matters is the pattern with depth.
- **The boreholes can tell it apart from the others.**

**What the existing parameters can do:**
- f_ref, f_ins and f_adv change the column with nearly the same depth pattern (correlation 0.9–0.98, measured in
  September; chapter 3). They act as one "warmer/colder" knob.
- That knob can absorb the even +0.5 °C near the surface, but not an error that grows with depth.

```
 depth      error with free air        one "warmer/colder" knob        past-climate offset
            (above 3800 m)             (f_ref, f_ins, f_adv)           (spin-up air colder)
 15–20 m    +0.5 °C   ■                 shifts  ■■■                     barely  ·
 30–50 m    +2.3      ■■■■              shifts  ■■■                     some    ■■
 50–100 m   +3.3      ■■■■■■            shifts  ■■■                     more    ■■■■
 100–200 m  +3.8      ■■■■■■■           shifts  ■■■                     most    ■■■■■■
```
*(expected shapes; the knob runs measure them)*

**The past-climate offset** is the spin-up's air shifted colder (`firnice_spinup_dT`):
- The deep, old ice feels it most; the top 20 m hardly.
- It stands for the colder climate before the run. The deep Mont Blanc holes record about 1.5 °C of warming since 1900
  (Gilbert & Vincent 2013).
- **It mirrors the flow model's spin-up** (Janosch). There, the 1961–90 mass balance and the deformation factor A are
  perturbed until each glacier matches its inventory volume and length. Here, the spin-up climate is made colder until
  the ice temperatures match.
- **One difference:** volume and length are known for every glacier, ice temperatures only at the 20 borehole sites.
  So the offset must be one value for the Alps, calibrated at the boreholes and carried to all glaciers, like f_ref.

**The guard against compensating again:**
- every parameter keeps a physically plausible range;
- a parameter at its bound is read as a missing process, not as a result;
- the values are reported as "effective".

## The runs

- **Five parts per band:** v1.0 uses five orientation parts (Janosch, 23 Oct), so the factors were rebuilt for the
  free air on five parts (orientation session, from the stored pixel factors: north faces slightly colder, south
  faces slightly warmer than with the 2 m air).
- **The base `fb1`:** free air, no cooling, f_ref at all 8 values (40 runs).
- **The knobs `fd2`** at f_ref 0.40, 0.55 and 0.70 (90 runs):
  - spin-up air 1 and 2 K colder;
  - f_adv 0.5 and 1.5;
  - f_ins 0.785 and 1.25.
- **The reading, fixed before the runs** (`scripts/icetemp_analysis/freeair/freeair_diag.py`, icetemp `b866c10`): a
  knob matches the deep error if both hold:
  - its depth pattern above 3800 m correlates with the base's error by depth at |r| ≥ 0.9;
  - it halves the 30–200 m error while 15–20 m stays within ±0.5 °C.
- The scorer's dry run reproduces the first free-air test (2.71 °C on left-out massifs, 3.28 above 3800 m).

## What the runs give

*Pending (chain `week3/chain_freeair_base.sh`, scores in `week3/freeair/score_fb1/` and `week3/freeair/diag_fb1/`).
On Fri 23 Oct at 15:40 the base (40 runs) and the second wave (45) had finished without a failure; the third wave,
the merge and the scoring were left to run over the weekend.*

## What is still open

- Janosch's choice of the parameter space for the calibration (Step 3 of the plan).
- The calibration (Step 5 plan of 20 Oct) on the free air with that parameter space.
- The mass balance on the free air.

## Sources

- Plan: `plans/2026-10-23_freeair_base_plan.md`; decision log, 23 Oct.
- GloGEM `freeair-thermal` `d783f84` (no change for this chapter). icetemp `b866c10` (scorer, plan, decision log).
- Factors: `week3/orient_parts5_fa` (five parts, free air), hole parts `week3/orient_parts5`.

## Glossary

| term | meaning |
|---|---|
| free air | ERA5's pressure-level temperature read at the band's height |
| standing input rule | input data are judged by their own observations, never by model skill |
| f_ref | the share of meltwater that refreezes and heats |
| f_ins | the scaling of the snow and firn conductivity (insulation) |
| f_adv | the scaling of ice advection |
| past-climate offset | the spin-up's air shifted colder (`firnice_spinup_dT`) |
| depth pattern | how a parameter changes the modelled temperature at each depth |
| effective parameter | a parameter whose value also absorbs missing processes |
| fb1, fd2 | the new base; the knob runs |
| o4, o5 | today's air with three and five orientation parts |
| five parts | gentle, north, east, south and west shares of each band, each with its own melt factor |
| error on left-out massifs | each massif predicted with the f_ref chosen on the others |
