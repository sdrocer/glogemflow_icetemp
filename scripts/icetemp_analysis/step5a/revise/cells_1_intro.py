md(r'''
# 06 — The calibration target: can the boreholes tell f_ref and f_ins apart? (Step 5a)

> **In short.** A calibration needs a *target*: the way model and boreholes are compared. Today's
> target (T1) compares one mean per borehole entity, so it sees only the **level** of the
> profiles. With it, f_ins cannot be calibrated: a fit is expected to leave it two thirds of its
> prior uncertainty (spread left 0.66), and more refreezing can still be undone by less insulation
> (correlation 0.88). Comparing every 1 m depth bin (T2), or every depth bin and month (T3), also
> uses the **shape** with depth and the **season**; both are expected to tell the two apart
> (correlation 0.58 and 0.34, f_ins spread left 0.19 and 0.12, by both methods).
>
> **But the level and the shape ask for different parameters** (5.4 and 6.9 standard deviations
> apart). The shape alone wants about half the refreezing, at which the high firn of Monte Rosa and
> Col du Dôme would be about 2 °C too cold; a whole T2 fit lands in between, with that firn
> 0.6–0.85 °C too cold (T3: 1.0–1.3 °C).
>
> **What the rule fixed beforehand says depends on the error model.** With the August error model,
> as planned, **T2 is outcome 2**: only "consistent" fails, and the held-out cost is inside the
> margin by every method (+0.07 °C expected; 90 % interval −0.12 to +0.25 °C). **T3 is not
> resolved:** +0.08 °C by the stand-ins, +0.23 °C by the real runs, and no real run lies near its
> optimum. With an error model
> fitted to the misfits, which describes them better, **both are outcome 3**: the cost is +0.19 to
> +0.37 °C, above your 0.10 °C margin. And with the emulator's own uncertainty varying with the
> parameters, as the pipeline computes it, neither target is told apart on this 29-run design: a
> risk the campaign must handle.
>
> **Recommendation: build T2,** with the August error model and f_adv fixed at 1.0, and document
> the level/shape conflict. The rule selects T2 over T3 because T3 does not come "at no higher RMSE
> cost"; T3 stays an option once runs near its optimum exist. **Your decisions, starting with the
> error model, are in section 7,** with the limits of this analysis (7 glaciers, 32 runs, 6 folds;
> Monte Rosa holds 55 % of the points; the separation numbers are expectations, not calibrations).
>
> *Status (2026-10-09, revised the same day after two independent reviews): analysis only
> (Step 5a); no calibration was run, nothing in `src/` changed, no IDL. Ensemble `alps_s2` (32 runs
> with the snow fix, Step 4) on the 7 Alpine test glaciers, 66 borehole entities; `alps_r7` (before
> the snow fix) for comparison. Data and model rows come from
> `scripts/icetemp_analysis/target_analysis.py` and its cache; the expensive Stage 2–3 results
> (emulator, bootstrap, checks) are read from
> `/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/` (section 9 lists them). Section
> 8 compares every number recomputed here with the stage reports, and every number quoted in the
> revised text with the values computed here.*

**Contents.** 1 What a target is · 2 Data and rows · 3 The month-to-month error term · 4 Can the
targets tell the parameters apart? (4.1 the emulator's own uncertainty) · 5 Held-out RMSE,
consistency, honest error model (5.4 the error model decides) · 6 The checks (a)–(g) · 7 The
decision rule and what you decide · 8 Reproduction check · 9 Where the read numbers come from

| word | meaning |
|---|---|
| entity | one glacier's boreholes in one 10 m elevation band, e.g. Grenzgletscher@4455 (66 on the 7 glaciers) |
| point | an entity's readings pooled into a 1 m depth bin; 771 points have a model value in every run |
| row | one number compared with the model; what it is depends on the target |
| run | one forward GloGEM run with one parameter set: 32 in the ensemble; the 29 with f_ref ≥ 0.1 train the smooth stand-ins (below 0.1, Tête Rousse leaves the melting point and the model jumps) |
| prior sd | the spread a parameter has before the fit: f_ref 0.289 (uniform 0–1), f_ins 0.326, f_adv 0.244. All sensitivities and spreads are in these units |
| spread left | uncertainty after the fit ÷ uncertainty before it: 1 = nothing learned, 0.5 = halved |
| correlation | how strongly f_ref and f_ins can still trade off after a fit: ±1 = one fully replaces the other |
| expected | estimated from how fast the misfit rises around one parameter set (its slopes), or from a fit on a stand-in; no calibration was run |
| best point | the parameter set a target's smooth fit prefers |
| 5 best runs | the 5 real runs (of the 29) with the smallest misfit for that target |
| held-out | scored on a glacier system the fit did not see (6 systems: Monte Rosa = 01225 + 01230 counts as one) |
| surface, emulator | two smooth stand-ins for GloGEM between the runs: a polynomial per row (quadratic in section 4, cubic in section 5), and the calibration pipeline's statistical emulator |
| E | a stand-in's own uncertainty: how far its prediction may be off; for the emulator it grows with the distance to the nearest run |
| coverage gate | the calibrator's limit on how far from the nearest run a parameter set may lie |
| error model | Σ, what misfit the comparison expects: measurement error plus model error (section 1) |
| misfit ratio | rᵀ Σ⁻¹ r / n: how large the misfits are compared with what the error model expects (about 1 if it is honest) |
''')

md(r'''
## 1. What a target is, and why it decides what we can calibrate

**The heater and the wall.** Think of the firn as a room that is heated in summer and cooled in
winter.
- **f_ref is the heater:** how much of the meltwater's heat stays in the firn when the water
  refreezes.
- **f_ins is the wall:** how easily heat moves through the firn and, since the snow fix of Step 2,
  through the winter snow on the tongues.

A warmer room can come from a stronger heater or from a better wall. So the *mean* temperature of a
borehole cannot tell the two apart: along a "ridge" of parameter pairs, a stronger heater behind a
leakier wall (more refreezing, more conductive firn) gives the same mean as a weaker heater behind a
better one. Their fingerprints with depth and season differ:

```
 depth     heater (f_ref)                            wall (f_ins)
  0-5 m    summer heat, added from the top           how far the seasonal swing reaches, and how late
  5-20 m   a warm bulge down to where the water      less conductive firn damps the swing within a
           runs out                                  shallower depth, and delays it more
  > 20 m   only the long-term mean                   only the long-term mean
 (on ice bands only the wall acts, through the winter snow)
```

**Three targets.** They differ only in what one compared number averages over:

```
 one borehole entity (glacier, 10 m band)          what each target compares with the model
   depth   Jun-98  Nov-98  Jun-99
     2 m     x       x       x          T1  entity mean      1 number per entity            the level
     5 m     x       x       x          T2  depth            1 number per 1 m depth bin     + the shape with depth
    10 m     x       x       x          T3  depth + month    1 number per bin and month     + the season
    20 m     x
```

- **T1** (today's target) averages each entity's points with their depth weights, as
  `calibrator.py` does: it sees the level only.
- **T2** keeps each point: the level and the shape with depth. A point still mixes months: 403 of
  the 771 points average readings from more than one month, 341 from more than one season.
- **T3** keeps each month of a point apart, so it also sees the seasonal swing. A profile dated
  only to a range of months stays one row, compared with the model averaged over that range
  (section 2 explains why).

**How a target judges a parameter set.** The misfits are weighed by what the error model expects:

  misfit(θ) = rᵀ Σ⁻¹ r,   r = model(θ) − boreholes,   Σ = D + S (+ K + V for T3) (+ E)

| symbol | meaning |
|---|---|
| θ | the parameters (f_ref, f_ins); f_adv = 1.0 unless stated |
| r | the misfits, one per row |
| D | measurement error: glenglat's σ², ×3 for "estimated" profiles |
| S | model error, the August choice: 3.5 °C² × exp(−depth difference ÷ 15 m) between rows of one entity, zero between entities |
| K | T3 only: model error between the months of one point (section 3) |
| V | T3 only: for a reading dated to a range of months, how much the model varies over that range |
| E | the stand-in's own uncertainty (sections 4 and 5) |
| Σ⁻¹ | weighs each misfit by the error expected there, and counts correlated rows about once |

For T1, D + S is pooled to the entity means exactly as `calibrator.py` pools it. The per-glacier
offset δ stays out of the fit (your decision); it appears only as a sensitivity.

**The decision rule, fixed before looking** (plan of 30 September):

| check | passes if |
|---|---|
| told apart | f_ins spread left ≤ 0.5 and correlation ≤ 0.8, by two methods (surface and emulator), at the target's best point and at 4 of its 5 best runs |
| consistent | fitting the level alone and the shape alone gives parameters within 3 standard deviations |
| RMSE | held-out RMSE ≤ T1's + 0.10 °C (your margin) |
| honest error model | misfit ratio between 0.5 and 2 at the best run |

f_adv is freed only if, with three parameters, its spread left is ≤ 0.5 and f_ref and f_ins still
pass "told apart". Outcomes: **(1)** all pass → build T2, or T3 if it cuts the f_ins spread by a
further 20 % at no higher RMSE cost; **(2)** only "consistent" fails → the parameters partly make
up for a model error; build if within your margin, document the conflict, your call; **(3)** "RMSE"
fails → you decide between two parameters at the measured cost and f_ref alone; **(4)** nothing
tells them apart → f_ref alone, and test a parameter with a fingerprint of its own.

*How this notebook reads two points the rule left open:* the f_adv condition is read like "told
apart" (both methods, best point and 4 of 5 runs; section 4 also gives the literal reading), and
the RMSE check, whose method was not fixed, passes only if every method agrees (the stand-ins, the
real runs, and the stand-in's error near the fit; section 7). The error model is the August one
(3.5 °C², 15 m), as the plan set; the fitted one is the sensitivity it asked for (section 5.4).
''')
