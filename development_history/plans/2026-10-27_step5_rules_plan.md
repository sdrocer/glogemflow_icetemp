# Step 5, step 3: the rules that decide the result, fixed before any fit

*27 Oct 2026. For Janosch's approval. This fills in step 3 of the Step 5 plan (2026-10-20_step5_plan.md) for the
physics you chose today (option A: the free air, f_sur 1.4 K fixed, five orientation parts; f_ref and f_ins
calibrated). The 40 design runs compute meanwhile (chain week3/chain_step5.sh, ends about 18:30). Nobody looks at
them, and nothing is fitted, until these rules are committed.*

## In short

| # | rule | what I propose | needs you? |
|---|---|---|---|
| 1 | error model | **C6b**, chosen by evidence: the plain form plus a small error per reading in the top 15 m, and sites that may miss much more than others | approve |
| 2 | sites | by a written rule (holes within 1 km and 150 m); the hand-made table as a check | approve |
| 3 | site checks | three checks per site; one named check for the 1999 campaign's shallow bands | approve |
| 4 | f_ref and f_ins apart? | reported, not a gate (as planned) | no |
| 5 | variants beside the main fit | deep readings only; no site above 20 % of the information; only equilibrium and published holes; a wider f_ins prior; without the 1999 shallow bands | approve |
| 6 | very deep, very old ice | **your choice**: keep today's Colle Gnifetti rule, or one rule for every low-accumulation summit | **decide** |
| 7 | readings above 0 °C | set to 0 °C and flagged, instead of dropped | approve |
| 8 | the benchmark | each massif predicted from the other six; Step 5 against the f_ref grid choice on the same physics and five simple rules; **45 extra runs tonight** | approve |
| 9 | adoption | as planned: at worst 0.1 °C worse than the grid choice (90 % confidence), overall and in the zone, with honest uncertainty ranges | no |
| 10 | the map's uncertainty | as planned: the larger of the calibration's own spread and the spread between left-out massifs | no |

## 1. The error model, chosen by evidence

**What it is.** The error model says how far the model may normally miss a reading, and which readings tend to miss
together. It decides how much each borehole counts in setting f_ref and f_ins, and how wide the uncertainty is. The
panel of 20 Oct found today's form wrong: the model's misses came out three times larger than the form expects.

**How I chose.** Eight forms, listed in the plan before I looked, were each fitted to the misses of the existing run
closest to the Step 5 settings (free air, f_sur 1.5, f_ref 0.85, f_ins 1.017; 2115 readings, 122 bands, 198 holes,
29 sites). The score is AIC: how well a form explains the misses, with a penalty for every extra term. A form must
also pass its own check: measured in its own units of "expected miss", the misses must average about 1 (0.5 to 2),
overall, for the shallow readings (top 15 m) and for the deep ones.

| form | what it adds | AIC (lower is better) | check overall | shallow | deep | passes |
|---|---|---|---|---|---|---|
| C0 | today's form | 8197 | 3.10 | 7.62 | 1.73 | no |
| C1 | + a small error of each reading on its own | 5148 | 1.06 | 1.74 | 0.40 | no (deep) |
| C2 | that error only in the top 15 m | 4196 | 0.97 | 1.01 | 1.03 | yes |
| C3 | different depth links above and below 15 m | 4970 | 1.07 | 1.82 | 0.46 | no (deep) |
| C4 | an error that grows with the years between visits | 5131 | 1.06 | 1.74 | 0.40 | no (deep) |
| C5 | own error sizes at the ice-core saddles | 4763 | 1.03 | 1.54 | 0.51 | yes |
| C6 | C1 + some sites miss much more than others | 3909 | 0.34 | 0.61 | 0.13 | no |
| **C6b** | **C2 + some sites miss much more than others** | **3611** | **1.28** | **1.45** | **1.19** | **yes** |

- **The winner is C6b**, by 585 AIC units over the next form that passes. The same ranking holds with the hand-made
  site table and on the free-air runs without f_sur (week3/step5_error_forms).
- **Why C2 beats C1.** Below 15 m the year's swings have died out, so a reading there has almost no error of its own.
  Near the surface, the exact depth and the month matter, and a monthly model cannot follow them.
- **Why C6b beats C2.** With one error size for all sites, C2 fails site by site in both directions:
  - Tortin, Jungfraufirn and Taconnaz miss 3–12 times more than it expects;
  - the temperate sites miss far less, e.g. Hintereisferner and Sex Rouge 0.1 times, because the model rarely misses
    ice at 0 °C.

  C6b lets each site have its own error size, drawn from one shared spread. With it, only 2 of 29 sites fall outside
  what the form expects: Tortin (too large, 4 readings) and the Grenz tongue at 2015 m (too small, temperate).
- **What C6b means for the calibration:**
  - a site the model cannot fit at any f_ref (Tortin, Jungfraufirn, Taconnaz) counts less in setting the parameters
    than it would with a single error size;
  - a site it fits well counts relatively more.

  That is the honest weighting when the model has errors of its own that the parameters cannot remove, like the
  massif split.

```
  a reading's miss = the band's shared error, smooth in depth (about 2.3 C^2, reaching about 120 m)
                   + the hole's shared error (1.4 C^2)        + the hole-and-year's (small)
                   + the site's shared error (0.4 C^2)
                   + its own small error, top 15 m only (0.4 C^2)
                   + the reading's measurement error
                   all of it scaled per site by a factor drawn from one spread (heavy tails, nu 2.6):
                   some sites miss much more than others
```

**What it costs.**
- The calibration's likelihood today assumes one error size for all sites, so it needs a change of code: a site-by-site
  heavy-tailed likelihood.
- The change is checked against a standard library and by recovering known parameters from simulated data.
- A max-effort reviewer checks it.
- About one day of work. I'll send the free-air window a line before touching the shared calibration code.

**When the sizes are fixed.** As planned:
- re-estimated once at the first result;
- re-estimated inside each left-out-massif fold;
- the form itself stays fixed from here.

## 2. Sites by a written rule

- **The rule.** Holes within 1 km of each other and less than 150 m apart in height belong to the same site, chained
  (a hole joins a site if it is close to any hole of it).
- **The result.** 29 sites from 198 holes, against 23 in the hand-made table.
  - **The rule joins** what lies close together under different glacier names:
    - Lys with the Monte Rosa summit saddles (3898–4485 m);
    - the upper Taconnaz holes with Bossons at the Dôme du Goûter;
    - Gorner with the Grenz tongue.
  - **It splits** what lies far apart in height under one name:
    - the Grenz tongue and Hintereisferner into three sites each;
    - Jungfraufirn, Taconnaz, Bossons and Grubengletscher into two each.
- **The check.** The hand-made table (23 sites) gives the same choice of error model, with almost the same sizes. It
  is kept as a check on the main fit.

## 3. The checks per site

For each site, after the fit:

1. **Are its misses within what the error model expects?** Inside the 99 % range of the form's own law for that
   site's number of readings.
2. **Is its own offset unusually large?** The offset is the site's shared miss. It counts as unusual if it is more than
   2.9 times its uncertainty, which is about a 1 in 300 chance for one site and so about 1 in 10 for any of 29.
3. **Does leaving it out move the result?** Its influence counts as decisive if f_ref or f_ins moves by more than the
   result's own uncertainty.

- A site that fails any check is reported and gets an analysis without it.
- **Only check 3 can stop adoption:** then the result and the result without that site are both reported, and the
  map takes the wider uncertainty.
- **A named check, set before the fit:** the fit without the 1999 campaign's bands that have readings at 15–20 m only
  (Lys 4035–4135 m and the saddle bands).
  - The free-air window found that across 50 high bands the best f_ref follows the measured 15–20 m temperature (rank
    correlation 0.91).
  - Monte Rosa's high f_ref comes from those warm bands.
  - So this check shows how much a single campaign decides f_ref.

## 4. Can f_ref and f_ins be told apart?

- **As planned, reported and not a gate.** We report:
  - how strongly the two trade off in the result;
  - whether the data narrow f_ins beyond its starting range.
- **If they cannot be told apart:** f_ins stays at 1.0, f_ref is calibrated alone, and the trade-off goes into the
  uncertainty.

## 5. Variants beside the main fit

Each is run on the same 40 design runs, with no new GloGEM runs:

- deep readings only (below 15 m);
- no site above 20 % of the information on f_ref (the share each site adds to how sharply f_ref is fixed; a site
  above 20 % is down-weighted to 20 %);
- only boreholes rated in equilibrium and published;
- a wider starting range for f_ins (spread doubled);
- without the 1999 campaign's shallow-only bands (point 3).

They are reported beside the main fit. None of them replaces it.

## 6. Very deep, very old ice: your choice

**Today.** One rule leaves out the Colle Gnifetti readings below 40 m, at the saddle bands at and above 4400 m. There
the ice at depth is decades to centuries old and holds a climate older than the model's spin-up (1942–71), so the
model cannot know it.

**The plan asked for one rule for all glaciers instead of a rule for one site.** Two options:

- **(a) Keep today's rule as it stands,** stated as: "below 40 m at the summit saddles with less than 0.6 m w.e. of
  accumulation per year (measured)". Today only Colle Gnifetti meets this.
- **(b) The same words, applied everywhere:** below 40 m wherever the measured accumulation is under 0.6 m w.e. per
  year. This adds Dôme du Goûter (about 0.5 m w.e. per year measured) and any other low-accumulation summit with deep
  readings. I would list them before the fit.

**My recommendation: (b).** It is the same physics everywhere and answers the panel's point. The cost is a few deep
readings at Dôme du Goûter.

## 7. Readings above 0 °C

- **Today:** they are dropped when the data are read.
- **The plan:** set them to 0 °C and flag them. Firn and ice cannot be warmer than 0 °C, so such a reading is a
  temperate reading with a measurement error. Dropping it removes evidence that the ice is temperate.
- **Cost:** one change in the data reading, off by default and switched on for Step 5; the count is reported.

## 8. The benchmark: each massif predicted from the other six

```
  7 massifs:  Mont Blanc | Monte Rosa | Mischabel-Saas | Valais west | Bernese Alps | Bernina | Ortles
  each in turn:  fit the error model, f_ref and f_ins on the other six  ──►  a real GloGEM run at that result
                 choose f_ref on the grid on the same six              ──►  the grid run at that f_ref
                 fit each simple rule on the same six                  ──►  predict the left-out massif
                 then compare, at depth, overall and in the zone 3000–3800 m
```

- **What Step 5 is compared with:**
  - **the grid choice on the same physics.** f_ref chosen from 0.03, 0.10, 0.25, 0.40, 0.55, 0.70, 0.85 and 1.00, f_ins
    at its starting value, all on the free air with f_sur 1.4 and five parts. That is **40 runs**; the old 2.51 °C grid
    choice was on the 2 m air, so it is no longer the right yardstick.
  - **the model at its starting values** (f_ref 0.55, f_ins 1.0): **5 runs**.
  - **five simple rules,** fitted on the same six massifs:
    - a straight line in elevation;
    - a straight line in the band's mean air temperature (the free air);
    - the 5 bands nearest in elevation;
    - an elevation threshold for the classes;
    - Suter's limit of cold firn (north and south faces).
- **Scored:**
  - the error at depth, overall and in the zone;
  - the regime classes counted by site;
  - frozen beds, those close to melting apart from the very cold ones;
  - whether the stated 90 % uncertainty ranges contain 90 % of the readings.
- **The uncertainty of every comparison:** resampling whole sites.
- **The 45 extra runs:** on vierzack03 and 07 tonight, after the design chain, about 1.5 hours. They are needed before
  the benchmark, not before the fit.

## 9. Adoption (as planned)

- **Step 5 is adopted** if, on the left-out massifs, it is at worst 0.1 °C worse than the grid choice (with 90 %
  confidence), overall and in the zone, and its uncertainty ranges are about right (the 90 % ranges hold 80–97 % of
  the readings).
- **If not:** v1.0 keeps the grid choice with its spread, and every table is published.

## 10. The map's uncertainty (as planned)

- **For a glacier without boreholes,** the uncertainty is the larger of two:
  - the calibration's own spread;
  - the spread between the left-out-massif results.
- **The map's probabilities** come from drawing the parameters and the error model together. With C6b that includes
  a site's own error size.
- **The map states** that a probability refers to the band average.

## Order of work

```
  now                  after your OK              tonight                   tomorrow
  this document  ─►    code: heavy-tailed    ─►   45 benchmark runs    ─►   commit all rules  ─►  the fit
                       likelihood (+ tests,       (vierzack03/07)            (one commit, its
                       max-effort review);                                   hash recorded)
                       site rule; site checks;
                       variants; benchmark;
                       readings above 0 C
```

## Glossary

| term | meaning |
|---|---|
| error model | how far the model may normally miss a reading, and which readings tend to miss together |
| AIC | Akaike's information criterion: how well a form explains the misses, with a penalty for every extra term; lower is better |
| check (misfit per reading) | the misses measured in the form's own units of "expected miss"; about 1 if the form is right |
| heavy tails (C6b, nu) | some sites miss much more than others; nu says how much (small nu = large differences between sites) |
| w2 | a reading's own small error (its exact depth and month), here only in the top 15 m |
| site | holes within 1 km and 150 m of each other, chained |
| offset of a site | the miss the whole site shares |
| fold | one round of the benchmark: one massif left out, everything fitted on the other six |
| grid choice | f_ref picked from eight fixed values by the same leave-one-massif-out rule, f_ins at its start |
| zone | 3000–3800 m, the firn–no-firn transition |
| m w.e. | metres of water equivalent |
| f_ref, f_ins | share of the meltwater that refreezes and heats; scale of the snow's insulation |
| f_sur | how much colder the snow surface is than the air (1.4 K, fixed) |
| FRHS | the central thermal variant: no T_off in the thermal inputs, ERA5's lapse rate, Herron & Langway meltwater, snow refreezing |
