# Step 5, step 3: the rules that decide the result, fixed before any fit

**Approved by Janosch on 28 Oct: all rules as written, and option (b) for very deep, very old ice (rule 6). The 40
design runs finished on 27 Oct at 18:05 and have not been looked at; the 45 benchmark runs started on 28 Oct at
11:48.**

**The rules commit is 3467851 (28 Oct, about 15:50, pushed).** It holds every script that decides the result. The fit
starts only from it, or from a later commit that changes no rules code (this line was added in such a commit).

## The result (28 Oct evening): not adopted by rule 9; two decisions for you

*Everything below ran automatically overnight, on the rules commit (the run itself used `cb50eee`, which differs from
`3467851` only outside the rules code). Step 5 will be redone on the free-air mass balance (decision log, 28 Oct), so
this run also fixes how that redo is read.*

```
 60 design runs ─► the fit ─► site checks ─► 7 left-out-massif ─► 8 real GloGEM runs ─► the benchmark
 (40 + 20 near    on all 7    and variants    fits (one per        (the 7 fold results    (each massif predicted
  the result)     massifs     (rules 3, 5)    massif left out)     and the main result)   from the other six)
```

**1. The fit (all 7 massifs).**
- **f_ref 0.93 ± 0.02, f_ins 1.18 ± 0.06.** The two are told apart (correlation +0.50; the data narrow f_ins to 0.18
  of its starting spread), so rule 4 keeps both.
- **The error model fits:** the misses average 1.15 expected misses (shallow 1.29, deep 1.07).
- **The gate passes.** The emulator check fails (reported only): with the emulator's own uncertainty, the result could
  move by about 1 sd in f_ref and 2 sd in f_ins. The real GloGEM runs at the left-out results agree with the emulator
  (2.60 against 2.58 °C), so the benchmark does not rest on it.

**2. The benchmark (rule 8): each massif predicted from the other six.**

| method | error at depth, all | zone 3000–3800 m | above 3800 m | bands in the right class (of 105) | polythermal bands right (of 18) | sites right (of 22) |
|---|---|---|---|---|---|---|
| **Step 5** | **2.60 °C** | **2.50** | **3.07** | 73 | 9 | 9 |
| grid choice | 2.83 | 2.51 | 3.43 | 72 | 9 | 8 |
| starting values (f_ref 0.55, f_ins 1.0) | 2.45 | 2.20 | 2.86 | 75 | 9 | 10 |
| the 5 bands nearest in elevation | 2.88 | 2.45 | 3.47 | 76 | 0 | 12 |
| a line in elevation | 3.68 | 2.25 | 4.46 | 72 | 0 | 9 |
| a line in the free air | 3.78 | 2.36 | 4.58 | 72 | 0 | 9 |
| elevation thresholds | – | – | – | 74 | 0 | 11 |
| Suter's cold-firn limit | – | – | – | 60 | 0 | 6 |

- **Step 5 beats the grid choice** by 0.23 °C overall (90 % range: 0.42 better to 0.02 worse) and equals it in the
  zone.
- **The starting values do 0.14 °C better than Step 5,** but that is within the uncertainty (0.22 better to 0.58
  worse, resampling whole sites). They win at Mont Blanc and the Bernese Alps and lose at Bernina.
- **Above 3800 m the model clearly beats the simple rules** (2.9–3.4 against 3.5–4.6 °C). In the zone a line in
  elevation beats Step 5 (2.25 against 2.50 °C) and matches the starting values, as found on 21 Oct.
- **No simple rule gets a single polythermal band right.** The 5-nearest rule gets the most bands right by calling
  almost everything cold.
- **Frozen beds, the same for every model version:**
  - 11 of the 14 measured frozen beds found, among them 11 of the 12 very cold ones;
  - neither of the 2 close to melting (Tête Rousse at 65 m, Corvatsch);
  - one temperate bed called frozen (Hohsaas at 21.5 m).

**3. Rule 9: not adopted.**

| check | needed | Step 5 | |
|---|---|---|---|
| error, all | at worst 0.1 °C worse than the grid choice (90 % confidence) | 0.03 °C *better* | passes |
| error, zone | the same | 0.004 °C better | passes |
| the 90 % ranges | hold 80–97 % of the readings | **75 %** | **fails** |

- **By the rule, v1.0 keeps the grid choice with its spread:** f_ref 0.85 at f_ins 1.017 (the best grid value on all
  7 massifs; 0.70 is nearly as good, 2.34 against 2.32 °C), spread 0.16 between the seven fold choices.
- **Why the ranges are too narrow:** the misses belong to whole sites, not to single readings.

  | site | the model is | readings inside their range |
  |---|---|---|
  | Taconnaz 3680–3830 m | 6.2 °C too warm | 3 % |
  | Bossons 3930–4304 m | 3.1 °C too warm | 28 % |
  | Fiescher 3889–3893 m | 2.7 °C too warm | 38 % |
  | Breithornplateau | 4.0 °C too cold | 12 % |

- **Each left-out fit claims to know f_ref within ±0.02–0.04, but the fits disagree far more.** Six give 0.92–0.94; the
  one without Monte Rosa gives 0.40. The fit follows whichever of the two big massifs it contains: the Monte Rosa
  saddles want f_ref near 1.0, Mont Blanc about 0.25. This is the massif split again.

**4. A check after the fact, not a rule** (`step5/coverage_spread.py`): the same ranges with the parameters spread as
for the map.

| spread of f_ref and f_ins | readings inside their 90 % range | model too warm / too cold |
|---|---|---|
| each fold's own (±0.02 / ±0.06): rule 9 as coded | 75 % | 18 % / 7 % |
| the sd of the 7 fold results (±0.20 / ±0.20) | 82 % | 13 % / 5 % |
| the map's spread of rule 10 (±0.20 / ±0.37) | 82 % | 14 % / 5 % |
| one massif's own spread, (n − 1) × sd (±1.2 / ±1.2) | 92 % | 4 % / 4 % |

The calibration's own uncertainty is too sure of itself for an unseen massif. Add the disagreement between the
massifs, as rule 10 already does for the map, and the ranges become about right. Only one massif's own spread gives
90 % ranges that hold about 90 %, missing equally often on both sides.

**5. The site checks (rule 3).**
- **Check 1** (misses within the error model): every site except Grenzgletscher at 2015 m (the tongue, 10 readings).
- **Check 2** (unusual offset): none.
- **Check 3** (leaving the site out moves the result by more than its uncertainty): three sites are decisive.
  - Breithornplateau: f_ref +0.05, f_ins +0.09.
  - Bossons 3930–4304 m: f_ins −0.07.
  - The Colle Gnifetti saddles (Grenzgletscher 3898–4485 m): f_ref −0.09, f_ins −0.37.

  They widen the map's f_ins uncertainty to ±0.37. They do not stop adoption.
- **The named check** (without the 1999 campaign's 14 shallow bands) moves f_ref −0.015 and f_ins +0.045. It is not
  decisive.

**6. The variants (rule 5): f_ref stays at 0.84–0.93 in every variant; f_ins moves between 0.91 and 1.49.**

| variant | f_ref | f_ins |
|---|---|---|
| main | 0.93 | 1.18 |
| deep readings only | 0.88 | 0.96 |
| no site above 20 % of the information | 0.84 | 0.91 |
| holes in equilibrium and published only | 0.91 | 1.49 |
| f_ins prior twice as wide | 0.93 | 1.18 |
| without the 1999 shallow bands | 0.91 | 1.22 |

**Decisions for you.**

| # | question | options | recommendation |
|---|---|---|---|
| A | Adoption | (a) **follow rule 9**: v1.0 takes the grid choice (f_ref 0.85, f_ins 1.017) with its spread; Step 5's result, its gain in error and the check after the fact are reported beside it. (b) adopt Step 5 anyway, because with the map's spread the ranges hold 82 %; a deviation from the rule, logged like S and orientation | **(a)**: the rule was fixed to catch exactly an over-confident calibration, and it did. The two values of f_ref lie close (0.85 against 0.93), well inside the map's spread of ±0.20, and the paper can say the rule set in advance was followed |
| B | Rule 10's spread for glaciers without boreholes | (a) **one massif's own spread**, (n − 1) × sd; (b) the sd of the fold results | **(a)**: the only spread whose 90 % ranges hold about 90 % (92 %), missing equally often on both sides. Its cost: f_ref for an unseen massif is then anywhere in 0.03–1.0, so the map's probabilities stay wide in the zone. That is the honest state of the data |

**Two changes on 28 Oct, decided by Janosch after the check of the data (before any fit):**
- **Rule 6 (b) reads "0.6 m w.e. per year or less".**
  - The check showed that Dôme du Goûter has no readings below 40 m, so (b) does not add it. I had claimed it would,
    and that was wrong.
  - One Colle Gnifetti band, Grenzgletscher 4465 m on outline 01225, measured exactly 0.6 (range 0.2–1.4), with
    readings to 83 m in ice about a century old. "Under 0.6" would have put its 15 deep points back into the fit.
  - With "or less", (b) leaves out the same points as today's rule (85 band depths, holding 107 pooled readings). It
    differs from it only for future data.
- **Rule 7 is withdrawn.** A reading above 0 °C measures air or water in the hole, not ice, so it cannot go into the
  calibration. Such readings stay dropped, as before.
  - In the calibration holes there are 148,773 of them, mostly at Ortles, Hohsaas, Alphubel and Tortin. Their median
    is +0.1 °C, but some reach +19 to +34 °C.

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

## How the rules are coded (28 Oct, fixed before the fit)

The approved rules leave a few details open. They are fixed here, in the code and in this section, before anything is
fitted. Scripts: `scripts/icetemp_analysis/step5/` (fit_step5.py, site_checks.py, variants.py, benchmark.py,
stage2_design.py); the error model in `src/icetemp/calibration/error_model.py`.

1. **Error sizes (rule 1).**
   - The fit starts from the C6b sizes and fits f_ref and f_ins.
   - At that first result the sizes are re-estimated, and the parameters fitted again.
   - The sizes are estimated as the likelihood uses them: with no free mean, and with the readings' own errors fixed.
   - A second re-estimate replaces them only if it explains the misses better by more than 2 log-units. Otherwise the
     second fit stands.
   - Each left-out-massif fold repeats the same steps on its six massifs.
2. **Unusual site offset (rule 3, check 2).** "More than 2.9 times its uncertainty" is coded as its stated meaning: a
   chance below 1 in 300, two-sided, under the error model's own law. For the heavy-tailed C6b that law is a t
   distribution, so the line lies further out than 2.9: at about 5.1 times its uncertainty for nu 2.6. Under C6b, 2.9
   times has a chance of 1 in 74 per site. For a Gaussian the two are the same. **Janosch chose the chance, 28 Oct.**
3. **Leaving a site out (rule 3, check 3).** All points of the site's bands are left out and the posterior recomputed.
   The site is decisive if f_ref or f_ins moves by more than the main result's posterior spread (one sd).
4. **f_ref and f_ins apart (rule 4).** By Step 5a's existing rule, they are told apart if they correlate at 0.8 or
   less and the data narrow f_ins to at most half of its starting spread. If not, f_ins is fixed at 1.0, and f_ref's
   spread from the joint fit, which contains the trade-off, is reported as its uncertainty.
5. **The variants (rule 5).**
   - **"No site above 20 % of the information":** each site's information on f_ref is computed at the result. Sites
     above 20 % get a weight that brings them down to 20 %.
   - **"Published":** a hole whose glenglat source is not a personal communication. "In equilibrium" means bands rated
     'true'.
   - **"Wider f_ins range":** the starting spread is doubled, from 0.35 to 0.7.
6. **The 1999 bands (rule 3 named check, rule 5).** Bands whose every reading is from 1999, none deeper than 25 m: 9
   Monte Rosa saddle bands and 5 Lys bands.
7. **The benchmark (rule 8).**
   - **What is scored:** the error at depth is M1, as for every earlier grid choice: the root-mean-square error below
     15 m, each outline counting equally.
   - **The simple rules:**
     - the elevation line and the air-temperature line are fitted by least squares, each outline counting equally;
     - the air temperature is the band's mean free air from October 1971 to September 2025, at the nearest grid
       point;
     - "5 nearest bands" means the mean of the five bands of the other massifs closest in elevation;
     - the elevation threshold is two heights (temperate below, polythermal between, cold above), chosen to classify
       the most bands right;
     - Suter's limit is 3400 m on north faces and 4150 m on south faces (Suter et al. 2001), and 3775 m on east,
       west and gentle parts.
   - **Classes per site:** a site counts as right if more than half its bands are.
   - **Frozen beds:** "close to melting" means measured at −2 °C or warmer.
   - **The 90 % ranges** for Step 5 come from each fold's posterior and its error model, for a site the fold has not
     seen.
   - **"With 90 % confidence":** the 90th percentile of the difference to the grid choice, over 4000 resamplings of
     whole sites, must be at most 0.1 °C.
8. **The 20 further runs (step 4.1).** The box is centred on the result. On each side it reaches the larger of 2.5
   posterior sd and the farthest left-out-massif result.
9. **The emulator's own uncertainty (the agreement check).** A per-site t cannot carry the coupling of the
   emulator's uncertainty between sites. Leaving the coupling out would understate it and make the check lean
   towards "agree" (the review of 28 Oct). So the check compares the Gaussian form with the same sizes, with and
   without the emulator's uncertainty, where that coupling is computed exactly. The fit itself uses the emulator mean
   only.
10. **The overall misfit check of the gate.** A correct heavy-tailed model gives a misfit per reading below 0.5 about
    a third of the time, so the old window of 0.5–2 would often fail it. For the t forms the window is therefore the
    99 % range of the form's own law (by simulation). The 0.5–2 window stays for the Gaussian forms.

**Details fixed after an independent review of the scripts (28 Oct afternoon, before the fit).** The review found
four errors in the code, all corrected and tested on synthetic data before the commit. In the order they mattered:
- the site resampling scored two sites of one outline as two outlines;
- the "published" variant also dropped published readings that shared a depth with an unpublished one;
- each left-out-massif fit started from error sizes that had seen the left-out massif;
- check 3 had no consequence in the code.

The details below were open in the words and are fixed here.

11. **f_ref and f_ins apart (rule 4) in detail.**
    - "Half of its starting spread" is half of the starting range's own sd, 0.326, so f_ins must be fixed to 0.163
      or better. (The nominal 0.35 is cut at 0.3 and 2.0.)
    - It is judged on the fit of item 1, which uses the emulator mean only, as the fit itself does.
    - Each left-out-massif fit applies the rule itself, on its six massifs.
    - When f_ins is fixed, the error sizes are not estimated again.
12. **One uncertainty everywhere.** "The result's own uncertainty" in check 3, in the box of the 20 further runs and
    in the map is the joint fit's sd. That sd holds the f_ref–f_ins trade-off, also when f_ins is fixed (rule 4).
13. **If rule 4 fixes f_ins:**
    - The 20 further runs still vary f_ins around 1.0, so that the refit on all 60 runs decides rule 4 again.
    - The wider-f_ins variant compares two fits with f_ins free, one per starting range, so it shows the starting
      range's effect alone.
    - Doubling the nominal spread (0.35 → 0.7) raises the range's own sd from 0.326 to 0.442.
14. **Check 3's consequence (rule 3).**
    - A decisive site does not stop adoption.
    - The map's uncertainty then covers both results: per parameter, at least the site's shift and the sd of the fit
      without it.
    - The named 1999 check is reported only.
15. **The left-out-massif fits:**
    - they start from the C6b sizes;
    - they train on all other bands, including the 8 of Ecrins, Oetztal and Titlis, which have no readings below 15 m
      and so are never scored.
16. **The site resampling.**
    - The second copy of a site goes with the second copies of the other sites on its outline. So drawing every site
      once gives M1 exactly, and copies count as separate outlines.
    - The grid choice and the left-out-massif results stay as they are in each resample; they are not chosen anew.
17. **The box of the 20 further runs** is the same on both sides: the farthest left-out-massif result, on either
    side, sets both.
18. **The yardsticks:**
    - the grid runs are at f_ins 1.017, the value of every earlier grid;
    - the starting-values run is at f_ins 1.0.
19. **"Published" is judged per hole.**
    - Two Grenzgletscher holes with a published source but profiles from a personal communication stay in.
    - The variant leaves out single readings, not whole depths that a published hole shares.
20. **What decides adoption.**
    - The gate (convergence, design runs near the result, the misfit window, the gap between the best point and the
      mean) and the emulator agreement are reported beside the decision.
    - The decision itself is rule 9. It scores real GloGEM runs at the left-out-massif results, so it does not rest
      on the emulator.
21. **Rule 10, "the spread between the left-out-massif results", is the sd of the 7 results, as approved.**
    - Each result comes from six of seven massifs, so for a mean, one massif's own value lies about 6 (n − 1) times
      further from the overall result than its left-out result does.
    - That sd therefore understates how much single massifs disagree.
    - The six-fold value is reported beside it. Which of the two the map uses is for Janosch to decide.
22. **The scripts refuse to run:**
    - on uncommitted rules code, data or config;
    - on a run whose parts did not run the expected point and settings;
    - on a missing prediction, or one at the wrong depth;
    - on a fit or left-out-massif fits computed on another set of design runs.

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
| M1 | the error at depth: root-mean-square miss below 15 m, each outline counting equally |
| outline | one glacier polygon of the inventory (RGI); several sites can lie on one outline |
| sd | standard deviation: the usual size of a spread |
| emulator | a fast statistical stand-in for GloGEM, trained on the design runs |
| gate | the fit's own health checks (convergence, design runs near the result, misfit window) |
| 90 % range, coverage | the range a method gives for a reading, which should contain the reading 9 times in 10; coverage is how often it does |
| decisive site | a site whose removal moves the result by more than the result's own uncertainty (check 3) |
| starting values | f_ref 0.55, f_ins 1.0: the values before any calibration |
| check after the fact | computed after the results were seen; reported, never used to decide |
| one massif's own spread | (n − 1) × the sd of the fold results: about how far one massif's own best value lies from the overall result |
| emulator check | whether the fit stays put when the emulator's own uncertainty is included |
