# Step 3b: the meltwater may go only as deep as the firn really is

*Plan of 2026-10-14 (evening), Step 3b of the approved plan of this afternoon (archived as
`glogemflow_icetemp/development_history/plans/2026-10-14_snow_firn_probabilistic_plan.md`). The
rest of that plan stands. Today's evidence is in book chapter 9 (icetemp `2100a8d`).*

## Context: what is wrong today, in plain words

Every month with melt, the thermal model lets the meltwater (times f_ref) warm the firn from the top
down until the water is used up. It may go as deep as a **fixed depth per band**: the Herron &
Langway pore close-off depth from the band's mean accumulation and temperature. That depth is
75–90 m at the zone sites, computed once and never changed.

Today's firn-model comparison shows that real firn is far thinner than that limit, and changes
with time:

| zone bands | firn in the firn model (GLAMOS) | thermal model: firn / water reaches |
|---|---|---|
| Corvatsch 3295–3365 m | none 1982–2013, ≤ 3 m otherwise | firn until 2023–25 / 13 m |
| Hohsaas 3275–3345 m | ≤ 17 m in the 1980s, gone by 1998–2003 | firn until about 2022 / 36 m |
| Jungfraufirn 3415–3535 m | about 48 m | firn / 68 m |

So the water heats ice that, in reality, has no pores: it is too deep, and for too many years.

```
 today                                          Step 3b
 surface ─ meltwater ─┐                         surface ─ meltwater ─┐
 snow  ░░░░░░░░░░░░░  │                         snow  ░░░░░░░░░░░░░  │
 firn  ▒▒▒▒▒▒▒▒▒▒▒▒▒  │  water warms             firn  ▒▒▒▒▒▒▒▒▒▒▒▒▒  │  water warms
 ice   ███████████████ │  through here            firn base ──────────┘  stops here: the rest runs off
 ice   ███████████████ │                          ice   ███████████████   stays as cold as conduction
 75-90 m (close-off)  ─┘  stops here              ice   ███████████████   and advection make it
```

## The idea: a firn ledger per band, kept by the mass balance

Each elevation band keeps a ledger of its annual firn layers. The rules are those of the firn model
of your first paper (Beer et al. 2026, *The Cryosphere* 20; Matthias's code).

```
 at the end of each hydrological year, per band (b = the year's band balance, m w.e.):
   b > 0   ->  a new layer of thickness b, stamped with the year
   b <= 0  ->  |b| is taken from the newest layers first (snow and firn melt from the top)
   layers older than 20 years become ice and leave the ledger (Ambach & Eisner 1966)
 firn thickness = sum of the layers (m w.e.), converted to metres with the model's firn density profile
 each month: meltwater reaches at most  min( firn thickness + this winter's snow ,  close-off depth )
```

**Consistency with GloGEM's mass balance (your question).** The thermal model and the mass balance
see the same snow and the same firn extent. The ledger only adds the one thing the mass balance
does not have: how thick that firn is.

```
 GloGEM mass balance (calibrated, unchanged)          thermal model
   band balance b each year  ──────────────────────►  the ledger (layers, ages) -> firn THICKNESS
   firn mask: 5-year mean b > 0  ──────────────────►  which bands are firn (as today)
   seasonal snow sno  ─────────────────────────────►  insulation, and the top of the water's path
   nothing comes back: the mass balance never reads the ledger
```

**What the thermal model takes from the mass balance today, and what the ledger changes.** All
of this is from the code, `firnice_temperature_model.pro`:

| what | today | with the ledger (3b) |
|---|---|---|
| **which bands are firn** | GloGEM's mask: 5-year mean band balance > 0 (`finalize_annual_massbalance.pro`) | same |
| **seasonal snow on ice bands: insulation** | GloGEM's snow `sno` fills the top node slices with snow conductivity (Calonne) × f_ins, so winter cold is damped | same |
| **seasonal snow on firn bands** | not separate: snow and firn share one density profile, so the snow adds no extra insulation there | same |
| **the firn column: density, conduction, heat capacity** | one fixed reference profile `fit_dens` for every firn band, reaching 830 kg/m³ at about 35–40 m and ice below; conduction × f_ins | same (see the remaining mismatch below) |
| **how much meltwater heat enters** | the month's melt + rain × f_ref, on firn bands only; on ice bands it runs off | same |
| **how deep the meltwater may go** | to the Herron–Langway close-off depth, fixed per band (75–90 m in the zone) | **to the ledger's firn thickness + this winter's snow, at most the close-off depth; it changes every year** |
| **surface temperature, advection, geothermal flux** | air temperature capped at 0 °C; advection; bed flux | same |

So only **one** line changes: how deep the meltwater may warm the column. Insulation, refreezing
amount, the firn mask and the snow stay exactly as they are.

**Why that line matters.** Today it is inconsistent even inside the thermal model:
- the conduction treats a firn band's column as firn down to about 35–40 m and as ice below;
- but the meltwater may warm it down to 75–90 m, inside what the model itself treats as ice.

The ledger brings the water limit back to the firn that the mass balance actually built.

**The mismatch that remains.** The firn density profile used for conduction is still one fixed
profile, whatever the ledger says:
- a band the ledger gives 3 m of firn still conducts as if it had 35–40 m;
- a band with 48 m of firn (Jungfraufirn) conducts as if it had 35–40 m.

Making the density profile follow the ledger is the natural next step. It is kept as a separate
switch (3b-ii), so that each change can be judged on its own.

- **Same input.** The ledger is driven only by GloGEM's own calibrated band balance (`bal`, the
  values stored in `baly`), the same numbers the mass balance writes. No observation enters at run
  time. The firn model on GLAMOS data is used only to check.
- **Same extent.** Which bands count as firn stays GloGEM's mask (5-year mean band balance > 0),
  and the snow is GloGEM's `sno`. The ledger can only make the water stop shallower, inside bands
  the mass balance already calls firn. It can never create firn where the mass balance has ice.
- **The mismatch removed is today's.** At present the depth comes from Herron & Langway: the steady
  firn column of a band's mean accumulation, 75–90 m. It ignores that the mass balance has had
  negative years there for decades. So today the thermal model already assumes far more firn than
  GloGEM's own balance has built. The ledger makes the depth follow that balance year by year.
- **No feedback.** The mass balance does not read the ledger, so its results stay byte-identical
  (gate check 2). That is also why the mask stays the mass balance's for now (decision 4). A
  ledger-based mask would have to be used by the mass balance too (the firn/ice melt factors),
  or the two would disagree.
- **Where the mass balance itself is wrong** (Corvatsch: too positive, so the firn is too
  extensive), both models are wrong together. That is fixed in the mass balance (Step 3a), and
  both models follow.

**Why this is defensible where the failed depth rule was not:**
- The ledger depends only on the mass balance, never on the temperatures it changes. So it cannot
  have the two stable states that sank the 'front' rule.
- It responds to the climate. Negative years eat the firn, so firn shrinks as the climate warms
  and grows back in good years. That is what a projection needs.
- No new free parameter: the 20-year ageing is the firn model's published constant.
- It runs wherever GloGEM runs: French and Italian glaciers, the projections. It is checked against
  your firn model on the 9 Swiss glaciers it covers.

**What it will and will not fix:**
- **It fixes the depth.** Water stops at the firn base instead of 75–90 m: Hohsaas and Jungfraufirn
  should keep cold ice below their firn.
- **It cannot fix the extent where the mass balance is wrong.** At Corvatsch the model's band
  balance is too positive (+0.3 against GLAMOS −0.6 m w.e. a year), so the ledger still holds firn,
  roughly 20 × 0.3 m w.e. That is Step 3a's job: the snowline recalibration, and possibly the
  targets of small glaciers. The two steps are complementary.

## Your decisions before I code (my recommendations)

1. **How a negative year removes firn.**
   - **Recommended: mass-conserving.** Exactly |b| is removed, newest layers first.
   - **Why not copy the firn model exactly:** while porting it, the agent found that the code
     removes more than |b| when several thin layers sit on top. Each layer thinner than |b| is wiped
     without reducing |b|, and then the full |b| is taken again from the next, thicker layer. That
     makes the firn disappear faster than the mass balance says. It may be worth a look for the
     published firn maps (your and Matthias's call).
   - I will quantify the difference between the two removal rules on the 9 glaciers offline.
2. **Ageing:** 20 years, as published. 15 and 30 years are run later as a sensitivity.
3. **Starting state in 1971:** the steady state of the spin-up's reference balance (20 layers of
   b_ref where b_ref > 0, none where it is negative), so the run does not start from zero firn.
   The firn model starts from zero in 1970, which biases its 1970s low.
4. **Which bands count as firn** stays the mass balance's mask (5-year mean balance > 0) for now.
   The ledger only limits the depth. A ledger-based mask is a later option, judged on its own.

## Implementation (GloGEM worktree from `GloGEMflow` `161dcdd`, branch `firn-ledger`; default unchanged)

- **`procedures/initialise/settings.pro`:** `firnice_perm_limit = 'firn'` as a new allowed value
  (the guard); `firnice_firn_age = 20` (years).
- **`procedures/initialise/initialise_firnicetemp_spinup.pro`:** the ledger arrays `fl_we[nb, 20]`
  and `fl_yr[nb, 20]`, filled with the steady state from `b_ref`.
  (`b_ref` is from `firnice_reference_forcing.pro`, already used for vertical advection.)
- **`procedures/processing/firnice_firn_ledger.pro` (new include):** the yearly update from `bal`
  (the band balance, as stored in `baly[ye,*]`). It is called from
  `procedures/processing/finalize_annual_massbalance.pro`, next to the existing firn-mask update
  (lines 30–40). Rules as above.
- **`procedures/processing/firnice_temperature_model.pro`:** where `perm_limit` is set from
  `z_perm_b` (the close-off depth × `perm_frac`), under 'firn' it becomes
  `z_perm_b < depth(ledger + sno)`. The depth comes from the thermal model's density profile
  `fit_dens`, by cumulative water equivalent. The spin-up uses the steady-state ledger.
- **A log `firnice_firnledger_<id>.dat`** per band and year: firn m w.e., depth in m, number of
  layers, oldest layer's year. Used for the checks and the book.
- **icetemp:**
  - the variant `FRHL` = FRH + `perm_limit='firn'` (`full_set_runs.py`); the template comment;
  - a Python reference of the ledger (`scripts/icetemp_analysis/firn_ledger.py`), reusing the
    agent's validated port of `firn_change.pro`.

## Checks before the run set (gate)

1. Default 'hl': byte-identical with the fair campaign (FRH at f_ref 0.85, 373 files).
2. Mass-balance-derived files identical.
3. **The IDL ledger equals the Python reference**, fed with the same band balances (the full-output
   mass-balance run `week2/mbout/regional`), to 1e-6 m w.e.
4. In every band and month, the wet depth is at most the ledger depth.
5. **Against the firn model on the 9 Swiss glaciers:** the ledger driven by GLAMOS balances must
   reproduce the firn model's firn thickness, up to the removal rule (decision 1). Driven by
   GloGEM's balances, it shows how the mass-balance errors carry over (reported).

## Run set and decision rule (fixed before scoring; scorer committed first)

- **Runs:** FRHL at the 7 f_ref values plus FRHL at 0.85 without deepening (8 runs, about 40 min),
  against the fair FRH runs.
- **'firn' replaces 'hl' if all hold:**
  - it gains, net, at least 2 robust zone entities (zone target D needs 2 against H);
  - no condition of rule A breaks against FRH (M1, robust class losses, warming);
  - no firn band that 'hl' keeps temperate at 10 m below the cold-firn limit is cooled by more
    than 0.2 °C (rule B's guard);
  - the frozen-bed table does not lose (12 of 14).
- **Reported:**
  - the wet depth against the firn model's thickness, and the firn line by decade;
  - Jungfraufirn, Lys, Fiescher and Breithorn, the temperate sites;
  - the snow-refreezing switch judged again on top (FRHS + 'firn').
- **Effort:** main session; one independent reviewer at max on the code before the gate.
- **Time:** about 2–3 days, before or around Alexandra's branch.

## When do the recalibrations happen?

There are two different recalibrations:

| what | when | why then |
|---|---|---|
| **mass balance with snow lines** (Step 3a: Alexandra's method) | from Friday: port her branch into a procedure on GloGEMflow (gate) → recalibrate the calibration glaciers (diagnostic) → decide whether the thermal runs use it. Roughly next week | it needs her code |
| **thermal Bayesian calibration** (f_ref, f_ins, f_adv; Step 5) | after 3a and 3b are decided: Alpine ensemble with the adopted model → T2 Bayesian calibration → Z2–Z4 measured once → v1.0 freeze | every structural change (forcing, snow, firn) moves the posterior, so it runs once, on the final model |

Step 3b's own run set is not a recalibration. It uses the prior-mean f_ins and f_adv, and f_ref
chosen out of sample, as in the fair test. After 3a, 3b is run once more on the recalibrated mass
balance.

## Verification

- The gate checks 1–5 above. `pytest` for icetemp (the training-config test extended to `'firn'`).
- A known-answer test of the Python ledger on a hand-made balance series: positive years, a
  negative year removing exactly |b|, ageing at 20 years.
- The scorer's known-answer check: FRH against itself gives no change.
- Book chapter 9, memory and commits as each step lands. No config files committed.

## Corrections after the code review (2026-10-14, before the run set)

- **The ledger keeps 21 layers, not 20:** layers stay while their year is at least the current
  year − 20, the convention of `firn_change.pro`. The code and the Python reference agree on this.
- **Corvatsch's ledger is thick, not thin.** GloGEM's band balance at 3295–3365 m averages +0.9 to
  +1.2 m w.e. a year, so the ledger holds 26–36 m w.e. until 1995 and 8–14 m w.e. in 2025. The
  column there is only 14–19 m deep, so the ledger changes nothing at Corvatsch. Its problem is the
  mass balance (Step 3a; the target is about 1 m w.e. a year too positive).
- **The snow term lags within the year.** In October–December, snow that survived the year is
  counted both in the new layer and in `sno`, until the snow store is reset. In negative years, the
  firn melted that summer stays in the path until the year ends. Both effects are small (little
  melt in those months; at most |b|) and are accepted. The alternative, the year-to-date balance,
  would need a monthly reference balance in the spin-up.
- **The floor.** With less than about 0.6 m w.e. of ledger plus snow, water still reaches node 1
  (1 m), the same floor as 'hl'.
- **The gate checks the logged limit.** The ledger log now carries the year's balance (full
  precision) and the deepest limit used (GloGEM `20ba803`). Check 4 compares the wet depth with
  that limit, not with the ledger alone, which flagged the snow.
- **Review result (f4935bf):**
  - no code error;
  - the IDL loop equals the Python reference on 905 915 random year-steps;
  - defaults byte-identical (373 files); the mass-balance logs identical;
  - 'firn' is never deeper than 'hl'. It is shallower in 28.5 % of band-years.
  - Gate-run wet depths, 'hl' → 'firn': Jungfraufirn 59–79 → 29–39 m by 2025; Hohsaas 29–44 → 14 m
    by 2022.
