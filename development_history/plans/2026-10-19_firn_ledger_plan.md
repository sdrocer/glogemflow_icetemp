# The firn ledger: would it pay off with the right firn lines?

*Plan of 2026-10-19 for your sign-off, revised after your comment. Terms are in the glossary at the
end.*

## Context

- **Your decision today:** the snowline recalibration is not adopted, and the work turns to the
  firn-layer tracking begun on 14 Oct, the **firn ledger**.
- **What the ledger is** (GloGEM branch `firn-ledger`):
  - each elevation band keeps its annual firn layers, built from GloGEM's own band balance;
  - layers older than 20 years become ice;
  - the meltwater may warm the column only as deep as that firn (`firnice_perm_limit = 'firn'`).
- **Where it stands:**
  - scored on 14 Oct, but **not adopted**:
    - only +1 zone entity improved (2 needed);
    - the cooling guard failed on 6 young firn-line bands;
  - on 15 Oct you approved a refined guard (only mature firn, at least 15 layers);
  - on 16 Oct the zone became 3000–3800 m;
  - both changes are written down but not yet in the scorer.
- **The suspicion:** the ledger fails because of its input. At Hohsaas it holds 8–15 m w.e. of firn in
  2020, while your firn model (Beer et al. 2026) has essentially none after 2003. GloGEM's too
  positive balance carries over.
- **Your question (comment on the plan):** if the model had the right snow and firn lines, would the
  ledger then improve the thermal model? If yes, it is worth investing in those lines. If not, it is
  not.
- **Intended outcome:** an answer to that question, with numbers. Then the ledger is rerun on the
  current code with the agreed amendments, and the next step (the firn column follows the ledger) is
  tested. Each is judged by a rule fixed before the runs.

## Step 0. Save what lives only in the session's `/tmp` folder (`/tmp` is wiped at reboot)

- **The scripts are committed** under `scripts/icetemp_analysis/firn/`:
  - the exact Python port of `firn_change.pro`;
  - the per-band comparison tables (firn model, GLAMOS band balance, thermal firn mask and wet depth,
    per glacier, year and RGI7 10 m band);
  - the ledger gate and review checks.
- **Their outputs are copied** to `…/week2/firnmodel/`.
- **Chapter 9's stale "What is next" bullets** of 14 Oct are cleaned up.

## Step 1. Check 5: the ledger against your firn model (rules, then input)

```
                          band balance, per year            rules                 firn per band and year
 firn model (your paper)  GLAMOS distributed grids  ──────► firn_change.pro ─────► reference firn
 5a: the rules            GLAMOS band balances      ──────► the ledger's rules ──► same as the reference?
 5b: the input            GloGEM's band balances    ──────► the ledger's rules ──► how far off, and where
```

**5a, the rules.** The ledger (`firn_ledger.py`) is driven by the GLAMOS band balances and compared
with the firn model band by band (the firn line by year; the thickness in the zone bands by decade).
Each known rule difference is switched on and off to see its size:
- removal: mass-conserving vs the original;
- firn only after 2 years in a row;
- the start: zero in 1970 vs the steady state;
- 20 vs 21 layers.

If 5a agrees, the rules are defensible with your paper. Any difference is a rule we change, or keep
on purpose and state.

**5b, the input.** The same with GloGEM's band balances (from the 14 Oct runs' ledger logs): per
glacier and decade, how far the firn line and thickness are off because of the mass balance.

**Output:** a table and figure (the firn line by year on the 9 glaciers: firn model, GLAMOS-driven
ledger, GloGEM-driven ledger); `firn/ledger_vs_firnmodel.py`.

## Step 2. The perfect-firn test: would the ledger pay off with the right firn lines? (your question)

**The idea:** the same thermal model and the same ledger rules, but on the glaciers where GLAMOS grids
exist, the ledger is fed the observed band balances instead of GloGEM's. Only the firn changes. The
energy (melt, rain, snow, air temperature) stays GloGEM's, so the comparison isolates the firn.

```
                                     the ledger is fed by          the firn mask (which bands are firn)
 FRH     ('hl', today)               -                              GloGEM
 FRHL    (the ledger as it is)       GloGEM's band balance          GloGEM
 FRHL-o  (observed firn thickness)   GLAMOS band balance            GloGEM
 FRHL-O  (observed firn line too)    GLAMOS band balance            GLAMOS (same rule: 5-year mean > 0)
```

**Where it can be tested:**
- 55 of the 105 scored borehole sites lie on the 9 GLAMOS glaciers, 15 of the 27 zone sites:
  Hohsaas 5, Corvatsch 4, Aletsch 4 (Jungfraufirn), Alphubel 2.
- Those four have GLAMOS grids over the whole run period (1970–2025; Aletsch from 1915).
- Fiescher, Gorner and Morteratsch start in 2005–2008; before that their bands keep GloGEM's
  balance (reported, not in the decision).

**Code (GloGEM, a diagnostic switch, never for production; default off, results byte-identical):**
- **The input:** `firnice_firn_obs_dir` points to files with the observed band balance per year and
  RGI7 band. They are written from the band tables by `firn/write_firn_obs.py`.
- **The ledger update** (`firnice_firn_ledger.pro`) takes the observed balance wherever one exists for
  that band and year, otherwise GloGEM's.
- **The start state** where observed: the steady state of the first ten observed years' mean balance.
- **`firnice_firn_obs_mask = 'y'`** (for FRHL-O) also sets the thermal firn mask from the observed
  5-year mean balance, the same rule as GloGEM's own mask. This happens in the reference forcing and
  in `finalize_annual_massbalance.pro`, thermal mask only.
- **The mass balance never reads any of it.**

**What the answer means (fixed before the runs).** All four variants are scored on the 15 zone sites
of the four fully observed glaciers: robust zone classes, and zone M1 on those outlines.

| | FRHL-O gains ≥ 2 robust zone sites over FRH | FRHL-O does not |
|---|---|---|
| **FRHL gains ≥ 2** | the ledger works already; better firn lines are a bonus | (unlikely) |
| **FRHL does not** | **the input is the bottleneck: investing in snow and firn lines pays off** | **the right firn lines would not fix these sites: don't invest there for the thermal model; look elsewhere** |

- **What else it shows:** FRHL-o against FRHL-O tells whether the firn thickness or the firn extent
  matters more.
- **Reported, not decided:** M1 on all 55 sites of the 9 glaciers.
- **Not tested:** the seasonal snow stays GloGEM's (there is no observed seasonal snow per band).
  This test is about the firn line; a snow-line test would need another input.

## Step 3. Rerun the ledger on the current code, with the agreed amendments

**Code (GloGEM):**
1. **The 360-day thermal year (your decision 6).** The thermal step comes from each month's real
   length (thermal model only).
   - One commit on local GloGEMflow.
   - Gate: mass balance identical; a thermal change of the expected size (about 1.4 %).
   - Recorded on the model-changes page (message the AMOC session first).
2. **Rebase `firn-ledger`** onto it (no conflict predicted), plus whatever check 5a decides about the
   rules. The Python reference is updated to match.

**Scorer (icetemp), committed before the runs:**
- **The zone:** 3000–3800 m (decided 16 Oct).
- **The cooling guard:** mature firn only, at least 15 layers (decided 15 Oct).
- **The rest of the 14 Oct rule unchanged.** 'firn' replaces 'hl' only if all hold:
  - a net robust zone gain of at least 2;
  - A's conditions against FRH;
  - the guard;
  - the frozen-bed table.
- **The book says** the 14 Oct runs failed under the old version.

## Step 4. 3b-ii: the firn column follows the ledger

**Today:** every firn band has the same fixed density profile, firn down to about 44 m and ice from
49 m. A band with 3 m of ledger firn still conducts as if it had 44 m.

```
 today (every firn band)              3b-ii (the column follows the ledger)
 0 m  ░ firn 250 kg/m³                0 m  ░ firn
      ▒ firn, denser with depth       8 m  ▒ firn base = ledger + this winter's snow
 35 m ▒ 830 (pores close)                  █ ice 900 kg/m³: conducts and stores heat like ice
 49 m █ ice                                █
```

**Code (new switch, off by default, needs `perm_limit = 'firn'`):**
- **The column split** (`firnice_temperature_model.pro:114-116`): firn density above the ledger's firn
  base, ice below, with the straddling node blended. The base is the same expression that already
  sets the water limit (`:189`), factored into one variable.
- **Firn bands with an empty ledger:** their snow goes through the ice-band snow branch.
- **Everything else follows,** because conduction, heat capacity, vertical velocity and the spin-up
  read the same density array.
- **Limitation:** firn thicker than the profile (above about 36 m w.e.) keeps the full profile.

**Gate:**
- default byte-identical;
- 'firn' without the switch identical to step 3's FRHL;
- a reviewer at max effort before the gate.

**Rule (your decision 2):** adopted unless worse against FRHL:
- M1 not worse by more than 0.1 °C;
- no net loss of robust zone sites;
- the guard holds;
- the frozen-bed table does not lose.

## The run set (one batch, all with the 360-day fix, 7 f_ref values each unless stated)

| set | why |
|---|---|
| H | the reference of zone target D |
| FRH | the baseline |
| FRHL | step 3: the ledger as it is |
| FRHL-o, FRHL-O | step 2: observed firn thickness; observed firn line too |
| FRHL + 3b-ii | step 4 (and FRHL-O + 3b-ii, reported, if the code is ready) |
| FRHL without deepening (1 run); ageing 15 and 30 years (2 runs at f_ref 0.85) | reported |

About 45 runs at about 30 minutes each, six at a time: about 4 hours.

## Not in this plan (later, each on its own)

- **Getting the snow and firn lines right in the mass balance,** if step 2 says it pays off: that is
  the next plan. Candidates:
  - the precipitation gradient with height;
  - the small-glacier targets;
  - a firn-line constraint in the calibration.
- **A firn mask from the ledger inside the mass balance:** it changes melt factors and refreezing, so
  the glaciers would need recalibrating.
- **The thermal Bayesian calibration** (Step 5): once, on the final model.

## Order and effort

1. **Step 0** (about 30 minutes).
2. **Step 1** (about half a day, main session).
3. **Code for steps 2–4,** in this order:
   - the 360-day fix and the rebase;
   - the observed-firn switch;
   - 3b-ii.
   A reviewer at max on the physics code before the gate.
4. **The scorer amendments,** committed before the runs.
5. **The batch** (about 4 hours, watched at low effort).
6. **Scoring** in the main session, then a new book chapter 10 "The firn ledger" (figures, glossary),
   memory, commits, and your decisions.

Time: about 3 days.

## Decisions for you

| # | question | recommendation |
|---|---|---|
| 1 | the 360-day fix goes in now, with this batch (so H and FRH are rerun too) | yes |
| 2 | 3b-ii judged as "adopted unless worse"; the ledger itself keeps the 14 Oct rule (+2 zone sites) | yes |
| 3 | where check 5a finds a rule difference to your firn model, the ledger follows your paper, except the removal (stays mass-conserving, decision of 14 Oct) | yes |
| 4 | the perfect-firn test is read on the 15 zone sites of Hohsaas, Corvatsch, Aletsch and Alphubel, with the threshold of the table above (≥ 2 robust sites) | yes |
| 5 | a new chapter 10, "The firn ledger" | yes |

## Verification

- **Step 0:** the copies are identical; the committed scripts regenerate the band tables.
- **Step 1:** with every rule switched to the firn model's, the GLAMOS-driven ledger matches the firn
  model to rounding. That is the test that the comparison itself is right.
- **Step 2:**
  - with the switch on but an empty directory, the run is identical to FRHL;
  - the logged ledger of an observed band equals the Python ledger fed the same GLAMOS balances;
  - the mass-balance files are identical to FRHL.
- **Step 3:**
  - the 360-day gate;
  - the rebased branch's defaults byte-identical to its own FRH;
  - the IDL ledger equals the Python reference;
  - the scorer amendments committed before the runs.
- **Step 4:** the gate above, and the reviewer's report.

## Glossary

| term | meaning |
|---|---|
| firn ledger | per elevation band, the annual firn layers built from the band balance; layers older than 20 years become ice |
| firn model | Matthias's `firn_change.pro` (Beer et al. 2026): the same bookkeeping, driven by the GLAMOS distributed mass balance on 9 Swiss glaciers |
| perfect-firn test | the thermal model with the ledger fed by observed (GLAMOS) instead of modelled band balances: what the ledger would do if the model's firn lines were right |
| FRHL-o, FRHL-O | the ledger fed with observed balances: o = firn thickness only (mask still GloGEM's), O = firn thickness and the firn line (mask from the same observations) |
| firn mask | which bands the thermal model treats as firn (today: 5-year mean band balance > 0) |
| firn line | the lowest elevation where firn survives the summer |
| m w.e. | metres of water equivalent |
| 'hl', 'firn' | the meltwater depth limit: Herron & Langway's close-off depth (75–90 m in the zone), or the ledger's firn thickness |
| 3b-ii | the firn density (so conduction and heat capacity) follows the ledger instead of one fixed profile |
| H, FRH, FRHL | thermal variants: H the earlier one with T_off; FRH the central one (no T_off, ERA5 lapse rate, Herron & Langway meltwater); FRHL = FRH with the ledger |
| f_ref | the refreezing factor: the share of melt and rain whose heat enters the column |
| zone, robust zone site | the 3000–3800 m band of the firn–no-firn transition; a borehole site whose temperature class is right across the f_ref values |
| zone target D | at least 2 more zone sites in the right temperature class than H |
| cooling guard | no firn band that 'hl' keeps temperate may be cooled by more than 0.2 °C; from now on only for mature firn (≥ 15 layers) |
| A's conditions | M1 not worse, no robust loss of temperature classes, warming error not worse |
| M1 | RMSE against borehole temperatures below 15 m, each outline weighted equally |
| frozen-bed table | how many measured frozen beds the model gets right (12 of 14 so far) |
| 360-day thermal year | the monthly thermal model counted every month as 30 days |
| GLAMOS | Glacier Monitoring Switzerland |
