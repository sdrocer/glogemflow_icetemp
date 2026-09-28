# Ice permeability, firn extent and the old measurements

*(The previous plan in this file — the thermal spin-up, Phases 0–3 — is finished and archived in
`development_history/plans/2026-09-24_thermal_spinup_plan.md`.)*

## Context — what prompted this

Since this week the model is compared with every borehole measurement at its own month, and the
run covers Oct 1971 – Sep 2025 so that the 1970s campaigns count too. Two questions followed:

1. **Old measurements:** is it right to spin up for, and fit, measurements from the 1970s?
2. **Ice permeability:** Hintereisferner (boreholes 1972–76) comes out ~5 °C too cold. Should
   meltwater be allowed to warm glacier *ice* (we switched that off on ice bands, `97dfbf2`)? What
   would be defensible in a paper, and would it help?

Three investigations (physics literature, the code, the borehole data) and an independent review
answered both. The short version: **ice permeability is not the problem the data point to. The
firn extent is.** Adding water to ice now would make most glaciers worse.

## What the evidence says, in plain words

### 1. What the boreholes show in ablation zones

```
  depth    measured (every ice site with a temperate part)     model today
    0 m  ┬  cold  ─┐                                          on some glaciers: temperate
         │         │ 15–25 m "cold skin" on small glaciers,   throughout (wrong firn), on
   20 m  │         │ 100–200 m cold core on the Gorner tongue others: cold throughout
         │  temperate below                                    (wrong ice)
         ▼
```

- Cold ice **above** temperate ice everywhere; never the reverse. This is the textbook structure
  of temperate ice formed in the firn area, carried down-glacier, with a winter cold wave and
  upward flow (emergence) building the cold skin (Blatter & Hutter 1991; Pettersson et al. 2003).
  The model already contains all those processes.
- Over all ice sites the model is **too warm**, not too cold (+0.5 to +1.1 °C at every depth).

### 2. Why Hintereisferner is too cold — and Tête Rousse and the Grenz tongue too warm

The model decides per band whether it is firn (meltwater percolates and warms it) or ice (no
water). It takes that from its own 1942–1971 mass balance: firn where the balance is positive.
GloGEM calibrates the mass balance for the *whole glacier*, not its shape with elevation, so the
firn line it implies is often far off:

| glacier (RGI7) | model firn line | reality | consequence |
|---|---|---|---|
| Hintereisferner (03116) | ~3185 m | ELA ~2950–3000 m in 1976 | upper boreholes modelled as cold bare ice, measured temperate firn: **−5 °C** |
| Tête Rousse (00777) | all bands firn | firn only near the top | whole glacier temperate, measured −2 °C: **+1.9 °C** |
| Gorner/Grenz (01225) | down to ~2600 m | ~3100–3200 m | tongue warm in 1979–80 (+2.6 °C), then a false cooling trend |

On top of that, after year 5 of the main run the firn mask is **replaced overnight** by one from the
modelled balances, and later bands **lose their firn instantly** once their balance turns negative.
Real firn does not vanish in a year — it is a store of old snow that has to be melted away first.
This creates the wrong 1980→2010 trend on the Grenz tongue (model cools, boreholes warm).

The 1970s geometry is a smaller effect: the Hintereisferner tongue has thinned up to ~100 m since
1976, but mapping the boreholes by position instead of elevation changes their temperature by only
0.1–0.6 °C.

### 3. Ways to let water into ice, and how defensible they are

| mechanism | idea | parameters | defensible? |
|---|---|---|---|
| **none (today)** | ice is impermeable at the grain scale; water runs off | — | **standard** (Cuffey & Paterson 2010 §9.4.3; Wilson & Flowers 2013) |
| **crevasse refreezing** (Gilbert et al. 2020, *TC* 14) | in crevassed areas a fraction of the melt enters the ice and refreezes top-down, as in firn | fraction *f* (~0.03), depth ≤ Nye crevasse depth (~20–30 m) | **high**: conserves energy, tested against boreholes (Himalaya) |
| cryo-hydrologic warming (Phillips et al. 2010/13) | ice warms towards 0 °C near water-filled cracks, time scale R²/κ | crack spacing R (~20 m) | medium: Greenland only, unlimited heat unless capped |
| water content in temperate ice (Aschwanden 2012) | temperate ice keeps ~1 % water as a heat reserve | ω | high, but needs an enthalpy formulation |
| melting surface (Gilbert et al. 2014b) | surface held at 0 °C in months with melt | none | medium; **already in GloGEM** (`firnice_melt_surface`, off) |

**Would they help?** Not now. All of them warm ice near the surface and would erase the measured
cold skins. They are the right tool only for a specific signature that we do not see yet: a correct
cold skin with temperate ice *missing below it*, in crevassed areas. The plan keeps the crevasse
mechanism as a backup and works on the firn first.

### 4. Two problems found on the way

- **~80 observations are silently skipped.** An observation is matched to the nearest model node
  only within half the node spacing. Nodes sit at 9 m and 14 m, so almost every **10 m reading**
  (52 of 53, the classic firn temperature depth) and 28 readings at 61–69 m get no model value.
- **The calibration cannot run.** Parts of the Python calibration still assume the retired
  parameter `z0`: one line fixes the third parameter at 15, which the emulator reads as
  `advection_scale = 15` (allowed range 0.3–1.7), others crash. The sensitivity ensembles were not
  affected; a real calibration would be.
- Also: the likelihood reduces each borehole entity to **one mean residual**, so the shape of the
  profile with depth never reaches the fit — one reason the parameters look like a single knob.

## The plan

Every step is tested one at a time, with the same parameters, against the boreholes and the 3000-
cycle checks where the spin-up is involved, as in Phase 2. One commit each.

### Step 0 — diagnostics (no physics change)
- Log the lowest thermal firn band each year and the number of bands that flip
  (`finalize_annual_massbalance.pro`); log the reference accumulation-area ratio (AAR)
  (`firnice_reference_forcing.pro`); write `firn_th` instead of the mass-balance mask to the
  validation file (`write_firnicetemp_validation.pro`).
- `data.py _classify_firn`: firn / ice / **unknown** instead of treating unknown as firn; resolve
  unknowns from the WGMS ELA of the measurement year where a series exists (Hintereisferner: from
  1953), else leave them out of zone statistics.

### Step 1 — make the calibration runnable (`z0` cleanup, Python)
`calibrator.py` (`_FULL_NAMES`, `priors.z0` at ~360/624), `pipeline.py` (`FIXED_PARAMS`),
`writeback.py` (third column → `delta_advection_scale`, which the IDL reader already expects),
`validation.py`, `baselines.py`, `plots.py`. Align the IDL clip of `insul_scale` ([0.3, 2.0]) with
its prior. **Test:** a calibration smoke run on an existing campaign without errors, with every
parameter inside its bounds.

### Step 2 — no more skipped observations (Python)
Replace "nearest node within tolerance" in `runner.match_model_to_obs` by **linear interpolation
between the two nodes around the observation**, month by month; NaN below the resolved column
(no extrapolation). The emulator then indexes rows by position (`emulator.assemble_matrix`).
Also save each design point's raw monthly profiles, so future matching changes need no IDL reruns.
```
   9 m node ●──────────── 14 m node ●        10.5 m reading: today skipped;
             ↑ 10.5 m = 70 % of the 9 m value + 30 % of the 14 m value
```
**Test:** a synthetic linear profile recovered exactly at 10.5 m and 65 m; +80 matched points on the
7-glacier run; for already-matched points the change stays small (< 0.2 °C) except at real
boundaries.

### Step 3 — the firn extent (the main lever)

**What we learned since the first draft.** In its calibration period (2000–2019) GloGEM puts the
ELA about right (Gorner 3271 m, Tête Rousse 3182 m, Hintereisferner 3201 m). The firn line is wrong
only in the *earlier* decades: under the 1961–1990 climate the model lowers it by ~750 m on
Gorner, where in reality the Alpine ELA rose by perhaps 150–300 m since the 1960s. So the model's
firn line **reacts roughly three times too strongly to climate**. A likely reason: all three
glaciers sit at the edges of their allowed calibration ranges (degree-day factor at its maximum 4.5,
precipitation factor at its minimum 1.3), the combination with the strongest melt response to
temperature. The glacier-wide balance alone cannot tell these parameter pairs apart; snowlines can.

```
  ELA ▲                 real (roughly)          model (calibrated 2000–2019)
      │ 3300 ─                    ●─────● ✓  ────────────────────●─────●  right here
      │ 3000 ─        ●──────●                                   ↑
      │ 2700 ─                              far too steep a rise │
      │ 2500 ─                             ●──────●              │
      └──────────── 1942–71 ─────── 1991–2000 ─────── 2000–19 ──────▶
```

**Why not anchor the thermal firn line separately (the first draft's 3b).** You are right that it
would decouple the thermal firn line from the mass balance, which still uses its own. And with a
model that reacts 3× too strongly, an anchored line moved by the model's changes would overshoot
by the 2000s. So the plan keeps **one firn line, shared with the mass-balance model**, and fixes
the cause where it sits:

**3a. Firn with memory, driven by the mass-balance model's own balance.** The thermal model keeps
a firn store per band (m w.e.) that grows with the model's balance and shrinks with it, up to the
firn mass above the pore close-off depth. A band turns to ice only once its firn has actually been
melted away. Same balance as the mass-balance model, so the two firn lines follow the same source;
only the overnight switch after year 5 disappears.
```
  firn store ▲        balance < 0 for years: the store empties slowly
             │███████▇▇▆▅▄▃▂▁                  → band becomes ice only here ↓
             └────────────────────────────────────────────────────▶ years
```

**3b. Measure the sensitivity problem before fixing it.** From the main run 1972–2025, log each
glacier's modelled ELA per year and compare with observed ELA series (WGMS: Hintereisferner since
1953; plus other Alpine WGMS glaciers in the calibration set). If the model's ELA indeed rises ~3×
too fast, that is a mass-balance calibration problem for all of GloGEM, not a thermal one.

**3c. The fix at the source: snowline calibration** (coordinate with Alexandra / Lander).
Transient snowlines add exactly the missing second constraint: where the snow line sits through
the summer depends on how accumulation and melt are split, not only on their difference. What
exists today: a standalone daily RGI6 script on the old `snow-calibration` branch (calibrates the
precipitation factor per glacier by fitting modelled to observed snowlines, 2000–2020), and Alps
snowline data (Sommer et al. 2026, RGI7 region 11) in Alexandra's directory, not readable for us.
Nothing on `main` uses it yet. Porting it into the current (monthly, RGI7) structure is a GloGEM-
wide change, so it needs your and their agreement. Once a snowline-calibrated parameter file
exists, the ice-temperature runs can read it without code changes.

**Until then:** the thermal model uses 3a with the current mass balance and we accept that the
firn extent in the early decades is too large on some glaciers — and report it. Optional interim,
only for the spin-up if you want it: widen the calibration bounds that are hit (precipitation factor
below 1.3, degree-day factor), which lets the geodetic fit choose a less sensitive combination.

**Expected, and how we check it:** 3a removes the false 1980→2010 cooling on the Grenz tongue;
Colle Gnifetti (always firn) unchanged — if it moves, something leaked; Tête Rousse and
Hintereisferner improve only once the mass balance's sensitivity is fixed (3c). Firn line per
year checked against WGMS ELA series.

### Step 4 — test the existing melting-surface switch
Run with `firnice_melt_surface = 'y'`. Keep it only if it reduces the near-surface ice bias without
warming the firn sites; if kept, weight it by the share of melt days in the month (one melt day
should not pin the whole month to 0 °C).

### Step 5 — fix strain heating
A units error makes it ~900× too small (`firnice_temperature_model.pro` ~358: divide by `cap_fit`
only), and the stress must not exceed the basal value below `thick`. Effect: only deep ice on thick
glaciers (Gorner's temperate bottom layer); off in calibration unless the basal rows improve.
**Test:** a column with constant heating warms by Q·Δt/(ρc).

### Backup (not scheduled) — crevasse refreezing
Kept in reserve; we first improve the firn (steps 3–4). Revisit only if, after the firn work,
observations at the melting point are still missed *below a correctly modelled cold skin*
(residual ≈ 0 at 0–10 m, > 1 °C at 15–40 m) and the misses sit in or just below stretching
(crevassed) bands. The mechanism would then be Gilbert et al. (2020): in bands where the Nye
crevasse depth exceeds ~10 m (from the velocity difference between bands), a fraction `f_crev`
(prior ~0.03, range 0–0.1) of the melt refreezes top-down to at most the Nye depth (≤ 30 m),
reusing the firn fill loop and flagging those layers for the drift spin-up.

### Old measurements — recommendation
- **Keep the 1971 start (A).** The 1970s readings are Hintereisferner's only ablation data and
  half of every trend pair; they are also closest to the spin-up equilibrium.
- **Use the old-vs-new pairs as a held-out trend test (D):** Colle Gnifetti 1976–83 vs 2003–21,
  Grenz tongue 1979–80 vs 2008–10, Col du Dôme 1973–80 vs 2005–17. Taking differences cancels level
  errors such as the accumulation profile, so it tests exactly what the level-based fit cannot.
- **Don't start in 1991 (B):** it throws those data away and hides the firn transient instead of
  fixing it.
- **Down-weighting (C) as "larger measurement error" does nothing** (0.1 °C against a model error of
  ~2 °C). If old entities are still clearly worse after step 3, add a model-error term for
  pre-1979 entities instead.

### Then: one new training campaign
Only after steps 1–3 and the step 4 decision, since every physics change invalidates the emulator.
The date-matched ensemble running now (`alps_dm`) still shows how the matching alone changes the
sensitivity; it will be superseded.

## Verification (per step)
- Step 1: calibration smoke run, parameters within bounds; unit test of the parameter expansion.
- Step 2: synthetic interpolation test; matched-point count and change distribution.
- Step 3: borehole bias by glacier and zone (same parameters), firn line vs WGMS ELA, the year-by-
  year firn log from step 0, drift spin-up convergence and a 3000-cycle check.
- Steps 4–5: the signatures named in each step, row-level residuals by depth.
- Book: a new chapter 6 ("Firn extent and water in ice") as the steps land.

## Decisions for you
| topic | recommendation |
|---|---|
| firn extent | **3a** firn memory on the mass-balance model's own balance (one shared firn line); **3b** measure the ELA sensitivity against WGMS series; **3c** snowline calibration as the real fix, to agree with Alexandra / Lander |
| ice permeability | **backup only**: crevasse refreezing (Gilbert 2020) kept in reserve; firn first |
| old measurements | **A + D**; C only as a model-error term if still needed |
| order | 0 → 1 → 2 → 3 → 4 → 5, then one training campaign |
