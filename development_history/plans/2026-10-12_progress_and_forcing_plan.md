# Where we stand, the high-elevation forcing, and the road to the Alps-wide runs

*Plan of 2026-10-12 (evening). Written for you to read when you are back. Nothing has been run or
changed since the thermal-forcing runs (icetemp `61cb4bf`). The plan rests on:*
- *read-only analyses of those runs;*
- *three searches: GloGEM's radiation code, the terrain and sunshine at the boreholes, and the
  literature;*
- *the calibration table of all 4073 CentralEurope glaciers;*
- *a critical review of the draft by an independent agent, whose corrections are included.*

*The plan of this morning is archived as `2026-10-12_thermal_forcing_plan.md`.*

## Context

Today's runs took the mass-balance temperature offset (T_off) out of the thermal model's energy
inputs (run F). That fixed Mont Blanc (bias +15.6 → −0.8 °C) and Colle Gnifetti. But it made sunny
high firn far too cold: Lys went from +1.8 to −6.8 °C, and Breithorn and Fiescher cooled too. So F
failed the check we had fixed before the runs. You asked three things:

1. Which way to go for the high sites, and are there other ways than removing T_off?
2. What else could lower the RMSE and raise the predictive power?
3. How far are we from the final aim, which is forecasting the thermal regime of all sorts of
   Alpine glaciers? When can we be "zufrieden" and start the real runs?

## The short answers

1. **High sites.** None of the four options is right as it stands. T_off does two jobs that are not
   its own:
   - At windy summits it invents melt, where in reality the snow leaves by wind.
   - On sunny high firn it stands in for melt that the degree-day model misses.

   The defensible path is to keep T_off out of the thermal energy inputs, but to **judge that
   fairly**. The failed check ran both forcings at one shared refreezing fraction f_ref. A single
   f_ref penalises the forcing that removes melt, because f_ref is exactly the knob that sets how
   much heat the meltwater brings.

   So we run one set of runs where each forcing gets its own f_ref, chosen out of sample. The
   forcing is then adopted by a rule fixed in advance. Two things are tested separately, on data
   independent of the thermal model:
   - the lapse rate, against measured mountain air temperatures;
   - sunshine (exposure), which the literature names as the control on high firn, in Python. We
     build it only if the test clearly passes.
2. **Other improvements.** Three sites carry 80 % of the error below 15 m. On the well-behaved
   glaciers we are already as close as a model with one forcing per elevation band can get. The gains
   lie in a few understood failures: Taconnaz, and the 3000–3500 m zone.

   The most important missing evidence is not an RMSE at all. It is **whether the model gets the
   warming over time right**. Repeat measurements exist at Colle Gnifetti, Col du Dôme and Tête
   Rousse.
3. **Progress.** The model works where its forcing is right: below 15 m, 1.3 °C in ice and 1.7 °C
   on the Step 5a glaciers. On other glaciers it does not yet beat a simple rule such as "cold above
   X m".

   Proposal: one gate (the forcing, decided by the fixed rule) inside a time box of 3 weeks, or 5
   with sunshine. Then freeze v1.0, calibrate once, evaluate once and run the Alps.

## 1. How we are doing

### What the model does today (prior-mean parameters, scored set: 23 glaciers, 110 entities)

| below 15 m | result | reading |
|---|---|---|
| ice below the firn | RMSE 1.3 °C | good |
| Step 5a glaciers (used to develop the model) | RMSE 1.7 °C | at the floor of a band model (see below) |
| the other glaciers | RMSE 3.9 °C (2.0 without Mont Blanc and Taconnaz) | a few large failures |
| warm or cold regime, all glaciers | regime skill 0.67 (Step 5a 0.79, others 0.51) | |

*Regime skill* is the average of two shares: warm sites found and cold sites kept. A coin toss scores
0.5, a perfect model 1.

**With three classes, which is what the study needs, it is weaker.** Polythermal sites are the
hazard cases and finding polythermal glaciers is a main goal of the study:

| below 15 m, today | cold found | **polythermal found** | temperate found | 3-class skill |
|---|---|---|---|---|
| all glaciers | 50 of 71 | **6 of 13** | 2 of 6 | 0.50 (chance: 0.33) |
| Step 5a glaciers | 35 of 45 | 5 of 8 | 1 of 2 | 0.63 |
| other glaciers | 15 of 26 | **1 of 5** | 1 of 4 | 0.34 |

- Run F finds 0 of the 5 polythermal sites on the other glaciers: it loses Lys at 4085 m.
- The polythermal profiles the model does catch are on Tête Rousse and the Gorner tongue.
- It misses those on Alphubel S, Hohsaas, Colle Gnifetti at 3985 m and Lys.

From here on, **three-class skill and polythermal recall are the main regime measures**, with the
two-class numbers only beside them.

### Two yardsticks we did not have

**The floor.** Boreholes on the same outline, within 60 m of elevation and 5 years of each other,
differ by 2.2–2.7 °C (RMS). At Colle Gnifetti and Dôme du Goûter they differ by up to 8–10 °C. A
model with one forcing per elevation band cannot come closer to a single borehole than about
1.5–1.9 °C. The Step 5a glaciers are at 1.7 °C, so tuning cannot buy more there.

**A simple rule.** Each rule below is fitted leaving one glacier out, and both columns use the same
entities and points:

| below 15 m | rule | GloGEM today |
|---|---|---|
| regime skill, all glaciers | 0.67 ("warm below X m" or "warm above air temperature Y") | 0.67 |
| RMSE, all (rule = straight line on air temperature) | 2.9 °C | 2.6 °C |
| Step 5a glaciers: skill / RMSE | 0.75 / 2.5 °C | 0.80 / 1.8 °C |
| other glaciers: skill / RMSE | 0.56 / 3.8 °C | 0.51 / 3.9 °C |

Resampling the 20 glaciers with points below 15 m gives a 90 % range of 0.53–0.77 for the regime
skill and 1.7–3.4 °C for the RMSE. Comparing two model versions on the same glaciers is sharper
than that. The fair-test script will report those paired differences.

### Where the error sits (share of the squared error below 15 m, today)

| glacier (RGI7) | share | what goes wrong |
|---|---|---|
| Taconnaz (00774), the 2008 holes at 3410–4080 m | 41 % | −6.5 °C measured at 3685–3785 m, 0 °C modelled |
| Monte Rosa: Colle Gnifetti, Grenzgletscher slopes and tongue (01225, 01230) | 27 % | Colle Gnifetti on 01230 +3.5 °C (offset) |
| Mont Blanc summit (00786) | 12 % | only 2 points, but +15.6 °C (offset) |
| Bossons (00773) | 6 % | |
| Fiescher (02611), Tête Rousse (00777), Lys (01205) | 3–4 % each | |

How often the regime is right, by elevation:

| elevation | regime right | |
|---|---|---|
| below 3000 m | 67 % | |
| 3000–3500 m | **42 %** | 10 of the 11 misses make slightly cold sites temperate: Tête Rousse, Corvatsch, Altels, Hohsaas |
| 3500–4000 m | 65 % | |
| above 4000 m | 85 % | 95 % with F |

**Verdict.** The model is good where its forcing is right. Its failures are few, large and high up.
But two things are missing:
- we cannot yet claim it beats a simple rule on new glaciers;
- it has never been tested on change over time, which is what a forecast does.

Those, not another 0.1 °C, stand between us and the production runs.

## 2. The high sites

### What the offset really is

```
 the mass-balance calibration, glacier by glacier (Huss & Hock 2015)
   c_prec  1.3 ── 2.3    ── budget closed? stop
   DDF     1.75 ── 4.5   ── budget closed? stop
   T_off   no bound      ── until the budget closes
```

- **It grows with height.**
  - Look at the 181 CentralEurope glaciers with their own measured target and a mean elevation
    above 3400 m. 93 % sit at both bounds (least snow, most melt) and still need a median T_off of
    +3.1 °C.
  - Of the 58 glaciers reaching above 4000 m, 42 carry 2 °C or more.
  - Below 3000 m the median is about 0.
  - It is slightly larger on south-facing glaciers (mean 0.78 against 0.40 °C).

  So at the top of the Alps the model keeps too much snow even at the bounds, and T_off melts it
  away. This is budget closure, and it grows steadily with height. In an Alps-wide map, T_off would
  put melt exactly on the high glaciers where cold firn lives.
- **Two opposite causes hide behind it.**
  - Mont Blanc: its target is a real measurement (−0.06 ± 0.18 m w.e./a). Its firn is as cold as its
    air, so the snow leaves by wind and sublimation, not by melt.
  - Lys: its firn is near 0 °C. That means real melt, more than the degree-day model gives.
- **T_off does not tell which high sites are wet.** We measure this with a ranking skill: how often a
  wet borehole comes out warmer than a dry one, where 0.5 is a coin.
  - T_off itself scores 0.5.
  - Melt computed with T_off ranks the sites no better than melt without it: 0.86 against 0.87,
    above 3300 m.
- **Without T_off the model orders the sites better.** Rank agreement with the measured temperatures
  below 15 m rises from 0.68 to 0.75, and on the other glaciers from 0.40 to 0.57.

  What F gets wrong is the **level**. At f_ref 0.55 it gives sunny high firn too little water: it
  finds 50 % of the wet sites at 10–20 m, instead of 62 %.
- **Literature.** The GloGEM papers treat T_off as an error in the air temperature (Huss & Hock 2015;
  Zekollari et al. 2019, 2024). We therefore argue the opposite with two tests that do not involve
  the thermal model's parameters:
  - the same site on two outlines;
  - dry firn as a thermometer.

  The equifinality literature supports us: the three factors trade off, and they "may also be
  accounting for aspects of the model physics that are not properly accounted for" (Rounce et al.
  2020; Schuster et al. 2023).

### What decides wet or dry firn up high: height and sunshine

```
 firn at 10-20 m, boreholes above 3300 m, against the band air temperature
    0 °C ┤ ● ●●  ●● ●  ●  ●●  ●            wet: meltwater refreezes, firn held near 0 °C
         │
   -5 °C ┤   ○    ○ ○○  ○ ○ ○
         │            ○  ○○  ○○ ○ ○        dry: firn ≈ air temperature + 0-3 °C
  -15 °C ┤                        ○ ○○ ○
         └──────────────────────────────── air temperature (colder →)
   neighbouring holes in the same band can sit in different states, up to 8-10 °C apart
```

- **The literature.**
  - Suter et al. (2001): cold firn is probable above about 3400 m on north faces, but only above
    4150 m on south faces. That is 750 m, or about 5 °C of air temperature. "Within a few hundred
    metres of altitude, firn temperatures are mainly determined by the direct incoming solar
    radiation rather than the MAAT."
  - Suter & Hoelzle (2002): height and radiation explain 82 % at Monte Rosa.
  - Gilbert et al. (2014): at Col du Dôme, melt follows the potential solar radiation.
  - Mattea et al. (2021): around Colle Gnifetti, melt ranges from under 1 cm w.e./a on steep shaded
    slopes to over 30 cm. 17–34 % of it falls at air temperatures below 0 °C, which a degree-day
    model cannot produce.
  - Hoelzle et al. (2011): Grenzgletscher firn warmed by 6.8 °C in 1991–2008, driven by meltwater.
- **Our boreholes.** Clear-sky sunshine over the year, as a share of flat ground, from swissALTI3D
  and RGE ALTI:
  - The sunniest sites are warm: Alphubel S 1.17, Lys 1.10, Jungfraufirn 1.03.
  - The shadiest are cold: Taconnaz 0.63–0.79, where the horizon to the south reaches 26°.
  - At Seserjoch and on the Grenzgletscher slope there are pairs at the same height facing opposite
    ways.
  - But summer differences are small (0.83–1.02), and Fiescher and Breithorn are warm with little
    sun (0.88–0.91).

  Sunshine explains part, not all.
- **GloGEM cannot do it yet.** Its radiation melt model ('3') has three problems:
  - it reads a radiation input that is never filled, so it would stop;
  - it uses only the noon sun;
  - read with its code, the RGI7 aspect classes come out one class too far north. I checked this
    against the DEMs; details under "Afterwards".

  Only the geometry idea is reusable. But the radiation input exists locally: ERA5-Land monthly
  downward solar radiation for 1950–2022, and the MeteoSwiss satellite product.

### The options, judged

| option | verdict |
|---|---|
| **F, each forcing with its own f_ref**: energy inputs without T_off | **first: one fair run set** |
| **lapse rate**: the June floor (−6.5 K/km) or ERA5's own rate (your option 3) | **decided by measured air temperatures, not by the runs.** The floor and f_ins both move dry firn, so the runs cannot tell them apart. The floor was added in June, while T_off was still in, against a warm bias. It acts on the surface but not on the melt, so whichever rate wins applies to both. Effect at 4400 m: about 1 °C. **Thermal model only:** the floor lives in `dtdz_icetemp`, which only the thermal files read; the mass balance keeps `dtdz` and its files stay identical. If the stations show ERA5's rate is off at height, that is a finding for the mass-balance side, not a change we make there. |
| sunshine for the water that refreezes | **test first** (Python). Build only if it passes with band sunshine. |
| your option 2: drop only the part inherited from the cell | fixes Colle Gnifetti 01230 and Jungfraufirn (both inherited), not Mont Blanc (own +5.8). c_prec and DDF were calibrated around the inherited part. **No.** |
| your option 4: today's forcing plus documentation | **fallback only.** T_off grows with height, so the Alps-wide map would be biased where cold firn lives. |
| a cap, or one offset for the region | arbitrary. **No.** |
| less snow on summits (wind erosion) | the real cause sits on the mass side. **Report it to the mass-balance side through you.** |

### Recommendation

```
 independent of the thermal model (Python, week 1):
   lapse rate: rebuilt air temperature vs measured mountain stations ─► floor or ERA5 rate
   sunshine: does it explain wet vs dry beyond the air temperature?  ─► build or not (you)
                                   │
                                   ▼
 one run set: today's forcing (H) | offset-free (F*H, with the chosen lapse rate)
              each at f_ref 0.10 ... 1.00; f_ins checked at two values
              the forcing and f_ref chosen leaving one MASSIF out
                                   │
                                   ▼
 rule fixed now: the offset-free forcing is adopted if it is NOT WORSE by more than
 0.3 °C (RMSE below 15 m, each outline weighted equally) and 0.05 (three-class skill),
 overall AND in the 3000-3500 m zone; otherwise today's forcing stays.
 Polythermal recall and the frozen-bed check are reported beside it (13 sites: too few for a margin).
 Ties go to the offset-free forcing on physical grounds.
                                   │
                                   ▼
              v1.0 forcing; the loser is kept as a projection sensitivity
```

Why this holds up:
- The two tests against T_off do not involve the thermal model's parameters.
- Each forcing gets its own f_ref, the parameter that scales the meltwater heat; f_ins is checked
  at two values.
- The forcing and f_ref are chosen without the massif they are scored on, so the scores stay out of
  sample.
- The rule and its margins are fixed before the runs.
- In the paper we also report the failed fixed-parameter check, and why it was not a fair test.
- H (the published Herron & Langway form, and latent heat per metre of water) is in both variants,
  because it fixes a bug.
- If the offset-free forcing wants f_ref at 1.0, the model is missing melt, and that counts as
  evidence for sunshine.
- Keeping the losing forcing as a projection sensitivity shows in the forecast how much the choice
  matters. With T_off kept, high sites would change regime earlier.

*What I expect:* the offset-free forcing chooses f_ref around 0.7–0.9, wins back part of Lys,
Fiescher and Breithorn, and loses a little at 3000–3500 m. That is within the margins, so it is
adopted. If not, today's forcing stays, with Mont Blanc and Colle Gnifetti 01230 documented.

*Note on sunshine.* A plain sunshine factor on the degree-day melt changes summer melt by only
±15 %. That is too little for Taconnaz or Lys. The strong effects reported in the literature come
from the melt threshold: melt starts only when the surface energy balance turns positive, and near
that threshold a little more sun can double it.

Capturing this needs a simplified energy balance for the water that refreezes (Oerlemans 2001;
Giesen & Oerlemans 2012):

    melt = (1 − albedo) × sunshine + C0 + C1 × air temperature, when positive

It must be computed with the sub-monthly spread, as the degree-day melt is. GloGEM's unused model
'3' already has this form. It brings 2–3 new constants (albedo, C0, C1), taken from the literature
or calibrated. That is new physics, hence the test first.

## 3. Other things that limit the model, ranked

**1. Change over time, never tested.**
- *Evidence:*
  - Same holes or close repeats: Tête Rousse holes 66 and 69 (2010→2015) and 423–425; CG05-1,
    CG18-1 and SJ98-A.
  - Grenzgletscher at 4250 m, 1991→2013. Hoelzle et al. (2011) report +6.8 °C there.
  - Col du Dôme at 4245 m, 1994→2017. Vincent et al. (2020) report +1.5 °C at 50 m.
- *What to do:* compare modelled and measured warming at the same depth and months. The monthly
  profile series are already saved, so this is Python only.
  - The list of pairs is fixed before computing anything.
  - The noise level comes from holes drilled in the same campaign: CG08-1/2/3 and CG13-1 to 4.
  - Caveats: all repeat sites lie on Step 5a glaciers, so this tests the *trend*, not the level;
    and the spin-up settles to the 1942–71 climate.
- *Buys:* the key evidence for forecasts.

**2. Taconnaz (00774).**
- *Evidence:*
  - 41 % of the error below 15 m: −6.5 °C measured, 0 °C modelled at 3685–3785 m.
  - It is not the offset (0.12) and not the meltwater depth (run P).
  - Meltwater reaches 14–34 m in most years and 80–100 m in some.
- *What to do:*
  - Read Gilbert et al. (2015), where these holes come from.
  - Runs on Taconnaz alone: f_ref 0 or 0.55 against f_adv 0.3 or 1.7, plus f_ref 0.1 and 0.2. Is it
    the water (shade, or runoff on steep slopes) or cold firn from the Dôme du Goûter not carried
    down?
- *Buys:* up to 41 % of the error. Steep, cold hanging glaciers are also the hazard cases.

**3. The 3000–3500 m zone: a priority, because of hazards.** Small, cold, steep glaciers sit
here, and this is where regimes change first under warming. So v1.0 must get *better* here, not just
not worse. Section 3a below has the work package.

**4. Data we lose, now audited.** I checked all 225 Alpine glenglat holes. Each was placed in
GloGEM's RGI7 glacier masks and compared with the GLIMS-id lookup the calibration uses.

| what | holes | details | fix |
|---|---|---|---|
| lost through the GLIMS id | **6** | Only Chessjengletscher. You were right: all 6 holes (2983–3046 m) lie inside **01161** (Feegletscher), whose GLIMS id is G007871E46080N. glenglat gives G007919E46068N, which is not in RGI7. Its 5 entities exist but never run. | assign holes to outlines **by position**, not by GLIMS id |
| on the wrong outline | 3 | CG95-1, CG08-3 and CG13-3 carry 01225's id but lie inside 01230 | same |
| on an outline edge | 17 | 12–55 m outside every mask (cells are 40–100 m). Mostly the Colle Gnifetti saddle, on the 01225/01230 divide. Also two Lisjoch holes, labelled Grenzgletscher (01225) but nearest to Lys (01205), and two Dôme du Goûter holes labelled Bossons or Taconnaz. | position, with the RGI7 polygons rather than the coarse masks; list the edge cases |
| mixed entity | — | each glacier-and-band entity takes the outline of its *first* hole, so mixed-outline entities exist at the saddle | split entities by outline |
| **dropped by the 10–30 m rule** | **7** | Entities need a reading at 10–30 m. That drops the deepest profiles we have: 6 Grenzgletscher tongue holes at 2545–2572 m (ryser2013, to 200–303 m) and Gornergletscher (eisen2009, to 359 m). This is exactly the deep ice the below-15 m target is about. | require a reading below 15 m **or** at 10–30 m |
| no usable readings | 18 | before the run starts (hughes1939 ×5, agassiz1847), equilibrium marked false (6, incl. St. Annafirn and 3 Corvatsch holes), only shallower than 10 m (5, incl. your Corvatsch and Tortin loggers), no readings (Hintergrat) | none (the "false" profiles would need a bias correction, a separate topic) |

- Chessjengletscher (01161) lies in the 3000–3500 m zone, in the ablation area. It becomes the
  **sealed** test glacier there: never used in a decision, scored once at the end.
- Its loggers give readings below 15 m only for August–September 2025, because the runs end in
  September 2025. The ERA5 file runs to December 2025, so the window could be extended later.
- Adding the 7 deep tongue profiles and moving holes to their real outlines changes the data set.
  So it happens **before** the fair test (Step 1c), and every number is recomputed once.

**5. The metrics, fixed before any run.**
- *Evidence:*
  - The pooled RMSE is dominated by a few outlines: 01225 has 191 points, Bossons 104, Tête Rousse
    93, Taconnaz 67.
  - The regime classification is tolerance-sensitive: half the glaciers lie within 0.3 °C of the
    melting point.
- *What to do:* write these into `forcing_fair_test.py` before the runs:
  - RMSE below 15 m with each outline weighted equally; the pooled RMSE beside it;
  - **regime with three classes** (cold, polythermal, temperate) at each hole's tolerance, using
    `thermal_structure.classify`:
    - three-class skill (the average of the three recalls);
    - the share of polythermal sites found, and the full table of measured against modelled class;
    - the two-class numbers only beside them;
  - **the bed**: at the 44 Alpine holes that reach the bed (15 glaciers, among them Tête Rousse,
    Gorner, Grenzgletscher, Bossons, Altels, Alphubel S, Hohsaas, Chessjen), measured bed
    temperature against the model's bottom node at that band. Is the bed frozen or at the melting
    point? The model column ends at the band's ice thickness, so this needs no new run;
  - 3000–3500 m reported separately;
  - the simple rules beside GloGEM, fitted the same way. For three classes, the rules get two
    thresholds. A rule from elevation or air temperature has no depth structure, so it cannot
    really produce a polythermal profile. This is where a physical model has to earn its place.
- *Buys:* honest numbers, and the first thing a reviewer asks.

Not worth doing now:
- where the meltwater stops (run P had no effect);
- strain heating and the melting-surface option (both tried, both worse);
- geothermal heat (the ice is already at 1.3 °C);
- more tuning on the Step 5a glaciers (they are at the floor).

## 3a. The 3000–3500 m zone (hazard priority)

A first look at the entities in the zone (today's run, below 15 m) shows two kinds of misses:

| kind | sites (RGI7) | measured / modelled | what the model does there |
|---|---|---|---|
| **near-melting-point flips**: temperature off by less than 1 °C, but the regime flips | Hohsaas (01112), Jungfraufirn at 3415/3535 m (02596), Altels (01265), upper Tête Rousse (00777) | e.g. −0.38 / −0.07, −0.75 / −0.02, −0.28 / −0.02 | right temperature to within a few tenths, wrong class |
| **clear misses**: measured cold, modelled 0 °C | Corvatsch (02103) at 3295–3325 m; lower Tête Rousse at 3135–3165 m; Taconnaz at 3415 m | −0.8 to −2.7 / −0.01; −1.5 to −2.5 / −0.04 to −0.8; −2.6 / −0.04 | meltwater reached these bands in most years while they had firn (Corvatsch to 14 m, Tête Rousse 29–54 m before its lower bands became ice, Taconnaz 59 m); thick winter snow (1 May: 1.2–1.8 m w.e.) |

```
 measured, at depth:   cold ice/firn ─────────── −1 to −3 °C   (winter cold wins)
 modelled:             summer water heat ─────── held at 0 °C  (meltwater wins every year)
 the balance that decides it:  winter cooling (depends on snow thickness)  vs  refreezing water
```

Three candidate causes, each with a literature anchor, tested in this order:

1. **Too little winter cooling: the snow is too thick.** Gilbert et al. (2012) show that snow
   thickness controls the cold ice of Tête Rousse. The model has 1.8 m w.e. there on 1 May, and the
   glacier sits at the lowest snow factor (c_prec 1.3). This is the same high-elevation over-
   accumulation that T_off hides.
   - *Test:* compare modelled winter snow with measured winter balances (GLAMOS/WGMS where they
     exist: Sex Rouge, Hintereisferner, glaciers near Hohsaas) and with Gilbert et al. (2012) for
     Tête Rousse.
   - *Sensitivity:* f_ins (less insulation) on the zone's glaciers.
2. **Too much meltwater, too deep.** The model lets water soak into cold firn layer by layer, which
   is an upper limit. On cold, thin firn, refrozen ice layers block deeper percolation and the rest
   runs off.
   - *Test:* the wet-depth log against the measured profile shapes.
   - *Sensitivity:* f_ref on the zone's glaciers.
   - If this is the cause, the fix is a percolation rule for cold firn, which is new physics with
     its own short plan. The literature on blocking ice layers and runoff (e.g. Machguth et al.
     2016, Mattea et al. 2021) will be read first.
3. **Heat left over from the firn years.** I had first suspected "firn where there is ice". The
   model's firn line, decade by decade, says otherwise. Firn covers almost all of these small
   glaciers in the 1970s–80s (you say that is right), then the firn line rises:

   | glacier (RGI7) | glacier spans | 1970s | 1990s | 2000s | 2010s | 2020s |
   |---|---|---|---|---|---|---|
   | Tête Rousse (00777) | 3035–3275 m | 3035 | 3047 | 3134 | 3187 | 3248 |
   | Corvatsch (02103) | 2885–3375 m | 2910 | 3009 | 3094 | 3159 | 3285 |
   | Altels (01265) | 2635–3265 m | 2715 | 2818 | 2895 | 2977 | 3083 |
   | Hohsaas (01112) | 2785–3965 m | 2960 | 3050 | 3183 | 3217 | 3317 |

   - At Tête Rousse in 2010 the model agrees with glenglat's labels: ice at 3135–3165 m, firn at
     3185–3205 m.
   - So the firn extent in the measurement years looks right there. The Corvatsch sites (3295–3325 m,
     1999–2001) lay above the model firn line, but their labels are mostly unknown.
   - The likelier problem is memory: bands that were wet, temperate firn until the 1990s turn into
     ice but keep that heat. The lower Tête Rousse ice is still near 0 °C in the model in 2010; it is
     measured at −1.5 to −2.5 °C.
   - *Test:* the modelled temperature at these bands through the year they lost their firn, against
     the repeat readings at Tête Rousse (2010 → 2015 → 2022).
   - If the firn heat of the 1970s–90s was too large, causes 1 and 2 above are the same problem seen
     earlier.

**What matters most for hazards: is part of the glacier frozen to its bed?** As you said, the
classic danger is a polythermal glacier with a frozen terminus and a temperate core higher up. That
is the setting for glacier avalanches and for water pockets held back by a cold margin; Tête Rousse
is the textbook case. So the zone, and every glacier, gets these measures, fixed before the runs:
- **At the bed** (44 holes reach it): frozen or at the melting point, measured against modelled.
- **Along the glacier:** which bands are frozen at the bed and which are temperate, giving a
  glacier type: temperate, cold, or polythermal with a frozen front or frozen upper part.
  - Checked against what is known for our glaciers, for example the cold margin of Tête Rousse, the
    cold core of the Gorner tongue, and the hanging glaciers Taconnaz and Altels.
- **In each profile:** whether the model has a cold layer where one is measured, and the
  temperature error near the melting point.
- The three-class profile classification stays the main measure. These sit beside it.
- For the Alps-wide product, the bed temperature and the glacier type become outputs per band and
  year.

**How it enters the plan.**
- The diagnosis (Python, existing runs) runs in week 1.
- Sensitivity runs on a catchment of the zone's glaciers (00777, 02103, 01265, 01112) take minutes
  each, in week 2.
- In the fair test the zone has its own margin: the adopted forcing must not make the zone worse.
- If a cause is found that a model change can fix, it goes into v1.0, and the time box grows by
  that work (your call).
- Chessjengletscher (01161) is the sealed test of whatever is changed here.

## 4. When are we "zufrieden"?

**One gate:** Z1, the forcing decided by the rule above, with the metrics and code frozen. When Z1
is passed, v1.0 is frozen.
- The time box is three weeks, or five if you choose to build sunshine, plus whatever fix the
  3000–3500 m zone needs.
- Whatever is open by then becomes a stated limitation.
- If we are faster, we move on earlier.

Then we calibrate once and evaluate once. The paper reports the rest as results, not as new rounds
of tuning:

| | measured once after calibration | today |
|---|---|---|
| Z2 | Glaciers not used for calibration, below 15 m. GloGEM against the two rules; the rules are fitted on the Step 5a glaciers only, the same information the calibration gets. | RMSE: GloGEM 3.9, rule 4.1 °C (GloGEM better). Two-class skill: GloGEM 0.51, rule 0.60 (rule better). Three-class: GloGEM 0.34, with 1 of 5 polythermal sites found; the rules' numbers are still to be computed. |
| Z3 | Warming over time at the pre-listed repeat sites. Pooled modelled-to-measured slope between 0.5 and 2, and the sign right at most sites. | not tested (can be run in week 1 on today's runs) |
| Z4 | Every remaining large failure has a stated cause and falls in a class the Alps-wide product can flag, e.g. shaded steep faces. | partly |

**One stop sign:** if Z3 fails, meaning the slope is outside 0.5–2 or the sign is wrong at most
sites, we take one diagnosis round of at most a week before the production runs. A forecast that
gets past warming wrong is not defensible.

Why stop here: our 20 glaciers cannot tell smaller improvements apart (section 1). A forecast is
defended by the right physics, the time test and an honest uncertainty, not by the last tenth of a
degree.

```
 week 1   Python: data fixes (outlines by position, deep profiles), lapse-rate test,
                  sunshine test (a), time test on today's runs, 3000-3500 m diagnosis
          GloGEM: lapse-rate switch;  icetemp: --fref/--fins, fair-test script (metrics frozen)
 week 2   fair run set + Taconnaz runs + 3000-3500 m sensitivity runs
          ──► you decide the forcing (Z1), sunshine, and any fix for the 3000-3500 m zone
          [+2 weeks if you choose to build sunshine]
 week 3   Alpine ensemble with that forcing ─► calibration ─► Z2, Z3, Z4 measured once
          ──► v1.0 frozen
 then     cost benchmark (one batch of ~200 Alpine glaciers, thermal module on)
          Alps-wide hindcast 1940-2025 ─► compare with Suter et al. (2001, 2002) cold-firn limits
          projections: GMIP4 climate, several GCMs x SSPs x posterior parameter sets
          (+ the losing forcing as a sensitivity)
          ─► regime per glacier and band to 2100, cold-firn maps, the years regimes change
```

## 5. Steps after your OK

**Step 0: housekeeping.**
- Archive this plan and this morning's plan (appendix) in
  `glogemflow_icetemp/development_history/plans/`.
- Update memory: the T_off–height finding holds for measured targets; Mont Blanc's target is
  measured, not a fill value; this plan.

**Step 1: tests independent of the thermal model (Python, week 1).**

*(a) Lapse rate.*
- Rebuild ERA5 air temperatures with `forcing_check.py`. Make its hard-coded −0.0065 floor
  (line 51) a setting first.
- Compare them, with the floor and with ERA5's own rate, against measured monthly air temperature
  at high stations: Jungfraujoch (3571 m), Corvatsch (about 3300 m), Weissfluhjoch (2691 m) and
  Säntis (2502 m).
  - The data come from MeteoSwiss open data (download), or else from the local MeteoSwiss gridded
    TabsM in `/itet-stor/jabeer/glazioarch/GlacioClimateData/MeteoSwiss/verified/`.
  - Add Aiguille du Midi or Capanna Margherita if a series can be found.
- Rule, fixed now: the version with the smaller mean annual bias at the stations above 3000 m wins.
  The lower stations only show how the bias changes with height. A difference below 0.3 °C keeps
  today's floor.
- The chosen rate then applies to the thermal surface **and** the thermal melt.

*(b) Sunshine, model-free part.*
- New `scripts/icetemp_analysis/exposure_check.py`.
- Terrain at each borehole:
  - slope and aspect from swissALTI3D 10 m (Switzerland, Monte Rosa) or RGE ALTI 5 m (France);
  - horizons from EU-DEM (all three under `/itet-stor/jabeer/glazioarch/GeoBaseData/HeightModels/`);
  - from these, clear-sky sunshine per month as a share of flat ground.
- Per band: the same, averaged over the band's area, using the RGI7 glacier DEMs in
  `/itet-stor/jabeer/glogem/geometricdata/rgiv7/dem/centraleurope/`.
- Fixed now: sites above 3300 m, and "wet" means a 10–20 m mean above −1 °C.
- Question: does adding summer sunshine to the air temperature improve the wet/dry ranking? Each
  massif is left out in turn, and the permutation shuffles whole glaciers.
- Pass: ranking skill +0.05 or more and p < 0.05, **with band sunshine**. Passing only with site
  sunshine means a band model cannot use it, so it becomes a stated limitation.
- Raster tools: we add rasterio (and pyproj) to the icetemp environment, as you decided.

*(c) Data fixes, before any run (icetemp `data.py`).*
- Assign each hole to the RGI7 outline that contains it. Use the RGI7 polygons (the masks are too
  coarse at the Colle Gnifetti divide); fall back to the GLIMS id only when no outline contains the
  hole.
- Split entities that mix outlines.
- Let entities in if they have a reading below 15 m *or* at 10–30 m.
- Print a table of every hole that is left out, and why, so that no data disappears silently again.
- Effect: Chessjengletscher (01161, sealed), 3 Colle Gnifetti holes and the 7 deep tongue profiles.
  The scored set and the Step 5a set are recomputed once; the old numbers stay in the book as they
  were.
- Regression check: with the new rules switched off, the loaded data set is identical to today's
  cache.

*(d) Sunshine, residual part.* Do the errors below 15 m of the adopted forcing, at its chosen f_ref,
follow the summer sunshine? This runs after Step 4.

**Step 2: GloGEM (branch GloGEMflow).** One switch, default = today.
- `firnice_thermal_lapse = 'floor' | 'reanalysis'`. The change is one line at
  `initialise_firnicetemp_spinup.pro:28` (`dtdz_icetemp = dtdz`). Every use goes through it: lines 35
  and 64, `firnice_reference_forcing.pro:113`, `firnice_temperature_model.pro:249` and the spin-up.
- If the station test picks the floor, the thermal melt in the 'n' branch of
  `firnice_thermal_forcing.pro` uses it too. That makes F consistent.
- Same settings guard and runner/config pattern as `firnice_forcing_toff` (GloGEM `8d10930`, icetemp
  `0f78c85`).
- Harness:
  - With the default, the output must be byte-identical.
  - With 'reanalysis', only the surface temperature in `firnice_snowclim_*` at or below the climate
    cell must stay identical. Horizontal advection couples the bands, and line 64 changes the
    diagnostic `firnice_temp_validation_*` everywhere, so those two checks do not apply.
  - The mass-balance files must stay identical.
- If you decide to build sunshine: a switch `firnice_thermal_sun`, which gets its own short plan.

**Step 3: icetemp.**
- `full_set_runs.py`:
  - add `--fref` and `--fins`;
  - build the run tag from variant, f_ref and f_ins, and refuse to run into an existing folder with
    other settings, because an old `.done` file would silently skip the run (runner.py:691);
  - write the real parameters to `design.csv`, which the scorer reads;
  - pass the new switch on through runner, config and pipeline, and extend the training-config test.
- New `forcing_fair_test.py`, reusing the loading code of `forcing_runs_score.py`:
  - the metrics from section 3, item 5, written in before the runs;
  - the massif list, fixed now: Mont Blanc, Monte Rosa, Mischabel–Saas, Bernese Alps, Valais west
    and Diablerets, Titlis, Écrins, Bernina, Ortles, Ötztal;
  - the forcing and f_ref chosen leaving one massif out;
  - paired differences per massif;
  - the Colle Gnifetti same-site gap and the dry-firn check;
  - the rule above.

**Step 4: the runs.**
- Variants H and F*H at f_ref {0.10, 0.25, 0.40, 0.55, 0.70, 0.85, 1.00}, on all glaciers. Then, at
  each variant's chosen f_ref, f_ins at 0.785 and 1.25.
- That is about 17 new runs of about 30 min each, run in batches as the IDL licences allow.
- Gate first. The new code at the defaults must reproduce the clean `f6c1453` runs (`ce25_R`,
  `ce25_H`, `ce25_FH`) on the three gate outlines (00786, 00774, 03116). That also proves that reusing
  H and FH at 0.55 is valid.

**Step 5: Taconnaz.**
- About six runs of a few minutes each, on the gate catchment, which reproduces the full set: the
  2×2 above plus f_ref 0.1 and 0.2.

**Step 6: time test.**
- On the existing runs (today, F, H, FH) in week 1, and again on v1.0.
- Pairs and noise level as in section 3, item 1, written in before computing.

**Step 7: the 3000–3500 m work package (section 3a).**
- Week 1, Python on the existing runs: the three causes, per entity.
  - Winter snow (`firnice_snowmay_<id>.dat`) against measured winter balances and Gilbert et al.
    (2012).
  - Wet depth (`firnice_wetdepth_<id>.dat`) against the measured profile shapes.
  - The firn line through time (`firnice_firnline_<id>.dat`) and the band temperatures across the
    year each band lost its firn (`temp_ID*_<id>.dat`), against the Tête Rousse repeats.
  - The frozen-bed check at the 44 holes that reach the bed.
- Week 2, on a catchment of 00777, 02103, 01265 and 01112 (minutes per run): f_ref, f_ins and f_adv
  sensitivity, to see which lever can make the clear misses cold without breaking the
  near-melting-point sites.
- Then a short plan for the physical fix, if one is identified.

**Step 8: your decision.** One table: the out-of-sample scores per variant, with the rules beside
them, the lapse-rate and sunshine tests, the time test and Taconnaz. You decide the v1.0 forcing,
and whether to build sunshine.

**Afterwards (its own plan):**
- Alpine ensemble with the adopted forcing → T2 calibration → Z2–Z4 → v1.0.
- Book: chapter 8, and a dated correction in chapter 6.
- For the mass-balance side, through you:
  - the model-changes text;
  - the T_off–height finding;
  - the aspect-class shift (see below);
  - a warning that `potential_solarradiation.pro` smooths the band slope in place, which would
    change the flow if model '3' were switched on.

*The aspect-class shift, checked.* For each 10 m band I computed the real aspect from the glacier
DEM and compared it with the class stored in GloGEM's band files:

| band files | class → measured angle from north (median per class) | coding | agrees with |
|---|---|---|---|
| RGI6 (8 glaciers, 1631 bands) | 1 → 9°, 2 → 29°, 3 → 104°, 4 → 148°, 5 → 162° | 1 = N … 5 = S | the code (66 % of bands) |
| RGI7 (12 glaciers, 1948 bands) | 0 → 16°, 1 → 32°, 2 → 84°, 3 → 145°, 4 → 162° | 0 = N … 4 = S | the code shifted by one (58 % of bands) |

- GloGEM's code expects 1 = N … 5 = S (`asp_class[aspect-1]`, class 0 set to 1).
- Read with that code, every RGI7 band except the north-facing ones is taken one class (45°)
  further north, and "south" never occurs.
- **This is not an RGI problem.** RGI gives one aspect per glacier. The per-band classes come from
  GloGEM's own geometry preprocessing for RGI7, which switched to counting from 0.
- So the report goes to whoever built the RGI7 band files (the GloGEM side), not to RGI.
- It has no effect on any run so far, because only the unused melt model '3' reads the class.
  It would matter as soon as anyone uses band aspect, including a sunshine factor built on the
  classes. Our sunshine test therefore computes aspect from the DEM.
- Production:
  - cost benchmark;
  - output design: per band and year, the regime, the temperature at 10, 15, 20 and 50 m and at the
    bed (frozen or not), and the wet depth; per glacier and year, the thermal type (temperate,
    cold, polythermal with frozen front or frozen upper part);
  - hindcast;
  - projections with posterior parameter sets, so the f_ref–f_ins trade-off shows up as forecast
    spread.

## How to run this plan: models, effort, ultracode

What is possible:
- **This session's model and effort** are yours to set (`/model`, the effort setting, `/fast`).
  I cannot switch them myself in the middle of the plan.
- **Work I hand to sub-agents** can get its own model (Opus, Sonnet, Haiku, Fable) and, when you
  ask for it, its own effort level. If you approve the table below, I follow it for every delegated
  step. That is the "automatic" part that works.
- **Ultracode / workflows** run many agents in parallel from a script. They pay off for wide,
  independent work and for adversarial checking. They do not pay off for sequential steps that
  wait on your decisions or that edit the one GloGEM tree, where parallel edits are a real risk:
  remember the concurrent-session incidents. They also cost several times the tokens.

My honest read of this plan:
- The bottlenecks are IDL wall time and licences, your decisions, and small code changes that must
  be exactly right. They are not breadth of thinking.
- So ultracode is worth it at **two points only**:
  - the pre-registration review before the fair test;
  - the v1.0 freeze review.

  A single reviewer already found six serious problems in this plan; a small panel with different
  angles (statistics, glaciology, code) is the cheapest insurance we can buy before weeks of runs.

| step | what it needs | recommendation |
|---|---|---|
| pre-registration: metrics, margins, massifs, time-test pairs, frozen-bed measure (before Step 4) | getting it right once; mistakes here cost weeks | **ultracode**: a small review panel of 4–6 agents (Opus, max) + one synthesis; you decide |
| data fixes in `data.py` (Step 1c) | exact, tested code | main session, high/xhigh; regression test; `/code-review` afterwards |
| GloGEM switch + harness (Step 2) | small, byte-identity checks | main session, high; no parallel edits |
| week-1 Python diagnostics: stations, sunshine, time test, 3000–3500 m, frozen bed | independent, read-only | 2–4 background sub-agents in parallel (Opus, high); I check each result. A workflow is not needed. |
| literature to read (Gilbert 2012 & 2015, Machguth 2016, Mattea 2021, Giesen & Oerlemans 2012) | careful extraction of numbers and caveats | one research sub-agent (Opus, high) |
| launching and watching IDL runs | mechanical | low effort is enough (Sonnet or Haiku sub-agent, or this session at low/medium effort on run days) |
| physics fixes: sunshine energy balance, cold-firn percolation, frozen-bed logic | new physics that must be defensible | main session, **max**, with its own plan and one reviewer |
| v1.0 freeze review and the paper's validation section | everything checked together | **ultracode** review panel (Opus, max) |
| book chapter, model-changes text, paper drafts | clear writing for others | main session, high |

- **Fable:** I have no reliable basis to judge how it compares on this kind of work. I would not
  use it for decision-critical steps without a small trial.
- **Day to day:** max effort in this session is right on design and decision days. On days that
  are mostly launching and watching runs, a lower effort saves time and tokens without risk.

## Decisions for you

| topic | recommendation |
|---|---|
| judge F by a fair run set: own f_ref, leave-one-massif-out, margins 0.3 °C and 0.05 | yes. The margins are yours to set. |
| decide the lapse rate with measured mountain air temperatures, not inside the run set | yes |
| sunshine: test first, build for v1.0 only if it passes with band sunshine | yes; then your call: v1.0 (+2 weeks) or v2 |
| one gate (Z1) and a time box of 3 weeks (5 with sunshine); Z2–Z4 measured once; stop sign on Z3 | yes |
| keep Chessjengletscher (01161) sealed as the one unseen glacier, in the 3000–3500 m zone | yes |
| data fixes before any run: holes to outlines by position, entities split by outline, deep profiles in (a reading below 15 m *or* at 10–30 m) | yes; all numbers recomputed once |
| 3000–3500 m zone: diagnosis in week 1, sensitivity runs in week 2, its own margin in the fair test; a physical fix goes into v1.0 if one is found | yes (hazard priority) |
| report the aspect-class coding to the GloGEM side (not RGI) | yes, through you |
| regime measures: three classes (polythermal recall as the headline) plus the frozen-bed check; two-class only beside them | yes |
| models and effort per step as in the table above; ultracode only for the two review panels | yes; your call |
| tell the mass-balance side: 93 % of the high glaciers with measured targets sit at both bounds, median T_off +3.1 °C | yes, through you |

## Verification

- **Step 1:**
  - The rebuild is `forcing_check.py`'s, validated to 0.005 °C.
  - The sunshine code reproduces flat-ground clear-sky values and the per-site numbers of the
    search.
  - The station comparison uses the same months and elevations as the station records.
- **Step 2:** the harness as specified above.
- **Step 3:**
  - The leave-one-massif-out code is checked on a case with a known answer.
  - A run tag clash refuses to run.
  - `tests/` pass.
- **Step 4:**
  - Gate equal to the clean `f6c1453` runs.
  - All entities have profiles.
  - Mass-balance files unchanged across variants.
  - Spin-ups converged.

## Housekeeping

- GloGEM edits only on GloGEMflow, in small commits, and nothing changes while runs are going.
- IDL licences are checked before every launch; no process is killed.
- No config file is committed. Your yaml, line 38, still says `explained_variance: 0.99` and needs
  `0.999`.
- Figures only if you ask, and then under `glogemflow_development/figures/<topic>/`.

