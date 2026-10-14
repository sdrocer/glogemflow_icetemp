# Design note: a probabilistic thermal-regime product (Step 4 of the plan of 2026-10-14)

*Design only; computed at production. Janosch's addition of 2026-10-14. Daniel's 2024 plan of
action proposed a "probability of being temperate"; its definitions are to be matched once the
document is found.*

## Why

The model is weakest exactly where beds and firn lie just below the melting point:
- half the sites are within 0.3 °C of it;
- class changes between model versions often flip with the tolerance.

A hard class (cold / polythermal / temperate) states more than we know there. A probability says
how sure the model is, and the check below shows whether that claim holds.

## What the product gives

Per glacier, per elevation band, per year (or 5-year mean):

| quantity | definition |
|---|---|
| P(cold), P(polythermal), P(temperate) | share of ensemble members whose profile below 15 m is classified so (`thermal_structure.classify_profile`, the same rule as for the boreholes) |
| P(frozen bed) | share of members with the bed node colder than the melting point minus the tolerance (bed log `firnice_bed_<id>.dat`, runs without deepening) |
| cold-ice elevation range | the bands where P(cold ice at 10–20 m, or at the bed) ≥ 0.5, with the 0.1 and 0.9 contours as its uncertainty |
| glacier type probabilities | share of members with a frozen front and a temperate upper part (the hazard type), the reverse, cold-based, or temperate-based (`bed_types` rules) |
| year of change | first year in which P(temperate) at a band passes 0.5, with the spread across members |

## Where the spread comes from

```
 member = (posterior parameter set) x (forcing) x (measurement tolerance)
   posterior sets: N draws of (f_ref, f_ins, f_adv) from the Bayesian calibration (e.g. N = 50)
   forcing:        FRH (central); H as a structural sensitivity, shown separately and as a mixture
   tolerance:      0.05 C / the hole's own / 0.5 C, as in the scoring (averaged)
   projections:    x GCMs x SSPs (kept apart: scenario is not uncertainty)
```

## Do the probabilities mean what they say? (the check, on held-out massifs only)

- For every borehole entity, the probabilities come from the ensemble calibrated without its
  massif (leave one massif out, as in the fair test).
- Scores against the measured classes:
  - the Brier score per class, and its skill against climatology (the class frequencies);
  - a reliability diagram (predicted probability against observed frequency, in bins);
  - the same for the frozen bed (the 14 modelled bed holes).
- Parameter spread alone usually underestimates structural error. If the probabilities are
  overconfident, they are widened in one physically motivated way: a temperature perturbation
  ε ~ N(0, σ) at classification, with σ = 1.5–1.9 °C, the neighbouring-borehole floor (chapter 8).
  This value is fixed now, not tuned to the scores.
- If still unreliable, the product reports the gap and keeps the hard classes beside the
  probabilities.

## Cost and output

- **Cost:** N × 2 runs per glacier and period. The cost benchmark of the plan (one batch of about
  200 Alpine glaciers) sets N and the time resolution of the output.
- **Output:** one netCDF per region, with variables (glacier, band, year), CF conventions, and
  provenance (GloGEM commit, posterior file, forcing).

## Open

- Daniel's 2024 definitions (temperate probability: of what, at which depth?).
- Whether P(polythermal) should be split into "cold over temperate" and "temperate over cold"
  (layer order, as in the fair test's reported table).
