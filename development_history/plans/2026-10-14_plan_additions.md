# Additions to the plan (2026-10-14, Janosch)

## 1. Snow and firn cover from data

The zone diagnosis of 2026-10-14 (chapter 8) traced the 3000–3500 m misses to too much winter snow
and too much firn. The same excess snow decides whether the snow-refreezing switch keeps Tête
Rousse's frozen bed.

- **Step 1, diagnostic.** Recalibrate the mass balance of the zone glaciers with the precipitation
  factor allowed below 1.3, against the same geodetic targets and checked against GLAMOS. Then
  rerun the thermal model on them, with and without snow refreezing. Write a plain-words note for
  the mass-balance side (Alexandra, Lander, Matthias).
- **Step 2, the firn model of Janosch's first paper** (with Matthias Huss). It is driven by GLAMOS
  mass balances and gives firn maps and firn depths for all Swiss glaciers. Uses:
  - check the thermal model's firn mask year by year;
  - replace the meltwater percolation depth (Herron & Langway close-off, 53–80 m, deeper than
    whole columns at the zone sites) by the real firn thickness. This is a criterion independent
    of the thermal state, which the failed depth rule lacked;
  - outside Switzerland and for projections: a firn-thickness bookkeeping inside GloGEM,
    calibrated and validated against it.
- **Step 3, the Alpine snowline product** (ESSD preprint essd-2026-35). Uses:
  - evaluate the model's snow line and firn line per glacier and year;
  - a snowline-constrained calibration of the snow (c_prec) and melt, with the mass-balance side.

## 2. A probabilistic product (deliverable)

The model is weakest exactly where beds lie just below the melting point, so hard classes say more
than we know.

- **What it gives,** per glacier and band and by year:
  - the probability of being cold, polythermal or temperate;
  - the probability of a frozen bed;
  - the elevation range where cold ice is likely.
- **Where the spread comes from:**
  - the posterior parameter sets of the Bayesian calibration;
  - the forcing (FRH central, H as structural sensitivity);
  - the measurement tolerance used for the classes (0.05 °C / the hole's own / 0.5 °C).
- **Check:** on held-out massifs, are the probabilities reliable (reliability diagram, Brier score
  against the boreholes)? Parameter spread alone usually underestimates structural error. If it
  does here, the probabilities are widened, or the gap is stated.
- Daniel's 2024 plan of action already proposed a "probability of being temperate". Its
  definitions are to be matched.
