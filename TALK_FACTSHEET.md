# Fact sheet — Tier-3 Bayesian calibration talk (Obergurgl, 2026-08-24)

**Give this to any assistant helping with the slides. Every number below was measured on the real
campaign artefacts. Do not add numbers that are not here, and do not round them "for clarity".**

---

## The two figures (already made — do not invent new ones)

| slide | file | what it shows |
|---|---|---|
| 1 | `figs/calibration_scheme_explained.png` / `.pdf` | A: the three-tier ladder. B: one real borehole — measured vs model, gap shaded. |
| 2 | `figs/paramspace_target_vs_posterior.png` / `.pdf` | Where the model meets the target (hatched) vs where the calibration posterior actually lands (orange cloud). |

**Hard cap: 2 slides, ideally 1.** Do not propose additional figures. If cut to one slide, keep
panel A of figure 1 and speak the three headline numbers below.

---

## The problem

- Goal: find where polythermal glaciers are in the Alps (glacial-hazard motivation, Tête
  Rousse-type frozen-front glaciers).
- Only **25** Alpine glaciers have borehole temperatures (glenglat). Parameters are needed for
  **~4000**.
- After splitting each glacier into 10 m elevation bands, the calibration set is **116 entities**
  / **111 with a thermal-structure label** / **107 scorable by the emulator**, holding **695**
  usable depth measurements.
- Class balance: **82 of 111 cold, 29 warm** (27 strictly polythermal, 2 temperate).

## The scheme

- **Tier 1** — grid-search the temperature model against each borehole glacier separately.
- **Tier 2** — a transfer model predicting parameters from climate + elevation, so unmeasured
  glaciers get values. (**Keep Tier-2 in the talk**: Tier-3's exported output is *defined* as a
  correction on top of it — the file the IDL model reads contains deltas relative to Tier-2.)
- **Tier 3** — Kennedy–O'Hagan Bayesian calibration: a fast emulator of the real IDL model, plus a
  spatial discrepancy term, plus observation error. Sampled with emcee.
- Two calibrated parameters: `perm_frac` (meltwater percolation depth scaling) and `dT_scale`
  (surface firn-insulation scaling). `z0` is held fixed — verified it cannot reach the physics.

## THE THREE HEADLINE NUMBERS

1. **The likelihood was over-confident by 13×.** A calibration check that should read ~1 read
   **13.0**. After the fix: **1.34**.
2. **Uncertainty was understated ~186×.** `perm_frac` went from **±0.0009** (pinned against its
   bound) to **±0.169**.
3. **Best temperature predictor, worst structure predictor.** Leave-one-out RMSE **1.99 °C**
   (best of three, adopted). But it finds only **6.9%** of warm glaciers where a simple method
   finds **55.2%**.

## Performance, leave-one-out, graded on the real model

| method | warm found | cold kept | balanced | plain accuracy |
|---|---|---|---|---|
| Tier-2 | 55.2% | 65.9% | **60.5%** | 55.9% |
| k-NN | 37.9% | 67.1% | 52.5% | 55.9% |
| **Tier-3 (KO)** | **6.9%** | 95.1% | 51.0% | 70.3% |
| do nothing ("all cold") | 0% | 100% | 50.0% | **73.9%** |

**The methodological point worth making:** plain accuracy is meaningless here. Calling everything
cold scores **73.9%**; our calibration scores **70.3%** — it was performing *worse than doing
nothing* while looking respectable. Use the balanced score (warm found + cold kept)/2.

## Adoption rule (gate, then bar)

Gate — all four must pass; campaign 8 passes all four, its predecessor failed:
R-hat ≤ 1.05 (**1.030**) · ≥5 real training points near the answer (**11**) · calibration check in
[0.5, 2.0] (**1.342**) · sampler within 5 log-units of an independent grid search (**−0.08**).
Bar — must clear an absolute floor (warm found **and** cold kept both > 50%) *and* beat both
baselines on the balanced score. **KO fails the bar.**

## Sensitivity — "how did you pick 3.5 for the model-error term?"

Re-ran the whole calibration at 2.3 and 5.0 (the span of three independent estimates that
converged on ~3.5). The conclusion does not move:

| model-error s² | calibration check | R-hat | posterior perm_frac | posterior dT_scale | warm found |
|---|---|---|---|---|---|
| 2.3 | 1.96 (edge of band) | **1.07 — fails** | 0.189 ± 0.172 | 1.303 ± 0.495 | gate failed |
| **3.5 (used)** | **1.34** | 1.03 | 0.192 ± 0.169 | 1.244 ± 0.535 | 6.9% |
| 5.0 | 0.97 | 1.01 | 0.198 ± 0.163 | 1.118 ± 0.559 | 6.9% |

Two things to say if asked:
- The answer is **insensitive**: at 5.0 the classification numbers are *identical* to 3.5
  (6.9% / 95.1% / 51.0%), and the posterior barely moves across all three.
- 2.3 is demonstrably **too small**: the calibration check sits at the edge of the acceptable band
  and the chains stop mixing (R-hat 1.07 > 1.05). That is the expected signature of a likelihood
  that is still too sharp — it makes the surface rough and hard to sample.

## Why it fails, and what is next

- 82 of 111 boreholes are cold, so the temperature fit is maximised by being cold everywhere.
  Making the statistics honest made the model *more* confident in the cold answer.
  **The remaining gap is physics, not statistics.**
- **The model can hit the target**: of 625 parameter settings swept, **113 clear both-above-50%**,
  best at **75.9% warm found / 73.1% cold kept**. It is a *selection* problem.
- Across the posterior, the probability of clearing both-above-50% is **3.3%**.
- Next steps: (a) the thermal spin-up — already written, currently switched off; (b) coupling the
  physical refreezing model into the ice column — it exists but is only 8 m deep and is entirely
  disconnected from the ice-temperature module.

---

## DO NOT SAY

- **"62.1% polythermal recall"** — does not reproduce. The honest number at that setting is 24%.
- **Campaign 7's posterior as a result** — it never passed its own gate.
- Any claim that Tier-3 currently beats the baselines at finding polythermal glaciers. It does not.
- Do not describe the emulator's uncertainty as "conservative" — it was anti-conservative,
  understating its leading direction by 242×, now fixed.
- Do not quote in-sample sweep numbers as validated performance; the 75.9/73.1 pair is in-sample
  and shows *capability*, not *skill*.
