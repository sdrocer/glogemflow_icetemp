# Prompt for the presentation assistant

Copy everything below the line into the Claude session preparing the slides.

---

I am preparing a conference talk for the Global Glacier Modelling Workshop in Obergurgl
(2026-08-24) about a Bayesian calibration scheme for englacial ice temperatures in GloGEM.

**Read `glogemflow_icetemp/TALK_FACTSHEET.md` first and treat it as the only source of numbers.**
Every figure in it was measured on real campaign output. Do not introduce any number that is not
in that file, do not round its numbers "for clarity", and obey its DO NOT SAY list. If you think a
number is missing, ask me — do not estimate one.

## Slides

There are **two figures already made**. Use them; do not design new ones.

1. `glogemflow_icetemp/figs/calibration_scheme_explained.png` (`.pdf` also available)
   Panel A: the three-tier ladder. Panel B: one real borehole — measured vs modelled profile with
   the gap shaded, over the line *measured = model + systematic bias + noise*.
2. `glogemflow_icetemp/figs/paramspace_target_vs_posterior.png` (`.pdf` also available)
   Hatched = parameter settings where the model meets the target; orange cloud = where the
   calibration posterior actually lands. They barely overlap.

**Slide budget: 2, absolute maximum 3.** The story below is longer than 2 slides of content —
that is intentional. Most of it is the **spoken** narrative; the slides are visual anchors, not
the script. Build speaker notes for the narrative and keep the slides sparse. If you believe a
third slide is genuinely needed, propose it and say what it earns.

## The story I want to tell, in this order

**1. Why this model is hard to calibrate.** Lead with the difficulty, not the method.
- Only 25 Alpine glaciers have borehole temperatures; ~4000 need parameters.
- The forward model is expensive IDL, so it cannot be run inside a sampler — an emulator stands
  in for it.
- The model is genuinely wrong in structured ways (not just noisy), and that error varies from
  place to place, so it has to be modelled rather than ignored.
- **The deepest problem, and the punchline of the whole talk:** what we can *measure* is
  continuous temperature, but what we *want* is a categorical answer — is this glacier
  polythermal or not. Those two are not the same target, and this work is largely the story of
  discovering how far apart they are.

**2. What we calibrate, and what those parameters physically mean.** Two parameters:
- `perm_frac` — scales how deep meltwater percolates before refreezing. It sets the depth limit
  of the latent-heat release in the ice column (it scales the Herron–Langway firn-ice transition
  depth, typically 10–20 m near the ELA and 50–80 m at high-accumulation sites). Physically: *how
  far down does refreezing warm the ice.*
- `dT_scale` — scales how much warmer the firn surface is than the air, because refreezing
  meltwater releases latent heat there. It multiplies an empirical per-elevation-band correction
  in the surface boundary condition. Ice-covered bands get 40% of the firn value, since seasonal
  snow insulates less than perennial firn. Physically: *how strongly does the surface decouple
  from air temperature.*
- A third parameter, `z0`, was dropped: we verified it cannot reach the physics at all in the
  current model, and its posterior never narrowed in any campaign. Worth one sentence — it is an
  honest negative result.

**3. What came out, and what I tried.** Chronological, honest:
- Early campaigns pooled all boreholes on a glacier into one entity, which compared a deep tongue
  borehole against the model column at the summit. Fixed by splitting into 10 m elevation bands.
- Several physics bugs in the temperature module were found and fixed along the way.
- An adversarial review of the scheme then found the real problem: **the likelihood was
  over-confident by a factor of 13**. It claimed to know the answer far better than it did, which
  pinned the posterior against its own bounds and excluded better answers with false confidence.
- The fix: an explicit model-error term, carrying the emulator's uncertainty properly instead of
  as a diagonal, and adding elevation to the spatial correction. The calibration check went from
  13.0 to 1.34 and the uncertainty widened ~186×.
- A stricter adoption rule was added: four gates, then an absolute performance bar.
- Result: **best temperature predictor of the three methods, worst structure predictor.** It fails
  its own bar for finding polythermal glaciers.
- Why: 82 of 111 boreholes are cold, so the temperature fit is maximised by being cold
  everywhere. Making the statistics honest made it *more* confident in the cold answer. **The
  remaining gap is physics, not statistics.**
- The model *can* hit the target — 113 of 625 parameter settings clear it — so it is a selection
  problem. That is figure 2.
- Next: a thermal spin-up that is already written but switched off, and coupling the physical
  refreezing model into the ice column (it exists, but is only 8 m deep and disconnected).

## Tone

This is a talk about catching and fixing my own method's flaw, not a results-victory talk. The
honest framing is stronger with this audience and should not be softened. Do not oversell; do not
add hedging either.

## What I want from you

1. A slide-by-slide outline (2 slides, or 3 with justification) placing the two figures.
2. Speaker notes carrying the narrative above — this is the main deliverable.
3. Three or four anticipated questions with short answers, drawn only from the fact sheet
   (the sensitivity table answers "why 3.5?"; the class-balance argument answers "why is your
   accuracy worse than doing nothing?").

Ask me before inventing anything you cannot source from the fact sheet.
