# Step 5: the Bayesian calibration on FRHS (draft for the review panel, 20 Oct)

*Draft. The review panel reads this together with the book (chapters 1-11), the plans and the code; its findings go into
the version Janosch reviews. Terms in the glossary.*

## Context

- **The project's goal** (Janosch, 19-20 Oct): a physically sound, calibrated and defensible englacial temperature model
  in GloGEM that gets the thermal regime (cold, polythermal, temperate) and frozen beds right wherever there are glaciers,
  with the firn-no-firn zone 3000-3800 m of particular interest.
  - **The deliverable:** a probability map, P(class) and P(frozen bed) per glacier, band and year, for the past and the
    projections.
  - **It is checked by three ways of leaving data out:** one massif (7 folds, transfer in space), one glacier system
    (outlines leak), one period (older readings fit, newer tested: transfer in time, closest to a projection).
  - **The paper** must defend the model, the calibration and the map.
- **The physics is final for v1.0** (decided 20 Oct):
  - the central thermal variant is FRHS: no T_off in the thermal inputs, ERA5's lapse rate, Herron & Langway meltwater,
    latent heat per m of water, and meltwater refreezing in the seasonal snow on ice bands (S);
  - the firn ledger, the station lapse rates and the full refreezing model were tested and not adopted (chapters 10-11).
- **Known structural errors,** to be stated, not fixed before v1.0:
  - **The ice geometry is fixed at 2003** for the spin-up (1942-71 climate) and the run (1971-2025). The model ice is a
    median half the drilled thickness at the holes that reach the bed, but this does not line up with the warm misfits
    (test 1a, rank correlation 0.07).
  - **GloGEM's mass balance before 2000 is 0.6-0.8 m w.e. a⁻¹ too positive** against the 7 long GLAMOS series, mostly
    winter snow. Correcting it (test 1c) leaves the zone and the sites unchanged and moves only the cold high firn,
    by 0.5-1 C in opposite directions on Mont Blanc and Monte Rosa; overall M1 2.51 -> 2.61 C.
  - **ERA5 is about 5 C too cold in winter at 3300-3600 m** (chapter 8).
  - **One Alps-wide f_ref fits the massifs differently:** in stage 2 the massifs chose 0.85, Mont Blanc 1.0, Monte Rosa
    0.55.
  - **The floor:** a band model cannot come closer to one borehole than about 1.5-1.9 C.
- **Where FRHS stands** (stage 2, out of sample with f_ref chosen leaving one massif out, f_ins fixed 1.017):
  - M1 2.51 C, zone M1 2.51 C;
  - 61 of 81 cold, 9 of 18 polythermal and 2 of 6 temperate bands right;
  - 11 of 14 measured frozen beds found.
- **Every Step-5 number so far** used the old forcing and 7 glaciers (chapter 7).

## What Step 5 does

```
  design runs (GloGEM, FRHS)  ──►  emulator (fast stand-in)  ──►  likelihood + priors  ──►  MCMC: posterior of f_ref, f_ins
   40 + 20 parameter sets           PCA + Gaussian process        readings per hole,           ──►  held-out checks, then the
   23 glaciers, 105 bands                                         year and 1 m depth bin            probability map
```

1. **Parameters:**
   - f_ref (share of the meltwater that refreezes and heats): uniform 0.1-1.0;
   - f_ins (snow insulation): truncated normal, mean 1, sd 0.35, range 0.3-2.0;
   - f_adv fixed at 1.0.
2. **Data:**
   - glenglat CentralEurope, holes on the RGI7 outline that contains them, CentralEurope_fair without the sealed
     Chessjengletscher (01161): 23 glaciers, 122 entities, 198 holes, 23 sites;
   - one point per borehole, calendar year and 1 m depth bin (2115 points at all depths, 1135 below 15 m);
   - the Colle Gnifetti points below 40 m are left out of the fit (scored apart).
3. **Target 'depth' (T2):**
   - every point is compared with the model at its depth, averaged over its readings' months;
   - the per-entity offset is estimated after the fit, not inside it.
4. **The error model** (new, chapter 11 and Step 3):
   - r = model - reading ~ N(0, D + C);
   - D = each reading's own variance (its hole's uncertainty; ×3 for estimated profiles);
   - C = band s2·exp(-|Δz|/ℓ) + borehole sh2 + borehole-and-year sy2 + site ss2.
   - **The sizes are fixed before the campaign,** estimated by REML from the stage-2 FRHS residuals (all depths, each
     massif at its stage-2 f_ref): s2 2.04 C² over ℓ 17 m; sh2 2.60 C²; sy2 0.05 C²; ss2 2.10 C²; free mean +0.42 C.
   - **Profiles (2 log-likelihood units):** s2 1.80-2.31; ℓ 13-22 m; sh2 2.30-3.34; sy2 0.03-0.07; ss2 0.99-4.46.
   - **Simulation check** (200 data sets simulated from the estimate on the real layout, refit): every size recovered
     without bias. 5-95 %: s2 1.7-2.4, ℓ 14-21 m, sh2 2.1-3.2, sy2 0.03-0.07, ss2 0.7-4.0. The one ridge is s2 with ℓ
     (correlation 0.9).
   - **By depth:** fit separately, the points at or above 15 m want ℓ at its lower bound (1 m; sh2 3.6) and those below
     15 m want ℓ ≈ 120 m (s2 1.7, sh2 3.4, ss2 3.5). The single exponential is a compromise between a short seasonal
     and a long deep correlation.
   - **Leaving out the site term** costs 26 log-likelihood units. Today's model error (3.5 C² over 15 m, nothing shared)
     fits 12 638 units worse.
   - **Without the two Grenz sites,** the shape changes (s2 3.2 over 56 m, sh2 0.8).
   - **Sensitivity arms,** each a refit of the same emulator, no new runs: A0 the estimate; A1 no site term; A2 site
     term at the top of its interval (4.46); A3 ℓ at the far end of its interval (13.3 m, s2 1.73); A4 estimated without
     Grenz; REF today's.
5. **Design and runs:**
   - stage 1: 40 Latin-hypercube runs over the prior box;
   - emulator, then MCMC;
   - stage 2: 20 runs within the stage-1 result ± 2.5 sd;
   - emulator and MCMC again;
   - each run takes 20-36 min, spread over the hosts with free cores (12 workers), so about 2 h per stage.
6. **Emulator:** PCA keeping 99.9 % of the variance, a Gaussian process per component, leave-one-out check of its
   variance. The fit uses the emulator mean; the fit with the emulator's uncertainty varying is an adoption check.
7. **MCMC:** emcee, 32 walkers, 5000 steps, 1000 burn-in.

## How it is judged (fixed before the runs)

- **The gate:**
  - R-hat ≤ 1.05;
  - at least 5 design runs within 0.5 standardised units of the posterior mean;
  - misfit ratio (Mahalanobis per point) 0.5-2;
  - mode gap ≤ 5 log units;
  - new: the standardised site offsets inside their expected range, and no site carrying more than 25 % of the misfit.
- **Adoption:** held out by glacier system (each system predicted from a fit to the others):
  - the RMSE beats the prior-mean run;
  - it is at most T1's + 0.10 C (T1 = the older one-number-per-entity target);
  - and the emulator check passes (the fit with emulator uncertainty moves by at most 1 posterior sd).
- **Then one real GloGEM run at the posterior mean,** scored like stage 2 (per band at depth, sites, frozen beds)
  against FRHS at f_ref 0.85, f_ins 1.017.
- **Reported, not gating:**
  - leave one massif out;
  - leave one period out (older readings fit, newer tested; needs new code);
  - the sensitivity arms A1-A4 and REF;
  - classes and frozen beds per band at depth;
  - every band's depth profile.

## After Step 5 (not in this plan)

- The probability map: P(class) and P(frozen bed) from posterior draws, plus the tolerance bracket and the
  structural-error spread; reliability (Brier score, reliability diagram) in each leave-out scheme.
- The projections; the v1.0 freeze; the paper.

## Open questions for the panel (from the author)

1. Is a single Alps-wide (f_ref, f_ins) defensible when massifs prefer different f_ref, or should the model error carry
   it (and how would a reviewer see that)?
2. The error sizes come from residuals of the same data the calibration then uses. Is that defensible (empirical Bayes),
   and how should it be reported?
3. Are 23 glaciers / 23 sites in one region enough to claim "wherever there are glaciers"? What can the paper claim, and
   what must be framed as a regional result?
4. Is the adoption rule (beat prior mean and T1 + 0.1 C held out by glacier system) the right benchmark, or should it be
   against simple baselines a reviewer would expect (e.g. elevation/temperature regression, k-NN, MAAT-based rule of
   thumb)?
5. Which of the known structural errors could sink the paper, and what sensitivity runs are needed to pre-empt that?

## Glossary

| term | meaning |
|---|---|
| FRHS | the central thermal variant: FRH (no T_off, ERA5 lapse rate, Herron & Langway meltwater) with S (refreezing in the snow on ice bands) |
| f_ref, f_ins, f_adv | refreezing share of meltwater; snow insulation scale; advection scale |
| entity, band, site | site × 10 m band × outline; site = one borehole location (23); borehole = one hole (198) |
| T1, T2 | calibration targets: one number per entity; one per 1 m depth bin |
| emulator | statistical stand-in for GloGEM trained on the design runs |
| REML | estimating the error sizes from residuals while allowing a free mean |
| M1 | RMSE over the readings below 15 m, each outline weighted equally |
| zone | 3000-3800 m, the firn-no-firn transition |
| floor | how close a band model can come to one borehole, about 1.5-1.9 C |
| glacier system | glaciers whose outlines share sites (e.g. Monte Rosa), held out together |
