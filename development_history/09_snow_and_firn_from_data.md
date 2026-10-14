# 9. Snow and firn from data: what the mass balance hands the thermal model

> **In short.** Chapter 8 traced the misses at 3000–3500 m to too much winter snow and too much
> firn on small high glaciers. The same excess snow decides whether meltwater refreezing in the
> winter snow on ice keeps Tête Rousse's frozen bed.
>
> Letting the precipitation factor fall below its 1.3 bound brings the winter snow to GLAMOS's
> values, and with less snow Tête Rousse's frozen bed comes back. But the ELA moves only 10–50 m
> instead of the 120–200 m GLAMOS shows. With only one glacier-wide target, less snow is simply
> balanced by less melt. A second observation must pin where the snow survives: snow lines
> (Alexandra von der Esch's snowline calibration) and firn maps (Matthias Huss's firn model,
> Janosch's first paper).
>
> *Status (2026-10-14): Step 1 of the plan of 2026-10-14 done; the firn model and the snow-line
> scoring are running; Alexandra's code arrives as a GloGEM branch by the end of the week.*

## What we believed

- The mass-balance calibration (Huss & Hock 2015, as in GloGEM): one target per glacier, the
  geodetic balance 2000–2020. It sets the precipitation factor (bounded 1.3–2.3 in CentralEurope),
  then the degree-day factor, then the temperature offset.
- The thermal model takes from it the snow (insulation), the firn mask (where meltwater enters) and
  the accumulation (burial).
- We assumed the firn mask was good enough. Chapter 6 compared the ELA with WGMS for a few
  glaciers, but never for small high ones.

## What we measured

### The precipitation bound (Step 1; diagnostic, 4 zone glaciers)

Two calibrations of Tête Rousse (00777), Hohsaas (01112), Altels (01265) and Corvatsch (02103)
were run with GloGEM's own calibration, unchanged.
- **'ref':** bound 1.3, as the regional calibration.
- **'low':** bound 0.8, from a copy of the regional parameter table. The shared file and the code
  stay untouched.

The 'ref' run reproduces the regional calibration: c_prec identical, winter balance within 1 %,
ELA within 10 m. DDF and T_off trade off slightly differently, because on 4 glaciers the cell
means that seed T_off differ.

| glacier | c_prec | winter balance, m w.e. (GLAMOS) | ELA 2010s / 2020s, m (GLAMOS) |
|---|---|---|---|
| Corvatsch | 1.30 → 0.95 | 1.20 → 0.87 (≈ 0.83, 2013–25) | 3179 → 3192 / 3272 → 3283 (≈ 3380, 2013–25) |
| Hohsaas | 1.30 → 0.80 | 1.66 → 1.03 (≈ 0.97) | 3222 → 3247 / 3333 → 3345 (≈ 3460, 2018–25) |
| Altels | 1.30 → 0.80 | 1.55 → 0.94 | 2986 → 3043 / 3095 → 3143 |
| Tête Rousse | 1.30 → 0.80 | 1.63 → 1.00 | 3190 → 3201 / 3185 → 3202 |

```
 one target (the geodetic balance), two knobs that trade off:
   less snow (c_prec down) -> less melt needed (DDF down) -> where they balance (the ELA) barely moves
 so the winter snow comes right, the snow line does not
```

The thermal model on both calibrations (FRH, f_ref 0.85; below 15 m, 18 entities, 163 points):

| | 'ref' | 'low' |
|---|---|---|
| RMSE / bias, FRH | 0.96 / +0.65 °C | 0.85 / +0.51 °C |
| RMSE, with snow refreezing on ice (FRHS) | 1.08 °C | 0.95 °C |
| Tête Rousse hole 64 (32.5 m, measured −2.34 °C): FRH / FRHS | −0.43 / −0.07 °C | −0.91 / −0.39 °C |
| Corvatsch, 4 entities | temperate (bias +0.8 to +2.7 °C) | identical |
| robust class changes | | none |

- **Tête Rousse improves.** With less snow, the snow-refreezing switch keeps its frozen bed: hole 64
  is at −0.39 °C, below the melting point. So the switch failed in chapter 8 because of the excess
  snow, as suspected.
- **Corvatsch does not move at all.** Its firn line rises by only 15 m, so the borehole bands stay
  in firn and meltwater keeps them at 0 °C.
- **Hohsaas gets slightly too cold.** The bias at 3275 m goes from −0.86 to −1.21 °C.

## What was wrong

- **Calibrating to one glacier-wide target cannot fix where snow survives.** The precipitation
  factor sets how much snow falls; the snow line also depends on how the melt and the snow change
  with height. Lowering c_prec trades snow for melt and leaves the ELA nearly in place.
- **Possible cause, not yet tested:** the precipitation gradient (1.5 % per 100 m above the ERA5
  cell, the regional table). On small glaciers high above their cell it adds much snow at the top.
  The GLAMOS band winter balances of Corvatsch are nearly flat with height.

## What is next

- **The firn model** (`glogemflow_icetemp/firn_model/firn_change.pro`, Beer et al. 2026), driven
  by the GLAMOS distributed mass balances, on our 9 Swiss calibration glaciers. Its firn extent and
  thickness by band and year are compared with the thermal firn mask and the meltwater depth.
- **The snow lines:**
  - today's calibration against the Alpine transient snow lines (Sommer et al., ESSD preprint
    essd-2026-35) and the GLAMOS ELAs;
  - then Alexandra's snowline calibration, ported from her GloGEM branch into a procedure on
    GloGEMflow.
- **Firn thickness for the meltwater depth:** a bookkeeping per band with the firn model's rules
  (Step 3b of the plan).

## Sources

- Plan: `plans/2026-10-14_snow_firn_probabilistic_plan.md`, `plans/2026-10-14_plan_additions.md`
- Diagnostic calibrations:
  `/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week2/zonecal/{ref,low}/`; the
  regional-table copy is in `zonecal/data_lowc/`
- Runs: `data/bayescal/zcal{ref,low}_{FRH,FRHS}_r0.85_i1.017_a1.00[_nodeep]_w0`
- GLAMOS: `/itet-stor/jabeer/glazioarch/GlacioBaseData/MassBalance/glacier/{corvatsch,hohsaas}_fix.dat`
- Literature: Huss & Hock (2015); von der Esch et al. (in review, snowline calibration of GloGEM);
  Sommer et al. (2026, ESSD preprint); Beer et al. (2026, *The Cryosphere* 20, 5453–5473)
