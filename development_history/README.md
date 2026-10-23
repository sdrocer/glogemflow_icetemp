# Development history of the GloGEM ice-temperature model

How the firn/ice temperature module got where it is: what we believed at each stage, what the
measurements showed, what was wrong, and why each change was made. Written so the reasoning can
be defended later, in a paper, a thesis or a talk.

## Chapters

| # | chapter | period | one-line summary |
|---|---|---|---|
| 1 | [Starting point](01_starting_point.md) | Jun – 1 Sep 2026 | glenglat, the first calibration framework, `perm_frac` / `dT_scale` |
| 2 | [Parameter redesign](02_parameter_redesign.md) | Aug–Sep 2026 | `f_ref`, `f_ins`, `f_adv`; energy and mass kept separate |
| 3 | [Why the parameters look alike](03_why_the_parameters_look_alike.md) | Sep 2026 | `f_ref` and `f_ins` confounded; the depth window is too narrow |
| 4 | [The depth window](04_depth_window.md) | Sep 2026 | bed node, implicit solver, conduction stencil, column to the borehole |
| 5 | [The thermal spin-up](05_thermal_spinup.md) | 24–28 Sep 2026 | a one-year climate, a slow spin-up, and the drift correction |
| 6 | [Firn extent and the mass-balance calibration](06_firn_extent_and_calibration.md) | 28 Sep 2026 – | not water in ice but the firn line; the runs used other glaciers' mass-balance parameters |
| 7 | [The snow, the comparison and the calibration target](07_snow_comparison_target.md) | 30 Sep 2026 – | each reading at its date; the winter snow on ice barely insulated (a solver bug); the fix lowers the RMSE in all 32 runs; target T2 (every depth bin) tells f_ref and f_ins apart, coded and tried |
| 8 | [The thermal forcing, the data we were losing, and the tests a forecast needs](08_thermal_forcing_and_the_tests_a_forecast_needs.md) | 12 Oct 2026 – | T_off closes budgets, it is not the air; a fair test needs own parameters; stations pick ERA5's lapse rate; sunshine fails at band level; the model warms half as fast as measured; a third of the deep data was being lost; the fair test makes the offset-free forcing central (not robust: H a full sensitivity), and without T_off the massifs want opposite meltwater heat |
| 9 | [Snow and firn from data](09_snow_and_firn_from_data.md) | 14–19 Oct 2026 | the precipitation bound fixes the winter snow but not the snow line; less snow restores Tête Rousse's frozen bed; snow lines and firn maps must pin the firn |
| 10 | [The firn ledger](10_the_firn_ledger.md) | 19–20 Oct 2026 | the ledger's rules reproduce the firn model; fed GloGEM's balance it holds 10–30 m w.e. of firn where observations leave none |
| 11 | [The missing heat](11_the_missing_heat.md) | 19–20 Oct 2026 – | meltwater refreezing in the snow on ice bands (S) repairs the small zone glaciers, but warms the already warm tongues and costs Tête Rousse's frozen bed; GloGEM's full refreezing model was broken and changes nothing; comparisons now per borehole and year |
| 12 | [Sun and shade, burial, and the air](12_sun_shade_burial_air.md) | 20–22 Oct 2026 | orientation fails its rule twice but is kept (three parts per band: error on left-out massifs 2.51 → 1.96 °C, the shade pattern inside sites nearly gone); burial 2–20× too fast and the air 2 °C too cold, neither explains the massif split; known limitations |
| 13 | [The air from above](13_the_air_from_above.md) | 22–23 Oct 2026 – | ERA5's 2 m air lapsed up is 5 °C too cold in winter at the summits; its own free air is right (+0.6 °C), no other product is; in the thermal model the free air with Gilbert's surface cooling fails its rule (2.55 against 1.86 °C): the shallow high firn gets right, the deep ice at the high Mont Blanc sites 2–4 °C too warm |

Decisions made on the data, with whether they followed a rule fixed in advance: [decision log](decision_log.md).

## How chapters are written

Each chapter follows the same five steps:

1. **What we believed** — the state of the model and our assumptions.
2. **What we measured** — the test, and its numbers.
3. **What was wrong** — explained in plain words, with a figure where it helps.
4. **What we changed and why** — with commit hashes (GloGEM and icetemp repos).
5. **What is still open.**

Each chapter starts with a short summary and ends with its sources: plans, notebooks, runs and
commits.

## Figures

Figures sit in the text where the concept is explained:

- **SVG schematics and plots** in `figures/`, made by one script per chapter
  (`figures/make_chNN.py`) or by the notebook named in the chapter. Rerun it to redraw the
  figures when numbers change; for a script, pass a
  folder as argument to also get PNG previews there.
- **Mermaid diagrams** written directly in the Markdown for step-by-step processes. GitHub and the
  VS Code preview draw them automatically.

## Plans

Approved plans are archived in `plans/`, dated, as they were approved.

## Making a PDF

```bash
cd development_history
pandoc README.md 0*.md -o history.pdf --resource-path=. --pdf-engine=xelatex
```

Mermaid blocks need the `mermaid-filter` pandoc filter to appear in the PDF.
