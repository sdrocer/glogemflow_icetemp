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
| 7 | [The snow, the comparison and the calibration target](07_snow_comparison_target.md) | 30 Sep 2026 – | each reading at its date; the winter snow on ice barely insulated (a solver bug); the fix lowers the RMSE in all 32 runs; calibration target next |

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
