# WAH-*i* MATLAB GUI (`wahi_gui`)

This folder contains the **interactive MATLAB graphical user interface** for **WAH-*i***  

(*Widefield Acoustics Heuristic – inverse, iterative*), An array-geometry optimisation tool for designing microphone arrays with customised localisation accuracy.

The GUI provides an end-to-end workflow for:

- Designing and editing microphone geometries
- Evaluating 3D fields of localisation accuracy
- Running forward WAH evaluation and WAH-*i* geometry optimisation
- Exporting optimisation runs for reproducible analysis

---

## Launching the GUI

From MATLAB, add `app/` and `src/` to the path, then run:

```matlab
addpath("app", "src")
wahi_gui
```

The GUI source depends on the core files in `src/`, particularly
`BatCallLocaliser.m`, `iwah_run_optimiser_core.m`, `plotCSV3D.m`, and
`formatLatex.m`.

![iwah_gui_start](../img/iwah_gui_start.png)

*The WAH-i GUI at launch*

## Files Produced by the GUI

- **Export Run** writes a MAT file containing the complete optimiser history
  `H` and `runInfo`. These files can be inspected with
  `../mat/iwah_singleRun_analyses.m`.
- **Export Geometry** writes the current microphone coordinates as CSV for
  reuse, construction, or CAD workflows.
- **Report** creates a summary figure and optionally saves it as PNG.

The GUI exposes signal, grid, accuracy, physical-constraint, and optimiser
settings. Numerical results are conditional on these choices; users should
rerun the optimisation with parameters representative of their species,
recording system, observation volume, and deployment constraints.

## Preferred Font

For the app to look its best, [install the Lato font](https://fonts.google.com/specimen/Lato) on your computer.
