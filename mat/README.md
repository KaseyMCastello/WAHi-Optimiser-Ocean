# MATLAB Analysis Scripts (`mat/`)

This folder contains the MATLAB scripts used for batch optimisation, post-hoc
analysis, table generation, and manuscript figures associated with the
**WAH-*i*** framework.

The scripts are intended for:

- Reproducing analyses reported in the accompanying paper
- Running batch optimisation experiments
- Analysing and visualising results exported from the WAH-*i* GUI

---

## Before Running

Add `src/` to the MATLAB path, then review the **USER OPTIONS** near the top of
the selected script. In particular, check the input run directory, output
directory, metric, and save switch. `iwah_batch_run_analyses.m` resolves its
default `data/runs/`, `fig/`, `results/`, and `src/` paths from the repository
location, so it can be launched from any MATLAB working directory.

Most figure scripts use:

```matlab
saveFigures = 0;  % inspect figures without writing files
saveFigures = 1;  % export final figures
```

For scripts with a `METRIC` option, run once with `"pos"` and once with `"ang"`
when both positional and angular outputs are required.

## Reproducing the Batch Analysis

1. Run `run_iwah_experiment_batch.m` to generate optimisation histories, or use
   the supplied histories in `../data/runs/`.
2. Run `iwah_batch_run_analyses.m` for each required metric.
3. Run `iwah_geometry_drift_plots_etc.m` for coordinate drift, compactness, and
   power-law diagnostics.
4. Run `economy_gain_table.m` after the positional and angular summary CSV files
   have been generated.
5. Run `plot_canonical_configurations.m` to regenerate the configuration-ID
   reference figure.

---

## Script Overview

### Core experiment control

- **`run_iwah_experiment_batch.m`**

  Runs five repeated optimisations by default for each canonical CSV geometry
  and for random starts with 4, 6, 8, 10, and 12 microphones. Each output MAT
  file contains the optimiser history `H` and run metadata.

### Batch-level analyses

- **`iwah_batch_run_analyses.m`**

  Aggregates repeated runs and produces pass-rate trajectories, average-run
  error distributions, accuracy-versus-distance summaries, a representative
  before/after optimisation map, improvement summaries, and pooled statistics.

- **`economy_gain_table.m`**

  Combines `summary_stats_pos.csv` and `summary_stats_ang.csv` to calculate
  performance gain per microphone and rank configurations by economy.

### Geometry and stability diagnostics

- **`iwah_geometry_drift_plots_etc.m`**

  Produces the configuration-ID drift figure and optimisation diagnostics, and
  exports geometry-compactness and power-law scaling tables.

- **`plot_canonical_configurations.m`**

  Plots the six canonical starting geometries with configuration IDs 1--6 and
  microphone indices. `figureSizePixels` controls the canvas and `saveFigure`
  controls PDF/PNG export.

### Single-run analysis (GUI-export friendly)

- **`iwah_singleRun_analyses.m`**

  Analyses a MAT file exported from the GUI. Set `RUN_FILE`, `AFTER`, and
  `METRIC` in the user-options block to inspect before/after spatial error,
  distance-binned distributions, pass-rate history, and geometry changes.

---

## Outputs

The scripts generate figures and tables that are typically written to:

- `fig/` – publication-quality figures
- `results/` – CSV and LaTeX-ready tables
- `reports/` – run-specific diagnostic plots

Principal figure outputs include:

- `canonical_configurations.pdf`
- `optimisation_demo.pdf`
- `pass_trajectories.pdf`
- `histograms_pos_error.pdf` and `histograms_ang_error.pdf`
- `accuracy_pos_distance.pdf` and `accuracy_ang_distance.pdf`
- `pos_improvements.pdf` and `ang_improvements.pdf`
- `geometry_drifts.pdf` and `geometry_diagnoses.pdf`

Exact paths are configurable in each script. Keep saving disabled while
adjusting styles to avoid overwriting publication outputs unintentionally.

---

## License

These scripts are released under **GPLv3**, consistent with the main WAH-*i*
software licence.
