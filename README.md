# WAH-*i* -  Microphone Array Geometry Optimisation Algorithm

**Companion repository for the research paper:**

**WAH-*i*: Optimising Microphone Array Geometry for Customised Localisation Accuracy**

[Read the Paper](https://doi.org/10.64898/2026.02.07.704547)

This repository contains the complete software, analyses, datasets, and application developed as part of the WAH-*i* research programme. The project introduces a principled optimisation framework for designing microphone array geometries that achieve a desired spatial field of localisation accuracy.

 

### [Quickstart and Usage Guide](https://raviumadi.github.io/WAHi-Optimiser/quickstart/)



---

## Overview

**WAH-*i*** (Widefield Acoustics Heuristic — inverse, iterative) is an optimisation framework that:

- Directly targets **volumetric localisation accuracy**, rather than heuristic geometric criteria
- Supports **task-specific accuracy thresholds** and pass-rate objectives
- Enables **systematic comparison of array geometries** under realistic constraints
- Scales from **small, portable arrays** to larger multichannel systems

The repository includes:

- A full-featured **graphical user interface**
- **Standalone installable application packages**
- The batch-run pipeline and analysis scripts used in the paper
- The canonical starting geometries and example optimised geometries
- The optimisation histories, figures, and exported summary tables reported in the paper

See `app/` for the MATLAB GUI source, `app_install/` for standalone installers,
and `mat/` for the reproducible analysis workflow.

## Requirements

- MATLAB R2023b or newer
  - Optimization Toolbox
  - Signal Processing Toolbox
  - Statistics and Machine Learning Toolbox
  - Communications Toolbox

The standalone applications do not require a MATLAB licence, but do require the
compatible MATLAB Runtime. See `app_install/README.md`.

---

## **Graphical User Interface (GUI)**

The interactive WAH-*i* GUI allows full control over optimisation parameters, constraints, and diagnostics.

📁 **Relevant folders**

- app/ — GUI source and documentation
- app_install/ — standalone, installable application packages

![iwah_gui_done_dense_grid](img/iwah_gui_done_dense_grid.png)

The GUI supports:

- Iterative optimisation with live diagnostics
- Manual geometry refinement and restart strategies
- Export to CSV (CAD workflows) and MAT (further analysis)
- Automated report generation

## Batch Runs and Analyses

Batch experiments and analysis scripts used in the paper are provided for full reproducibility.

📁 **MATLAB scripts**

- `mat/run_iwah_experiment_batch.m` — generate repeated canonical and random-start runs
- `mat/iwah_batch_run_analyses.m` — batch summaries and the principal error figures
- `mat/iwah_geometry_drift_plots_etc.m` — geometry drift, compactness, and scaling diagnostics
- `mat/economy_gain_table.m` — economy and responsiveness summary
- `mat/plot_canonical_configurations.m` — canonical-geometry reference figure
- `mat/iwah_singleRun_analyses.m` — inspect one GUI-exported run

📁 **Data**

- data/runs/ — batch optimisation runs used in the manuscript
- data/single_runs/ — example GUI-driven optimisation runs

📁 **Outputs**

- fig/ — publication figures
- reports/ — GUI-generated summary reports
- results/ — processed tables and statistics

## Array Geometries

📁 **configs/**

Contains CSV files defining microphone array geometries used throughout the study, including:

- Canonical starting configurations: pyramid, tetrahedron, octahedron, double
  tetrahedron, pentagonal antiprism, and icosahedron
- Example high-pass-rate four-microphone geometries exported from WAH-*i*

Random starting configurations are generated programmatically during batch runs and
are therefore not stored as fixed CSV files. CSV geometries can be loaded directly
into the GUI or used as starting conditions for batch experiments.

The canonical configuration IDs used in the manuscript are:

| ID | Starting geometry | Microphones |
|---:|---|---:|
| 1 | Pyramid | 4 |
| 2 | Tetrahedron | 4 |
| 3 | Octahedron | 6 |
| 4 | Double tetrahedron | 8 |
| 5 | Pentagonal antiprism | 10 |
| 6 | Icosahedron | 12 |

## Source Code

📁 **src/**

Core algorithmic and utility components, including:

- Optimisation engine
- Localisation model
- Plotting and export utilities

## License

- **Software:** GNU General Public License v3 (GPLv3)
- **Documentation:** CC BY-NC-ND 4.0

See LICENSE.md for details.

## Disclaimer

This software and associated analyses are provided **as-is**, without warranty of any kind.

No guarantee is made regarding suitability for any particular application or experimental setup.

Users are responsible for validating results in their own deployment contexts.

## Support & Acknowledgement

This work is self-funded and developed as part of open, independent research.

If you find WAH-*i* useful and would like to support further development:

☕ **Buy me a coffee:** https://buymeacoffee.com/raviumadi

You can find my other research and open tools at:

🌍 **https://biosonix.io**

---

Whenever possible, consider supporting or promoting **science education in developing countries**.

<hr>
<h2>☕ Support my open science projects</h2>


<p>
This project is developed and maintained independently as part of my open research work.
If you find it useful and would like to support continued development, documentation,
and free public releases, consider buying me a coffee.
</p>
<p>
<a href="https://buymeacoffee.com/raviumadi"
   target="_blank"
   style="
     display: inline-block;
     padding: 10px 16px;
     background-color: #FFDD00;
     color: #000;
     font-weight: 600;
     border-radius: 6px;
     text-decoration: none;
     border: 1px solid #e6c800;
   ">
  ☕ Buy me a coffee
</a>
</p>

<p><em>All tools remain free for academic and research use.</em></p>
