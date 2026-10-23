# WAH-*i* Ocean GUI (`wahi_gui`)

Interactive MATLAB app for designing and optimising hydrophone arrays for whale localisation.

```matlab
addpath("app", "src")
wahi_gui
```

`wahi_gui` adds `../src` to the path automatically when run from the repository.

Workflow:

1. **Array Design** — generate a pattern (dx/dy spacing), edit coordinates in the table, export the CSV.
2. **Environment** — optional bathymetry, sound-speed model (constant / 1-D / NCOM), noise and detection
   settings, clock and survey uncertainty. Check the SNR-vs-range plot.
3. **Calls** — species preset or imported template (vector of points or .wav).
4. **Field & Optimise** — set the field (radius up to 100 km, spacing, depths, sector), *Evaluate field*,
   *Optimise geometry*, *Array position sweep*, report and exports.

*Menu → Save session / Load session* stores every setting together with the array, bathymetry,
sound-speed data and call template. *Menu → Tips* has workflow and interpretation notes.

See the top-level README for file formats and model assumptions.
