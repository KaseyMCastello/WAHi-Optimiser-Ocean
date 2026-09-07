---
layout: default
title: WAH-i
permalink: /
nav_order: 1
include_slider: true
---

# WAH-*i*-Optimiser

WAH-*i* (Widefield Acoustics Heuristic — *inverse iterative*; pronounced ***waa-hee*** [**/ˈwɑː.hiː/**]) is a MATLAB-based framework and a software tool for designing and optimising 3D microphone-array geometries for TDOA localisation. The GUI supports rapid iteration, evaluation over a user-defined Field of Accuracy (FoA), and iterative optimisation under real-world constraints (arm length, minimum spacing, etc.).

<a class="button-pill button-pill--primary" href="https://doi.org/10.1111/2041-210X.70407" target="_blank" rel="noreferrer">Read the published paper</a>

<details class="citation-card">
  <summary>
    <span>
      <strong>Cite WAH-<em>i</em></strong>
      <small>Published in <em>Methods in Ecology and Evolution</em></small>
    </span>
    <span class="citation-card__chevron" aria-hidden="true"></span>
  </summary>

  <div class="citation-card__body">
    <div class="citation-card__doi">
      <div>
        <span class="citation-card__label">DOI</span>
        <a href="https://doi.org/10.1111/2041-210X.70407" target="_blank" rel="noreferrer">10.1111/2041-210X.70407</a>
      </div>
      <button class="copy-button" type="button" data-copy-value="10.1111/2041-210X.70407">Copy DOI</button>
    </div>

    <div class="citation-card__bibtex">
      <div class="citation-card__heading">
        <span class="citation-card__label">BibTeX</span>
        <button class="copy-button" type="button" data-copy-target="wahi-bibtex">Copy BibTeX</button>
      </div>
      <pre><code id="wahi-bibtex">@article{umadiWAHOptimisingMicrophone2026,
  title   = {WAH-i: Optimising Microphone Array Geometry for Customised Localisation Accuracy},
  author  = {Umadi, Ravi},
  journal = {Methods in Ecology and Evolution},
  year    = {2026},
  doi     = {10.1111/2041-210X.70407},
  url     = {https://doi.org/10.1111/2041-210X.70407}
}</code></pre>
    </div>
  </div>
</details>

<section class="wahi-gallery">
  <h2>WAH-<em>i</em> in action</h2>

  <div class="wahi-slider">
    <button class="wahi-slider-btn prev" aria-label="Previous slide">&#10094;</button>

    <div class="wahi-slider-track">

      <figure class="wahi-slide">
        <img src="{{ '/assets/img/home/opti_run_3.png' | relative_url }}"
             alt="Field of Accuracy example 1">
        <figcaption>WAH<em>i</em> Run: With diagonal left-right octants; 6 mics at 0.75 m max arm-length. After 20 iterations, the optimiser converged to 81% pass rate, with 95% grid points below 30 cm inaccuracy within a field of 4m radius, with a inner exlucsion zone up to 1.5 m. Report >>  </figcaption>
      </figure>

      <figure class="wahi-slide">
        <img src="{{ '/assets/img/home/iwah_report_cross_octants.png' | relative_url }}"
             alt="Field of Accuracy example 2">
        <figcaption>The optimiser generates a report output, showing all comparaitive metrics, before and after optimisation run. </figcaption>
      </figure>

      <figure class="wahi-slide">
        <img src="{{ '/assets/img/home/mic_drift_cross_ocatnat_run.png' | relative_url }}"
             alt="Optimisation run">
        <figcaption>Also, a microphone position log is produced. The user may also export the final microphone configuration as a CSV file.</figcaption>
      </figure>

      <figure class="wahi-slide">
        <img src="{{ '/assets/img/home/optimisation_run.png' | relative_url }}"
             alt="Optimisation run">
        <figcaption>Another optimisation run, with a thinner shell, achieving a 98% pass rate.</figcaption>
      </figure>

    </div>

    <button class="wahi-slider-btn next" aria-label="Next slide">&#10095;</button>
  </div>

</section>

## What WAH-*i* is for?

WAH-*i* is a practical tool for designing, evaluating, and optimising microphone array geometries under realistic localisation constraints. 

At its core, WAH-*i* answers a simple but hard question:

**"Given a limited number of microphones and physical placement constraints, how should they be arranged to maximise reliable localisation performance over a region of space I actually care about?"**

## Designed for real-world array design

WAH-*i* is built for situations where arrays are not idealised and not free to grow arbitrarily:
* limited arm length or mounting radius
* minimum microphone spacing
* near to mid-field accuracy
* asymmetric or task-specific regions of interest

Instead of optimising abstract metrics, WAH-*i* evaluates performance over a user-defined Field of Accuracy (FoA) and explicitly quantifies where localisation succeeds or fails.

<section class="wahi-gallery">
  <h2>Field View</h2>

  <div class="wahi-slider">
    <button class="wahi-slider-btn prev" aria-label="Previous slide">&#10094;</button>

    <div class="wahi-slider-track">

     <figure class="wahi-slide">
        <img src="{{ '/assets/img/home/3d_side_1.png' | relative_url }}"
              alt="Report output">
        <figcaption>The chosen octants, post-optimisation view. Explore the field of achieved accuracy via the in-built 3D view tools.</figcaption>
        </figure>

      <figure class="wahi-slide">
      <img src="{{ '/assets/img/home/3d_side_2.png' | relative_url }}"
            alt="Report output">
      <figcaption>Any combination of octants can be chosen for a complete customised field of accuracy (FoA). A combination of inner and outer radius values further allows narrowing of the depth of the FoA, allowing a shell-like area of any depth.</figcaption>
      </figure>

      <figure class="wahi-slide">
      <img src="{{ '/assets/img/home/3d_side_3.png' | relative_url }}"
            alt="Report output">
      <figcaption>Top-view of the chosen octants.</figcaption>
      </figure>

      </div>

    <button class="wahi-slider-btn next" aria-label="Next slide">&#10095;</button>
  </div>
</section>

## What makes WAH-*i* different

Unlike closed-form array designs or single-shot optimisation tools, WAH-*i*:
* treats array design as an iterative, constraint-driven process
* exposes the geometry–performance trade-off explicitly
* supports staged optimisation (coarse → fine)
* prioritises robust pass rate rather than best-case accuracy

The optimiser perturbs microphone coordinates directly and evaluates performance using the same localisation pipeline that would be used in practice. This ensures that improvements are meaningful, interpretable, and reproducible.

## Typical use cases

WAH-*i* is particularly useful when you need to:
* design a portable or field-deployable microphone array
* compare candidate geometries under identical signal and noise conditions
* refine an existing array rather than start from scratch
* explore how performance changes with geometry, coverage, and constraints
* generate reproducible optimisation trajectories for research or validation

Field researchers can often stop once a geometry meets their target pass rate.
Method developers can go further—modifying the localisation backend, spatial sampling strategy, or optimisation logic directly in code.

## From exploration to deployment

The GUI supports rapid, single-geometry exploration, while the accompanying batch scripts enable large-scale comparative analyses across many geometries and optimisation runs. Together, they form a complete workflow from concept → optimisation → evaluation → deployment.

## What you’ll find here

- [**Install**]({{ '/install/' | relative_url }}): MATLAB Runtime + packaged app installation, and development setup.
- [**Download**]({{ '/download/' | relative_url }}): get the installable packages for your OS.
- [**Quickstart**]({{ '/quickstart/' | relative_url }}): the shortest path from “installed” to “first optimisation run”.
- [**GUI walkthrough**]({{ '/gui/' | relative_url }}): what each panel does, and a suggested workflow.
- [**Optimisation tips**]({{ '/optimisation-tips/' | relative_url }}): staged optimisation, interpreting metrics, and avoiding common traps.
- [**FAQ**]({{ '/faq/' | relative_url }}): deployment, paths, performance, and troubleshooting.

## Recommended first steps
1. Read [**Install**]({{ '/install/' | relative_url }})
2. Follow [**Quickstart**]({{ '/quickstart/' | relative_url }})
3. Use [**Optimisation tips**]({{ '/optimisation-tips/' | relative_url }}) when pushing for high pass rates

<!-- <hr>
<h2>Support my open science projects</h2>

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
     padding: 5px 8px;
     background-color: #00c351cf;
     color: #000;
     font-weight: 600;
     border-radius: 6px;
     text-decoration: none;
     border: 1px solid #00e6b4;
   ">
  ☕ Buy me a coffee
</a>
</p> -->

