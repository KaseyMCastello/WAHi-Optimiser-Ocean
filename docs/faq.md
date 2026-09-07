---
layout: default
title: FAQ
permalink: /faq/
nav_order: 6
---

# FAQ

## What does WAH-*i* mean?

The acronym stands for *Widefield Acoustics Heuristic inverse-iterative*. It is a continuation of the [Array WAH project](https://github.com/raviumadi/Array_WAH), where a method is developed to estimate the spatial localisation accuracy of a given microphone array geometry using the time-difference-of-arrival (TDOA) method. This work extends it to answer the question: *If I need a certain level of accuracy in a given spatial region, what is the best geometric configuration?*

## Do I need MATLAB to run WAH-*i*?

**No.**
The deployed application runs using **MATLAB Runtime**, which is **free** and does **not** require a MATLAB licence.
The installer will automatically prompt you to install the correct MATLAB Runtime version if it is not already present on your system.

If you intend to **modify the code or develop new methods**, a full MATLAB installation is required.

---

## Where can I find the research behind WAH-*i*?

WAH-i is directly associated with peer-reviewed research on microphone-array optimisation and acoustic localisation.

A curated list of related publications is provided here:

[Widefield acoustics heuristic: advancing microphone array design for accurate spatial tracking of echolocating bats](https://doi.org/10.1186/s12862-025-02441-4)

[WAH-*i*: Optimising Microphone Array Geometry for Customised Localisation Accuracy](https://doi.org/10.1111/2041-210X.70407)

[BATSY4-PRO: An Open-Source Multichannel Ultrasound Recorder for Field Bioacoustics](https://doi.org/10.1101/2025.08.11.669530)

[ESPERDYNE: A dual-band heterodyne monitor and ultrasound recorder for bioacoustic field surveys](https://doi.org/10.1111/2041-210x.70241)

These publications describe the theoretical background, optimisation heuristics, and evaluation methodology implemented in the software.

---

## Is there a portable recorder or microphone array that can use WAH-*i*-optimised geometries?

**Yes!**
WAH-*i* was developed alongside **BATSY4-PRO**, an open-source, portable multichannel ultrasound recorder designed for field bioacoustics.

Optimised geometries exported from WAH-i (CSV) can be directly adopted or adapted for:
- Custom microphone arms
- Rigid or semi-rigid array frames
- Field-deployable recording systems

More information on Batsy4-Pro can be found here:
[BATSY4-PRO Project](https://github.com/raviumadi/Embedded_Ultrasonics/tree/main/Batsy4-Pro)

---

## Do I need to pay to use WAH-i?
**No.**
WAH-*i* is **free of charge** for all users.

The software is released under the **GNU GPL v3** licence:
- Academic and research use is fully permitted
- Commercial use must comply with GPLv3 terms

There are **no hidden fees, subscriptions, or feature locks**.

---

## How can I support this work?
WAH-*i* is developed and maintained as a self-funded research project.

If you find the tool useful and would like to support continued development of this and all of my other open-source projects, you may:
- [**Buy me a coffee**](https://buymeacoffee.com/raviumadi) ☕
- Share how you used WAH-*i* in your research
- Cite the relevant paper(s) if appropriate

Hearing how the tool is used in real research contexts is genuinely motivating and helps guide future development.

---

## Where can I find more of your research and tools?
A broader collection of open-source tools, hardware projects, and research outputs is available at: [**biosonix**](https://biosonix.io)

This includes:
- Bioacoustics instrumentation
- Localisation and tracking frameworks
- Field-deployable recording systems
- Educational and outreach resources

---

## I have a question that isn’t covered here
If something is unclear, or if you are adapting WAH-*i* for a new application, feel free to:
- Open a **GitHub Issue**
- Reach out via the contact details on [**biosonix**](https://biosonix.io)

Constructive questions and feedback are always welcome.
