# **WAH-*i* App — Installation**

This folder contains standalone installers for the **WAH-*i*** (Widefield
Acoustics Heuristic -- inverse, iterative) MATLAB GUI.

## **Folder structure**

```
app_install/
├── winx64/wahi_gui_web.exe.zip   → Windows 64-bit web installer
└── mac/wahi_gui_web.zip          → macOS web installer
```

## **Installation notes**

- The application is distributed as a **standalone compiled MATLAB app**.
- **MATLAB Runtime (MCR)** is required to run the app.
- These are web installers. If MATLAB Runtime is **not already installed**, an
  internet connection is required so the installer can download the compatible
  runtime during setup.
- No full MATLAB licence is required to run the installed application.
- For the app to look its best, [install the Lato font](https://fonts.google.com/specimen/Lato) on your computer. 

## **Platform notes**

- **Windows:** Extract `winx64/wahi_gui_web.exe.zip`, run
  `wahi_gui_web.exe`, and follow the on-screen instructions.
- **macOS:** Extract `mac/wahi_gui_web.zip`, open `wahi_gui_web.app`, and
  follow the installer. You may need to approve it in *System Settings →
  Privacy & Security* on first launch.

## Installed Apps

### Windows 10 Installation

![Win10InstallationExample](../img/win_done.PNG)

### MacOS Installation

![macosInstallationExample](../img/mac_done.png)

## Reports

View a run summary by clicking **Report**. A typical report looks like:

![iwah_gui_run__4mics__20260203_100305](../reports/iwah_gui_run__4mics__20260203_100305.png)

## Export Run

Use **Export Run** to save the complete run history as a MAT file. The
`mat/iwah_singleRun_analyses.m` script provides a starting point for detailed
analysis and custom downstream pipelines.

## Flexibility

The app is designed to support a wide range of optimisation strategies, recognising that multiple combinations of parameter choices, manual adjustments, and staged runs can lead to the desired accuracy within a specified field. Rather than enforcing a single workflow, WAH-*i* exposes the full set of algorithmic controls described in the accompanying research, enabling users to explore and refine solutions that best match their experimental constraints and goals. Users are encouraged to develop their own optimisation workflows and downstream analysis pipelines tailored to their specific applications.

## **Caution on Computational Load**

Localisation accuracy improves with increasing grid resolution, larger array aperture, and higher microphone counts; however, these gains come at a substantial computational cost. Finer grids, extended fields, and high-order arrays can significantly increase memory usage and processing time, particularly during optimisation runs. Users should ensure that adequate computational resources are available before attempting large-scale or high-resolution configurations.

All development and primary testing were performed on Apple M1/M2 systems. Limited testing was conducted on an Intel i5–class processor, where performance may be reduced for demanding configurations.

## Disclaimer

This software is provided **as is**, without any warranties of any kind.

All standard academic and software disclaimers apply.

Use at your own risk.
